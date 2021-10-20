#' @importFrom IRanges IntegerList
#' @importFrom S4Vectors DataFrame List
#' @importFrom futile.logger flog.info
#' @importFrom iterators nextElem idiv
#' @importFrom methods callNextMethod new validObject
#' @importFrom stats complete.cases
NULL

#' Transcription start site (TSS) hidden Markov model (HMM).
#'
#' @description
#'
#' Find promoter regions in nascent RNA data using Andre Martin's 2014 hidden
#' Markov model for PRO-cap reads, and Pariksheet Nanda's model extension to
#' replace the background reference with a full gene transcription assay.
#'
#' @details
#'
#' The model takes inputs of raw nascent RNA counts, aggregates 10 basepair
#' tiles by maximum value, and categorizes them into 3 types of observations:
#' \enumerate{
#' \item no signal (TAP+ == 0)
#' \item enriched  (TAP+ >  TAP-)
#' \item depleted  (TAP- >= TAP+ > 0)
#' }
#' ...and searches for 3 groups of hidden states:
#' \enumerate{
#' \item background,
#' \item peaked TSS regions, and
#' \item non-peaked TSS regions
#' }
#'
#' Both the peaked and non-peaked TSS regions each require 3 states to describe
#' them:
#' \itemize{
#' \item Non-peaked regions flanked by single low intensity transitions and a
#' moderately long intensity center.
#' \item Peaked regions flanked by one-or-more low intensity transitions and a
#' short intense center.
#' }
#' Therefore, in total the HMM has 7 states for PRO-cap and an additional state
#' for PRO-seq:
#' \describe{
#' \item{B}{Background}
#' \item{N1}{Non-peaked TSS transition state}
#' \item{N2}{Non-peaked TSS repeating state}
#' \item{N3}{Non-peaked TSS transition state}
#' \item{P1}{Peaked TSS moderate signal}
#' \item{P2}{Peaked TSS high signal}
#' \item{P3}{Peaked TSS moderate signal}
#' \item{GB}{Gene-body signal if background reference data provides it}
#' }
#'
#' Finally, after the hidden states are obtained from Viterbi decoding, the
#' hidden states are converted back into stranded GRanges.
#'
#' @section Constructor:
#'
#' model <- TSSHMM()
#' model <- TSSHMM(bg_genebody = TRUE)
#'
#' `TSSHMM()` returns a default model object to be trained and then used
#' for decoding.  An optional argument `bg_genebody` uses PRO-seq for
#' background instead of PRO-cap; the default is `bg_genebody = FALSE`.
#'
#' @name TSSHMM-class
#' @aliases TSSHMM
#' @exportClass TSSHMM
setClass(
    "TSSHMM",
    slots = c(
        bg_genebody = "logical",
        trans = "matrix",
        emis = "matrix",
        emis_tied = "integer"
    )
)
setMethod(
    "initialize",
    signature = "TSSHMM",
    definition = function(.Object, ..., bg_genebody = FALSE,
                          trans = matrix(nrow = 0, ncol = 0),
                          emis = matrix(nrow = 0, ncol = 0),
                          emis_tied = vector("integer")) {
        ## Initialize the model with default values.
        bg_genebody(.Object) <- bg_genebody
        ## Apply any user parameters.
        parameters(.Object) <-
            list(trans = trans, emis = emis, emis_tied = emis_tied)
        validObject(.Object)
        ## Boilerplate per ?initialize to pass '...' arguments to future
        ## subclasses.
        .Object <-
            callNextMethod(.Object, ..., trans = transitions(.Object),
                           emis = emissions(.Object),
                           emis_tied = emissions_tied(.Object))
        .Object
    }
)

#' @export
TSSHMM <- function(bg_genebody = FALSE, trans = NULL, emis = NULL,
                   emis_tied = NULL) {
    ## Apply any user parameters.  new() requires us to match the inputs to the
    ## slot types.  Setting an empty value is treated as a noop.
    if (is.null(trans))
        trans <- matrix(nrow = 0, ncol = 0)
    if (is.null(emis))
        emis <- matrix(nrow = 0, ncol = 0)
    if (is.null(emis_tied))
        emis_tied <- vector("integer")
    new("TSSHMM",
        bg_genebody = bg_genebody,
        trans = trans,
        emis = emis,
        emis_tied = emis_tied)
}

setValidity(
    "TSSHMM",
    function(object) {
        is_valid <-
            .Call(C_is_model_valid,
                  PACKAGE = "tsshmm",
                  transitions(object),
                  emissions(object),
                  emissions_tied(object),
                  start(object))
        if (! is_valid) {
            return("Error: Could not generate TSSHMM() due to C errors above.")
        }
        TRUE
    })

#' @rdname TSSHMM-class
#' @section Accessors:
#'
#' `dim(x)` returns the number of states and emissions.
#'
#' `parameters(model)`, `parameters(model) <- list(trans = ..., emis = ...,
#' emis_tied = ..., bg_genebody = ...)` gets or sets 2 matrices: the transition
#' state probability matrix, and the discrete observable emission probability
#' matrix.
#'
#' After training a model, you may wish to save the parameters to later reload
#' to eliminiate retraining the model in the future.
#'
#' @param x S4 object of a pre-designed hidden Markov model.
#' @exportMethod dim
setMethod(
    "dim",
    signature = "TSSHMM",
    definition = function(x) {
        dim <- dim(emissions(x))
        names(dim) <- c("states", "emissions")
        dim
    }
)

## Internal method.
setGeneric("bg_genebody",
           function(x) standardGeneric("bg_genebody"))
setMethod(
    "bg_genebody",
    signature = "TSSHMM",
    definition = function(x) x@bg_genebody
)

## Internal method.
setGeneric("bg_genebody<-",
           function(x, value) standardGeneric("bg_genebody<-"))
setMethod(
    "bg_genebody<-",
    signature = c("TSSHMM", "logical"),
    definition = function(x, value) {
        stopifnot(length(value) == 1L)
        ## These states must match the enum order in src/models.c
        names_trans <- c("B", "N1", "N2", "N3", "P1", "P2", "P3", "GB")
        names_emis <- c("no_signal", "enriched", "depleted")
        n_states <- length(names_trans)
        n_emis <- length(names_emis)
        emis <- ##   no signal  enriched  depleted
            matrix(c(0.90,      0.05,     0.05, # B
                     0.09,      0.90,     0.01, # N1
                     0.09,      0.90,     0.01, # N2 - tied to N1
                     0.09,      0.90,     0.01, # N3 - tied to N1
                     0.10,      0.45,     0.45, # P1
                     0.10,      0.45,     0.45, # P2 - tied to P1
                     0.10,      0.45,     0.45, # P3 - tied to P1
                     0.05,      0.05,     0.90  # GB
                     ),
                   nrow = n_states, byrow = TRUE,
                   dimnames = list(names_trans, names_emis))
        stopifnot(all(zapsmall(rowSums(emis)) == 1L))
        ## 0 is the untied sentinel value.  Tied values are indices of the
        ## state they are tied to.
        N1 <- 2L
        P1 <- 5L
        ##             B   N1  N2  N3  P2  P2  P3  GB
        emis_tied <- c(0L, N1, N1, N1, P1, P1, P1, 0L)
        trans <- ##  B      N1     N2   N3   P1     P2   P3    GB
            matrix(c(0.990, 0.005, NA , NA , 0.005, NA , NA  , NA   , # B
                     NA   , NA   , 1.0, NA , NA   , NA , NA  , NA   , # N1
                     NA   , NA   , 0.5, 0.5, NA   , NA , NA  , NA   , # N2
                     0.250, NA   , NA , NA , NA   , NA , NA  , 0.750, # N3
                     NA   , NA   , NA , NA , 0.500, 0.5, NA  , NA   , # P1
                     NA   , NA   , NA , NA , 0.450, 0.1, 0.45, NA   , # P2
                     0.125, NA   , NA , NA , NA   , NA , 0.50, 0.375, # P3
                     0.010, NA   , NA , NA , NA   , NA , NA  , 0.990  # GB
                     ),
                   nrow = n_states, byrow = TRUE,
                   dimnames = list(names_trans, names_trans))
        stopifnot(all(zapsmall(rowSums(trans, na.rm = TRUE)) == 1L))
        if (! value) { # Remove gene body state.
            emis_tied <- emis_tied[-nrow(emis)]
            emis <- emis[-nrow(emis), ]
            trans <- trans[-nrow(trans), -ncol(trans)]
            n_states <- nrow(trans)
            B <- 1L
            N3 <- 4L
            P3 <- 7L
            trans[N3, B] <- 1.0
            trans[P3, c(B, P3)] <- c(0.5, 0.5)
            stopifnot(all(zapsmall(rowSums(trans, na.rm = TRUE)) == 1L))
        }
        ## Populate the model.
        transitions(x) <- trans
        emissions(x) <- emis
        emissions_tied(x) <- emis_tied
        x@bg_genebody <- value
        x
    }
)

## Internal method.
setGeneric("transitions",
           function(x) standardGeneric("transitions"))
setMethod(
    "transitions",
    signature = "TSSHMM",
    definition = function(x) x@trans
)

## Internal method.
setGeneric("transitions<-",
           function(x, value) standardGeneric("transitions<-"))
setMethod(
    "transitions<-",
    signature = c("TSSHMM", "matrix"),
    definition = function(x, value) {
        ## Sanity check.
        stopifnot(all(zapsmall(rowSums(value, na.rm = TRUE)) == 1L))
        ## When first instantiating the object, accept the value as is.
        if (identical(dim(x@trans), c(0L, 0L))) {
            x@trans <- value
            return(x)
        }
        ## Empty input is a noop for skipped user parameter instantiation.
        if (identical(value, matrix(nrow = 0, ncol = 0))) {
            return(x)
        }
        ## Otherwise make sure the dimensions and NA values are consistent.
        stopifnot(all(dim(x@trans) == dim(value)))
        na_idx <- function(x) which(is.na(x))
        stopifnot(all(na_idx(x@trans) == na_idx(value)))
        ## Preserve dimnames.
        dimnames <- dimnames(x@trans)
        x@trans <- value
        dimnames(x@trans) <- dimnames
        x
    }
)

## Internal method.
setGeneric("emissions",
           function(x) standardGeneric("emissions"))
setMethod(
    "emissions",
    signature = "TSSHMM",
    definition = function(x) x@emis
)

## Internal method.
setGeneric("emissions<-",
           function(x, value) standardGeneric("emissions<-"))
setMethod(
    "emissions<-",
    signature = c("TSSHMM", "matrix"),
    definition = function(x, value) {
        ## Sanity check.
        stopifnot(all(zapsmall(rowSums(value)) == 1L))
        ## ## When first instantiating the object, accept the value as is.
        if (identical(dim(x@emis), c(0L, 0L))) {
            x@emis <- value
            return(x)
        }
        ## Empty input is a noop for skipped user parameter instantiation.
        if (identical(value, matrix(nrow = 0, ncol = 0))) {
            return(x)
        }
        ## Otherwise make sure the dimensions are consistent and match tied
        ## emissions.
        stopifnot(all(dim(x@emis) == dim(value)))
        stopifnot(! identical(emissions_tied(x), vector("integer")))
        emis_tied <- emissions_tied(x)
        emis_tied_to <- unique(emis_tied)
        emis_tied_to <- emis_tied_to[emis_tied_to > 0L]
        for (i in seq_along(emis_tied_to)) {
            expected <- value[emis_tied_to[i], ]
            stopifnot(all(apply(value[emis_tied == emis_tied_to[i], ], 1,
                                function(row) all(row == expected))))
        }
        ## Preserve dimnames.
        dimnames <- dimnames(x@emis)
        x@emis <- value
        dimnames(x@emis) <- dimnames
        x
    }
)

## Internal method.
setGeneric("emissions_tied",
           function(x) standardGeneric("emissions_tied"))
setMethod(
    "emissions_tied",
    signature = "TSSHMM",
    definition = function(x) x@emis_tied
)

## Internal method.  Requires transitions to be instantiated!
setGeneric("emissions_tied<-",
           function(x, value) standardGeneric("emissions_tied<-"))
setMethod(
    "emissions_tied<-",
    signature = c("TSSHMM", "integer"),
    definition = function(x, value) {
        ## Empty input is a noop for skipped user parameter instantiation.
        if (identical(value, vector("integer"))) {
            return(x)
        }
        ## Otherwise make sure the dimensions are consistent and the domain
        ## matches the number of transition states.
        stopifnot(length(value) == nrow(transitions(x)))
        stopifnot(all(value >= 0L))
        stopifnot(all(value <= length(value)))
        x@emis_tied <- value
        x
    }
)

## Internal method.  Requires transitions to be instantiated!
setMethod(
    "start",
    signature = "TSSHMM",
    definition = function(x) {
        trans <- transitions(x)
        start <- vector("numeric", nrow(trans))
        B <- 1
        N1 <- 2
        P1 <- 5
        vec <- c(trans[B, N1], trans[B, P1])
        start[N1] <- vec[1] / sum(vec)
        start[P1] <- vec[2] / sum(vec)
        stopifnot(zapsmall(sum(start)) == 1L)
        start
    }
)

#' @rdname TSSHMM-class
setGeneric("parameters", function(x) standardGeneric("parameters"))
#' @rdname TSSHMM-class
#' @exportMethod parameters
setMethod(
    "parameters",
    signature = "TSSHMM",
    definition = function(x) {
        list(trans = transitions(x),
             emis = emissions(x),
             emis_tied = emissions_tied(x),
             bg_genebody = bg_genebody(x))
    }
)

#' @rdname TSSHMM-class
setGeneric("parameters<-",
           function(x, value) standardGeneric("parameters<-"))
#' @rdname TSSHMM-class
#' @param value named list containing 2 model matrices: the transition state
#'     probability matrix named "trans", and the discrete observable emission
#'     probability matrix named "emis".
#' @exportMethod "parameters<-"
setMethod(
    "parameters<-",
    signature = "TSSHMM",
    definition = function(x, value) {
        transitions(x) <- value[["trans"]]
        emissions_tied(x) <- value[["emis_tied"]]
        emissions(x) <- value[["emis"]]
        x
    }
)

#' @rdname TSSHMM-class
#' @section Displaying:
#'
#' `show(model)` summarizes `parameters(model)` and the parameter dimesions,
#' namely, the number of transition states and number of discrete observable
#' emissions.
#'
#' @param object S4 object to display.
setMethod(
    "show",
    signature = "TSSHMM",
    definition = function(object) {
        params <- parameters(object)
        dim <- dim(object)
        cat(paste(as.character(class(object)), "object with",
                  dim[1], "hidden states and",
                  dim[2], "emissions:\n"))
        cat("Transition matrix:\n")
        print(params$trans)
        cat("Emission matrix:\n")
        print(params$emis)
        invisible() # Documented return for show().
    }
)

params_idx <- function(x) {
    params <- parameters(x)
    emis_tied <- emissions_tied(x)
    emis_tied[emis_tied == 0] <- which(emis_tied == 0)
    dim <- dim(x)
    n_emis <- prod(dim)
    list(trans = which(! is.na(params$trans)),
         emis = matrix(1:n_emis,
                       nrow = dim[1],
                       ncol = dim[2])[unique(emis_tied), ])
}

params_names <- function(mat, idx) {
    from <-
        matrix(rownames(mat),
               nrow = nrow(mat),
               ncol = ncol(mat))[idx]
    to <-
        matrix(colnames(mat),
               byrow = TRUE,
               nrow = nrow(mat),
               ncol = ncol(mat))[idx]
    ## Use a dot as the separator, because DataFrame() column names don't allow
    ## symbols like - or : unfortunately.
    paste0(from, ".", to)
}

df_updates <- function(model, n) {
    params <- parameters(model)
    idx <- params_idx(model)
    names_trans <- params_names(params$trans, idx$trans)
    names_emis <- params_names(params$emis, idx$emis)
    cols <- c("samples", names_trans, names_emis)
    ## Allocation with -1.0 uninitialized probability values.
    matrix(-1.0, nrow = n, ncol = length(cols), dimnames = list(NULL, cols))
}

#' @rdname TSSHMM-class
#' @section Model Training:
#'
#' obs <- encode_obs(signal, background)
#' obs <- sample(obs) # shuffle
#' model <- train(model, obs)
#'
#' `encode_obs(signal, background)` returns an `IntegerList` of encoded
#' observation states for training.  The arguments, `signal` and `bg`, are
#' stranded, single base `GRanges` with integer scores.
#'
#' The maximum length for each element in `obs` is not configurable.  The
#' maximum length has been set to the upper limit of the Baum-Welch matrix
#' maximum size.  Using a lower maximum length reduces training accuracy,
#' because hidden states such as the background (B) and gene body (G) with
#' their >= 0.99 probability require long contiguous sequences for more
#' accurate training.
#'
#' `train(model, obs)` returns the `model` with updated parameters of
#' transition and emission probabilities.  The argument, `obs` produced by the
#' `encode_obs()` function is an `IntegerList` of encoded observation states.
#' One must shuffle the observations before running training.
#'
#' After training, you may wish to save the model parameters so that they can
#' be reloaded in a later R session as explained in the examples.  Note that
#' the time used for encoding the observations typically takes a magnitude
#' longer than the actual training.
#'
#' To train the model, the sparse reads of the signal and background need to be
#' exploded into dense encoded windows categorized as either enriched, depleted
#' or no-read observations, which are then processed by the Baum-Welch EM
#' algorithm to update the model transition and emission probabilities.
#'
#' The training data are divided into obs large enough to provide sufficient
#' samples for the calculating the background state transitions, but small
#' enough to work within numerical precision limits and to make the training
#' process more observable.  On each batch, the input training data is
#' transformed into dense training observations and then the Baum-Welch
#' algorithm is run.  After each batch, the model state is shown alongside a
#' time estimate to complete training if the logging level has not been reduced
#' from the INFO level.
#'
#' @param signal Stranded, single base \code{GRanges} with integer score.
#' @param bg Stranded, single base \code{GRanges} with integer score.
#' @param nrow Integer, rows to convert at a time.  Primarily used for unit
#'     tests.  Lowering or raising would use would use less and more memory,
#'     respectively, because the encoding process suboptimally requires 70x
#'     more memory than the input data.
#' @export
encode_obs <- function(signal, bg, nrow = 1e3) {
    check_valid_hmm_reads(signal)
    check_valid_hmm_reads(bg)
    ## Train using both strands.  Use range() to estimate number of bases
    ## to feed for training.
    range <- range(c(signal, bg))
    width <- 1e5
    flog.info(sprintf("Split %g bases into %d windows",
                      sum(lengths(range)), width))
    ## Choose an initial batch size of most 100kb width windows and a small
    ## number of sequence rows and measure memory use for setting up this
    ## training dataset.
    regions <- unlist(tile(range, width = width), use.names = FALSE)
    n_batches <- ceiling(length(regions) / nrow)
    flog.info(sprintf(paste(
        "Creating %d obs with up to %d rows of windows each"),
        n_batches, nrow))
    obs <- IntegerList()
    completed <- 0
    t_diff <- 0
    i <- 0
    iterator <- idiv(length(regions), chunkSize = nrow)

    while (TRUE) {
        i <- i + 1
        finished <- FALSE
        tryCatch(chunk <- nextElem(iterator),
                 error = function(e) {
                     if (e$message == "StopIteration") {
                         finished <<- TRUE
                     } else { ## Unhandled error.
                         stop(e) # nocov
                     }
                 })
        if (finished) {
            flog.info("Finished creating obs!")
            break
        }
        ## Begin measure time used for generating this batch of training data.
        t_start <- Sys.time()
        gr <- regions[(completed+1):(completed+chunk)]
        completed <- completed + chunk
        flog.info(
            sprintf("%4d: Tiling and encoding %d regions for training",
                    i, length(gr)))
        windows <- mapply(tile_with_rev,
                          x = as(gr, "GRangesList"),
                          rev = as.vector(strand(gr) == "-"),
                          MoreArgs = list(width = 10))
        ## This encode() line uses a massive 70x peak memory than final
        ## result according to peakRAM::peakRAM() and could stand to be
        ## optimized.  It's not clear how much of the large memory use is
        ## from running the unlist(List(...)).  In any case, an optimal way
        ## to address this would be to change the preceding function
        ## tile_with_rev() to natively generate GRanges using vectorized
        ## rev input to eliminate wrapping the function with mapply() which
        ## produces the undesirable list output.
        obs <- append(obs, encode(signal, bg, unlist(List(windows))))
        ## End measure time used for generating this batch of data.
        t_end <- Sys.time()
        t_diff <- difftime(t_end, t_start, units = "secs")
        t_diff <- t_diff + t_diff
        ## Report batch generation speed in windows per second, the batch
        ## number and remaining total obs, time for processing this
        ## batch, ETA to process remaining obs, total time elapsed, and
        ## total time of since beginning of training.
        elapsed_mins <- as.numeric(t_diff, units = "mins")
        rate_mins <- completed / elapsed_mins
        remaining <- length(regions) - completed
        eta_mins <- remaining / rate_mins
        flog.info(sprintf(paste("%4d: This batch: %.0f secs",
                                "Elapsed: %.1f mins",
                                "Completed: %d/%d regions",
                                "Rate: %.1f regions/min",
                                "ETA: %.0f mins",
                                sep = ", "),
                          i, as.numeric(t_diff, units = "secs"),
                          elapsed_mins,
                          completed, length(regions),
                          rate_mins,
                          eta_mins))
        ## Repeat until completing a single pass of all data.
    }
    obs
}

#' @rdname TSSHMM-class
#' @param model S4 object of a pre-designed hidden Markov model.
#' @param obs IntegerList of encoded observation windows.
setGeneric("train", function(model, obs) standardGeneric("train"))
#' @rdname TSSHMM-class
#'
#' @exportMethod train
setMethod(
    "train",
    signature = c("TSSHMM", "IntegerList"),
    definition = function(model, obs) {
        if (! length(obs)) {
            return(model)
        }

        flog.info(
            sprintf("%4d: Model initial transition and emission matrices:", 0))
        flog.info(sprintf("%s", as(model, "character")))

        ## Begin measure time used for generating this batch of training
        ## data.
        t_start <- Sys.time()
        flog.info("Running Baum-Welch")
        ret <-
            .Call(C_train,
                  PACKAGE = "tsshmm",
                  unlist(obs, use.names = FALSE),
                  lengths(obs),
                  transitions(model),
                  emissions(model),
                  emissions_tied(model),
                  start(model))
        converged <- ret[[1]]
        transitions(model) <- ret[[2]]
        emissions(model) <- ret[[3]]
        ## End measure time used for training the model on this batch of
        ## data.
        t_end <- Sys.time()
        t_diff <- difftime(t_end, t_start, units = "secs")

        ## Model change.
        flog.info(sprintf("Converged? %d", converged))
        flog.info("Model transition and emission matrices:")
        flog.info(sprintf("%s", as(model, "character")))

        ## Report training speed in windows per second, and total time elapsed.
        elapsed_mins <- as.numeric(t_diff, units = "mins")
        rate_mins <- length(obs) / elapsed_mins
        flog.info(sprintf(paste("Elapsed: %.1f mins",
                                "Contigs: %d",
                                "Rate: %.1f contigs/min",
                                sep = ", "),
                          as.numeric(t_diff, units = "mins"),
                          length(obs),
                          rate_mins))
        flog.info("Finished training!")
        model
    }
)

stranded <- function(x) split(x, strand(x))[-3] # -3 removes "*" strand.

#' @rdname TSSHMM-class
setGeneric("viterbi",
           function(model, signal, bg, ...) standardGeneric("viterbi"))
#' @rdname TSSHMM-class
#' @param ... Optional arguments passed on from generics to methods.
#' @param tol Relative change to HMM steady state used to calculate
#'     intra-promoter distance.  Preference should be given to
#'     adjusting the `n` windows parameter before adjusting this.
#' @param n Maximum number of background windows after a promoter,
#'     beyond which one is certain that the promoter has ended.  See
#'     the appendix in the vignette to understand the theory behind
#'     using `n` and `tol` in this calculation.  `n` is the upper
#'     limit; the actual number of background windows is calculated
#'     from the model.
#' @section Model Evaluation:
#'
#' `viterbi(model, signal, background, tol = 1e-3, n = 200)` returns
#' `GRanges` of active promoter or enhancer regions along with the
#' decoded hidden states for each window.  The `tol` and `n`
#' parameters are passed on to the intra-promoter distance estimation
#' function that add flanking windows to ensure that no cluster of proximal
#' promoters is broken up.
#'
#' Running `flog.threshold(DEBUG)` before running `viterbi` logs additional
#' information of what the function is doing.
#' @exportMethod viterbi
setMethod(
    "viterbi",
    signature = c("TSSHMM", "GRanges", "GRanges"),
    definition = function(model, signal, bg, tol = 1e-3, n = 200) {
        check_valid_hmm_reads(signal)
        check_valid_hmm_reads(bg)
        ranges <- reduce(
            c(signal, bg, flank(c(signal, bg),
                                ## Convert prom_dist() windows to basepairs by
                                ## multiplying by 10.
                                width = 10 * prom_dist(model, tol = tol, n = n),
                                start = FALSE,
                                use.names = FALSE)))
        flog.info(sprintf("Decode %g bases across %d regions",
                          sum(lengths(ranges)), length(ranges)))
        unlist(
            mendoapply(
                viterbi_by_strand, signal = stranded(signal), bg = stranded(bg),
                MoreArgs = list(model = model, ranges = ranges)),
            use.names = FALSE)
    }
)

#' @name TSSHMM-class
#' @rdname TSSHMM-class
#' @section Coercion:
#'
#' `as(model, "character")` or `as.character(model)` compactly pastes
#' transition and emission matrices in the same line for logging.
#'
#' @importFrom Rdpack reprompt
#' @references \insertRef{core_analysis_2014}{tsshmm}
#' @seealso \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4254663/}
#' @examples
#' # Model evaluation.
#' signal <- GRanges(paste0("chr1:", c(100, 110, 200, 300), ":+"))
#' score(signal) <- rep(5L, 4)
#' signal
#' (bg <- GRanges())
#' (model <- TSSHMM())
#' (promoters_peaked <- viterbi(model, signal, bg))
#' # There is only 1 promoter in thise region because the first two signal
#' # values filling 2x 20 bp windows are captured by the HMM as a promoter.
#' # The remaining 2 signal peaks with no surrounding signal are ignored.
#' stopifnot(length(promoters) == 1)
#'
#' # Model training.
#' #train(model, signal, bg)
#' (params <- parameters(model))
#'
#' # Save and reload parameters.
#' \dontrun{save(params, file = "model_trained_on_foo_dataset.RData")
#' # ... in a later R session.
#' load(file = "model_trained_on_foo_dataset.RData")
#' }
#' (model <- TSSHMM())
#' parameters(model) <- params
#' model
NULL

## https://stackoverflow.com/a/5173906
precision_single <- function(x) {
    if (abs(x - round(x)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', as.character(x)),
                       ".", fixed = TRUE)[[1]][[2]])
    } else {
        return(0)
    }
}

precision <- function(x) max(sapply(x[complete.cases(x)], precision_single))

#' Convert matrix to string vectors, substituting NA values with spaces.
#'
#' @param x Matrix which may have NA values.
#' @param collapse characters between matrix columns.
#' @return Character vectors of matrix rows with maximum inferred per column.
sprint_precision_or_na_spaces <- function(x, collapse = "  ") {
    precisions <- apply(x, 2, precision)
    ## Cap precisions to keep the matrices readable.
    precisions <- sapply(precisions, function(x) min(x, 5))
    n <- nrow(x)
    str <- vector("character", length = n)
    for (i in seq(n)) {
        fmt <- sprintf("%%.%sf", precisions)
        nas <- ! complete.cases(x[i, ])
        if (sum(nas)) {
            spaces <- sapply(precisions[nas] + 2,
                             function(x) paste0(rep(" ", x),
                                                collapse = ""))
            fmt[nas] <- spaces
        }
        fmt <- paste(fmt, collapse = collapse)
        str[i] <- do.call(sprintf,
                          c(fmt, lapply(c(x[i, !nas]), mean)))
    }
    str
}

## Reexport the coerce methods below as described in Bioconductor S4 Objects
## lab exercise:
## https://www.bioconductor.org/help/course-materials/2011/AdvancedRFeb2011Seattle/ImplementingS4Objects-lab.pdf # nolint
## Here we use ROxygen tags to manage the NAMESPACE file for us
## instead of hand editing the file as suggested by the above aged exercise.
#' @name coerce
#' @importFrom methods coerce
#' @export
methods::coerce

#' Allow TSSHMM to be converted to character vector for logging.
#'
#' @name coerce
#' @aliases coerce,TSSHMM,character-method
#' @docType methods
setAs(
    "TSSHMM", "character",
    function(from) {
        params <- parameters(from)
        str_trans <- sprint_precision_or_na_spaces(params$trans)
        str_emis <- sprint_precision_or_na_spaces(params$emis)
        n <- nrow(params$trans)
        to <- vector("character", length = n)
        for (i in seq(n)) {
            to[i] <- paste0(str_trans[i], " | ", str_emis[i])
        }
        to
    }
)

#' Allow TSSHMM to be converted to character vector for logging.
#'
#' @param x TSSHMM model to convert to character vector for logging.
#' @return Character vector of column concatenated transition and emission
#'     matrices.
#'
#' @aliases as.character,TSSHMM-method
#' @exportMethod as.character
#' @docType methods
setMethod(
    "as.character",
    "TSSHMM",
    function(x) as(x, "character"))
