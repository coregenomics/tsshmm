#' @importFrom S4Vectors List
#' @importFrom futile.logger flog.info
#' @importFrom iterators nextElem idiv
#' @importFrom methods callNextMethod
#' @importFrom stats complete.cases
NULL

#' Transcription start site (TSS) hidden Markov model (HMM).
#'
#' @description
#'
#' Initialize, train, and decode nascent RNA reads to find transcription start
#' sites using Andre Martin's 2014 hidden Markov model.
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
#' Therefore, in total the HMM has 7 states:
#' \describe{
#' \item{B}{Background}
#' \item{N1}{Non-peaked TSS transition state}
#' \item{N2}{Non-peaked TSS repeating state}
#' \item{N3}{Non-peaked TSS transition state}
#' \item{P1}{Peaked TSS moderate signal}
#' \item{P2}{Peaked TSS high signal}
#' \item{P3}{Peaked TSS moderate signal}
#' }
#'
#' Finally, after the hidden states are obtained from Viterbi decoding, the
#' hidden states are converted back into stranded GRanges.
#'
#' @section Constructor:
#'
#' model <- new("TSSHMM")
#'
#' `new("TSSHMM")` returns a default model object to be trained and then
#' used for decoding.
#'
#' The TSS HMM model is implemented using the General Hidden Markov Model
#' (GHMM) C library.  Due to the complexity of the GHMM's C-interface, the R
#' wrapper provides no setters to modify the number of states or the number
#' of observations of the model.
#'
#' @name TSSHMM-class
#' @aliases TSSHMM
#' @exportClass TSSHMM
setClass("TSSHMM", slots = c(external_pointer = "externalptr"))
setMethod(
    "initialize",
    signature = "TSSHMM",
    definition = function(.Object, ...) {
        ## Boilerplate per ?initialize to pass '...' arguments to future
        ## subclasses.
        .Object <- callNextMethod()
        ## Initialize the model and track the C struct pointer.
        .Call(C_model_tsshmm, PACKAGE = "tsshmm", .Object@external_pointer)
        .Object
    }
)

#' @rdname TSSHMM-class
setGeneric("parameters", function(model) standardGeneric("parameters"))
#' @rdname TSSHMM-class
#' @section Accessors:
#'
#' `parameters(model)`, `parameters(model) <- list(trans = ..., emis = ...)`
#' gets or sets 2 matrices: the transition state probability matrix, and the
#' discrete observable emission probability matrix.
#'
#' After training a model, you may wish to save the parameters to later reload
#' to eliminiate retraining the model in the future.
#' @exportMethod parameters
setMethod(
    "parameters",
    signature = "TSSHMM",
    definition = function(model) {
        n_states <- vector("integer", length = 1);
        n_emis <- vector("integer", length = 1);
        .Call(C_model_sizes, PACKAGE = "tsshmm", n_states, n_emis,
              model@external_pointer)
        trans <- vector("numeric", length = n_states * n_states);
        emis <- vector("numeric", length = n_states * n_emis);
        .Call(C_model_matrices, PACKAGE = "tsshmm", trans, emis,
              model@external_pointer)
        list(trans = t(matrix(trans, nrow = n_states)),
             emis = t(matrix(emis, ncol = n_states)))
    }
)

#' @rdname TSSHMM-class
setGeneric("parameters<-",
           function(model, value) standardGeneric("parameters<-"))
#' @rdname TSSHMM-class
#' @param value named list containing 2 model matrices: the transition state
#'     probability matrix named "trans", and the discrete observable emission
#'     probability matrix named "emis".
#' @exportMethod "parameters<-"
setMethod(
    "parameters<-",
    signature = "TSSHMM",
    definition = function(model, value) {
        .Call(C_model_set_matrices, PACKAGE = "tsshmm", c(t(value$trans)),
              c(t(value$emis)), model@external_pointer)
        model
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
        dim <- dim(params$emis)
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

#' @rdname TSSHMM-class
#' @param model S4 object of a pre-designed hidden Markov model.
#' @param signal Stranded, single base \code{GRanges} with integer score.
#' @param bg Stranded, single base \code{GRanges} with integer score.
setGeneric("train", function(model, signal, bg) standardGeneric("train"))
#' @rdname TSSHMM-class
#' @section Model Training:
#'
#' `train(model, signal, background)` returns the model with trained transition
#' and emission probabilities.  The arguments, `signal` and `bg`, are stranded,
#' single base `GRanges` with integer scores.
#'
#' After training, you may wish to save the parameters so that they can be
#' reloaded in a later session as explained in the examples.
#'
#' To train the model, the sparse reads of the signal and background need to be
#' exploded into dense encoded windows categorized as either enriched, depleted
#' or no-read observations, which are then processed by the Baum-Welch EM
#' algorithm to update the model transition and emission probabilities.
#'
#' The training data are randomized and divided into batches large enough to
#' provide sufficient samples for the calculating the background state
#' transitions, but small enough to work within numerical precision limits and
#' to make the training process more observable.  On each batch, the input
#' training data is transformed into dense training observations and then the
#' Baum-Welch algorithm is run.  After each batch, the model state is shown
#' alongside a time estimate to complete training if the logging level has not
#' been reduced from the INFO level.
#'
#' @exportMethod train
setMethod(
    "train",
    signature = c("TSSHMM", "GRanges", "GRanges"),
    definition = function(model, signal, bg) {
        check_valid_hmm_reads(signal)
        check_valid_hmm_reads(bg)
        ## Train using both strands.  Use range() to estimate number of bases
        ## to feed for training.
        range <- range(c(signal, bg))
        flog.info(sprintf("Train using %g bases", sum(lengths(range))))

        ## Start monitoring memory and time use.
        NULL

        ## Choose an initial batch size of least 10k windows and a small number
        ## of sequence rows and measure memory use for setting up this training
        ## dataset.  Select regions in random order to be encoded.
        width <- 1e5
        set.seed(123)
        regions <- sample(unlist(tile(range, width = width), use.names = FALSE))
        rows <- 1e3
        flog.info(sprintf(paste(
            "Creating %d batches with up to %d rows each"),
            ceiling(length(regions) / rows), rows))
        completed <- 0
        t_elapsed <- 0
        i <- 0
        iterator <- idiv(length(regions), chunkSize = rows)

        flog.info(
            sprintf("%4d: Model initial transition and emission matrices:", 0))
        flog.info(sprintf("%s", as(model, "character")))

        while (TRUE) {
            i <- i + 1
            completed <- FALSE
            tryCatch(chunk <- nextElem(iterator),
                     error = function(e) {
                         if (e$message == "StopIteration") {
                             completed <<- TRUE
                         } else { ## Unhandled error.
                             stop(e)
                         }
                     })
            if (completed) {
                flog.info("Training complete!")
                break
            }
            ## Begin measure time used for generating this batch of training data.
            t_start <- Sys.time()
            gr <- regions[(completed+1):(completed+chunk)]
            completed <- completed + chunk
            flog.info(
                sprintf("%4d: Tiling and encoding %d regions for training", i, length(gr)))
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
            obs <- encode(signal, bg, unlist(List(windows)))
            flog.info(sprintf("%4d: Running Baum-Welch", i))
            converged <- NA
            .Call(C_train, PACKAGE = "tsshmm", converged,
                  model@external_pointer, unlist(obs, use.names = FALSE),
                  lengths(obs))
            ## End measure time used for training the model on this batch of data.
            t_end <- Sys.time()
            t_diff <- difftime(t_end, t_start, units = "secs")
            t_elapsed <- t_elapsed + t_diff

            ## Model change.
            flog.info(sprintf("%4d: Converged? %d", i, converged))
            flog.info(
                sprintf("%4d: Model transition and emission matrices:", i))
            flog.info(sprintf("%s", as(model, "character")))

            ## Adjust the batch size towards, say 80% of free memory, also
            ## reporting the ratio of R training data and C-level Baum-Welch.
            NULL

            ## Report training speed in windows per second, the batch number and
            ## remaining total batches, time for processing this batch, ETA to
            ## process remaining batches, total time elapsed, and total time of
            ## since beginning of training.
            elapsed_mins <- as.numeric(t_elapsed, units = "mins")
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
    }
)

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
#' (model <- new("TSSHMM"))
#' (promoters_peaked <- viterbi(model, signal, bg))
#' # There is only 1 promoter in thise region because the first two signal
#' # values filling 2x 20 bp windows are captured by the HMM as a promoter.
#' # The remaining 2 signal peaks with no surrounding signal are ignored.
#' stopifnot(length(promoters) == 1)
#'
#' # Model training with saving and reloading parameters.
#' #train(model, signal, bg)
#' (params <- parameters(model))
#' \dontrun{
#' save(params, file = "model_trained_on_foo_dataset.RData")
#' # ... in a later R session.
#' load(file = "model_trained_on_foo_dataset.RData")
#' }
#' (model <- new("TSSHMM"))
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
