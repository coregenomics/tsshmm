#' @importFrom S4Vectors List
#' @importFrom iterators nextElem idiv
#' @importFrom methods callNextMethod
#' @importFrom stats complete.cases
NULL

#' Transcription start site (TSS) hidden Markov model (HMM).
#'
#' @description
#'
#' Initialize, train, and decode nascent RNAreads using a pre-designed hidden
#' Markov model.
#'
#' @details
#'
#' The TSS HMM model is implemented using the General Hidden Markov Model
#' (GHMM) C library.  Due to the complexity of the GHMM's C-interface, the R
#' layer provides no setters to modify the model aside from running the train
#' function.
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
#' @section Constructor:
#'
#' model <- new("TSSHMM")
#'
#' `new("TSSHMM")` returns a default model object to be trained and then
#' used for decoding.
#'
#' @section Accessor:
#'
#' `parameters(model)` returns 2 matrices: the transition state probability
#' matrix, and the discrete observable emission probability matrix.
#'
#' @section Displaying:
#'
#' `show(model)` summarizes `parameters(model)` and the parameter dimesions,
#' namely, the number of transition states and number of discrete observable
#' emissions.
#'
#' @section Model Training:
#'
#' `train(model, signal, background)` returns the model with trained transition
#' and emission probabilities.  The arguments, `signal` and `bg`, are stranded,
#' single base `GRanges` with integer scores.
#'
#' @section Model Evaluation:
#'
#' `decode(model, signal, background)` returns `GRanges` of active promoter or
#' enhancer regions along with the decoded hidden states for each window.
#'
#' @section Coercion:
#'
#' `as(model, "character")` compactly pastes transition and emission matrices
#' in the same line for logging.
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

setGeneric("parameters", function(model) standardGeneric("parameters"))
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

setGeneric("train", function(model, signal, bg) standardGeneric("train"))
setMethod(
    "train",
    signature = c("TSSHMM", "GRanges", "GRanges"),
    definition = function(model, signal, bg) {
        ## Train using both strands.  Use range() to estimate number of bases
        ## to feed for training.
        range <- range(c(signal, bg))
        flog.debug(sprintf("Train using %g bases", sum(lengths(range))))

        ## Start monitoring memory and time use.
        NULL

        ## Choose an initial batch size of least 10k windows and a small number
        ## of sequence rows and measure memory use for setting up this training
        ## dataset.  Select regions in random order to be encoded.
        width <- 1e5
        set.seed(123)
        regions <- sample(unlist(tile(range, width = width), use.names = FALSE))
        rows <- 1e3
        flog.debug(sprintf(paste(
            "Creating %d batches with %d rows each",
            "to then maximize optimal rows / batch size to available memory"),
            ceiling(length(regions) / rows), rows))
        completed <- 0
        t_elapsed <- 0
        i <- 0
        iterator <- idiv(length(regions), chunkSize = rows)

        flog.debug(
            sprintf("%4d: Model initial transition and emission matrices:", 0))
        flog.debug(sprintf("%s", as(model, "character")))

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
                flog.debug("Training complete!")
                break
            }
            ## Begin measure time used for generating this batch of training data.
            t_start <- Sys.time()
            gr <- regions[(completed+1):(completed+chunk)]
            completed <- completed + chunk
            flog.debug(
                sprintf("%4d: Tiling and encoding %d regions for training", i, length(gr)))
            windows <- mapply(tile_with_rev,
                              x = as(gr, "GRangesList"),
                              rev = as.vector(strand(gr) == "-"),
                              MoreArgs = list(width = 10))
            ## This encode() line uses a massive 70x peak memory than final
            ## result according to peakRAM::peakRAM() and could stand to be
            ## optimized.  It's not clear how much of the large memory use is
            ## from running the unlist(List(...))
            obs <- encode(signal, bg, unlist(List(windows)))
            flog.debug(sprintf("%4d: Running Baum-Welch", i))
            converged <- NA
            .Call(C_train, PACKAGE = "tsshmm", converged,
                  model@external_pointer, unlist(obs, use.names = FALSE),
                  lengths(obs))
            ## End measure time used for training the model on this batch of data.
            t_end <- Sys.time()
            t_diff <- difftime(t_end, t_start, units = "secs")
            t_elapsed <- t_elapsed + t_diff

            ## Model change.
            flog.debug(sprintf("%4d: Converged? %d", i, converged))
            flog.debug(
                sprintf("%4d: Model transition and emission matrices:", i))
            flog.debug(sprintf("%s", as(model, "character")))

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
            flog.debug(sprintf(paste("%4d: This batch: %.0f secs",
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

## Reexport the `setAs` coercion methods below as described in Bioconductor S4
## Objects lab exercise:
## https://www.bioconductor.org/help/course-materials/2011/AdvancedRFeb2011Seattle/ImplementingS4Objects-lab.pdf # nolint
## Here we use ROxygen tags to manage the NAMESPACE file for us instead of hand
## editing the file as suggested by the above aged exercise.  The `setAs`
## functions themselves are intentionally not documented so as to not clobber
## the base R `as` help.  This creates the following check warning:
##
## W  checking for missing documentation entries (...s)
##    Undocumented S4 methods:
##      generic 'coerce' and siglist 'TSSHMM,character'
##    All user level objects in a package (including S4 classes and methods)
##    should have documentation entries.
##
## Bioconductor core packages document `as` alongside the class descriptions or
## constructors, so these coercion functions should be upstreamed.
#' @importFrom methods coerce
#' @exportMethod coerce
methods::coerce
## Allow TSSHMM to be converted to character vector for logging.
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
