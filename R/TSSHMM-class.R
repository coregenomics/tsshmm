#' @importFrom iterators nextElem idiv
#' @importFrom methods callNextMethod
NULL

#' Transcription start site (TSS) hidden Markov model (HMM).
#'
#' The TSS HMM model is implemented using the General Hidden Markov Model
#' (GHMM) C library.  Initialize generates a pre-designed hidden Markov model;
#' due to the complexity of the C-interface of the model, the model is managed
#' at the C-layer and the R layer has no control over setting up the model: R
#' merely manages the external reference pointer and provides the train,
#' decoding and print functions; there are no setters to arbitrarily modify the
#' model besides running the train function.
#'
#' @usage model <- new("TSSHMM")
#'
#' @slot external_pointer Memory address of C GHMM model.
#'
#' @return new("TSSHMM") returns a default model object to be trained and then
#'     used for decoding.
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

#' Train the transcription start site (TSS) hidden Markov model (HMM).
#'
#' To train the model, the sparse reads of the signal and background need to be
#' exploded into dense encoded windows categorized as either enriched, depleted
#' or no-read observations, and processed by the Baum-Welch EM algorithm to
#' update the model transition and emission probabilities.  The larger memory
#' required for the 2 steps of generating the dense training data batch and for
#' the large matrix to run the Baum-Welch algorithm are empirically determined
#' and adjusted to use 80% of the free system memory by adjusting the batch
#' size; similarly the time estimate to complete training is also calculated
#' after processing each batch.
#'
#' @rdname TSSHMM-class
#'
#' @param signal Stranded, single base \code{GRanges} with integer score.
#' @param bg Stranded, single base \code{GRanges} with integer score.
#'
#' @return train() returns the model with trained transition and emission
#'     probabilities.
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

        while (TRUE) {
            i <- i + 1
            tryCatch(chunk <- nextElem(iterator),
                     error = function(e) {
                         if (e == "StopIteration") {
                             break
                         } else { ## Unhandled error.
                             stop(e)
                         }
                     })
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

## setMethod(
##     "show",
##     signature = "TSSHMM",
##     definition = function(object) {
##         invisible(NULL)
##     }
## )
