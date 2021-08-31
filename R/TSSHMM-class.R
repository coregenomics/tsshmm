#' @importFrom methods callNextMethod
NULL

#' Transcription start site (TSS) hidden Markov model (HMM).
#'
#' Implementation of the TSS HMM model using the General Hidden Markov Model
#' (GHMM) library.  Initialize a pre-designed hidden Markov model.  Due to the
#' complexity of the C-interface of the model, the model is managed in the
#' C-layer and the R layer merely manages the external reference pointer and
#' provides the print function; there are no setters to arbitrarily modify the
#' model besides running the training and decoding functions.
#'
#' @return A default model object to be trained and then used for decoding.
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

## setMethod(
##     "show",
##     signature = "TSSHMM",
##     definition = function(object) {
##         invisible(NULL)
##     }
## )
