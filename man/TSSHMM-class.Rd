% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TSSHMM-class.R
\docType{class}
\name{TSSHMM-class}
\alias{TSSHMM-class}
\alias{TSSHMM}
\alias{dim,TSSHMM-method}
\alias{parameters}
\alias{parameters,TSSHMM-method}
\alias{parameters<-}
\alias{parameters<-,TSSHMM-method}
\alias{show,TSSHMM-method}
\alias{encode_obs}
\alias{train}
\alias{train,TSSHMM,RleList-method}
\alias{viterbi}
\alias{viterbi,TSSHMM,GRanges,GRanges-method}
\title{Transcription start site (TSS) hidden Markov model (HMM).}
\usage{
\S4method{dim}{TSSHMM}(x)

parameters(x)

\S4method{parameters}{TSSHMM}(x)

parameters(x) <- value

\S4method{parameters}{TSSHMM}(x) <- value

\S4method{show}{TSSHMM}(object)

encode_obs(signal, bg, nrow = 1000)

train(model, obs)

\S4method{train}{TSSHMM,RleList}(model, obs)

viterbi(model, signal, bg, ...)

\S4method{viterbi}{TSSHMM,GRanges,GRanges}(model, signal, bg, tol = 0.001, n = 1000)
}
\arguments{
\item{x}{S4 object of a pre-designed hidden Markov model.}

\item{value}{named list containing 2 model matrices: the transition state
probability matrix named "trans", and the discrete observable emission
probability matrix named "emis".}

\item{object}{S4 object to display.}

\item{signal}{Stranded, single base \code{GRanges} with integer score.}

\item{bg}{Stranded, single base \code{GRanges} with integer score.}

\item{nrow}{Integer, rows to convert at a time.  Primarily used for unit
tests.  Lowering or raising would use would use less and more memory,
respectively, because the encoding process suboptimally requires 70x
more memory than the input data.}

\item{model}{S4 object of a pre-designed hidden Markov model.}

\item{obs}{RleList of encoded observation windows.}

\item{...}{Optional arguments passed on from generics to methods.}

\item{tol}{Relative change to HMM steady state used to calculate
intra-promoter distance.  Preference should be given to
adjusting the \code{n} windows parameter before adjusting this.}

\item{n}{Maximum number of background windows after a promoter,
beyond which one is certain that the promoter has ended.  See
the appendix in the vignette to understand the theory behind
using \code{n} and \code{tol} in this calculation.  \code{n} is the upper
limit; the actual number of background windows is calculated
from the model.}
}
\description{
Find promoter regions in nascent RNA data using Andre Martin's 2014 hidden
Markov model for PRO-cap reads, and Pariksheet Nanda's model extension to
replace the background reference with a full gene transcription assay.
}
\details{
The model takes inputs of raw nascent RNA counts, aggregates 10 basepair
tiles by maximum value, and categorizes them into 3 types of observations:
\enumerate{
\item no signal (TAP+ == 0)
\item enriched  (TAP+ >  TAP-)
\item depleted  (TAP- >= TAP+ > 0)
}
...and searches for 3 groups of hidden states:
\enumerate{
\item background,
\item peaked TSS regions, and
\item non-peaked TSS regions
}

Both the peaked and non-peaked TSS regions each require 3 states to describe
them:
\itemize{
\item Non-peaked regions flanked by single low intensity transitions and a
moderately long intensity center.
\item Peaked regions flanked by one-or-more low intensity transitions and a
short intense center.
}
Therefore, in total the HMM has 7 states for PRO-cap and an additional state
for PRO-seq:
\describe{
\item{B}{Background}
\item{N1}{Non-peaked TSS transition state}
\item{N2}{Non-peaked TSS repeating state}
\item{N3}{Non-peaked TSS transition state}
\item{P1}{Peaked TSS moderate signal}
\item{P2}{Peaked TSS high signal}
\item{P3}{Peaked TSS moderate signal}
\item{GB}{Gene-body signal if background reference data provides it}
}

Finally, after the hidden states are obtained from Viterbi decoding, the
hidden states are converted back into stranded GRanges.
}
\section{Constructor}{


model <- TSSHMM()
model <- TSSHMM(bg_genebody = TRUE)

\code{TSSHMM()} returns a default model object to be trained and then used
for decoding.  An optional argument \code{bg_genebody} uses PRO-seq for
background instead of PRO-cap; the default is \code{bg_genebody = FALSE}.
}

\section{Accessors}{


\code{dim(x)} returns the number of states and emissions.

\code{parameters(model)}, \code{parameters(model) <- list(trans = ..., emis = ..., emis_tied = ..., bg_genebody = ...)} gets or sets 2 matrices: the transition
state probability matrix, and the discrete observable emission probability
matrix.

After training a model, you may wish to save the parameters to later reload
to eliminiate retraining the model in the future.
}

\section{Displaying}{


\code{show(model)} summarizes \code{parameters(model)} and the parameter dimesions,
namely, the number of transition states and number of discrete observable
emissions.
}

\section{Model Training}{


obs <- encode_obs(signal, background)
model <- train(model, obs)

\code{encode_obs(signal, background)} returns an \code{RleList} of encoded
observation states for training.  The arguments, \code{signal} and \code{bg}, are
stranded, single base \code{GRanges} with integer scores.

The maximum length for each element in \code{obs} is not configurable.  The
maximum length has been set to the upper limit of the Baum-Welch matrix
maximum size.  Using a lower maximum length reduces training accuracy,
because hidden states such as the background (B) and gene body (G) with
their >= 0.99 probability require long contiguous sequences for more
accurate training.

\code{train(model, obs)} returns the \code{model} with updated parameters of
transition and emission probabilities.  The argument, \code{obs} produced by the
\code{encode_obs()} function is an \code{RleList} of encoded observation states.  One
may set a seed value the observations before running \code{train()} to control
the observation shuffling

After training, you may wish to save the model parameters so that they can
be reloaded in a later R session as explained in the examples.  Note that
the time used for encoding the observations typically takes a magnitude
longer than the actual training.

To train the model, the sparse reads of the signal and background need to be
exploded into dense encoded windows categorized as either enriched, depleted
or no-read observations, which are then processed by the Baum-Welch EM
algorithm to update the model transition and emission probabilities.

The training data are divided into obs large enough to provide sufficient
samples for the calculating the background state transitions, but small
enough to work within numerical precision limits and to make the training
process more observable.  On each batch, the input training data is
transformed into dense training observations and then the Baum-Welch
algorithm is run.  After each batch, the model state is shown alongside a
time estimate to complete training if the logging level has not been reduced
from the INFO level.
}

\section{Model Evaluation}{


\code{viterbi(model, signal, background, tol = 1e-3, n = 200)} returns
\code{GRanges} of active promoter or enhancer regions along with the
decoded hidden states for each window.  The \code{tol} and \code{n}
parameters are passed on to the intra-promoter distance estimation
function that add flanking windows to ensure that no cluster of proximal
promoters is broken up.

Running \code{flog.threshold(DEBUG)} before running \code{viterbi} logs additional
information of what the function is doing.
}

\section{Coercion}{


\code{as(model, "character")} or \code{as.character(model)} compactly pastes
transition and emission matrices in the same line for logging.
}

\examples{
# Model evaluation.
signal <- GRanges(paste0("chr1:", c(100, 110, 200, 300), ":+"))
score(signal) <- rep(5L, 4)
signal
(bg <- GRanges())
(model <- TSSHMM())
(promoters_peaked <- viterbi(model, signal, bg))
# There is only 1 promoter in thise region because the first two signal
# values filling 2x 20 bp windows are captured by the HMM as a promoter.
# The remaining 2 signal peaks with no surrounding signal are ignored.
stopifnot(length(promoters) == 1)

# Model training.
#train(model, signal, bg)
(params <- parameters(model))

# Save and reload parameters.
\dontrun{save(params, file = "model_trained_on_foo_dataset.RData")
# ... in a later R session.
load(file = "model_trained_on_foo_dataset.RData")
}
(model <- TSSHMM())
parameters(model) <- params
model
}
\references{
\insertRef{core_analysis_2014}{tsshmm}
}
\seealso{
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4254663/}
}
