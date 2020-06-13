# tsshmm 0.4.0 (2020-06-13)

## New features

- Support both strands for `hmm()` and `tss()`.

## Significant user-visible changes

- Negative signal and background scores are no longer supported.
- Unstranded signal and background ranges are no longer supported.

## Bug fixes and improvements

- Remove the fragile starting state assumption of the published HMM by instead
  considering all starting states as equally likely.


# tsshmm 0.3.0 (2020-06-02)

## New features

- Added `tss()` peak detection accelerated by C API.


# tsshmm 0.2.0 (2020-05-31)

## New features

- Added `hmm()` high level GRanges interface to `viterbi()`.

## Significant user-visible changes

- Added package vignette illustrating the model.

## Bug fixes and improvements

- Package passes BiocCheck without any errors.


# tsshmm 0.1.0 (2020-05-30)

## New features

- Added `viterbi()` decoding accelerated by C API.

## Bug fixes and improvements

- Added GitLab continuous integration to regularly install and run `R CMD
  check` and using scripts in `tools/` directory.

- C-level tests for `viterbi()` using multiple models are in a separate
  repository [gitlab.com/omsai/viterbi/](https://gitlab.com/omsai/viterbi/)
