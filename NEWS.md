# tsshmm 0.6.0 (2020-07-21)

## New features

- `hmm()` now returns a metadata column containing all hidden states as an
  `IntegerList()` to inform whether the region is a peaked or non-peaked
  promoter, with the state of each window to aid visualization of the model
  behavior output track against the input sequencing data track.
- The timing and intermediate steps of `hmm()` can be inspected by importing
  the `futile.flogger` library and setting `flog.threshold(DEBUG)`.

## Significant user-visible changes

- `hmm()` on negative strand data no longer suffers from a large speed penalty.
  Previously, calculations performed on negative strand were extremely
  inefficient due to a single call `endoapply(..., rev)` to reverse windows.
  This call has been eliminated by instead patching `tile()` methods used by
  `GRanges` and `IRanges` to directly produce reversed windows.

## Bug fixes and improvements

- `viterbi()` has been vectorized at the R layer by allowing for a list of
  regions to be provided resulting in about a 4x speedup.
- `viterbi()` has further sped up using OpenMP multithreading at the C layer.
  The speed of `viterbi()` is now comparable to the speed of `tss()`.
- Cleaned up `devtools::check()` warnings and notes.


# tsshmm 0.5.0 (2020-06-15)

## Significant user-visible changes

- `hmm()` and `tss()` now have documentation and examples.
- `hmm()` now returns `GRanges` instead of strand split `GRangesList`.

## Bug fixes and improvements

- Require sorted signal input to `tss()` required by underlying C algorithm.


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
