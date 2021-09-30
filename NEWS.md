# tsshmm 0.7.0 (2020-09-17)

## New features

- Support Baum-Welch training.  The trained model can be loaded and saved using
  the `params()` accessor and setter.
- Support genome-wide Viterbi inference and remove the limiting regions input.

## Significant user-visible changes

- New models are created by the S4 `TSSHMM-class` for training and inference.
- `hmm()` has been renamed to `viterbi()`
- Added `train()` for model training.
- Added `params()` accessor and setter to load and save the model transition
  and emission matrices.
- `tss()` ties are now broken by strand.  In other words, if identical peaks
  are found on the positive strand, the left-most peak is chosen, but on the
  negative strand, the right-most peak is chosen.  This was feature was added
  to alleviate the weaker signal seen when plotting the CA transcription
  initiation motif of the negative strand compared to stronger signal seen on
  the positive strand.
- Documentation in the vignette appendix now describes the derivation and
  calculation of the minimum distance between neighboring promoter regions.
  This distance is necessary to flank reads so that Viterbi can be run
  genome-wide without exhausting RAM with the dense genomic windows.
- Documentation in the vignette appendix now includes `sessionInfo()`

## Bug fixes and improvements

- Generating large GRanges of reversed windows has been significantly sped up
  from 52 minutes down to 2 seconds by hacking the GRanges and IRanges `tile()`
  core with a drop-in replacement called `tile_with_rev()`.  This is
  effectively a one-line code change to the core Bioconductor functions which
  will be upstreamed.
- All C functions are now documented using doxygen markup, and C documentation
  consistency is also checked by continuous integration.
- HMM computation is now handled using the published GHMM library.  This was
  necessary for complex requirements of Baum-Welch training of the model such
  as tying emission states.  A system installed GHMM library is automatically
  preferred; if no system GHMM is detected and no GHMM_ROOT to a prefix
  installation is supplied, the fallback bundled GHMM dependency is instead
  patched, compiled and installed.
- Build system has migrated from `Makevars` to autotools `configure.ac` and
  `src/Makefile.am` to manage the GHMM dependency, and run recursive `make`
  when using the bundled GHMM library.


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
