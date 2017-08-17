## Resubmission
This is a resubmission. In this version I have implemented the following changes:

## Changes of existing code
* Progress bars are implemented using `pbapply` to estimate the run time of algorithms.
* More efficient code allows to speed up computations by an order of magnitude for some functions. 
* Within `align_chromatograms`, retention times are not rounded anymore prior to executing the algorithm. Instead, calculations are internally based on rounded values. The output contains the input retention time to allow cross-reference between raw data and aligned data.
* Axis labels in `plot.GCaling` and `gc_heatmap` were changed to enhance clearity. 
* Helpfile were rewritten.

## Extensions
* `draw_chromatogram` is added and allows to represent peak lists graphically as chromatograms
* A second vignette `GCalignR How does the Algorithm works?` gives a more detailed tutorial on the concepts.
* A CITATION file was added to link the package to the preprint of our manuscript. 

## Release summary

## Test environments
* devtools::build_win()
* local Mac OS X El Capitan 10.11
* local Windows 10 Home
* Travis-CI

## R CMD check results
There were no ERRORs or WARNINGs or Notes

## win_build check results
There were no ERRORs or WARNINGs and 1 Note:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Meinolf Ottensmann <meinolf.ottensmann@web.de>'
* Possibly mis-spelled words in DESCRIPTION
The flaged words are spelled correctly

## travis-ci results
There were no ERRORs or WARNINGs or Notes
