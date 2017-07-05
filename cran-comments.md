## Resubmission
This is a resubmission. In this version I have implemented the following changes:

## Changes of existing code
* Within `align_chromatograms`, retention times are not rounded to a precision of two decimals anymore prior to executing the algorithm. Instead, calculations are internally based on the the value rounded to two decimals. After alignment, retention times are returned as they appear in the input.
* Axis labels in `plot.GCaling` and `gc_heatmap` were changed to enhance clearity. 

## Extensions
* `draw_chromatogram` is added and allows to represent peak lists graphically as chromatograms
* A second vignette `GCalignR How does the Algorithm works?` gives a more detailed tutorial on the concepts.

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
