## Resubmission
This is a resubmission. In this version I have implemented the following changes:

* Within align_chromatograms, retention times are not rounded to a precision of two decimals anymore prior to executing the algorithm. Instead, calculations are internally based on the the value rounded to two decimals.
* Axis labels in plot functions were changed minorly. 

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

## travis-ci results
There were no ERRORs or WARNINGs or Notes
