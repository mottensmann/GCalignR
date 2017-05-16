# GCalignR 0.1.0.9000

___

## Changes since release 0.1.0

### Algorithm
* Retention times are not rounded to two decimals anymore. Calculations still capture a precision of two decimals for purely computational reasons. 
* Within the aligned results, retention times correspond to the input values. Linear adjustments are only used internally and are documented within the Logfile found in the output.

### warning messages
* Warnings addressing formatting issues are now more explicit and partly rephrased to avoid ambiguity

### Plots 
* Added horitontal axis to barplots summarising peak numbers in `plot.GCalign`
* Changed to more prominent colours in binary heatmaps with `gc_heatmap`

### Vignette
* Extended workflow and more comprehensive explainations
___



