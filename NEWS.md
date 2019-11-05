# GCalignR 1.0.2.9
## Recent changes

* Added a new boolean parameter 'remove_empty' for the function 'align_chromatograms'. If samples are empty  (ie. no peak) this parameter allows to remove those cases from the dataset to avoid problems in post-hoc analyses. By default FALSE, ie. all but the blank samples are kept.
___
# GCalignR 0.1.0.9000
## Changes since release 0.1.0

### Algorithm
* Retention times are not rounded to two decimals anymore. Calculations still capture a precision of two decimals for purely computational reasons. 
* Within the aligned results, retention times correspond to the input values. Linear adjustments are only used internally and are documented within the Logfile found in the output.

### Dependencies
The bioconductor package (MassSpecWavelet)[http://bioconductor.org/packages/MassSpecWavelet/] is added as *in R* solution for picking peaks from GC data. 


### warning messages
* Warnings addressing formatting issues are now more explicit and partly rephrased to avoid ambiguity

### Plots 
* Added horitontal axis to barplots summarising peak numbers in `plot.GCalign`
* Changed to more prominent colours in binary heatmaps with `gc_heatmap`

### Vignette
* Extended workflow and more comprehensive explainations
___



