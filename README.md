
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GCalignR [<img src="vignettes/GCalignRLogo.png" height="200" align="right"/>](https://github.com/mottensmann/GCalignR)

![Build
Status](https://travis-ci.org/mottensmann/GCalignR.svg?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/GCalignR)](https://cran.r-project.org/package=GCalignR)
[![](http://cranlogs.r-pkg.org/badges/grand-total/GCalignR)](https://cran.r-project.org/package=GCalignR)
[![](https://img.shields.io/badge/doi-10.1371/journal.pone.0198311-Darkorange.svg)](https://doi.org/10.1371/journal.pone.0198311)
[![](https://img.shields.io/badge/Altmetric-14-Darkorange.svg)](https://www.altmetric.com/details/43624695)

`GCalignR` provides simple functions to align peak lists obtained from
Gas Chromatography Flame Ionization Detectors (GC-FID) based on
retention times and plots to evaluate the quality of the alignment. The
package supports any other one-dimensional chromatograpy technique that
enables the user to create a peak list with at least one column
specifying retention times as illustrated below.

<img src="vignettes/Two_Chromas_Peak_List.png" style="display: block; margin: auto;" />

As with other software you need to get used to the input format which is
shown in the illustration:

-   Row 1: Sample names
-   Row 2: Variable names
-   Row 3-N: GC data
    -   Each block belongs to a sample as shown for sample A (green) and
        sample B (orange) above

### Installing GCalignR:

The latest release v1.0.3 is on `CRAN`. [Click
here](https://github.com/mottensmann/GCalignR/releases) for an overview
of past releases and a brief description of applied changes.

``` r
install.packages("GCalignR", dependencies = T)
```

*The developmental is always available on GitHub*

``` r
    if (!("devtools" %in% rownames(installed.packages()))) { 
    install.packages("devtools")
    } else if (packageVersion("devtools") < 1.6) {
    install.packages("devtools")
    }
    devtools::install_github("mottensmann/GCalignR", build_vignettes = TRUE)
```

### Get started with GCalignR

To get started read the vignettes:

``` r
browseVignettes("GCalignR")
```

Basic usage of the main function to align peaks:

-   `peak_data`: List of data frames, each corresponding to a sample
-   `rt_col_name`: column name of retention time values
-   `max_linear_shift`: Here, no adjustment of systematic linear drift
-   `max_diff_peak2mean`: Here, sort all peaks strictly by retention
    time
-   `min_diff_peak2peak`: Here, try to merge peaks when rt differs by
    less than 0.1

``` r
library(GCalignR)
aligned <- align_chromatograms(data = peak_data[1:10], # list of data frame 
                               rt_col_name = "time", # retention time
                               max_linear_shift = 0, #
                               max_diff_peak2mean = 0, 
                               min_diff_peak2peak = 0.08) 
#> Run GCalignR
#> Start: 2022-02-09 15:45:43
#> 
#> Data for 10 samples loaded.
#> No reference was specified. Hence, a reference will be selected automatically ...
#>  
#> 'C3' was selected on the basis of highest average similarity to all samples (score = 0.12).
#> Start correcting linear shifts with "C3" as a reference ...
#> 
#> Start aligning peaks ...  this might take a while!
#> 
#> Merge redundant rows ...
#>  
#> Alignment completed!
#> Time: 2022-02-09 15:46:06
```

*The parameter values above differ from the defaults shown in the paper
and the package vignette. In a nutshell, we now suggest in most cases to
set `max_diff_peak2mean = 0`. This way peaks are first simply sorted
based on the given retention time value and then purely
`min_diff_peak2peak` specifies which peaks will be evaluated for a
merge. Additionally, this enables the possibility for a considerable
boost in computation speed!*

If you encounter bugs or if you have any suggestions for improvement
(for instance on how to speed up the algorithm!), just contact
meinolf.ottensmann\[at\]web.de

*Also I´m happy to provide help if you can´t get it to work. Usually it
is easy to solve small problems. However, in order to simplify this
process please send a short description of the problem along with the
code you have been using as a script file (.R) together with a minimal
example input file (.txt).*

### Reference

[Ottensmann M, Stoffel MA, Nichols HJ, Hoffman JI (2018) GCalignR: An R
package for aligning gas-chromatography data for ecological and
evolutionary studies. PLoS ONE 13(6): e0198311.
https://doi.org/10.1371/journal.pone.0198311](https://doi.org/10.1371/journal.pone.0198311)
