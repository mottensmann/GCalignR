
GCalignR
========

![Build Status](https://travis-ci.org/mastoffel/GCalignR.svg?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/GCalignR)](https://cran.r-project.org/package=GCalignR) [![](http://cranlogs.r-pkg.org/badges/grand-total/GCalignR)](https://cran.r-project.org/package=GCalignR)

`GCalignR` provides simple functions to align gas-chromatography data based on retention times and plots to evaluate the quality of the alignment.

Installing GCalignR:

-   Get the latest development version from github with

``` r
    if (!("devtools" %in% rownames(installed.packages()))) { 
    install.packages("devtools")
    } else if (packageVersion("devtools") < 1.6) {
    install.packages("devtools")
    }
    devtools::install_github("mastoffel/GCalignR", build_vignettes = TRUE)
```

-   If the installation fails try installing the packages devtools, ggplot2 and vegan with the commmands below (tested with a new installation of R 3.2.5)

``` r
    if (!("devtools" %in% rownames(installed.packages()))) {
    install.packages("devtools")
    } else if (packageVersion("devtools") < 1.6) {
    install.packages("devtools")
    }
    install.packages("ggplot2",dependencies = TRUE)
    install.packages("vegan",dependencies = TRUE)
    devtools::install_github("mastoffel/GCalignR", build_vignettes = TRUE)    
```

### Get started with GCalignR

To get started read the vignette:

``` r
    vignette("GCalignR_step_by_step", package = "GCalignR")
```

If you encounter bugs or if you have any suggestions for improvement, just contact meinolf.ottensmann\[at\]web.de

### Reference

[Ottensmann, M., Stoffel, M.A., Hoffman, J.I. GCalignR: An R package for aligning Gas-Chromatography data. *In review*.](https://doi.org/10.1101/110494)
