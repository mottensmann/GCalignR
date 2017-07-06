
GCalignR
========

![Build Status](https://travis-ci.org/mastoffel/GCalignR.svg?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/GCalignR)](https://cran.r-project.org/package=GCalignR) [![](http://cranlogs.r-pkg.org/badges/grand-total/GCalignR)](https://cran.r-project.org/package=GCalignR)

`GCalignR` provides simple functions to align peak lists obtained from Gas Chromatography Flame Ionization Detectors (GC-FID) based on retention times and plots to evaluate the quality of the alignment. The package supports any other one-dimensional chromatograpy technique that enables the user to create a peak list with at least one column specifying retention times as illustrated below.

<img src="internal/Two_Chromas_Peak_List.png" width="576" style="display: block; margin: auto;" />

Installing GCalignR:

-   The current release 0.1.0 is on CRAN

``` r
install.packages("GCalignR", dependencies = T)
```

-   Get the latest development version from GitHub with

``` r
    if (!("devtools" %in% rownames(installed.packages()))) { 
    install.packages("devtools")
    } else if (packageVersion("devtools") < 1.6) {
    install.packages("devtools")
    }
    devtools::install_github("mastoffel/GCalignR", build_vignettes = TRUE)
```

-   If the installation fails try installing the packages devtools, ggplot2 and vegan with the commands below (tested with a new installation of R 3.2.5)

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

To get started read the vignettes:

``` r
browseVignettes("GCalignR")
```

If you encounter bugs or if you have any suggestions for improvement, just contact meinolf.ottensmann\[at\]web.de

### Reference

[Ottensmann, M., Stoffel, M.A., Hoffman, J.I. GCalignR: An R package for aligning Gas-Chromatography data. *In review*.](https://doi.org/10.1101/110494)
