
GCalignR
========

![Build Status](https://travis-ci.org/mottensmann/GCalignR.svg?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/GCalignR)](https://cran.r-project.org/package=GCalignR) [![](http://cranlogs.r-pkg.org/badges/grand-total/GCalignR)](https://cran.r-project.org/package=GCalignR)

`GCalignR` provides simple functions to align peak lists obtained from Gas Chromatography Flame Ionization Detectors (GC-FID) based on retention times and plots to evaluate the quality of the alignment. The package supports any other one-dimensional chromatograpy technique that enables the user to create a peak list with at least one column specifying retention times as illustrated below.

<img src="vignettes/Two_Chromas_Peak_List.png" width="864" style="display: block; margin: auto;" />

### Installing GCalignR:

-   The current release 1.0.1 is on CRAN.

``` r
install.packages("GCalignR", dependencies = T)
```

-   Get the latest developmental version

### Get started with GCalignR

To get started read the vignettes:

``` r
browseVignettes("GCalignR")
```

If you encounter bugs or if you have any suggestions for improvement, just contact meinolf.ottensmann\[at\]web.de

### Reference

[Ottensmann M, Stoffel MA, Nichols HJ, Hoffman JI (2018) GCalignR: An R package for aligning gas-chromatography data for ecological and evolutionary studies. PLoS ONE 13(6): e0198311. https://doi.org/10.1371/journal.pone.0198311](https://doi.org/10.1371/journal.pone.0198311)
