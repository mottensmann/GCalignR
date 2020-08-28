
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GCalignR [<img src="vignettes/GCalignRLogo.png" height="200" align="right"/>](https://github.com/mottensmann/GCalignR)

![Build
Status](https://travis-ci.org/mottensmann/GCalignR.svg?branch=master)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/GCalignR)](https://cran.r-project.org/package=GCalignR)
[![](http://cranlogs.r-pkg.org/badges/grand-total/GCalignR)](https://cran.r-project.org/package=GCalignR)
[![](https://img.shields.io/badge/doi-10.1371/journal.pone.0198311-Darkorange.svg)](https://doi.org/10.1371/journal.pone.0198311)
[![](https://img.shields.io/badge/Altmetric-13-Darkorange.svg)](https://www.altmetric.com/details/43624695)
[![](https://img.shields.io/badge/cited%20in%20Web%20of%20Science%20Core%20Collection--blue.svg)](http://apps.webofknowledge.com/InboundService.do?customersID=LinksAMR&mode=FullRecord&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=PARTNER_APP&SrcAuth=LinksAMR&SID=F37s1YuMNyRR8cqMmcR&UT=WOS%3A000434384900030)

`GCalignR` provides simple functions to align peak lists obtained from
Gas Chromatography Flame Ionization Detectors (GC-FID) based on
retention times and plots to evaluate the quality of the alignment. The
package supports any other one-dimensional chromatograpy technique that
enables the user to create a peak list with at least one column
specifying retention times as illustrated below.

<img src="vignettes/Two_Chromas_Peak_List.png" style="display: block; margin: auto;" />

As with other software you need to get used to the input format which is
shown in the illustration:

  - Row 1: Sample names
  - Row 2: Variable names
  - Row 3-N: GC data
      - Each block belongs to a sample as shown for sample A (green) and
        sample B (orange) above

### Installing GCalignR:

The latest release v1.0.3 is on `CRAN`. [Click
here](https://github.com/mottensmann/GCalignR/releases) for an overview
of past releases and a brief description of applied changes.

``` r
install.packages("GCalignR", dependencies = T)
```

*The developmental (currently identical to the CRAN release) is always
available on GitHub*

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

If you encounter bugs or if you have any suggestions for improvement
(for instance on how to speed up the algorithm\!), just contact
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
