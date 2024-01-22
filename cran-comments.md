## Submission GCalignR version 1.0.6

This is a minor release, removing an unused function argument (gc_peak_df) from the function align_peaks. 

## Test environments
* devtools::check_mac_release()
* devtools::check_win_devel()
* devtools::check_win_release()
* devtools::check_rhub()
* local Windows 11 Home

## R CMD check results

There were no ERRORs or WARNINGs

There are 2 NOTEs:

Note 1:

```
checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''
```
According to [R-hub issue #560](https://github.com/r-hub/rhub/issues/560), this seems to be an Rhub issue and so can likely be ignored. 

Note 2:

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
According to [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.
