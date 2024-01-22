
# GcalignR 1.0.6

------------------------------------------------------------------------

- Removing unused argument `gc_peak_df` from `align_peaks`

# GCalignR 1.0.5

------------------------------------------------------------------------

- Bugfix in `choose_optimal_reference` that always selected the first
  sample as a reference. Thanks to Heberto del Rio who pointed this out
  on <https://github.com/mottensmann/GCalignR/issues/27>

# GCalignR 1.0.3.9

- **Speedboost** when setting `max_diff_peak2mean = 0`: In this special
  case there is no need to use a time-consuming iterative approach but
  peaks can be sorted simply based on absolute values. This is
  implemented in two steps. (1) Across all samples, unique retention
  times are extracted, sorted in increasing temporal order and written
  to a template data frame. (2) For each sample, peaks are matched to
  the corresponding row of the template data frame.
- Small bug fixed that caused problems when plotting x-axis labels in
  `gc_heatmap`.
- Added a test for detecting inconsistently ordered retention times
  within samples. Retention times are expected in increasing order,
  starting with the lowest number. If this assumption is violated,
  retention times are reordered and a warning is shown.

# GCalignR 1.0.3

------------------------------------------------------------------------

- Added `fill = TRUE` as a parameter in `utils::read.table` when reading
  data from text within internal functions. *Loading GC data with
  utils::read.table failed in cases of missing values in a column
  (i.e. empty). This is the correct behaviour as missing data should
  always be coded explicitly by ‘NA’*
- Tibbles are now coerced to data frames
- Added a new boolean parameter `remove_empty` for the main function
  `align_chromatograms`. If samples are empty (i.e.. no peak) this
  parameter allows to remove those cases from the dataset to avoid
  problems in post-hoc analyses. By default `FALSE`, i.e.. all but the
  blank samples are kept.
- Added a new boolean parameter `permute` for the functions
  `align_chromatograms` and `align_peaks`. This allows to change the
  default behaviour of random permutation of samples during the
  alignment and might be useful if exact replication is needed.

# GCalignR 1.0.2

------------------------------------------------------------------------

- The accompanying manuscript is published
  <https://doi.org/10.1371/journal.pone.0198311> and the citation has
  been added
- The function *beta* `read_empower2` allows to import HPLC data that
  has been generated using the EMPOWER 2 software

# GCalignR 1.0.1

------------------------------------------------------------------------

**Bugfixes**

- A bugfix was applied for handling multiple blanks correctly.
- Progressbars are removed in non-interactive R sessions

------------------------------------------------------------------------

# GCalignR 1.0.0

**New functions implemented**

- `choose_optimal_reference` offers an automatism to pick suitable
  references.
- `draw_chromatograms` allows to represent a peak list in form of
  chromatogram.
- `remove_blanks`allows to get rid of peaks that represent contamination
  after aligning a dataset
- `remove_singletons` allows to remove single peaks from the dataset
  after aligning
- `merge_redundant_rows` allows to merge rows that were not recognised
  as redundant during the alignment by increasing the threshold value
  for the evaluation of similarity

**Algorithm**

- Using `pbapply`, we implemented progress bars to inform the user about
  the progress and the estimated running time of intermediate steps in
  the alignment of peak lists.
- By implementing more efficient code, we were able to speed up the
  processing, especially picking references is faster by an order of
  magnitude.
- Retention times are not rounded to two decimals anymore. Calculations
  still capture a precision of two decimals for computational reasons.
- Within the aligned results, retention times correspond to the input
  values. Linear adjustments are only used internally and are documented
  within the Logfile accessible in the output.
- Reference samples that are used for the coarse alignment of retention
  times can be picked using a novel algorithm that determines the
  average similarity across the dataset.

**warning messages**

- Warnings addressing formatting issues are now more explicit and partly
  rephrased to avoid ambiguity.

**Plots**

- Added horizontal axis to barplots summarising peak numbers in
  `plot.GCalign`.
- Changed to more prominent colours in binary heatmaps with
  `gc_heatmap`.
- The function `draw_chromatograms` was added as another visualisation
  tool.

**Vignettes**

- We included a second vignette that explains the algorithm and the
  supported data in detail.

**Documentation**

- Helpfiles were rewritten to enhance clarity.

------------------------------------------------------------------------
