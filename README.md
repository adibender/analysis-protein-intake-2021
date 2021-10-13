# Association of protein intake and time-to-event outcomes in critically ill patients

This repository contains the code used for the analyses presented in the manuscript

"Protein intake and outcome of critically ill patients – analysis of a large international database using piece-wise exponential additive mixed models (PAMMs)"

To reproduce the results presented in the maunscript, run the file in `rerun-analysis.R` (it is not necessary to run the first script `preproc-data.R`, as the resulting data is already contained within this repository. The original, uprocessed SAS data bases could not be shared.)

Note that this could take some time (and computational ressources). Optionally,
you can run the individual lines of code in `rerun-analysis.R` sequentially and restart your R session in between.


The session info is given below:


```r
> devtools::session_info()
─ Session info ───────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.0.3 (2020-10-10)
 os       Ubuntu 20.04.3 LTS
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Berlin
 date     2021-10-13

─ Packages ───────────────────────────────────────────────────────────────────
 package     * version    date       lib source
 assertthat    0.2.1      2019-03-21 [1] CRAN (R 4.0.2)
 backports     1.2.1      2020-12-09 [1] CRAN (R 4.0.3)
 cachem        1.0.6      2021-08-19 [1] CRAN (R 4.0.3)
 callr         3.7.0      2021-04-20 [1] CRAN (R 4.0.3)
 checkmate     2.0.0      2020-02-06 [1] CRAN (R 4.0.2)
 cli           3.0.1      2021-07-17 [1] CRAN (R 4.0.3)
 codetools     0.2-16     2018-12-24 [4] CRAN (R 4.0.0)
 colorspace    2.0-2      2021-06-24 [1] CRAN (R 4.0.3)
 crayon        1.4.1      2021-02-08 [1] CRAN (R 4.0.3)
 DBI           1.1.1      2021-01-15 [1] CRAN (R 4.0.3)
 desc          1.3.0      2021-03-05 [1] CRAN (R 4.0.3)
 devtools      2.4.2      2021-06-07 [1] CRAN (R 4.0.3)
 dplyr       * 1.0.7      2021-06-18 [1] CRAN (R 4.0.3)
 ellipsis      0.3.2      2021-04-29 [1] CRAN (R 4.0.3)
 fansi         0.5.0      2021-05-25 [1] CRAN (R 4.0.3)
 fastmap       1.1.0      2021-01-25 [1] CRAN (R 4.0.3)
 foreach       1.5.1      2020-10-15 [1] CRAN (R 4.0.3)
 Formula       1.2-4      2020-10-16 [1] CRAN (R 4.0.3)
 fs            1.5.0      2020-07-31 [1] CRAN (R 4.0.3)
 generics      0.1.0      2020-10-31 [1] CRAN (R 4.0.3)
 ggplot2     * 3.3.5      2021-06-25 [1] CRAN (R 4.0.3)
 glue          1.4.2      2020-08-27 [1] CRAN (R 4.0.2)
 gtable        0.3.0      2019-03-25 [1] CRAN (R 4.0.2)
 iterators     1.0.13     2020-10-15 [1] CRAN (R 4.0.3)
 lattice       0.20-41    2020-04-02 [4] CRAN (R 4.0.0)
 lava          1.6.9      2021-03-11 [1] CRAN (R 4.0.3)
 lazyeval      0.2.2      2019-03-15 [1] CRAN (R 4.0.2)
 lifecycle     1.0.1      2021-09-24 [1] CRAN (R 4.0.3)
 magrittr      2.0.1      2020-11-17 [1] CRAN (R 4.0.3)
 Matrix        1.2-18     2019-11-27 [4] CRAN (R 4.0.0)
 memoise       2.0.0      2021-01-26 [1] CRAN (R 4.0.3)
 mgcv        * 1.8-33     2020-08-27 [4] CRAN (R 4.0.2)
 munsell       0.5.0      2018-06-12 [1] CRAN (R 4.0.2)
 mvtnorm       1.1-2      2021-06-07 [1] CRAN (R 4.0.3)
 nlme        * 3.1-149    2020-08-23 [4] CRAN (R 4.0.2)
 numDeriv      2016.8-1.1 2019-06-06 [1] CRAN (R 4.0.2)
 pammtools   * 0.5.7      2021-06-21 [1] CRAN (R 4.0.3)
 patchwork   * 1.1.1      2020-12-17 [1] CRAN (R 4.0.3)
 pec           2020.11.17 2020-11-16 [1] CRAN (R 4.0.3)
 pillar        1.6.3      2021-09-26 [1] CRAN (R 4.0.3)
 pkgbuild      1.2.0      2020-12-15 [1] CRAN (R 4.0.3)
 pkgconfig     2.0.3      2019-09-22 [1] CRAN (R 4.0.2)
 pkgload       1.2.1      2021-04-06 [1] CRAN (R 4.0.3)
 prettyunits   1.1.1      2020-01-24 [1] CRAN (R 4.0.2)
 processx      3.5.2      2021-04-30 [1] CRAN (R 4.0.3)
 prodlim       2019.11.13 2019-11-17 [1] CRAN (R 4.0.2)
 ps            1.6.0      2021-02-28 [1] CRAN (R 4.0.3)
 purrr       * 0.3.4      2020-04-17 [1] CRAN (R 4.0.2)
 R6            2.5.1      2021-08-19 [1] CRAN (R 4.0.3)
 Rcpp          1.0.7      2021-07-07 [1] CRAN (R 4.0.3)
 remotes       2.4.0      2021-06-02 [1] CRAN (R 4.0.3)
 rlang         0.4.11     2021-04-30 [1] CRAN (R 4.0.3)
 rprojroot     2.0.2      2020-11-15 [1] CRAN (R 4.0.3)
 scales        1.1.1      2020-05-11 [1] CRAN (R 4.0.2)
 sessioninfo   1.1.1      2018-11-05 [1] CRAN (R 4.0.2)
 survival      3.2-13     2021-08-24 [1] CRAN (R 4.0.3)
 testthat      3.0.4      2021-07-01 [1] CRAN (R 4.0.3)
 tibble        3.1.5      2021-09-30 [1] CRAN (R 4.0.3)
 tidyr       * 1.1.4      2021-09-27 [1] CRAN (R 4.0.3)
 tidyselect    1.1.1      2021-04-30 [1] CRAN (R 4.0.3)
 timereg       2.0.0      2021-05-20 [1] CRAN (R 4.0.3)
 usethis       2.0.1      2021-02-10 [1] CRAN (R 4.0.3)
 utf8          1.2.2      2021-07-24 [1] CRAN (R 4.0.3)
 vctrs         0.3.8      2021-04-29 [1] CRAN (R 4.0.3)
 withr         2.4.2      2021-04-18 [1] CRAN (R 4.0.3)
```
