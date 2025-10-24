
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pppower

<!-- badges: start -->

[![R-CMD-check](https://github.com/yiqunchen/pppower/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yiqunchen/pppower/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of `pppower` is to conduct power analysis based on the
Prediction-Powered Inference (PPI) framework. It is a convenient tools
for planning and evaluating studies with PPI estimators. The package
provides analytical formulas, Monte Carlo estimators, plotting helpers,
and sample-size solvers for both mean estimation and linear regression
(PPI-OLS) workflows.

## Installation

You can install the development version of `pppower` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yiqunchen/pppower")
```

## Vignettes

- **Prediction-Powered Power Planning**

Open in R:

``` r
vignette("intro-ppi", package = "pppower")
```

    #> ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    #>  setting  value
    #>  version  R version 4.3.1 Patched (2023-07-19 r84711)
    #>  os       Rocky Linux 9.4 (Blue Onyx)
    #>  system   x86_64, linux-gnu
    #>  ui       X11
    #>  language (EN)
    #>  collate  en_US.UTF-8
    #>  ctype    en_US.UTF-8
    #>  tz       US/Eastern
    #>  date     2025-10-21
    #>  pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/ (via rmarkdown)
    #>  quarto   NA
    #> 
    #> ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    #>  ! package     * version    date (UTC) lib source
    #>    brio          1.1.3      2021-11-30 [2] CRAN (R 4.3.1)
    #>    bslib         0.5.1      2023-08-11 [2] CRAN (R 4.3.1)
    #>    cachem        1.0.8      2023-05-01 [2] CRAN (R 4.3.1)
    #>    callr         3.7.3      2022-11-02 [2] CRAN (R 4.3.1)
    #>    cli           3.6.1      2023-03-23 [2] CRAN (R 4.3.1)
    #>    commonmark    1.9.0      2023-03-17 [2] CRAN (R 4.3.1)
    #>    crayon        1.5.2      2022-09-29 [2] CRAN (R 4.3.1)
    #>    curl          7.0.0      2025-08-19 [1] CRAN (R 4.3.1)
    #>    desc          1.4.3      2023-12-10 [1] CRAN (R 4.3.1)
    #>    devtools      2.4.5      2022-10-11 [2] CRAN (R 4.3.1)
    #>    diffobj       0.3.5      2021-10-05 [2] CRAN (R 4.3.1)
    #>    digest        0.6.33     2023-07-07 [2] CRAN (R 4.3.1)
    #>    ellipsis      0.3.2      2021-04-29 [2] CRAN (R 4.3.1)
    #>    evaluate      1.0.5      2025-08-27 [1] CRAN (R 4.3.1)
    #>    fansi         1.0.4      2023-01-22 [2] CRAN (R 4.3.1)
    #>    fastmap       1.1.1      2023-02-24 [2] CRAN (R 4.3.1)
    #>    fs            1.6.3      2023-07-20 [2] CRAN (R 4.3.1)
    #>    glue          1.6.2      2022-02-24 [2] CRAN (R 4.3.1)
    #>    htmltools     0.5.6      2023-08-10 [2] CRAN (R 4.3.1)
    #>    htmlwidgets   1.6.2      2023-03-17 [2] CRAN (R 4.3.1)
    #>    httpuv        1.6.11     2023-05-11 [2] CRAN (R 4.3.1)
    #>    jquerylib     0.1.4      2021-04-26 [2] CRAN (R 4.3.1)
    #>    jsonlite      1.8.7      2023-06-29 [2] CRAN (R 4.3.1)
    #>    knitr         1.44       2023-09-11 [2] CRAN (R 4.3.1)
    #>    later         1.3.1      2023-05-02 [2] CRAN (R 4.3.1)
    #>    lifecycle     1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
    #>    magrittr      2.0.3      2022-03-30 [2] CRAN (R 4.3.1)
    #>    memoise       2.0.1      2021-11-26 [2] CRAN (R 4.3.1)
    #>    mime          0.12       2021-09-28 [2] CRAN (R 4.3.1)
    #>    miniUI        0.1.2      2025-04-17 [1] CRAN (R 4.3.1)
    #>    pillar        1.9.0      2023-03-22 [2] CRAN (R 4.3.1)
    #>    pkgbuild      1.4.8      2025-05-26 [1] CRAN (R 4.3.1)
    #>    pkgconfig     2.0.3      2019-09-22 [2] CRAN (R 4.3.1)
    #>    pkgload       1.4.1      2025-09-23 [1] CRAN (R 4.3.1)
    #>  P pppower     * 0.0.0.9000 2025-10-08 [?] load_all()
    #>    prettyunits   1.1.1      2020-01-24 [2] CRAN (R 4.3.1)
    #>    processx      3.8.2      2023-06-30 [2] CRAN (R 4.3.1)
    #>    profvis       0.3.8      2023-05-02 [2] CRAN (R 4.3.1)
    #>    promises      1.2.1      2023-08-10 [2] CRAN (R 4.3.1)
    #>    ps            1.7.5      2023-04-18 [2] CRAN (R 4.3.1)
    #>    purrr         1.0.2      2023-08-10 [2] CRAN (R 4.3.1)
    #>    R6            2.5.1      2021-08-19 [2] CRAN (R 4.3.1)
    #>    rcmdcheck     1.4.0      2021-09-27 [2] CRAN (R 4.3.1)
    #>    Rcpp          1.0.11     2023-07-06 [2] CRAN (R 4.3.1)
    #>    remotes       2.5.0      2024-03-17 [1] CRAN (R 4.3.1)
    #>    rlang         1.1.1      2023-04-28 [2] CRAN (R 4.3.1)
    #>    rmarkdown     2.30       2025-09-28 [1] CRAN (R 4.3.1)
    #>    roxygen2      7.2.3      2022-12-08 [2] CRAN (R 4.3.1)
    #>    rprojroot     2.1.1      2025-08-26 [1] CRAN (R 4.3.1)
    #>    rstudioapi    0.15.0     2023-07-07 [2] CRAN (R 4.3.1)
    #>    sass          0.4.7      2023-07-15 [2] CRAN (R 4.3.1)
    #>    sessioninfo   1.2.3      2025-02-05 [1] CRAN (R 4.3.1)
    #>    shiny         1.7.5      2023-08-12 [2] CRAN (R 4.3.1)
    #>    stringi       1.7.12     2023-01-11 [2] CRAN (R 4.3.1)
    #>    stringr       1.5.0      2022-12-02 [2] CRAN (R 4.3.1)
    #>    testthat    * 3.1.10     2023-07-06 [2] CRAN (R 4.3.1)
    #>    tibble        3.2.1      2023-03-20 [2] CRAN (R 4.3.1)
    #>    urlchecker    1.0.1      2021-11-30 [2] CRAN (R 4.3.1)
    #>    usethis       3.2.1      2025-09-06 [1] CRAN (R 4.3.1)
    #>    utf8          1.2.3      2023-01-31 [2] CRAN (R 4.3.1)
    #>    vctrs         0.6.3      2023-06-14 [2] CRAN (R 4.3.1)
    #>    waldo         0.6.2      2025-07-11 [1] CRAN (R 4.3.1)
    #>    withr         3.0.2      2024-10-28 [1] CRAN (R 4.3.1)
    #>    xfun          0.40       2023-08-09 [2] CRAN (R 4.3.1)
    #>    xml2          1.3.5      2023-07-06 [2] CRAN (R 4.3.1)
    #>    xopen         1.0.0      2018-09-17 [2] CRAN (R 4.3.1)
    #>    xtable        1.8-4      2019-04-21 [2] CRAN (R 4.3.1)
    #>    yaml          2.3.7      2023-01-23 [2] CRAN (R 4.3.1)
    #> 
    #>  [1] /users/mguo/R/4.3
    #>  [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
    #>  [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
    #> 
    #>  * ── Packages attached to the search path.
    #>  P ── Loaded and on-disk path mismatch.
    #> 
    #> ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
