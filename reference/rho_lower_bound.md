# Lower confidence bound for the prediction-outcome correlation

Planning with the point estimate of \\\rho\_{Yf}\\ from a small pilot
gives an underpowered study about half of the time, because the estimate
is as likely to be too high as too low. This function returns the
one-sided Fisher \\z\\ lower bound
\$\$\tanh\\\operatorname{atanh}(\|\hat\rho\|) - z\_{c}/\sqrt{m -
3}\\,\$\$ truncated at 0 and carrying the sign of \\\hat\rho\\, where
\\m\\ is the number of labeled pairs in the pilot and \\c\\ is the
confidence level. Square it to get a lower bound for \\R^2\\.

## Usage

``` r
rho_lower_bound(rho_hat, pilot_n, confidence = 0.9)
```

## Arguments

- rho_hat:

  Pilot estimate of the correlation between the outcome and the
  predictions.

- pilot_n:

  Number of labeled pairs in the pilot (more than 3).

- confidence:

  One-sided confidence level (default 0.9).

## Value

The lower bound, with the sign of `rho_hat`.

## Details

The power functions
[`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md),
[`power_ppi_paired()`](https://yiqunchen.github.io/pppower/reference/power_ppi_paired.md),
and
[`power_ppi_ttest()`](https://yiqunchen.github.io/pppower/reference/power_ppi_ttest.md)
apply this bound when a pilot size is given.

## Examples

``` r
rho_lower_bound(0.7, pilot_n = 50)
#> [1] 0.5917579
rho_lower_bound(0.7, pilot_n = 50)^2   # lower bound for R^2
#> [1] 0.3501775
```
