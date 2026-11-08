# Plan a prediction-powered study in one call (agent-friendly interface)

A single entry point to the main design calculations, with inputs in
planning terms. It is used by the MCP tools
([`pppower_mcp_tools()`](https://yiqunchen.github.io/pppower/reference/pppower_mcp_tools.md)),
by the agent skill in the `agents/` folder of the source repository, and
by the command-line script
`agents/skills/ppi-power-analysis/scripts/ppi_design.R`, so all three
give the same answers. Leave the sample size argument empty to solve for
it, or give it to compute power.

## Usage

``` r
ppi_design(
  design = c("mean", "two_sample", "two_proportion", "paired", "2x2", "inputs"),
  ...
)
```

## Arguments

- design:

  One of `"mean"` (one-sample mean or prevalence), `"two_sample"`,
  `"two_proportion"`, `"paired"`, `"2x2"`, or `"inputs"` (convert
  reported metrics into planning inputs).

- ...:

  Arguments for the chosen design:

  `"mean"`

  :   `delta`, `N`, and either `r2` (with optional `sigma_y2`) or
      `prevalence`, `sensitivity`, `specificity`; optional `n`, `power`,
      `alpha`, and `pilot_n` with `confidence` to plan at the lower
      confidence bound for a pilot `r2`
      ([`rho_lower_bound()`](https://yiqunchen.github.io/pppower/reference/rho_lower_bound.md)).

  `"two_sample"`

  :   `delta`, `N_per_group`, `r2_A`, optional `r2_B`, `sigma_y2_A`,
      `sigma_y2_B`, `n_per_group`, `power`, `alpha`, `pilot_n`,
      `confidence`.

  `"two_proportion"`

  :   `p_A`, `p_B`, `N_per_group`, `sensitivity`, `specificity`,
      optional `n_per_group`, `power`, `alpha`.

  `"paired"`

  :   `delta`, `N`, `sigma_d2`, `r2_d`, optional `n`, `power`, `alpha`,
      `pilot_n`, `confidence`.

  `"2x2"`

  :   `p_exposed`, `p_control`, `N`, `sensitivity`, `specificity`,
      optional `effect_measure` (`"OR"` or `"RR"`), `prev_exposed`, `n`,
      `power`, `alpha`.

  `"inputs"`

  :   any of `mse` with `sigma_y2`; `prevalence` with `sensitivity` and
      `specificity`; `pilot_r` with `pilot_n` (and `confidence`).

## Value

A named list with the labeled sample size or power, the classical
comparison, the percent reduction in labels, and a `notes` vector of
assumptions.

## Examples

``` r
ppi_design("mean", delta = 0.2, N = 1000, r2 = 0.49)
#> $design
#> [1] "one-sample mean"
#> 
#> $outcome
#> [1] "continuous"
#> 
#> $target_power
#> [1] 0.8
#> 
#> $alpha
#> [1] 0.05
#> 
#> $labeled_n_required
#> [1] 110
#> 
#> $classical_n_required
#> [1] 197
#> 
#> $label_reduction_percent
#> [1] 44.2
#> 
#> $r2
#> [1] 0.49
#> 
#> $unlabeled_N
#> [1] 1000
#> 
#> $notes
#> [1] "Assumes the prediction model is trained before the study and that labeled units are chosen completely at random from the study units."                                                  
#> [2] "Uses the asymptotic Wald test with the optimal (oracle) tuning parameter. With an estimated tuning parameter, add about one label for Gaussian data, and more for rare binary outcomes."
#> [3] "Plan with a conservative (lower) R^2 if prediction quality is uncertain; overstating R^2 gives an underpowered study, while the test stays valid."                                      
#> 
ppi_design("mean", delta = 0.2, N = 1000, r2 = 0.49, pilot_n = 50)
#> $design
#> [1] "one-sample mean"
#> 
#> $outcome
#> [1] "continuous"
#> 
#> $target_power
#> [1] 0.8
#> 
#> $alpha
#> [1] 0.05
#> 
#> $labeled_n_required
#> [1] 136
#> 
#> $classical_n_required
#> [1] 197
#> 
#> $label_reduction_percent
#> [1] 31
#> 
#> $r2
#> [1] 0.3502
#> 
#> $unlabeled_N
#> [1] 1000
#> 
#> $notes
#> [1] "Assumes the prediction model is trained before the study and that labeled units are chosen completely at random from the study units."                                                  
#> [2] "Uses the asymptotic Wald test with the optimal (oracle) tuning parameter. With an estimated tuning parameter, add about one label for Gaussian data, and more for rare binary outcomes."
#> [3] "Plan with a conservative (lower) R^2 if prediction quality is uncertain; overstating R^2 gives an underpowered study, while the test stays valid."                                      
#> 
ppi_design("mean", delta = 0.05, N = 1000, prevalence = 0.3,
           sensitivity = 0.85, specificity = 0.85)
#> $design
#> [1] "one-sample mean"
#> 
#> $outcome
#> [1] "binary"
#> 
#> $target_power
#> [1] 0.8
#> 
#> $alpha
#> [1] 0.05
#> 
#> $labeled_n_required
#> [1] 458
#> 
#> $classical_n_required
#> [1] 660
#> 
#> $label_reduction_percent
#> [1] 30.6
#> 
#> $r2
#> [1] 0.4466
#> 
#> $unlabeled_N
#> [1] 1000
#> 
#> $notes
#> [1] "Assumes the prediction model is trained before the study and that labeled units are chosen completely at random from the study units."                                                  
#> [2] "Uses the asymptotic Wald test with the optimal (oracle) tuning parameter. With an estimated tuning parameter, add about one label for Gaussian data, and more for rare binary outcomes."
#> [3] "Plan with a conservative (lower) R^2 if prediction quality is uncertain; overstating R^2 gives an underpowered study, while the test stays valid."                                      
#> 
ppi_design("2x2", p_exposed = 0.4, p_control = 0.2, N = 1000,
           sensitivity = 0.8, specificity = 0.8, effect_measure = "RR")
#> $design
#> [1] "2x2 table, RR"
#> 
#> $target_power
#> [1] 0.8
#> 
#> $alpha
#> [1] 0.05
#> 
#> $labeled_n_required_total
#> [1] 135
#> 
#> $unlabeled_N
#> [1] 1000
#> 
#> $notes
#> [1] "Assumes the prediction model is trained before the study and that labeled units are chosen completely at random from the study units."                                                  
#> [2] "Uses the asymptotic Wald test with the optimal (oracle) tuning parameter. With an estimated tuning parameter, add about one label for Gaussian data, and more for rare binary outcomes."
#> [3] "Plan with a conservative (lower) R^2 if prediction quality is uncertain; overstating R^2 gives an underpowered study, while the test stays valid."                                      
#> 
ppi_design("inputs", pilot_r = 0.7, pilot_n = 60)
#> $r2_point
#> [1] 0.49
#> 
#> $r2_lower_bound
#> [1] 0.3633841
#> 
#> $confidence
#> [1] 0.9
#> 
#> $pilot_note
#> [1] "Plan with r2_lower_bound (one-sided Fisher z bound) to protect against an optimistic pilot."
#> 
```
