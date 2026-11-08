# Tools for LLM agents

Returns the main pppower design calculations as tool definitions created
with
[`ellmer::tool()`](https://ellmer.tidyverse.org/reference/tool.html), so
that an LLM agent (for example through
[`ellmer::chat_anthropic()`](https://ellmer.tidyverse.org/reference/chat_anthropic.html)
or an MCP client) can plan prediction-powered studies. Each tool takes
inputs in planning terms (effect size, number of unlabeled units,
\\R^2\\, or sensitivity and specificity) and returns the required
labeled sample size or the power, together with the classical comparison
and the main assumptions.

## Usage

``` r
pppower_mcp_tools()
```

## Value

A list of `ellmer` tool objects.

## Details

The tools are:

- `ppi_mean_design`:

  One-sample mean or prevalence.

- `ppi_two_sample_design`:

  Difference of two means.

- `ppi_two_proportion_design`:

  Difference of two proportions.

- `ppi_paired_design`:

  Paired mean difference.

- `ppi_2x2_design`:

  Odds ratio or relative risk in a 2x2 table.

- `ppi_planning_inputs`:

  Converts reported metrics (MSE, sensitivity and specificity, a pilot
  correlation) into conservative planning inputs.

Regression contrasts need pilot covariance blocks and are not exposed as
a tool; use
[`power_ppi_regression()`](https://yiqunchen.github.io/pppower/reference/power_ppi_regression.md)
directly.

Each tool returns a JSON object (as a character string) with named
fields.

## See also

[`pppower_mcp_server()`](https://yiqunchen.github.io/pppower/reference/pppower_mcp_server.md)

## Examples

``` r
tools <- pppower_mcp_tools()
vapply(tools, function(t) t@name, character(1))
#> [1] "ppi_mean_design"           "ppi_two_sample_design"    
#> [3] "ppi_two_proportion_design" "ppi_paired_design"        
#> [5] "ppi_2x2_design"            "ppi_planning_inputs"      
```
