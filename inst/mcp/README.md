# pppower for LLM agents (MCP)

`pppower` can serve its design calculations as tools for LLM agents, through the
Model Context Protocol (MCP). An agent such as Claude Code or Claude Desktop can
then answer questions like "how many images do I need to label to detect a 1%
change in melanoma prevalence, given a classifier with 85% sensitivity and
specificity and 10,000 unlabeled images?" by calling the package directly.

## Tools

| Tool | What it computes |
| --- | --- |
| `ppi_mean_design` | One-sample mean or prevalence. Continuous outcome via `r2`, binary outcome via prevalence, sensitivity, and specificity. |
| `ppi_two_sample_design` | Difference of two means (per-group labeled size or power). |
| `ppi_two_proportion_design` | Difference of two proportions with a binary classifier. |
| `ppi_paired_design` | Paired mean difference. |
| `ppi_2x2_design` | Odds ratio or relative risk in a 2x2 table. |
| `ppi_planning_inputs` | Converts reported metrics into planning inputs: R^2 from a mean squared error, R^2 from prevalence, sensitivity, and specificity, and a conservative R^2 (lower confidence bound) from a pilot correlation. |

Each design tool solves for the labeled sample size when `n` is left empty and
returns the power when `n` is given. Every result includes the classical
sample size or power for comparison, the percent reduction in labels, and the
main assumptions, as a JSON object. Regression contrasts need pilot covariance
blocks and are not exposed as a tool; use `power_ppi_regression()` in R.

## Setup

```r
install.packages(c("ellmer", "mcptools"))
remotes::install_github("yiqunchen/pppower")
```

**Claude Code**

```bash
claude mcp add pppower -- Rscript -e "pppower::pppower_mcp_server()"
```

**Claude Desktop** (add to `claude_desktop_config.json`)

```json
{
  "mcpServers": {
    "pppower": {
      "command": "Rscript",
      "args": ["-e", "pppower::pppower_mcp_server()"]
    }
  }
}
```

**Inside R with ellmer** (no MCP client needed)

```r
library(pppower)
chat <- ellmer::chat_anthropic()
chat$register_tools(pppower_mcp_tools())
chat$chat("I have 1000 unlabeled patients and a model with R^2 = 0.5 for systolic
           blood pressure. How many patients must I measure to detect a 3 mmHg
           change from 125 mmHg if the SD is 20 mmHg?")
```

## Suggested agent workflow

1. Identify the design (one-sample, two-sample, paired, 2x2) and the effect size
   on the outcome scale.
2. Ask for the number of unlabeled units that will have predictions.
3. Ask how good the predictions are. If the user knows R^2, use it. Otherwise
   call `ppi_planning_inputs` with what they know (mean squared error,
   sensitivity and specificity, or a pilot correlation), and prefer the
   conservative value.
4. Call the design tool. Report the labeled sample size, the classical sample
   size, and the assumptions in the `notes` field.
5. If prediction quality is uncertain, repeat the call over a range of R^2 and
   report the range of sample sizes.

## Test the server by hand

```bash
Rscript -e "pppower::pppower_mcp_server()"
```

The server speaks JSON-RPC over standard input and output. MCP clients start it
and stop it automatically.
