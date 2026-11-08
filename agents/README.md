# pppower for AI agents

This folder makes `pppower` easy to use from AI agents (Claude Code, Claude Desktop, and other clients that support Agent Skills or the Model Context Protocol). There are three ways in, and they all call the same function, `pppower::ppi_design()`, so they give the same answers.

| Entry point | Where | Best for |
| --- | --- | --- |
| Agent Skill | `skills/ppi-power-analysis/` | Any agent that reads `SKILL.md` skills. Gives the agent the workflow, the formulas, the pitfalls, and a reporting template. |
| MCP server | `pppower::pppower_mcp_server()` (see `inst/mcp/README.md`) | Agents that call tools directly and get JSON back. |
| Command line | `skills/ppi-power-analysis/scripts/ppi_design.R` | Agents with a shell but no MCP connection. |

## Install the skill

Claude Code looks for skills in `~/.claude/skills/` (all projects) or `.claude/skills/` (one project):

```bash
mkdir -p ~/.claude/skills
cp -r agents/skills/ppi-power-analysis ~/.claude/skills/
```

## Connect the MCP server

```r
install.packages(c("ellmer", "mcptools"))
remotes::install_github("yiqunchen/pppower")
```

```bash
claude mcp add pppower -- Rscript -e "pppower::pppower_mcp_server()"
```

## Try it from the shell

```bash
Rscript agents/skills/ppi-power-analysis/scripts/ppi_design.R mean --delta 0.2 --N 1000 --r2 0.49
```

```json
{
  "design": "one-sample mean",
  "labeled_n_required": 110,
  "classical_n_required": 197,
  "label_reduction_percent": 44.2,
  ...
}
```

## Skill layout

```
skills/ppi-power-analysis/
  SKILL.md                    workflow, decision table, pitfalls
  references/formulas.md      variances, sample sizes, equivalences with classical methods
  references/related-tools.md ppi_py, dsl, survey, PROCOVA, ipd: coverage and agreement
  references/reporting.md     sample size paragraph template and checklist
  scripts/ppi_design.R        command-line front end, prints JSON
```
