# Serve pppower as an MCP server for LLM agents

Starts a Model Context Protocol (MCP) server, through the mcptools
package, that exposes the tools from
[`pppower_mcp_tools()`](https://yiqunchen.github.io/pppower/reference/pppower_mcp_tools.md).
MCP clients such as Claude Code or Claude Desktop can then call the
design calculations directly. The function blocks while the server runs,
so it is meant to be launched by the client, e.g.
`claude mcp add pppower -- Rscript -e "pppower::pppower_mcp_server()"`.

## Usage

``` r
pppower_mcp_server(type = c("stdio", "http"), ...)
```

## Arguments

- type:

  Transport, `"stdio"` (default) or `"http"`.

- ...:

  Passed to
  [`mcptools::mcp_server()`](https://posit-dev.github.io/mcptools/reference/server.html),
  e.g. `host` and `port` for the HTTP transport.

## Value

Does not return while the server is running.

## See also

[`pppower_mcp_tools()`](https://yiqunchen.github.io/pppower/reference/pppower_mcp_tools.md)
