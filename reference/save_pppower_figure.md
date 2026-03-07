# Save a ggplot figure in both PDF and PNG formats

Follows figures4papers conventions: vector PDF with TrueType font
embedding, plus raster PNG at 400 DPI. White background, tight crop.

## Usage

``` r
save_pppower_figure(plot, path, width = 10, height = 5, dpi = 400)
```

## Arguments

- plot:

  A ggplot2 object.

- path:

  Base path without extension (e.g., "figures/power_curve").

- width:

  Figure width in inches (default 10).

- height:

  Figure height in inches (default 5).

- dpi:

  DPI for PNG output (default 400).

## Value

Invisible NULL. Side effect: writes `path.pdf` and `path.png`.
