# ============================================================================
# Centralized color palette and ggplot2 theme for pppower
# ============================================================================
# All plotting code in the package, vignettes, and simulation scripts should
# reference these definitions instead of hardcoding colors or themes.

# ---------------------------------------------------------------------------
# Color palette (colorblind-friendly, Okabe-Ito inspired)
# ---------------------------------------------------------------------------

#' @keywords internal
pppower_colors <- list(
  # --- Method colors ---
  ppi_pp      = "#D55E00",
  classical   = "#0072B2",


  # --- Lambda-mode colors ---
  vanilla     = "#1B9E77",
  oracle      = "#D95F02",
  user        = "#7570B3",

  # --- Prediction quality levels ---
  poor        = "#999999",
  fair        = "#E69F00",
  good        = "#56B4E9",
  excellent   = "#009E73",

  # --- Outcome type colors ---
  continuous  = "#56B4E9",
  binary      = "#E69F00",

  # --- Reference / annotation ---
  reference   = "#999999",
  power_target = "gray50",

  # --- Empirical vs. analytical (base R) ---
  empirical   = "#0072B2",
  analytical  = "#D55E00"
)

# Named vectors for convenient scale_*_manual() usage
pppower_method_colors <- c(
  "PPI++"     = pppower_colors$ppi_pp,
  "Classical" = pppower_colors$classical,
  "Vanilla PPI" = pppower_colors$vanilla
)

pppower_method_colors_full <- c(
 "PPI++ (Theoretical)"     = pppower_colors$ppi_pp,
 "PPI++ (Empirical)"       = pppower_colors$ppi_pp,
 "Classical (Theoretical)" = pppower_colors$classical,
 "Classical (Empirical)"   = pppower_colors$classical,
 "Vanilla PPI (Theoretical)" = pppower_colors$vanilla,
 "Vanilla PPI (Empirical)"   = pppower_colors$vanilla
)

pppower_quality_colors <- c(
  "Poor"      = pppower_colors$poor,
  "Fair"      = pppower_colors$fair,
  "Good"      = pppower_colors$good,
  "Excellent" = pppower_colors$excellent
)

pppower_lambda_colors <- c(
  "vanilla" = pppower_colors$vanilla,
  "oracle"  = pppower_colors$oracle,
  "user"    = pppower_colors$user
)

pppower_outcome_colors <- c(
  "Continuous" = pppower_colors$continuous,
  "Binary"     = pppower_colors$binary
)

# ---------------------------------------------------------------------------
# Line type and point shape conventions
# ---------------------------------------------------------------------------
# Theoretical curves:  solid   (linewidth = 1.8)
# Classical curves:    dashed  (linewidth = 1.8)
# Empirical points:    PPI++ uses filled circle  (shape = 16, size = 3.5)
#                      Classical uses triangle    (shape = 17, size = 3.5)
# Reference lines:     dotted  (linewidth = 0.8)
# Base R:              lwd = 2.5

pppower_method_linetypes <- c(
  "PPI++"       = "solid",
  "Classical"   = "dashed",
  "Vanilla PPI" = "dotdash"
)

# ---------------------------------------------------------------------------
# ggplot2 theme
# ---------------------------------------------------------------------------

#' Publication-ready ggplot2 theme for pppower
#'
#' Clean theme with no top/right panel borders, bold titles, bottom legend,
#' and no minor gridlines. Designed for figures4papers-style output.
#'
#' @param base_size Base font size (default 14).
#' @return A ggplot2 theme object.
#' @export
theme_pppower <- function(base_size = 14) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for theme_pppower().", call. = FALSE)
  }
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      # Titles
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle    = ggplot2::element_text(size = base_size - 1, color = "gray30"),
      # Axes
      axis.title       = ggplot2::element_text(face = "bold", size = base_size),
      axis.text        = ggplot2::element_text(size = base_size - 2),
      axis.ticks.length = ggplot2::unit(4, "pt"),
      axis.line        = ggplot2::element_line(linewidth = 0.8),
      # Remove top & right panel borders
      panel.border     = ggplot2::element_blank(),
      axis.line.x      = ggplot2::element_line(linewidth = 0.8),
      axis.line.y      = ggplot2::element_line(linewidth = 0.8),
      # Grid
      panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      # Legend
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      legend.background = ggplot2::element_blank(),
      legend.key       = ggplot2::element_blank(),
      # Facets
      strip.text       = ggplot2::element_text(face = "bold", size = base_size - 1),
      strip.background = ggplot2::element_rect(fill = "gray95", color = NA),
      # Margins
      plot.margin      = ggplot2::margin(8, 12, 8, 8)
    )
}

# ---------------------------------------------------------------------------
# Flexible figure saver
# ---------------------------------------------------------------------------

#' Save a ggplot figure in PDF and/or PNG formats
#'
#' Follows figures4papers conventions: vector PDF with TrueType font
#' embedding and optional raster PNG at 400 DPI. White background,
#' tight crop.
#'
#' @param plot A ggplot2 object.
#' @param path Base path without extension (e.g., "figures/power_curve").
#' @param width Figure width in inches (default 10).
#' @param height Figure height in inches (default 5).
#' @param dpi DPI for PNG output (default 400).
#' @param formats Character vector containing any of \code{"pdf"} and
#'   \code{"png"}.
#' @return Invisible NULL. Side effect: writes the requested files.
#' @export
save_pppower_figure <- function(plot, path, width = 10, height = 5, dpi = 400,
                                formats = c("pdf", "png")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for save_pppower_figure().", call. = FALSE)
  }
  formats <- unique(formats)
  if ("pdf" %in% formats) {
    # Try cairo_pdf for TrueType embedding; fall back to default pdf device
    pdf_ok <- tryCatch({ grDevices::cairo_pdf(tempfile()); grDevices::dev.off(); TRUE },
                       error = function(e) FALSE, warning = function(w) FALSE)
    if (pdf_ok) {
      ggplot2::ggsave(
        paste0(path, ".pdf"), plot, width = width, height = height,
        device = grDevices::cairo_pdf
      )
    } else {
      ggplot2::ggsave(
        paste0(path, ".pdf"), plot, width = width, height = height
      )
    }
  }
  if ("png" %in% formats) {
    ggplot2::ggsave(
      paste0(path, ".png"), plot, width = width, height = height, dpi = dpi,
      bg = "white"
    )
  }
  invisible(NULL)
}
