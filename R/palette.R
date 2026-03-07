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
  "Classical" = pppower_colors$classical
)

pppower_method_colors_full <- c(
 "PPI++ (Theoretical)"     = pppower_colors$ppi_pp,
 "PPI++ (Empirical)"       = pppower_colors$ppi_pp,
 "Classical (Theoretical)" = pppower_colors$classical,
 "Classical (Empirical)"   = pppower_colors$classical
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
# Theoretical curves:  solid   (linetype = 1, linewidth = 1)
# Classical curves:    dashed  (linetype = 2, linewidth = 1)
# Empirical points:    PPI++ uses circle (shape = 16, size = 2.5)
#                      Classical uses triangle (shape = 17, size = 2.5)
# Reference lines:     dotted  (linetype = 3, linewidth = 0.5)
# Base R:              lwd = 2

pppower_method_linetypes <- c(
  "PPI++"     = "solid",
  "Classical" = "dashed"
)

# ---------------------------------------------------------------------------
# ggplot2 theme
# ---------------------------------------------------------------------------

#' Publication-ready ggplot2 theme for pppower
#'
#' A minimal theme with bold titles, bottom legend, and no minor gridlines.
#' Merges the conventions of `theme_pub()` and `theme_validation()` used
#' in earlier simulation scripts.
#'
#' @param base_size Base font size (default 11).
#' @return A ggplot2 theme object.
#' @export
theme_pppower <- function(base_size = 11) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for theme_pppower().", call. = FALSE)
  }
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "gray40"),
      axis.title    = ggplot2::element_text(face = "bold"),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size),
      plot.margin      = ggplot2::margin(10, 10, 10, 10)
    )
}
