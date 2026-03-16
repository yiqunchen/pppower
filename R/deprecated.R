#' Deprecated functions in pppower
#'
#' @description
#' These functions have been renamed or merged into unified interfaces following
#' the \pkg{pwr}-style convention of "leave one NULL" (compute whichever of
#' power or sample size is left `NULL`).
#'
#' \describe{
#'   \item{`power_ppi_pp_mean`}{Use [power_ppi_mean()] instead.}
#'   \item{`power_ppi_pp_paired`}{Use [power_ppi_paired()] instead.}
#'   \item{`power_ppi_pp_paired_binary`}{Use [power_ppi_paired_binary()] instead.}
#'   \item{`power_ppi_pp_ttest`}{Use [power_ppi_ttest()] instead.}
#'   \item{`power_ppi_pp_ttest_binary`}{Use [power_ppi_ttest_binary()] instead.}
#'   \item{`n_required_ppi_pp`}{Use [power_ppi_mean()] or [power_ppi_regression()] with `n = NULL`.}
#'   \item{`n_required_ppi_pp_paired`}{Use [power_ppi_paired()] with `n = NULL`.}
#'   \item{`simulate_power_ppiplus_mean`}{Use [simulate_ppi_mean()] instead.}
#'   \item{`simulate_power_ppi_pp_ttest_binary`}{Use [simulate_ppi_ttest_binary()] instead.}
#'   \item{`simulate_power_ppi_mean`}{Use `simulate_ppi_vanilla_mean()` or [simulate_ppi_mean()].}
#' }
#'
#' @name deprecated
#' @keywords internal
NULL
