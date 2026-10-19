# ============================================================================
# Agent layer: pppower as tools for LLM agents (MCP server)
# ============================================================================
# The tools below wrap the main design functions with inputs that an agent can
# collect from a user in plain terms (effect size, number of unlabeled units,
# R^2 or sensitivity/specificity), and they return named lists that serialize
# cleanly to JSON. pppower_mcp_server() serves them over the Model Context
# Protocol through the mcptools package.

.mcp_notes <- function(extra = NULL) {
  c(
    "Assumes the prediction model is trained before the study and that labeled units are chosen completely at random from the study units.",
    "Uses the asymptotic Wald test with the optimal (oracle) tuning parameter. With an estimated tuning parameter, add about one label for Gaussian data, and more for rare binary outcomes.",
    "Plan with a conservative (lower) R^2 if prediction quality is uncertain; overstating R^2 gives an underpowered study, while the test stays valid.",
    extra
  )
}

# smallest integer n in [n_lo, n_hi] with power_fun(n) >= target
.mcp_solve_n <- function(power_fun, target, n_lo = 2, n_hi = 1e6) {
  if (power_fun(n_hi) < target) return(NA_integer_)
  if (power_fun(n_lo) >= target) return(as.integer(n_lo))
  root <- stats::uniroot(function(n) power_fun(n) - target, c(n_lo, n_hi))$root
  n <- ceiling(root)
  while (n > n_lo && power_fun(n - 1) >= target) n <- n - 1
  while (power_fun(n) < target) n <- n + 1
  as.integer(n)
}

.mcp_reduction <- function(n_ppi, n_cl) {
  if (is.na(n_ppi) || is.na(n_cl)) return(NA_real_)
  round(100 * (1 - n_ppi / n_cl), 1)
}

.mcp_check_r2 <- function(r2) {
  if (!is.numeric(r2) || length(r2) != 1 || r2 < 0 || r2 >= 1) {
    stop("r2 must be a single number in [0, 1).", call. = FALSE)
  }
}

# R^2 to plan with: the pilot value, or its lower confidence bound if pilot_n is given
.mcp_plan_r2 <- function(r2, pilot_n = NULL, confidence = 0.9) {
  .mcp_check_r2(r2)
  if (is.null(pilot_n)) return(r2)
  rho_lower_bound(sqrt(r2), pilot_n, confidence)^2
}

# ---------------------------------------------------------------------------
# Tool functions
# ---------------------------------------------------------------------------

.mcp_mean <- function(delta, N, r2 = NULL, sigma_y2 = 1, prevalence = NULL,
                      sensitivity = NULL, specificity = NULL, n = NULL,
                      power = 0.8, alpha = 0.05, pilot_n = NULL, confidence = 0.9) {
  if (!is.null(prevalence)) {
    if (is.null(sensitivity) || is.null(specificity)) {
      stop("For a binary outcome, give prevalence, sensitivity, and specificity.", call. = FALSE)
    }
    m <- binary_moments_from_sens_spec(p = prevalence, sens = sensitivity, spec = specificity)
    s_y2 <- m$sigma_y2; s_f2 <- m$sigma_f2; cv <- m$cov_y_f
    outcome <- "binary"
  } else {
    if (is.null(r2)) stop("Give r2 (continuous outcome) or prevalence with sensitivity and specificity (binary outcome).", call. = FALSE)
    r2 <- .mcp_plan_r2(r2, pilot_n, confidence)
    s_y2 <- sigma_y2; s_f2 <- sigma_y2; cv <- sqrt(r2) * sigma_y2
    outcome <- "continuous"
  }
  r2_used <- cv^2 / (s_y2 * s_f2)
  ppi_pow <- function(nn) suppressWarnings(power_ppi_mean(delta = delta, N = N, n = nn, alpha = alpha,
                                                         sigma_y2 = s_y2, sigma_f2 = s_f2, cov_y_f = cv))
  cl_pow <- function(nn) {
    se <- sqrt(s_y2 / nn); z <- stats::qnorm(1 - alpha / 2)
    stats::pnorm(-z + abs(delta) / se) + stats::pnorm(-z - abs(delta) / se)
  }
  if (is.null(n)) {
    n_ppi <- .mcp_solve_n(ppi_pow, power)
    n_cl <- .mcp_solve_n(cl_pow, power)
    list(design = "one-sample mean", outcome = outcome, target_power = power, alpha = alpha,
         labeled_n_required = n_ppi, classical_n_required = n_cl,
         label_reduction_percent = .mcp_reduction(n_ppi, n_cl),
         r2 = round(r2_used, 4), unlabeled_N = N, notes = .mcp_notes())
  } else {
    list(design = "one-sample mean", outcome = outcome, labeled_n = n, unlabeled_N = N,
         power_ppi = round(ppi_pow(n), 4), power_classical = round(cl_pow(n), 4),
         r2 = round(r2_used, 4), alpha = alpha, notes = .mcp_notes())
  }
}

.mcp_two_sample <- function(delta, N_per_group, r2_A, r2_B = NULL, sigma_y2_A = 1,
                            sigma_y2_B = NULL, n_per_group = NULL, power = 0.8, alpha = 0.05,
                            pilot_n = NULL, confidence = 0.9) {
  if (is.null(r2_B)) r2_B <- r2_A
  if (is.null(sigma_y2_B)) sigma_y2_B <- sigma_y2_A
  r2_A <- .mcp_plan_r2(r2_A, pilot_n, confidence)
  r2_B <- .mcp_plan_r2(r2_B, pilot_n, confidence)
  pw <- function(nn, rA, rB) power_ppi_ttest(
    delta = delta, n_A = nn, n_B = nn, N_A = N_per_group, N_B = N_per_group, alpha = alpha,
    sigma_y2_A = sigma_y2_A, sigma_f2_A = sigma_y2_A, cov_yf_A = sqrt(rA) * sigma_y2_A,
    sigma_y2_B = sigma_y2_B, sigma_f2_B = sigma_y2_B, cov_yf_B = sqrt(rB) * sigma_y2_B)
  ppi_pow <- function(nn) pw(nn, r2_A, r2_B)
  cl_pow <- function(nn) pw(nn, 0, 0)
  if (is.null(n_per_group)) {
    n_ppi <- .mcp_solve_n(ppi_pow, power); n_cl <- .mcp_solve_n(cl_pow, power)
    list(design = "two-sample difference in means", target_power = power, alpha = alpha,
         labeled_n_per_group = n_ppi, classical_n_per_group = n_cl,
         label_reduction_percent = .mcp_reduction(n_ppi, n_cl),
         r2_A = round(r2_A, 4), r2_B = round(r2_B, 4),
         unlabeled_N_per_group = N_per_group, notes = .mcp_notes())
  } else {
    list(design = "two-sample difference in means", labeled_n_per_group = n_per_group,
         unlabeled_N_per_group = N_per_group, power_ppi = round(ppi_pow(n_per_group), 4),
         power_classical = round(cl_pow(n_per_group), 4), alpha = alpha, notes = .mcp_notes())
  }
}

.mcp_two_proportion <- function(p_A, p_B, N_per_group, sensitivity, specificity,
                                n_per_group = NULL, power = 0.8, alpha = 0.05) {
  ppi_pow <- function(nn) power_ppi_ttest_binary(
    p_A = p_A, p_B = p_B, n_A = nn, n_B = nn, N_A = N_per_group, N_B = N_per_group,
    sens_A = sensitivity, spec_A = specificity, delta = p_A - p_B, alpha = alpha)
  cl_pow <- function(nn) {
    se <- sqrt(p_A * (1 - p_A) / nn + p_B * (1 - p_B) / nn); z <- stats::qnorm(1 - alpha / 2)
    d <- abs(p_A - p_B)
    stats::pnorm(-z + d / se) + stats::pnorm(-z - d / se)
  }
  if (is.null(n_per_group)) {
    n_ppi <- .mcp_solve_n(ppi_pow, power); n_cl <- .mcp_solve_n(cl_pow, power)
    list(design = "two-sample difference in proportions", target_power = power, alpha = alpha,
         labeled_n_per_group = n_ppi, classical_n_per_group = n_cl,
         label_reduction_percent = .mcp_reduction(n_ppi, n_cl),
         unlabeled_N_per_group = N_per_group, notes = .mcp_notes())
  } else {
    list(design = "two-sample difference in proportions", labeled_n_per_group = n_per_group,
         unlabeled_N_per_group = N_per_group, power_ppi = round(ppi_pow(n_per_group), 4),
         power_classical = round(cl_pow(n_per_group), 4), alpha = alpha, notes = .mcp_notes())
  }
}

.mcp_paired <- function(delta, N, sigma_d2, r2_d, n = NULL, power = 0.8, alpha = 0.05,
                        pilot_n = NULL, confidence = 0.9) {
  r2_d <- .mcp_plan_r2(r2_d, pilot_n, confidence)
  ppi_pow <- function(nn) suppressWarnings(power_ppi_paired(delta = delta, N = N, n = nn, alpha = alpha,
                                                           sigma_D2 = sigma_d2, rho_D = sqrt(r2_d)))
  cl_pow <- function(nn) {
    se <- sqrt(sigma_d2 / nn); z <- stats::qnorm(1 - alpha / 2)
    stats::pnorm(-z + abs(delta) / se) + stats::pnorm(-z - abs(delta) / se)
  }
  if (is.null(n)) {
    n_ppi <- .mcp_solve_n(ppi_pow, power); n_cl <- .mcp_solve_n(cl_pow, power)
    list(design = "paired mean difference", target_power = power, alpha = alpha,
         labeled_pairs_required = n_ppi, classical_pairs_required = n_cl,
         label_reduction_percent = .mcp_reduction(n_ppi, n_cl), r2_d = round(r2_d, 4), unlabeled_N = N,
         notes = .mcp_notes())
  } else {
    list(design = "paired mean difference", labeled_pairs = n, unlabeled_N = N,
         power_ppi = round(ppi_pow(n), 4), power_classical = round(cl_pow(n), 4),
         alpha = alpha, notes = .mcp_notes())
  }
}

.mcp_2x2 <- function(p_exposed, p_control, N, sensitivity, specificity,
                     effect_measure = "OR", prev_exposed = 0.5, n = NULL,
                     power = 0.8, alpha = 0.05) {
  em <- toupper(effect_measure)
  if (!em %in% c("OR", "RR")) stop("effect_measure must be 'OR' or 'RR'.", call. = FALSE)
  call_2x2 <- function(...) suppressWarnings(power_ppi_2x2(
    p_exp = p_exposed, p_ctrl = p_control, N = N, alpha = alpha, prev_exp = prev_exposed,
    sens = sensitivity, spec = specificity, effect_measure = em, ...))
  notes <- .mcp_notes(if (em == "OR")
    "The odds ratio uses one tuning parameter for the logistic model; group-specific tuning can save a few labels.")
  if (is.null(n)) {
    n_ppi <- call_2x2(n = NULL, power = power)
    list(design = paste("2x2 table,", em), target_power = power, alpha = alpha,
         labeled_n_required_total = n_ppi, unlabeled_N = N, notes = notes)
  } else {
    list(design = paste("2x2 table,", em), labeled_n_total = n, unlabeled_N = N,
         power_ppi = round(call_2x2(n = n, power = NULL), 4), alpha = alpha, notes = notes)
  }
}

.mcp_inputs <- function(sigma_y2 = NULL, mse = NULL, prevalence = NULL, sensitivity = NULL,
                        specificity = NULL, pilot_r = NULL, pilot_n = NULL,
                        confidence = 0.9) {
  out <- list()
  if (!is.null(mse)) {
    if (is.null(sigma_y2)) stop("Give sigma_y2 together with mse.", call. = FALSE)
    out$r2_from_mse <- max(0, 1 - mse / sigma_y2)
    out$mse_note <- "R^2 = 1 - MSE / Var(Y) is a lower bound on the squared correlation, so it is a conservative planning input."
  }
  if (!is.null(prevalence)) {
    if (is.null(sensitivity) || is.null(specificity)) {
      stop("Give prevalence, sensitivity, and specificity together.", call. = FALSE)
    }
    m <- binary_moments_from_sens_spec(p = prevalence, sens = sensitivity, spec = specificity)
    out$r2_binary <- m$cov_y_f^2 / (m$sigma_y2 * m$sigma_f2)
    out$prediction_positive_rate <- sensitivity * prevalence + (1 - specificity) * (1 - prevalence)
  }
  if (!is.null(pilot_r)) {
    if (is.null(pilot_n) || pilot_n <= 3) stop("Give pilot_n > 3 together with pilot_r.", call. = FALSE)
    out$r2_point <- pilot_r^2
    out$r2_lower_bound <- rho_lower_bound(pilot_r, pilot_n, confidence)^2
    out$confidence <- confidence
    out$pilot_note <- "Plan with r2_lower_bound (one-sided Fisher z bound) to protect against an optimistic pilot."
  }
  if (length(out) == 0) {
    stop("Give mse with sigma_y2, or prevalence with sensitivity and specificity, or pilot_r with pilot_n.", call. = FALSE)
  }
  out
}

#' Plan a prediction-powered study in one call (agent-friendly interface)
#'
#' A single entry point to the main design calculations, with inputs in
#' planning terms. It is used by the MCP tools ([pppower_mcp_tools()]), by the
#' agent skill in the `agents/` folder of the source repository, and by the
#' command-line script `agents/skills/ppi-power-analysis/scripts/ppi_design.R`,
#' so all three give the same answers. Leave the sample size argument empty to
#' solve for it, or give it to compute power.
#'
#' @param design One of `"mean"` (one-sample mean or prevalence),
#'   `"two_sample"`, `"two_proportion"`, `"paired"`, `"2x2"`, or `"inputs"`
#'   (convert reported metrics into planning inputs).
#' @param ... Arguments for the chosen design:
#'   \describe{
#'     \item{`"mean"`}{`delta`, `N`, and either `r2` (with optional `sigma_y2`)
#'       or `prevalence`, `sensitivity`, `specificity`; optional `n`, `power`,
#'       `alpha`, and `pilot_n` with `confidence` to plan at the lower
#'       confidence bound for a pilot `r2` ([rho_lower_bound()]).}
#'     \item{`"two_sample"`}{`delta`, `N_per_group`, `r2_A`, optional `r2_B`,
#'       `sigma_y2_A`, `sigma_y2_B`, `n_per_group`, `power`, `alpha`,
#'       `pilot_n`, `confidence`.}
#'     \item{`"two_proportion"`}{`p_A`, `p_B`, `N_per_group`, `sensitivity`,
#'       `specificity`, optional `n_per_group`, `power`, `alpha`.}
#'     \item{`"paired"`}{`delta`, `N`, `sigma_d2`, `r2_d`, optional `n`,
#'       `power`, `alpha`, `pilot_n`, `confidence`.}
#'     \item{`"2x2"`}{`p_exposed`, `p_control`, `N`, `sensitivity`,
#'       `specificity`, optional `effect_measure` (`"OR"` or `"RR"`),
#'       `prev_exposed`, `n`, `power`, `alpha`.}
#'     \item{`"inputs"`}{any of `mse` with `sigma_y2`; `prevalence` with
#'       `sensitivity` and `specificity`; `pilot_r` with `pilot_n` (and
#'       `confidence`).}
#'   }
#' @return A named list with the labeled sample size or power, the classical
#'   comparison, the percent reduction in labels, and a `notes` vector of
#'   assumptions.
#' @examples
#' ppi_design("mean", delta = 0.2, N = 1000, r2 = 0.49)
#' ppi_design("mean", delta = 0.2, N = 1000, r2 = 0.49, pilot_n = 50)
#' ppi_design("mean", delta = 0.05, N = 1000, prevalence = 0.3,
#'            sensitivity = 0.85, specificity = 0.85)
#' ppi_design("2x2", p_exposed = 0.4, p_control = 0.2, N = 1000,
#'            sensitivity = 0.8, specificity = 0.8, effect_measure = "RR")
#' ppi_design("inputs", pilot_r = 0.7, pilot_n = 60)
#' @export
ppi_design <- function(design = c("mean", "two_sample", "two_proportion",
                                  "paired", "2x2", "inputs"), ...) {
  design <- match.arg(design)
  f <- switch(design,
    mean = .mcp_mean, two_sample = .mcp_two_sample,
    two_proportion = .mcp_two_proportion, paired = .mcp_paired,
    "2x2" = .mcp_2x2, inputs = .mcp_inputs)
  f(...)
}

# Wrap a tool function so that it returns a JSON object with named fields,
# which MCP clients pass to the model as text. The wrapper keeps the formals of
# the original function, which ellmer uses to match the tool arguments.
.mcp_json_wrap <- function(f) {
  g <- function() {
    args <- as.list(environment())
    supplied <- !vapply(args, function(a) identical(a, quote(expr = )), logical(1))
    res <- do.call(f, args[supplied])
    as.character(jsonlite::toJSON(res, auto_unbox = TRUE, pretty = TRUE,
                                  digits = NA, null = "null", na = "null"))
  }
  formals(g) <- formals(f)
  g
}

# ---------------------------------------------------------------------------
# Exported: tool definitions and server
# ---------------------------------------------------------------------------

#' Tools for LLM agents
#'
#' Returns the main \pkg{pppower} design calculations as tool definitions
#' created with [ellmer::tool()], so that an LLM agent (for example through
#' [ellmer::chat_anthropic()] or an MCP client) can plan prediction-powered
#' studies. Each tool takes inputs in planning terms (effect size, number of
#' unlabeled units, \eqn{R^2}, or sensitivity and specificity) and returns the
#' required labeled sample size or the power, together with the classical
#' comparison and the main assumptions.
#'
#' The tools are:
#' \describe{
#'   \item{`ppi_mean_design`}{One-sample mean or prevalence.}
#'   \item{`ppi_two_sample_design`}{Difference of two means.}
#'   \item{`ppi_two_proportion_design`}{Difference of two proportions.}
#'   \item{`ppi_paired_design`}{Paired mean difference.}
#'   \item{`ppi_2x2_design`}{Odds ratio or relative risk in a 2x2 table.}
#'   \item{`ppi_planning_inputs`}{Converts reported metrics (MSE, sensitivity
#'     and specificity, a pilot correlation) into conservative planning inputs.}
#' }
#' Regression contrasts need pilot covariance blocks and are not exposed as a
#' tool; use [power_ppi_regression()] directly.
#'
#' Each tool returns a JSON object (as a character string) with named fields.
#'
#' @return A list of `ellmer` tool objects.
#' @seealso [pppower_mcp_server()]
#' @examplesIf requireNamespace("ellmer", quietly = TRUE)
#' tools <- pppower_mcp_tools()
#' vapply(tools, function(t) t@name, character(1))
#' @export
pppower_mcp_tools <- function() {
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    stop("The ellmer package is needed: install.packages(\"ellmer\").", call. = FALSE)
  }
  num <- function(desc, required = TRUE) ellmer::type_number(desc, required = required)
  alpha_arg <- num("Two-sided significance level. Default 0.05.", required = FALSE)
  power_arg <- num("Target power, used when solving for the sample size. Default 0.8.", required = FALSE)

  list(
    ellmer::tool(
      .mcp_json_wrap(.mcp_mean),
      name = "ppi_mean_design",
      description = paste(
        "Labeled sample size (or power) for testing a one-sample mean or prevalence with",
        "prediction-powered inference (PPI++), where a pretrained model predicts the outcome for",
        "all N unlabeled units and n units get gold-standard labels. For a continuous outcome give r2",
        "(the R^2 between predictions and outcome) and sigma_y2; for a binary outcome give prevalence,",
        "sensitivity, and specificity. Leave n empty to solve for the required n; give n to get power.",
        "Returns the classical sample size for comparison."),
      arguments = list(
        delta = num("Effect size: difference between the true mean and the null value, on the outcome scale."),
        N = num("Number of unlabeled units that will have predictions."),
        r2 = num("Continuous outcome: R^2 (squared correlation) between predictions and outcome.", required = FALSE),
        sigma_y2 = num("Continuous outcome: variance of the outcome. Default 1.", required = FALSE),
        prevalence = num("Binary outcome: expected prevalence in the study population.", required = FALSE),
        sensitivity = num("Binary outcome: sensitivity of the classifier.", required = FALSE),
        specificity = num("Binary outcome: specificity of the classifier.", required = FALSE),
        n = num("Labeled sample size. Give it to compute power; leave empty to solve for n.", required = FALSE),
        pilot_n = num(paste("Number of labeled pairs in the pilot that estimated R^2. If given, the plan uses",
                            "the one-sided lower confidence bound for R^2 instead of the pilot value."), required = FALSE),
        confidence = num("One-sided confidence level for the bound. Default 0.9.", required = FALSE),
        power = power_arg, alpha = alpha_arg
      ),
      annotations = ellmer::tool_annotations(read_only_hint = TRUE, idempotent_hint = TRUE)
    ),
    ellmer::tool(
      .mcp_json_wrap(.mcp_two_sample),
      name = "ppi_two_sample_design",
      description = paste(
        "Labeled sample size per group (or power) for comparing two group means with PPI++, with",
        "equal labeled and unlabeled sizes in the two groups. Give the R^2 of the predictions in each",
        "group. Leave n_per_group empty to solve for it. Returns the classical sample size for comparison."),
      arguments = list(
        delta = num("Difference in group means to detect."),
        N_per_group = num("Number of unlabeled units per group."),
        r2_A = num("R^2 between predictions and outcome in group A."),
        r2_B = num("R^2 in group B. Defaults to r2_A.", required = FALSE),
        sigma_y2_A = num("Outcome variance in group A. Default 1.", required = FALSE),
        sigma_y2_B = num("Outcome variance in group B. Defaults to sigma_y2_A.", required = FALSE),
        n_per_group = num("Labeled units per group. Give it to compute power.", required = FALSE),
        pilot_n = num(paste("Number of labeled pairs in the pilot that estimated R^2. If given, the plan uses",
                            "the one-sided lower confidence bound for R^2 instead of the pilot value."), required = FALSE),
        confidence = num("One-sided confidence level for the bound. Default 0.9.", required = FALSE),
        power = power_arg, alpha = alpha_arg
      ),
      annotations = ellmer::tool_annotations(read_only_hint = TRUE, idempotent_hint = TRUE)
    ),
    ellmer::tool(
      .mcp_json_wrap(.mcp_two_proportion),
      name = "ppi_two_proportion_design",
      description = paste(
        "Labeled sample size per group (or power) for comparing two proportions with PPI++, when a",
        "binary classifier with known sensitivity and specificity labels every unit. Leave n_per_group",
        "empty to solve for it."),
      arguments = list(
        p_A = num("Proportion in group A."),
        p_B = num("Proportion in group B."),
        N_per_group = num("Number of unlabeled units per group."),
        sensitivity = num("Sensitivity of the classifier."),
        specificity = num("Specificity of the classifier."),
        n_per_group = num("Labeled units per group. Give it to compute power.", required = FALSE),
        power = power_arg, alpha = alpha_arg
      ),
      annotations = ellmer::tool_annotations(read_only_hint = TRUE, idempotent_hint = TRUE)
    ),
    ellmer::tool(
      .mcp_json_wrap(.mcp_paired),
      name = "ppi_paired_design",
      description = paste(
        "Number of labeled pairs (or power) for a paired mean difference with PPI++, where predictions",
        "of the within-pair difference are available for N unlabeled pairs. Leave n empty to solve for it."),
      arguments = list(
        delta = num("Mean within-pair difference to detect."),
        N = num("Number of unlabeled pairs."),
        sigma_d2 = num("Variance of the within-pair difference."),
        r2_d = num("R^2 between the predicted and true within-pair differences."),
        n = num("Labeled pairs. Give it to compute power.", required = FALSE),
        pilot_n = num(paste("Number of labeled pairs in the pilot that estimated R^2. If given, the plan uses",
                            "the one-sided lower confidence bound for R^2 instead of the pilot value."), required = FALSE),
        confidence = num("One-sided confidence level for the bound. Default 0.9.", required = FALSE),
        power = power_arg, alpha = alpha_arg
      ),
      annotations = ellmer::tool_annotations(read_only_hint = TRUE, idempotent_hint = TRUE)
    ),
    ellmer::tool(
      .mcp_json_wrap(.mcp_2x2),
      name = "ppi_2x2_design",
      description = paste(
        "Total labeled sample size (or power) for testing an odds ratio or relative risk in a 2x2 table",
        "with PPI++, when exposure is known for everyone and a classifier predicts the binary outcome.",
        "Leave n empty to solve for it."),
      arguments = list(
        p_exposed = num("Outcome probability in the exposed group."),
        p_control = num("Outcome probability in the control group."),
        N = num("Total number of unlabeled units."),
        sensitivity = num("Sensitivity of the outcome classifier."),
        specificity = num("Specificity of the outcome classifier."),
        effect_measure = ellmer::type_enum(c("OR", "RR"), "Odds ratio (OR) or relative risk (RR). Default OR.",
                                           required = FALSE),
        prev_exposed = num("Fraction of units that are exposed. Default 0.5.", required = FALSE),
        n = num("Total labeled sample size. Give it to compute power.", required = FALSE),
        power = power_arg, alpha = alpha_arg
      ),
      annotations = ellmer::tool_annotations(read_only_hint = TRUE, idempotent_hint = TRUE)
    ),
    ellmer::tool(
      .mcp_json_wrap(.mcp_inputs),
      name = "ppi_planning_inputs",
      description = paste(
        "Converts what is usually reported about a prediction model into planning inputs for the other",
        "tools: R^2 from a mean squared error, R^2 for a binary classifier from prevalence, sensitivity,",
        "and specificity, and a conservative R^2 (lower confidence bound) from a pilot correlation.",
        "Use it before the design tools when the user does not know R^2 directly."),
      arguments = list(
        sigma_y2 = num("Outcome variance, needed with mse.", required = FALSE),
        mse = num("Mean squared error of the predictions.", required = FALSE),
        prevalence = num("Binary outcome prevalence.", required = FALSE),
        sensitivity = num("Classifier sensitivity.", required = FALSE),
        specificity = num("Classifier specificity.", required = FALSE),
        pilot_r = num("Correlation between predictions and outcome in a pilot sample.", required = FALSE),
        pilot_n = num("Size of the pilot sample.", required = FALSE),
        confidence = num("One-sided confidence level for the lower bound. Default 0.9.", required = FALSE)
      ),
      annotations = ellmer::tool_annotations(read_only_hint = TRUE, idempotent_hint = TRUE)
    )
  )
}

#' Serve pppower as an MCP server for LLM agents
#'
#' Starts a Model Context Protocol (MCP) server, through the \pkg{mcptools}
#' package, that exposes the tools from [pppower_mcp_tools()]. MCP clients such
#' as Claude Code or Claude Desktop can then call the design calculations
#' directly. The function blocks while the server runs, so it is meant to be
#' launched by the client, e.g.
#' `claude mcp add pppower -- Rscript -e "pppower::pppower_mcp_server()"`.
#'
#' @param type Transport, `"stdio"` (default) or `"http"`.
#' @param ... Passed to [mcptools::mcp_server()], e.g. `host` and `port` for
#'   the HTTP transport.
#' @return Does not return while the server is running.
#' @seealso [pppower_mcp_tools()]
#' @export
pppower_mcp_server <- function(type = c("stdio", "http"), ...) {
  type <- match.arg(type)
  if (!requireNamespace("mcptools", quietly = TRUE)) {
    stop("The mcptools package is needed: install.packages(\"mcptools\").", call. = FALSE)
  }
  mcptools::mcp_server(tools = pppower_mcp_tools(), type = type, session_tools = FALSE, ...)
}
