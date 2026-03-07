\name{power_ppi_pp_paired_binary}
\alias{power_ppi_pp_paired_binary}
\title{Power for PPI++ paired test with binary outcomes}
\description{
Computes power for a paired test when outcomes are binary and predictions are
derived from a classifier with sensitivity and specificity. The paired
differences are handled via a CLT approximation using population moments
from a joint Bernoulli model with within-pair correlation.
}
\usage{
power_ppi_pp_paired_binary(delta, N, n, alpha = 0.05, p_A, p_B, rho_within, sens, spec)
}
\arguments{
  \item{delta}{Effect size: mean difference \eqn{E[Y^A - Y^B]}. If NULL, uses \eqn{p_A - p_B}.}
  \item{N}{Number of unlabeled pairs with predictions.}
  \item{n}{Number of labeled pairs.}
  \item{alpha}{Two-sided significance level (default 0.05).}
  \item{p_A}{Marginal probability \eqn{P(Y^A=1)}.}
  \item{p_B}{Marginal probability \eqn{P(Y^B=1)}.}
  \item{rho_within}{Pearson correlation between \eqn{Y^A} and \eqn{Y^B}.}
  \item{sens}{Classifier sensitivity.}
  \item{spec}{Classifier specificity.}
}
\value{
Scalar power in \eqn{[0,1]}.
}
\details{
The joint probability \eqn{P(Y^A=1, Y^B=1)} is set to
\deqn{p_{11} = p_A p_B + \rho_{\text{within}} \sqrt{p_A(1-p_A)p_B(1-p_B)}}
and clipped to the feasible range.
}
\examples{
power_ppi_pp_paired_binary(
  delta = 0.08,
  N = 5000,
  n = 500,
  p_A = 0.34,
  p_B = 0.26,
  rho_within = 0.3,
  sens = 0.85,
  spec = 0.85
)
}
