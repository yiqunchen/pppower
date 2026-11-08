# Interactive PPI Sample Size Calculator

This calculator uses **prediction-powered inference (PPI/`PPI++`)**
formulas from this package for quick study planning examples. It is
designed for planning and intuition; for final protocol numbers, confirm
assumptions with your statistician.

Interactive PPI Tool

## Prediction-powered sample size calculator

Browser-based planner aligned with pppower formulas.

**Assumptions:** labeled outcomes are sampled completely at random from
the target population, and labeled/unlabeled sets are from the same
population.

**Non-response:** expected share of sampled units without a usable
labeled outcome. Recruitment-adjusted n inflates by
`1 / (1 - non-response)`.

**Methods:** details and derivations are on [Variance &
Methods](https://yiqunchen.github.io/pppower/articles/deep-dive-math.md).

PPI one-sample mean / prevalence

PPI difference in prevalence

PPI difference in means

PPI paired mean difference

PPI paired proportions

### Inputs

Outcome type Prevalence (binary) Mean (continuous)

Lambda mode Oracle (PPI++) Vanilla (lambda = 1)

Effect size delta

Unlabeled N

Alpha (two-sided)

Target power (%)

Binary prevalence inputs

Prevalence p = P(Y=1) (%)

Expected non-response (%)

Prediction quality

Sensitivity

Specificity

Continuous mean moments

Var(Y) = sigma_y2

Var(f) = sigma_f2

Cov(Y,f) = cov_y_f

Expected non-response (%)

Pilot size m (optional) If the correlation implied by these moments came
from a pilot of m labeled units, plan at its one-sided 90% lower
confidence bound.

### Results

PPI analysis n (base) **-**

PPI recruitment n (adjusted) **-**

Classical n (no prediction) **-**

Labeled saving vs classical **-**

Classic

-

PPI

-

PPI adj.

-

Formula note will appear here.

### Inputs

Group A prevalence pA (%)

Group B prevalence pB (%)

Unlabeled N_A

Unlabeled N_B

Alpha (two-sided)

Target power (%)

Prediction quality by group

Sensitivity A

Specificity A

Sensitivity B

Specificity B

Expected non-response A (%)

Expected non-response B (%)

### Results

PPI n per group (base) **-**

PPI recruitment n (adj A / B) **-**

Classical n per group **-**

Labeled saving vs classical **-**

Classic

-

PPI

-

PPI adj.

-

Formula note will appear here.

### Inputs

Difference in means to detect

Expected non-response (%)

Outcome SD, group A

Outcome SD, group B

Prediction R² in group A

Prediction R² in group B

Pilot size m (optional) If the R² in each group came from a pilot of m
labeled units, plan at its one-sided 90% lower confidence bound.

Unlabeled N_A

Unlabeled N_B

Alpha (two-sided)

Target power (%)

### Results

PPI n per group (base) **-**

PPI recruitment n per group (adj.) **-**

Classical n per group **-**

Labeled saving vs classical **-**

Classic

-

PPI

-

PPI adj.

-

Formula note will appear here.

### Inputs

Mean within-pair difference to detect

SD of the within-pair difference

R² of predicted vs true differences

Unlabeled pairs N

Alpha (two-sided)

Target power (%)

Expected non-response (%)

Pilot size m (optional) If the R² came from a pilot of m labeled units,
plan at its one-sided 90% lower confidence bound.

### Results

PPI pairs (base) **-**

PPI recruitment pairs (adj.) **-**

Classical pairs **-**

Labeled saving vs classical **-**

Classic

-

PPI

-

PPI adj.

-

Formula note will appear here.

### Inputs

Prevalence under condition A (%)

Prevalence under condition B (%)

Within-pair correlation of true outcomes

Unlabeled pairs N

Classifier sensitivity

Classifier specificity

Alpha (two-sided)

Target power (%)

Expected non-response (%)

### Results

PPI pairs (base) **-**

PPI recruitment pairs (adj.) **-**

Classical pairs **-**

Labeled saving vs classical **-**

Classic

-

PPI

-

PPI adj.

-

Formula note will appear here.

Formula details intentionally hidden here for readability. See [Variance
&
Methods](https://yiqunchen.github.io/pppower/articles/deep-dive-math.md)
for full equations.
