# Interactive PPI Sample Size Calculator

This calculator uses **prediction-powered inference (PPI/PPI++)**
formulas from this package for quick study planning examples.

Interactive PPI Tool

## Prediction-powered sample size calculator

Uses pppower-style moments and lambda logic. Outputs are labeled sample
sizes.

PPI one-sample mean / prevalence

PPI difference in prevalence

### Inputs

Outcome type Prevalence (binary) Mean (continuous)

Lambda mode Oracle (PPI++) Vanilla (lambda = 1)

Effect size delta

Unlabeled N

Alpha (two-sided)

Target power (%)

Binary prevalence inputs

Prevalence p = P(Y=1) (%)

Sensitivity

Specificity

Non-response (%)

Continuous mean moments

Var(Y) = sigma_y2

Var(f) = sigma_f2

Cov(Y,f) = cov_y_f

Non-response (%)

### Results

PPI labeled n (base) **-**

PPI labeled n (adjusted) **-**

Classical n (no prediction) **-**

Labeled saving vs classical **-**

Classic

-

PPI

-

PPI adj

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

Non-response (%)

### Results

PPI n per group (base) **-**

PPI n per group (adjusted) **-**

Classical n per group **-**

Labeled saving vs classical **-**

Classic

-

PPI

-

PPI adj

-

Formula note will appear here.

### PPI formulas used in this page

**One-sample oracle variance:** Var = sigma_y2 / n - (cov_yf^2 /
sigma_f2) \* N / (n \* (n + N))

**One-sample vanilla variance:** Var = sigma_y2 / n + sigma_f2 / N +
sigma_f2 / n - 2 \* cov_yf / n

**Binary moments from sens/spec:** cov_yf = p \* (1-p) \* (sens + spec -
1)

**Two-group prevalence:** power from the PPI two-sample variance with
group-specific oracle lambda, solved for n by search.

### Notes

- This page is now **PPI-specific** (not a generic textbook calculator).
- It is designed for planning and intuition; for final protocol numbers,
  confirm assumptions with your statistician.
