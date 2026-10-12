## ---------------------------
##
## Script name: 05-ipw.R
##
## Purpose of script: Inverse probability of censoring weights (IPCW) for
##                    (1) selection into the 2005 DNAm sub-sample and
##                    (2) attrition from the 2005 DNAm sample to the 2012
##                        analytic sample (complete 2012 outcomes).
##                    Each is estimated with weightit() + .cens(); the final weight
##                    is the product of the two (sequential censoring).
##
## Author: Elijah Watson
##
## Notes:
##  - Method: covariate balancing propensity scores (CBPS).
##  - Stage 1 (censoring = no 2005 DNAm clock data or no genetic PC1): baseline (1983) covariates; all women.
##  - Stage 2 (censoring = clock data but incomplete 2012 outcomes): 2005 covariates only;
##    women with clock data only (the at-risk set).
##  - Stabilization: each stage's weight is multiplied by the marginal probability of being
##    uncensored in that stage (numerator of the stabilized IPCW).
##  - The analytic sample IS restricted to women with genetic PC1 (pc1). Having 2005 clock
##    data AND pc1 defines being "uncensored" at stage 1, so the weights account for it.
##  - Requires a WeightIt version that supports .cens() in weightit().
##
## ---------------------------

options(scipen = 6, digits = 4)

library(tidyverse)
library(WeightIt)
library(cobalt)

source("04-transform-outcomes.R")

clock_vars <- c("pchorvath1", "pchannum", "pcphenoage", "pcgrimage", "dunedinpace")

# Define sample indicators and censoring indicators ---------------------------------------
mom_analysis_data <- mom_analysis_data %>%
  mutate(
    # 2005 age (years)
    age = age05,

    has_clocks = as.integer(if_all(all_of(clock_vars), ~ !is.na(.x))),

    # uncensored at stage 1 = all five clocks AND genetic PC1 observed
    has_dnam = as.integer(has_clocks == 1 & !is.na(pc1)),

    outcomes_complete = as.integer(
      !is.na(mmsesum) & !is.na(grip_mean) & !is.na(mean_tug) & !is.na(withiadl12) &
        !is.na(no_conditions) & !is.na(one_condition) & !is.na(two_or_more_conditions) &
        !is.na(hba1c_avg) & !is.na(bpsysm) & !is.na(log_crp)
    ),

    analytic_sample = as.integer(has_dnam == 1 & outcomes_complete == 1),

    # censoring indicators (1 = censored)
    cens_dnam = 1L - has_dnam,                                     # time 1
    cens_2012 = if_else(has_dnam == 1, 1L - outcomes_complete,     # time 2 (only defined if at risk)
                        NA_integer_)
  )

# Weight models ----------------------------------------------------------------------------
# Time 1: selection into the 2005 DNAm sub-sample (baseline covariates)
f_t1 <- .cens(cens_dnam) ~ wealth_pc_1983 + mage83_years + meduc83 + urbanic83

# Time 2: attrition to the 2012 analytic sample (2005 covariates only)
f_t2 <- .cens(cens_2012) ~ wealth_pc_2005 + urbanic05 + dhinc05 + mbmi05 + sbp05 +
  pchorvath1 + pchannum + pcphenoage + pcgrimage + dunedinpace

# Stage 1: all women in the baseline file. Censored (no DNAm) women get weight 0.
W1 <- weightit(f_t1, data = mom_analysis_data, method = "cbps",
               missing = "ind")   # missing-indicator approach for incomplete covariates

# Stage 2: only women with clock data (the at-risk set). Censored (incomplete 2012) get weight 0.
# Explicit row key so the (smaller) stage 2 weights merge back to the right women
mom_analysis_data <- mom_analysis_data %>% mutate(.row_id = row_number())

data_stage2 <- mom_analysis_data %>% filter(has_dnam == 1)

W2 <- weightit(f_t2, data = data_stage2, method = "cbps",
               missing = "ind")

# Stabilize: multiply by marginal probability of remaining uncensored in each stage
p_uncens1 <- mean(mom_analysis_data$has_dnam == 1)
p_uncens2 <- mean(data_stage2$outcomes_complete == 1)

# Stage 1 weights are on the full data (same row order as mom_analysis_data)
stopifnot(length(W1$weights) == nrow(mom_analysis_data))
stage1_w <- tibble(.row_id = mom_analysis_data$.row_id, ipw_stage1 = W1$weights * p_uncens1)

# Stage 2 weights are on the at-risk subset: attach by row key, not position
stopifnot(length(W2$weights) == nrow(data_stage2))
stage2_w <- tibble(.row_id = data_stage2$.row_id, ipw_stage2 = W2$weights * p_uncens2)

# Attach weights ---------------------------------------------------------------------------
# Women without clock data have no stage 2 weight (not at risk) -> 0.
# Women censored at either stage receive weight 0, so only analytic-sample women are positive.
mom_analysis_data <- mom_analysis_data %>%
  left_join(stage1_w, by = ".row_id") %>%
  left_join(stage2_w, by = ".row_id") %>%
  mutate(ipw_stage2 = coalesce(ipw_stage2, 0),
         ipw        = ipw_stage1 * ipw_stage2)

# Merge checks
stopifnot(
  "Stage 2 weights must be 0 for women without clock/PC1 data" =
    all(mom_analysis_data$ipw_stage2[mom_analysis_data$has_dnam == 0] == 0),
  "No weights should be NA" = !anyNA(mom_analysis_data$ipw)
)
stopifnot(
  "Positive weights should equal analytic sample" =
    sum(mom_analysis_data$ipw > 0, na.rm = TRUE) == sum(mom_analysis_data$analytic_sample)
)

# Weight diagnostics (weights are not trimmed) ---------------------------------------------
analytic_w <- mom_analysis_data %>% filter(analytic_sample == 1)

ess <- function(w) sum(w)^2 / sum(w^2)

weight_diag_tbl <- tibble(
  n      = nrow(analytic_w),
  mean   = mean(analytic_w$ipw),
  sd     = sd(analytic_w$ipw),
  min    = min(analytic_w$ipw),
  p01    = quantile(analytic_w$ipw, .01),
  p25    = quantile(analytic_w$ipw, .25),
  median = median(analytic_w$ipw),
  p75    = quantile(analytic_w$ipw, .75),
  p99    = quantile(analytic_w$ipw, .99),
  max    = max(analytic_w$ipw),
  ESS    = ess(analytic_w$ipw)
)

# Analytic datasets -------------------------------------------------------------------------
# Analytic sample: clock data + genetic PC1 + complete 2012 outcomes
analytic_data <- mom_analysis_data %>% filter(analytic_sample == 1)
