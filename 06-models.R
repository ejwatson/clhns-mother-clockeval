## ---------------------------
##
## Script name: 06-models.R
##
## Purpose of script: All models, figures and objects for the manuscript.
##   - Builds the analytic data (via 05-ipw.R), fits every model, builds Figures 1-3,
##     and exports journal-ready figure files (output/figures/) plus
##     output/analysis.RData, which mom-manuscript.qmd, mom-tables.qmd and
##     mom-supplement.qmd load. Run this script before rendering them.
##
## Notes:
##   - Models: weighted OLS / g-computation risk differences (marginaleffects), SEs
##     clustered on uncmomid; weights = stabilised IPCW (05-ipw.R), not trimmed.
##   - figures exported at 300 dpi, 300 mm wide (TIFF + PNG).
##
## ---------------------------

options(scipen = 6, digits = 4)

library(tidyverse)
library(sandwich)
library(lmtest)
library(marginaleffects)
library(gtsummary)

source("05-ipw.R")   # builds mom_analysis_data, analytic_data, weights, diagnostics

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

weight_var <- "ipw"   # stabilised IPCW from 05-ipw.R (not trimmed)

# ---- Variable labels / constants --------------------------------------------------------
clock_vars <- c("pchorvath1", "pchannum", "pcphenoage", "pcgrimage", "dunedinpace")
clock_lab  <- c(pchorvath1 = "PCHorvath", pchannum = "PCHannum", pcphenoage = "PCPhenoAge",
                pcgrimage = "PCGrimAge", dunedinpace = "DunedinPACE")
clock_gen  <- c(pchorvath1 = "First generation", pchannum = "First generation", pcphenoage = "Second generation",
                pcgrimage = "Second generation", dunedinpace = "Pace of ageing")
gen_cols   <- c("First generation" = "#420A68", "Second generation" = "#A52059", "Pace of ageing" = "#F37819")
clocks_z   <- paste0("z_", clock_vars)

# ---- Analysis data prep -----------------------------------------------------------------
zs <- function(x) as.numeric(scale(as.numeric(x)))

# Clocks are pre-processed as age acceleration: residuals from regressing each age-scale clock on chronological
# age in 2005 (analytic sample, unweighted). ra_* = residual in the clock's own units; z_* = standardised.
# DunedinPACE is a pace-of-ageing measure and is NOT residualised (ra_dunedinpace = raw value).
resid_age <- function(x, a) as.numeric(resid(lm(x ~ a, na.action = na.exclude)))
age_clocks <- c("pchorvath1", "pchannum", "pcphenoage", "pcgrimage")

prep <- function(df) {
  df %>%
    haven::zap_labels() %>%
    mutate(
      across(all_of(age_clocks), ~ resid_age(.x, age05), .names = "ra_{.col}"),
      ra_dunedinpace = dunedinpace,
      across(all_of(paste0("ra_", clock_vars)), zs, .names = "z_{sub('ra_', '', .col)}"),
      mmse_z = zs(mmsesum), grip_z = zs(grip_mean), tug_z = zs(mean_tug),
      hba1c_z = zs(hba1c_avg), hba1c_mmol = (hba1c_avg - 2.15) * 10.929,   # IFCC conversion: % -> mmol/mol
       sbp_z = zs(bpsysm), logcrp_z = zs(log_crp),
      iadl   = as.numeric(as.character(withiadl12) == "Yes"),
      crp_clean = if_else(crp <= 10, crp, NA_real_),
      bmi_z  = zs(mbmi05),
      # coded 1 = yes (see the table() checks below)
      smoke  = as.numeric(masmoke == 1),
      drink  = as.numeric(madrink == 1),
      edu_mom_z = zs(gradema), edu_own_z = zs(meduc83), wealth_z = zs(wealth_pc_2005),
      mm_cat4 = factor(case_when(multimorb_count == 0 ~ "No conditions",
                                 multimorb_count == 1 ~ "1 condition",
                                 multimorb_count == 2 ~ "2 conditions",
                                 multimorb_count >= 3 ~ "≥3 conditions"),
                       levels = c("No conditions", "1 condition", "2 conditions", "≥3 conditions"))
    )
}

d     <- prep(analytic_data)

n_clocks   <- sum(mom_analysis_data$has_clocks == 1)
n_dnam     <- sum(mom_analysis_data$has_dnam == 1)   # clocks + genetic PC1
n_analytic <- nrow(d)

# ---- Outcome specification --------------------------------------------------------------
outcome_spec <- tribble(
  ~outcome,    ~label,                    ~type,     ~digits, ~axis,
  "mmse_z",    "Cognition (z)",           "linear",  2, "Beta (SD)",
  "grip_z",    "Grip strength (z)",       "linear",  2, "Beta (SD)",
  "tug_z",     "Timed get up and go (z)",   "linear",  2, "Beta (SD)",
  "iadl",      "Difficulty with ≥1 IADL", "rd", 3, "Risk difference",
  "hba1c_mmol", "HbA1c (mmol/mol)",        "linear",  2, "Beta (mmol/mol)",
  "bpsysm",    "Systolic blood pressure (mmHg)",      "linear",  2, "Beta (mmHg)",
  "log_crp",   "C-reactive protein (mg/L)",              "log_pct", 1, "Percent change (%)",
  "mm_ge1",    "At least 1 condition (vs. <1)", "rd", 3, "Risk difference",
  "mm_ge2",    "At least 2 conditions (vs. <2)", "rd", 3, "Risk difference",
  "mm_ge3",    "At least 3 conditions (vs. <3)", "rd", 3, "Risk difference",
  # z-scored versions used in Figure 1
  "hba1c_z",   "HbA1c (z)",               "linear",  2, "Beta (SD)",
  "sbp_z",     "Systolic blood pressure (z)",         "linear",  2, "Beta (SD)",
  "logcrp_z",  "Log C-reactive protein (z)",             "linear",  2, "Beta (SD)"
)
functional_o <- c("mmse_z", "grip_z", "tug_z", "iadl")
biomarker_o  <- c("hba1c_mmol", "bpsysm", "log_crp")
multimorb_o  <- c("mm_ge1", "mm_ge2", "mm_ge3")
binary_o     <- c("iadl", multimorb_o)                                   # Figure 2
fig_cont_o   <- c("mmse_z", "grip_z", "tug_z", "hba1c_z", "sbp_z", "logcrp_z")  # Figure 1 (all z-scored)

get_spec <- function(o, col) {
  i <- match(o, outcome_spec$outcome)
  if (is.na(i)) c(type = "linear", digits = 2)[[col]] else outcome_spec[[col]][i]
}

# ---- Estimation -------------------------------------------------------------------------
# Weighted regression with SEs clustered on uncmomid (weights induce pseudo-replication).
#  linear  : weighted OLS; beta
#  log_pct : weighted OLS on log outcome; 100*(exp(beta)-1)
#  rd      : weighted logistic (quasibinomial) model + g-computation (marginaleffects):
#            average risk difference (probability scale) for a 1-unit (1-SD) increase in exposure, averaged
#            over the weighted analytic sample; cluster-robust variance.
fit_one <- function(data, outcome, exposure, covars, wcol, id = "uncmomid") {
  type <- get_spec(outcome, "type")
  data$.w <- data[[wcol]]
  fml <- reformulate(c(exposure, covars), response = outcome)
  cl_for <- function(fit) {
    cl <- data[[id]]
    if (!is.null(fit$na.action)) cl <- cl[-fit$na.action]
    cl
  }

  if (type == "rd") {
    fit <- glm(fml, data = data, weights = .w, family = quasibinomial())
    V   <- sandwich::vcovCL(fit, cluster = cl_for(fit), type = "HC1")
    ac  <- avg_comparisons(fit, variables = setNames(list(1), exposure),
                           newdata = fit$model, wts = fit$prior.weights, vcov = V)
    out <- tibble(est = ac$estimate, lo = ac$conf.low, hi = ac$conf.high,
                  p = ac$p.value, n = nobs(fit))
  } else {
    fit <- lm(fml, data = data, weights = .w)
    V   <- sandwich::vcovCL(fit, cluster = cl_for(fit), type = "HC1")
    ct  <- lmtest::coeftest(fit, vcov. = V)[exposure, ]
    b <- ct[[1]]; se <- ct[[2]]
    if (type == "log_pct") {
      out <- tibble(est = 100 * (exp(b) - 1), lo = 100 * (exp(b - 1.96 * se) - 1),
                    hi = 100 * (exp(b + 1.96 * se) - 1), p = ct[[4]], n = nobs(fit))
    } else {
      out <- tibble(est = b, lo = b - 1.96 * se, hi = b + 1.96 * se, p = ct[[4]], n = nobs(fit))
    }
  }
  bind_cols(tibble(exposure = exposure, outcome = outcome,
                   digits = as.integer(get_spec(outcome, "digits"))), out)
}

run_grid <- function(data, exposures, outcomes, covars, wcol) {
  expand_grid(exposure = exposures, outcome = outcomes) %>%
    pmap_dfr(function(exposure, outcome) fit_one(data, outcome, exposure, covars, wcol))
}

# Clocks are already residualised on 2005 age. Adjustment sets are cumulative:
#   model 1: age; model 2: + genetic PC1; model 3: + estimated cell proportions (CD8T, CD4T, NK, B cell, monocyte; granulocytes are the omitted reference because the six sum to 1)
# Clock -> 2012 outcome models adjust for age in 2012 (age_2012); models with a clock as the OUTCOME
# (risk factors -> clock, SES -> clock) adjust for age in 2005 (age05).
cell_pcs <- c("cd8t", "cd4t", "nk", "bcell", "mono")   # estimated white-blood-cell proportions (granulocytes omitted: the six sum to 1)
mk_cov <- function(a) list(base = a, cell = c(a, cell_pcs), pc = c(a, "pc1"), cellpc1 = c(a, "pc1", cell_pcs))
cov_o <- mk_cov("age_2012")   # clock -> outcome
cov_c <- mk_cov("age05")      # exposure -> clock

risk_vars <- c("bmi_z", "smoke", "drink")
ses_vars  <- c("edu_mom_z", "edu_own_z", "wealth_z")
all_o     <- c(functional_o, biomarker_o, multimorb_o)   # the ten 2012 outcomes in original units/binary

# Primary results (original units for biomarkers; binary outcomes as risk differences)
res_main <- run_grid(d, clocks_z, all_o, cov_o$base, weight_var)

# Figure 1: all continuous outcomes z-scored (clock -> outcome)
res_fig <- run_grid(d, clocks_z, fig_cont_o, cov_o$base, weight_var)

# (Figure 1, model 1 = res_fig)

res_xs  <- run_grid(d, risk_vars, clocks_z, cov_c$base, weight_var)
res_ses <- run_grid(d, ses_vars,  clocks_z, cov_c$base, weight_var)

# Full cumulative grids for the per-clock supplementary tables (S2-S6):
#   rows = 2012 outcomes (z-scored, binary, and biomarkers in original units), risk factors -> clock, SES -> clock
out_rows <- unique(c(fig_cont_o, binary_o, biomarker_o))
res_o <- list(base    = run_grid(d, clocks_z, out_rows, cov_o$base,    weight_var),
              pc1     = run_grid(d, clocks_z, out_rows, cov_o$pc,      weight_var),
              cellpc1 = run_grid(d, clocks_z, out_rows, cov_o$cellpc1, weight_var))
res_xs_pc1     <- run_grid(d, risk_vars, clocks_z, cov_c$pc,      weight_var)
res_xs_cellpc1 <- run_grid(d, risk_vars, clocks_z, cov_c$cellpc1, weight_var)

# Models 2-3 in Figures 1-3 (cumulative; same weights, same sample)
res_fig_pc1     <- run_grid(d, clocks_z, fig_cont_o, cov_o$pc,      weight_var)   # + genetic PC1
res_bin_pc1     <- run_grid(d, clocks_z, binary_o,   cov_o$pc,      weight_var)
res_ses_pc1     <- run_grid(d, ses_vars,  clocks_z,  cov_c$pc,      weight_var)
res_fig_cellpc1 <- run_grid(d, clocks_z, fig_cont_o, cov_o$cellpc1, weight_var)   # + cell composition
res_bin_cellpc1 <- run_grid(d, clocks_z, binary_o,   cov_o$cellpc1, weight_var)
res_ses_cellpc1 <- run_grid(d, ses_vars,  clocks_z,  cov_c$cellpc1, weight_var)

# ---- Formatting helpers ----------------------------------------------------------------
fmt_est <- function(est, lo, hi, d = 2) {
  z <- function(x) round(x, d) + 0   # avoids "-0.00"
  gsub("-", "\u2212", sprintf("%.*f (%.*f, %.*f)", d, z(est), d, z(lo), d, z(hi)), fixed = TRUE)   # true minus sign
}
fmt_est <- Vectorize(fmt_est, USE.NAMES = FALSE)

# Inline text helpers (estimates pulled directly from result tables)
ms  <- function(x, dg = 1) sprintf("%.*f", dg, mean(x, na.rm = TRUE))   # mean
sdv <- function(x, dg = 1) sprintf("%.*f", dg, sd(x, na.rm = TRUE))      # SD
pc  <- function(x, dg = 1) sprintf("%.*f%%", dg, 100 * mean(x, na.rm = TRUE))
wd  <- weight_diag_tbl
# pull a formatted estimate from a figure data frame (df_fig1/2/3) by clock, item and adjustment key
Ea <- function(df, clock, item, adj, dg = 2) {
  r <- df[df$clock == clock & df$item == item & df$adj == adj_lab[[adj]], ]
  fmt_est(r$est, r$lo, r$hi, dg)
}

# Per-unit estimates (Moqri et al. 2024 recommend reporting per SD AND per absolute unit):
# per 1 year of epigenetic age for the four age-scale clocks; per 0.1 for DunedinPACE.
# Age acceleration (residual) in original units: years for the four age-scale clocks, 0.1 PACE units for DunedinPACE
d <- d %>% mutate(ra_dunedinpace_01 = ra_dunedinpace * 10)
unit_clocks <- c("ra_pchorvath1", "ra_pchannum", "ra_pcphenoage", "ra_pcgrimage", "ra_dunedinpace_01")
res_unit <- run_grid(d, unit_clocks, c(fig_cont_o, binary_o), cov_o$base, weight_var)

# ---- Descriptives used in the text -------------------------------------------------------
cor_age <- map_dbl(setdiff(clock_vars, "dunedinpace"), ~ cor(d[[.x]], d$age05, use = "complete.obs"))
names(cor_age) <- clock_lab[setdiff(clock_vars, "dunedinpace")]
resid_df <- map_dfc(setNames(paste0("ra_", clock_vars), clock_lab[clock_vars]), ~ d[[.x]])
resid_cor <- cor(resid_df)

# ---- Figures --------------------------------------------------------------------------------
item_lab <- c(mmse_z = "Cognition", grip_z = "Grip strength", tug_z = "Timed get up and go",
              hba1c_z = "HbA1c", sbp_z = "Systolic blood pressure", logcrp_z = "Log C-reactive protein",
              iadl = "Difficulty with \u22651 IADL", mm_ge1 = "\u22651 condition",
              mm_ge2 = "\u22652 conditions", mm_ge3 = "\u22653 conditions",
              edu_mom_z = "Mother's education", edu_own_z = "Own education",
              wealth_z = "Household wealth")
adj_lab <- c(base = "Age", pc1 = "+ Genetic PC1", cellpc1 = "+ Cell composition")

# Stack the three adjustment sets for one figure. clock_col = "exposure" when the clock is the
# predictor (Figures 1-2) and "outcome" when the clock is the outcome (Figure 3).
stack_adj <- function(rs, clock_col) {
  imap_dfr(rs, function(r, nm) {
    r %>% transmute(clock = sub("^z_", "", .data[[clock_col]]),
                    item  = if (clock_col == "exposure") outcome else exposure,
                    est, lo, hi, adj = adj_lab[[nm]])
  })
}

# Figure 1: clock -> continuous 2012 outcomes (z-scored; SD scale)
df_fig1 <- stack_adj(list(base = res_fig,
                          pc1 = res_fig_pc1,
                          cellpc1 = res_fig_cellpc1), "exposure")
# Figure 2: clock -> binary 2012 outcomes (marginal risk difference, probability scale)
df_fig2 <- stack_adj(list(base = res_main %>% filter(outcome %in% binary_o),
                          pc1 = res_bin_pc1,
                          cellpc1 = res_bin_cellpc1), "exposure")
# Figure 3: SES -> clock
df_fig3 <- stack_adj(list(base = res_ses, pc1 = res_ses_pc1, cellpc1 = res_ses_cellpc1), "outcome")

# One row per measure; five clock panels sharing a common x-axis. Each estimate is drawn
# three times (dodged): age (circle), + genetic PC1 (triangle), + cell composition (square);
# adjustments are cumulative; all but the first are semi-transparent.
# Sizes are set for a 300 mm-wide export; text stays legible when scaled to a journal page.
adj_plot <- function(df, items, xlab, breaks = waiver(), base = 15, pt = 3.2, lw = 0.7,
                     strip = 14, ytxt = 13, xtxt = 12, xtitle = 14, leg = 13,
                     leg_tight = FALSE, us = FALSE) {
  df <- df %>%
    filter(item %in% items) %>%
    mutate(clock_f = factor(paste0(clock_lab[clock], " (z)"),
                            levels = paste0(clock_lab[clock_vars], " (z)")),
           gen     = factor(clock_gen[clock], levels = names(gen_cols)),
           item_f  = factor(item_lab[item], levels = rev(unname(item_lab[items]))),
           adj     = factor(adj, levels = unname(adj_lab)))
  dodge <- position_dodge(width = 0.85)
  ggplot(df, aes(x = est, y = item_f, xmin = lo, xmax = hi, colour = gen,
                 shape = adj, alpha = adj, group = adj)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_linerange(orientation = "y", linewidth = lw, position = dodge) +
    geom_point(size = pt, position = dodge) +
    facet_wrap(~ clock_f, nrow = 1) +                       # shared x-axis across the five panels
    scale_colour_manual(values = gen_cols, name = "Clock generation:",
                        labels = if (us) function(x) sub("ageing", "aging", x) else waiver()) +
    scale_shape_manual(values = setNames(c(16, 17, 15), unname(adj_lab)), name = "Adjustment:") +
    scale_alpha_manual(values = setNames(c(1, 0.55, 0.55), unname(adj_lab)), name = "Adjustment:") +
    guides(shape = guide_legend(override.aes = list(colour = "grey30")),
           alpha = guide_legend(override.aes = list(colour = "grey30"))) +
    scale_x_continuous(breaks = breaks) +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = base) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = strip),
          axis.text.y = element_text(face = "bold", size = ytxt),
          axis.text.x = element_text(angle = 45, hjust = 1, size = xtxt),
          axis.title.x = element_text(size = xtitle),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.text = element_text(size = leg),
          legend.title = element_text(size = leg),
          legend.spacing.y = if (leg_tight) unit(0, "pt"),
          legend.margin = if (leg_tight) margin(0, 0, 0, 0),
          legend.box.spacing = if (leg_tight) unit(3, "pt"))
}

figure1 <- adj_plot(df_fig1, fig_cont_o, "Standardised estimate (SD) \u00b1 95% CI")
figure2 <- adj_plot(df_fig2, binary_o,   "Risk difference \u00b1 95% CI", breaks = c(0, 0.05, 0.10))
figure3 <- adj_plot(df_fig3, ses_vars,   "Standardised estimate (SD) \u00b1 95% CI")

# ---- Export figures for the journal ---------------------------------------------------------
# IJE: separate image files (TIFF/PNG/EPS/JPG), >= 300 dpi at >= 300 mm wide (3600 px).
save_fig <- function(p, name, h_mm) {
  ggsave(file.path("output/figures", paste0(name, ".tiff")), p, width = 305, height = h_mm,
         units = "mm", dpi = 300, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path("output/figures", paste0(name, ".png")), p, width = 305, height = h_mm,
         units = "mm", dpi = 300, bg = "white")
}
save_fig(figure1, "Figure1", 210)
save_fig(figure2, "Figure2", 150)
save_fig(figure3, "Figure3", 140)

# Journal of Gerontology Series A: figures at final size (170 mm = double column; text >= 8 pt; 300 dpi)
dir.create("output/figures/seriesA", recursive = TRUE, showWarnings = FALSE)
sa_args <- list(base = 8, pt = 1.6, lw = 0.4, strip = 8, ytxt = 8, xtxt = 8, xtitle = 8, leg = 8, leg_tight = TRUE, us = TRUE)
save_fig_sa <- function(p, name, h_mm) {
  ggsave(file.path("output/figures/seriesA", paste0(name, ".tiff")), p, width = 170, height = h_mm,
         units = "mm", dpi = 300, device = "tiff", compression = "lzw", bg = "white")
}
save_fig_sa(do.call(adj_plot, c(list(df_fig1, fig_cont_o, "Standardized estimate (SD) \u00b1 95% CI"), sa_args)), "Figure1", 150)
save_fig_sa(do.call(adj_plot, c(list(df_fig2, binary_o, "Risk difference \u00b1 95% CI", breaks = c(0, 0.05, 0.10)), sa_args)), "Figure2", 110)
save_fig_sa(do.call(adj_plot, c(list(df_fig3, ses_vars, "Standardized estimate (SD) \u00b1 95% CI"), sa_args)), "Figure3", 100)

# Supplementary figures (also exported; embedded in the supplement)
supp_fig1 <- as.data.frame(as.table(resid_cor)) %>%
  rename(c1 = Var1, c2 = Var2, r = Freq) %>%
  mutate(c1 = factor(c1, levels = clock_lab), c2 = factor(c2, levels = rev(clock_lab))) %>%
  ggplot(aes(c1, c2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4) +
  scale_fill_gradient2(low = "#420A68", mid = "white", high = "#F37819", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/figures/SupplementaryFigureS1.png", supp_fig1, width = 160, height = 140, units = "mm", dpi = 300, bg = "white")

# ---- Save everything for the Word documents ------------------------------------------------
dir.create("output", showWarnings = FALSE)
save(list = ls(all.names = FALSE), file = "output/analysis.RData")
message("Saved output/analysis.RData and output/figures/*")
