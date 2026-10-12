# Longitudinal validation of epigenetic clocks among middle-aged Filipino women

Analysis code for "A longitudinal validation of epigenetic clocks among middle-aged Filipino women" (Cebu Longitudinal Health and Nutrition Survey, CLHNS). Epigenetic clocks measured in 2005 are related to ageing outcomes measured in 2012, with inverse probability of censoring weights for selection and attrition.

Code: https://github.com/ejwatson/clhns-mother-clockeval

## Data

The CLHNS data are not publicly available (privacy and ethical restrictions) and are not included. Scripts read the source files from a sibling `../Data/` folder; paths are set in `01-import.R`, which also sets the working directory (edit `setwd()` for your machine).

## Scripts

Run in order; each script sources the one before it, so running `06-models.R` runs the whole pipeline.

| Script | Purpose |
|---|---|
| `01-import.R` | Imports the source data files |
| `02-transform-ses.R` | Age, early-life and household socioeconomic variables, health behaviours |
| `03-transform-dnamclocks-2005.R` | Merges 2005 DNA methylation clocks, cell proportions, genetic PC1, 2005 age and blood pressure |
| `04-transform-outcomes.R` | Constructs the 2012 outcomes (cognition, grip strength, Timed Up and Go, HbA1c, blood pressure, CRP, IADL, multimorbidity) |
| `05-ipw.R` | Stabilised inverse probability of censoring weights (CBPS, two stages: DNAm sub-sample selection, then 2012 attrition); weight diagnostics |
| `06-models.R` | Age-acceleration residuals, all models, figures and the objects used by the documents; writes `output/analysis.RData` and `output/figures/` |

## Documents

Render after `06-models.R` (each document sources it if `output/analysis.RData` is missing):

```
quarto render mom-manuscript.qmd --to docx   # main text
quarto render mom-tables.qmd --to docx       # Table 1
quarto render mom-supplement.qmd --to docx   # Supplementary tables and figure
```

All three use `ije-reference.docx` for Word styles. Main-text word counts on the title page are read from `output/wordcount.csv` if present.

`06-models.R` also writes figures at the Journal of Gerontology Series A size (170 mm wide, 8 pt text) to `output/figures/seriesA/`.

## Methods in brief

- Clocks: PCHorvath, PCHannum, PCPhenoAge, PCGrimAge and DunedinPACE. The four age-scale clocks are regressed on age in 2005 and the residual is standardised; DunedinPACE is standardised only.
- Models: weighted linear regression for continuous outcomes; weighted logistic regression with g-computation (`marginaleffects`) for risk differences in binary outcomes. Standard errors are clustered on participant.
- Adjustment: models of 2012 outcomes adjust for age in 2012; models with a clock as the outcome adjust for age in 2005. Estimates are cumulative: age; plus genetic PC1; plus estimated cell proportions.

## Requirements

R with `tidyverse`, `haven`, `WeightIt` (version supporting `.cens()` in `weightit()`), `cobalt`, `sandwich`, `lmtest`, `marginaleffects`, `gtsummary` and `flextable`; Quarto for the documents.
