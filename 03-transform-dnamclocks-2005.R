## ---------------------------
##
## Script name: 03-transform-dnamclocks-2005.R
##
## Purpose of script: Cleaning and constructing 2005 DNAm clock data
##
## Author: Elijah Watson
##
## Date Created: June 17, 2024
##
##
## ---------------------------
##
## Notes: 
##
## ---------------------------

options(scipen = 6, digits = 4)

## ---------------------------

source("02-transform-ses.R")

# merging epigenetic clock data  ------------------------------------------
# Replicate samples are removed; non-replicate samples are used
mom_clock_2005 %>%
  count(basebrgy,basewman) %>% 
  filter(n>1)

mom_clock_2005 %>%
  filter((basebrgy == 1 & basewman==59) |
           (basebrgy == 3 & basewman==139) |
           (basebrgy == 3 & basewman==243) |
           (basebrgy == 3 & basewman==248) |
           (basebrgy == 8 & basewman==69) |
           (basebrgy == 9 & basewman==79) |
           (basebrgy == 10 & basewman==88) |
           (basebrgy == 11 & basewman==137) |
           (basebrgy == 13 & basewman==430) |
           (basebrgy == 13 & basewman==510) |
           (basebrgy == 17 & basewman==215) |
           (basebrgy == 29 & basewman==15))

#merging non-replicate mother clock data
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_clock_2005 %>%
                                   filter(!((basebrgy == 1 & basewman==59) |
                                              (basebrgy == 3 & basewman==139) |
                                              (basebrgy == 3 & basewman==243) |
                                              (basebrgy == 3 & basewman==248) |
                                              (basebrgy == 8 & basewman==69) |
                                              (basebrgy == 9 & basewman==79) |
                                              (basebrgy == 10 & basewman==88) |
                                              (basebrgy == 11 & basewman==137) |
                                              (basebrgy == 13 & basewman==430) |
                                              (basebrgy == 13 & basewman==510) |
                                              (basebrgy == 17 & basewman==215) |
                                              (basebrgy == 29 & basewman==15))) %>% 
                                   select(cd8t, cd4t, nk, bcell, mono, gran,
                                          pchorvath1, pchannum, pcphenoage, pcgrimage, dunedinpace, 
                                          uncmomid, basebrgy, basewman),
                                 by = c("basebrgy", "basewman"))


#Time of blood sample collection
mom_analysis_data <- left_join(mom_analysis_data,
                               mom_blood_2005 %>%
                                 select(timeblood, uncid) %>%
                                 rename(uncmomid = uncid),
                               by=c("uncmomid"))

# Age at the 2005 blood draw / DNAm assay, taken from the clock file or the blood file
# (measured for every woman with clock data); the icclongifile-based age (mage05/12) is used only
# where neither is available. If the age column is not detected correctly, set
# `age_var_clock` / `age_var_blood` below.
age_var_clock <- NULL   # e.g. "age"
age_var_blood <- NULL   # e.g. "age"
pick_age <- function(df, override) {
  if (!is.null(override)) return(override)
  nm <- grep("^age|age$", names(df), ignore.case = TRUE, value = TRUE)
  nm <- setdiff(nm, c("age05", "age_2012"))
  nm <- nm[!grepl("accel|pace|dnam|grim|pheno|horvath|hannum|dunedin|resid", nm, ignore.case = TRUE)]  # exclude clock columns
  if (length(nm)) nm[1] else NA_character_
}
age_var_clock <- pick_age(mom_clock_2005, age_var_clock)
age_var_blood <- pick_age(mom_blood_2005, age_var_blood)

if (!is.na(age_var_clock)) {
  mom_analysis_data <- left_join(mom_analysis_data,
                                 mom_clock_2005 %>%
                                   filter(!((basebrgy == 1 & basewman==59) | (basebrgy == 3 & basewman==139) |
                                            (basebrgy == 3 & basewman==243) | (basebrgy == 3 & basewman==248) |
                                            (basebrgy == 8 & basewman==69) | (basebrgy == 9 & basewman==79) |
                                            (basebrgy == 10 & basewman==88) | (basebrgy == 11 & basewman==137) |
                                            (basebrgy == 13 & basewman==430) | (basebrgy == 13 & basewman==510) |
                                            (basebrgy == 17 & basewman==215) | (basebrgy == 29 & basewman==15))) %>%
                                   transmute(basebrgy, basewman, age05_clockfile = as.numeric(.data[[age_var_clock]])),
                                 by = c("basebrgy", "basewman"))
} else mom_analysis_data$age05_clockfile <- NA_real_

if (!is.na(age_var_blood)) {
  mom_analysis_data <- left_join(mom_analysis_data,
                                 mom_blood_2005 %>% transmute(uncmomid = uncid, age05_bloodfile = as.numeric(.data[[age_var_blood]])) %>%
                                   distinct(uncmomid, .keep_all = TRUE),
                                 by = "uncmomid")
} else mom_analysis_data$age05_bloodfile <- NA_real_

# If a file stores age in months, convert (mean age in 2005 is ~48 years)
to_years <- function(x) if (isTRUE(mean(x, na.rm = TRUE) > 150)) x / 12 else x
mom_analysis_data <- mom_analysis_data %>%
  mutate(age05_icc = age05,
         age05 = coalesce(to_years(age05_clockfile), to_years(age05_bloodfile), age05_icc))

#merging genetic pcs
mom_analysis_data <-  left_join(mom_analysis_data,
                                mom_geneticpc %>% select(pc1, uncmomid) %>% drop_na(),
                                by = c("uncmomid"))



# 2005 systolic BP (mean of up to 3 readings; negative codes = NR/NA) -- used in 2012 attrition weights
mom_analysis_data <- left_join(mom_analysis_data,
                               mom_anthdiet_2005 %>% select(-mom_anthdiet_2005_observed),
                               by = "uncmomid")

mom_analysis_data <- mom_analysis_data %>%
  mutate(across(c(BPSYST1, BPSYST2, BPSYST3), ~ if_else(.x < 0, NA_real_, as.numeric(.x)))) %>%
  mutate(sbp05 = rowMeans(across(c(BPSYST1, BPSYST2, BPSYST3)), na.rm = TRUE),
         sbp05 = if_else(is.nan(sbp05), NA_real_, sbp05))
