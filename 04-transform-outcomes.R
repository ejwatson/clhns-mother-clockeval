## ---------------------------
##
## Script name: 04-transform-outcomes.R
##
## Purpose of script: Cleaning and constructing 2012 aging outcome data
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

source("03-transform-dnamclocks-2005.R")

#merging mom's age in 2012
mom_analysis_data <-  left_join(mom_analysis_data,
                                mom_person_2012 %>% 
                                  select(age, uncmomid) %>%
                                  rename(age_2012 = age) %>%
                                  mutate(age_2012 = as.double(age_2012)),
                                by = c("uncmomid")) 

#merging mom hba1c data (2012)
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_hba1c_2012 %>% select(nycoread1, nycoread2, nycoread3, uncmomid),
                                 by = c("uncmomid")) 

mom_analysis_data <- mom_analysis_data %>%
  mutate(hba1c_avg = ((nycoread1 + nycoread2 + nycoread3) / 3))

#merging mom MMSE data 
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mmse_corrected_2012 %>% select(uncmomid, mmsesum),
                                 by = c("uncmomid"))

#mothers BP in 2012
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_bp_2012 %>% select(uncmomid, bpsysm),
                                 by = c("uncmomid"))

# merging mom getupgo data (2012)
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_getupgo_2012 %>% select(timetest1, timetest2, timetest3, uncmomid),
                                 by = c("uncmomid"))

mom_analysis_data <- mom_analysis_data %>%
  mutate(mean_tug = (timetest1 + timetest2 + timetest3) / 3)


# merging mom gripstrength data (2012)
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_gripstrength_2012 %>% select(handed, uncmomid, rightgrip1, rightgrip2, rightgrip3, leftgrip1, leftgrip2, leftgrip3),
                                 by = c("uncmomid"))


mom_analysis_data <- mom_analysis_data %>%
  mutate(grip_mean = case_when(handed == 1 ~ ((rightgrip1 + rightgrip2 + rightgrip3) / 3),
                               handed == 2 ~ ((leftgrip1 + leftgrip2 + leftgrip3) / 3),
                               handed == 3 ~ ((rightgrip1 + rightgrip2 + rightgrip3 + leftgrip1 + leftgrip2 + leftgrip3) / 6)))   # ambidextrous: mean of all six trials

# merging mom hearing data (2012)
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_hearing_2012 %>% select(-basebrgy, -basewman, -basehhno),
                                 by = c("uncmomid"))


# merging mom lipid data (2012)
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_lipid_2012 %>% select(-basebrgy, -basewman, -basehhno),
                                 by = c("uncmomid"))

# merging mom iadl data (2012)
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_iadl_2012 %>% select(-basebrgy, -basewman),
                                 by = c("uncmomid"))

mom_analysis_data <- mom_analysis_data %>%
  mutate(withiadl12 = factor(withiadl12, 
                                levels = c(0, 1), 
                                labels = c("No", "Yes")))


# merging mom adl data (2012)
mom_analysis_data <-  left_join(mom_analysis_data,
                                mom_adl_2012 %>% select(-basebrgy, -basewman),
                                by = c("uncmomid"))

# merging mom morbidity 1 data (2012)
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_morbidity1_2012 %>% select(-basebrgy, -basewman, -basehhno),
                                 by = c("uncmomid"))

mom_analysis_data <- mom_analysis_data %>%
  mutate(multimorb_count = diabete + heartdis + cancer + hblood + anemia + arthriti,
         multimorb_cat = case_when(
           multimorb_count == 0 ~ 0,
           multimorb_count == 1 ~ 1,
           multimorb_count > 1 ~ 2
         ))

# Create binary variables for each multimorbidity category
mom_analysis_data <- mom_analysis_data %>%
  mutate(
    no_conditions = ifelse(multimorb_cat == 0, "Yes", "No"),
    one_condition = ifelse(multimorb_cat == 1, "Yes", "No"),
    two_or_more_conditions = ifelse(multimorb_cat == 2, "Yes", "No")
  )


# Cumulative multimorbidity thresholds (Figure 3): at least 1 / 2 / 3 conditions
mom_analysis_data <- mom_analysis_data %>%
  mutate(
    mm_ge1 = as.numeric(multimorb_count >= 1),
    mm_ge2 = as.numeric(multimorb_count >= 2),
    mm_ge3 = as.numeric(multimorb_count >= 3)
  )

# Convert binary and categorical variables to factors
mom_analysis_data <- mom_analysis_data %>%
  mutate(
    withiadl12 = as.factor(withiadl12),
    multimorb_cat = as.factor(multimorb_cat),
  )


# continuous mother aging outcomes
mom_analysis_data <- mom_analysis_data %>%
  mutate(mean_tug = case_when(mean_tug < 0 ~ NA,
                              mean_tug > 0 ~ mean_tug),
         grip_mean = case_when(grip_mean < 0 ~ NA,
                               grip_mean > 0 ~ grip_mean * 0.45359237),   # recorded in lb -> kg (SI)
         hba1c_avg = case_when(hba1c_avg < 0 ~ NA,
                               hba1c_avg > 0 ~ hba1c_avg))

# negative values are missing-value codes
mom_analysis_data <- mom_analysis_data %>% mutate(bpsysm = if_else(bpsysm > 0, bpsysm, NA_real_))

# merging mom crp data (2012)
mom_analysis_data <-  left_join(mom_analysis_data,
                                mom_crp_2012 %>% select(uncmomid, crppreflu, crpnonflu),
                                by = c("uncmomid"))

mom_analysis_data <- mom_analysis_data %>%
  mutate(crp = if_else(!is.na(crppreflu), crppreflu, crpnonflu)) %>%
  mutate(log_crp = if_else(crp <= 10, log(crp), NA_real_))
