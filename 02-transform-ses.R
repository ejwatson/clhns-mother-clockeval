## ---------------------------
##
## Script name: 02-transform-ses.R
##
## Purpose of script: Cleaning and constructing SES variables
##
## Author: Elijah Watson
##
## Date Created: May 13, 2024
##
##
## ---------------------------
##
## Notes:
##
## ---------------------------

options(scipen = 6, digits = 4)

## ---------------------------

library(tidyverse)

source("01-import.R")

#selecting relevant variables 
mom_analysis_data <- icclongifile %>%
  select(mage83, dhinc83, dhinc05, urbanic83, basewman, basebrgy, meduc83, mbmi05, urbanic05, mage83, mage05, mhght83, icclongifile_observed, uncchdid) %>%
  mutate(mage83_years = mage83 / 12,
         # age in 2005 (years); mage05 is recorded in months
         age05 = mage05 / 12)

# Check that mage05 is in months (mean age in 2005 should be about 48 years)
if (!isTRUE(dplyr::between(mean(mom_analysis_data$age05, na.rm = TRUE), 45, 52))) {
  warning("age05 mean = ", round(mean(mom_analysis_data$age05, na.rm = TRUE), 1),
          "; check units of mage05 (expected months; mean age in 2005 should be ~48).")
}

mom_analysis_data <- left_join(mom_analysis_data,
                               mother_person_wave0 %>% select(yearbrth, basebrgy, basewman),
                               by = c("basebrgy","basewman"))

##transform early life SES ----
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_2002_decision %>% select(gradepa,
                                                              occupa,
                                                              gradema,
                                                              occuma,
                                                              gradegpa,
                                                              occugpa,
                                                              gradegma,
                                                              occugma,
                                                              basebrgy,
                                                              basewman,
                                                              mom_2002_decision_observed) %>%
                                   rename(gradepa_02 = gradepa,
                                          occupa_02 = occupa,
                                          gradema_02 = gradema,
                                          occuma_02 = occuma,
                                          gradegpa_02 = gradegpa,
                                          occugpa_02 = occugpa,
                                          gradegma_02 = gradegma,
                                          occugma_02 = occugma),
                                 by = c("basebrgy","basewman")) 

#those who didnt answer in 2002 were asked in 2005:
mom_analysis_data <-  left_join(mom_analysis_data,
                                 mom_2005_decision %>% select(gradepa,
                                                              occupa,
                                                              gradema,
                                                              occuma,
                                                              gradegpa,
                                                              occugpa,
                                                              gradegma,
                                                              occugma,
                                                              basebrgy,
                                                              basewman, 
                                                              mom_2005_decision_observed) %>%
                                   rename(gradepa_05 = gradepa,
                                          occupa_05 = occupa,
                                          gradema_05 = gradema,
                                          occuma_05 = occuma,
                                          gradegpa_05 = gradegpa,
                                          occugpa_05 = occugpa,
                                          gradegma_05 = gradegma,
                                          occugma_05 = occugma),
                                 by = c("basebrgy","basewman")) 


#create variable for mother's mother's years of education 
mom_analysis_data <- mom_analysis_data %>%
  mutate(gradema_02 = case_when(gradema_02 < 0 ~ NA,
                                gradema_02 >= 0 ~ gradema_02),
         gradema_05 = case_when(gradema_05 < 0 ~ NA,
                                gradema_05 >= 0 ~ gradema_05),
         gradepa_02 = case_when(gradepa_02 < 0 ~ NA,
                                gradepa_02 >= 0 ~ gradepa_02),
         gradepa_05 = case_when(gradepa_05 < 0 ~ NA,
                                gradepa_05 >= 0 ~ gradepa_05),
         gradema = case_when(!is.na(gradema_02) ~ gradema_02,
                             is.na(gradema_02) ~ gradema_05,
                             is.na(gradema_02) & is.na(gradema_05) ~ NA),
         gradepa = case_when(!is.na(gradepa_02) ~ gradepa_02,
                             is.na(gradepa_02) ~ gradepa_05,
                             is.na(gradepa_02) & is.na(gradepa_05) ~ NA))


#wealth
mom_analysis_data <- left_join(mom_analysis_data,
                               wealth_pc_1983 %>% 
                                 rename(wealth_pc_1983 = pc) %>%
                                 select(uncchdid, wealth_pc_1983),
                               by = c("uncchdid"))

mom_analysis_data <- left_join(mom_analysis_data,
                               wealth_pc_2005 %>%
                                 rename(wealth_pc_2005 = pc) %>%
                                 select(uncchdid, wealth_pc_2005),
                               by = c("uncchdid"))

mom_analysis_data <- left_join(mom_analysis_data,
                               wealth_pc_2002 %>%
                                 rename(wealth_pc_2002 = pc) %>%
                                 select(uncchdid, wealth_pc_2002),
                               by = c("uncchdid"))


#also importing health behaviors here
mom_analysis_data <- left_join(mom_analysis_data,
                               mom_2005_morb %>%
                                 select(masmoke, madrink, basebrgy, basewman),
                               by = c("basebrgy","basewman"))


