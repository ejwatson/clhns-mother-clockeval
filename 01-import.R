## ---------------------------
##
## Script name: 01-import
##
## Purpose of script: Importing CLHNS data
##
## Author: Elijah Watson
##
## Date Created: October 13, 2023
##
##
## ---------------------------
##
## Notes: Only importing data in this script.
##        Data are tidied and transformed in scripts 02-04.
## ---------------------------

options(scipen = 6, digits = 4)

## ---------------------------

setwd("~/Library/CloudStorage/Dropbox/clhns-dnamclockeval-moms/Analysis")

library(tidyverse)
library(haven)

# Importing baseline and postpartum wave data ------------------------------------------
#Importing mother data
mbase2 <-  as_tibble(read_dta("../Data/1983-1986/mbase2.dta")) %>%
  mutate(mbase2_observed = 1)

#Importing mother person data from baseline
mother_person_wave0 <-  as_tibble(read_dta("../Data/1983-1986/person.dta")) %>%
  filter(wave==0 & RELNPRW1==30) %>%
  mutate(person_observed = 1)

# Importing longitudinal data with previously constructed variables  --------
icclongifile <- as_tibble(read_dta("../Data/icclongifile.dta")) %>%
  mutate(icclongifile_observed = 1)

# Importing 2005 survey Mother data -------------------------------------------
# 2005 mother morbidity data
mom_2005_morb <- as_tibble(read_dta("../Data/2005/Survey/morbidty.dta")) %>%
  mutate(mom_2005_morb_observed = 1)

# 2005 and 2002 mother decision data (data on her early life SES)
mom_2005_decision <- as_tibble(read_dta("../Data/2005/Survey/decision.dta")) %>%
  mutate(mom_2005_decision_observed = 1)

mom_2002_decision <- as_tibble(read_dta("../Data/2002/decision.dta")) %>%
  mutate(mom_2002_decision_observed = 1)


# 2005 mother hh data data
mom_2005_hhold <- as_tibble(read_dta("../Data/2005/Survey/hhold.dta")) %>%
  mutate(mom_2005_hhold_observed = 1)


# Importing 2005 Epigenetic clock data for Mothers -----------
mom_clock_2005 <- as_tibble(read_dta("../Data/2005/Epigenetic Clocks/Mothers clock data.dta"))  %>%
  mutate(mom_clock_2005_observed = 1)

mom_blood_2005 <- as_tibble(read_dta("../Data/2005/Epigenetic Clocks/BLOOD.dta"))  %>%
  mutate(mom_blood_2005_observed = 1)

# 2005 mother anthropometry (used for 2005 systolic BP; BPSYST1-3, negative codes = missing)
mom_anthdiet_2005 <- as_tibble(read_dta("../Data/2005/Survey/anthdiet.dta")) %>%
  select(uncmomid, BPSYST1, BPSYST2, BPSYST3) %>%
  mutate(mom_anthdiet_2005_observed = 1)

mom_geneticpc <- as_tibble(read_dta("../Data/2021-01-24 Calen mom 2016 and 2015 blood TL and date- fake TL.dta"))   %>%
                             mutate(mom_geneticpc_observed = 1)

# Importing post-2005 Mother data for clock evaluation analysis -------------------------------------------
# mom hba1c data (2012)
mom_hba1c_2012 <- as_tibble(read_dta("../Data/2012/Biomarker/HBA1C.dta")) %>%
  mutate(mom_hba1c_2012_observed = 1)

# mom MMSE data 
mmse_corrected_2012 <- as_tibble(read_dta("../Data/2012/x_core and created variables/MMSE and Clock Drawing Test/mmseclock_corrected 20140709.dta")) %>%
  mutate(mmse_corrected_2012_observed = 1)

# BP in 2012
mom_bp_2012 <- as_tibble(read_dta("../Data/2012/x_core and created variables/mainbp2012.dta")) %>%
  mutate(mom_bp_2012_observed = 1)

# mom getupgo data (2012)
mom_getupgo_2012 <- as_tibble(read_dta("../Data/2012/Biomarker/GETUPGO.dta")) %>%
  mutate(mom_getupgo_2012_observed = 1)

# mom gripstrength data (2012)
mom_gripstrength_2012 <- as_tibble(read_dta("../Data/2012/Biomarker/GRIPSTRGTH.dta")) %>%
  mutate(mom_gripstrength_2012_observed = 1)

# mom hearing data (2012)
mom_hearing_2012 <- as_tibble(read_dta("../Data/2012/Biomarker/HEARING.dta")) %>%
  mutate(mom_hearing_2012_observed = 1)

# mom lipid data (2012)
mom_lipid_2012 <- as_tibble(read_dta("../Data/2012/Biomarker/LIPID.dta")) %>%
  mutate(mom_lipid_2012_observed = 1)

# mom iadl (2012)
mom_iadl_2012 <- as_tibble(read_dta("../Data/2012/x_core and created variables/disability/iadl disability/2012 iadl disability vars.dta")) %>%
  mutate(mom_iadl_2012_observed = 1)

# mom adl (2012)
mom_adl_2012 <- as_tibble(read_dta("../Data/2012/x_core and created variables/disability/adl disability/2012 adl disability vars.dta")) %>%
  mutate(mom_adl_2012_observed = 1)

# mom morbidity1 (2012)
mom_morbidity1_2012 <- as_tibble(read_dta("../Data/2012/Survey/MORBIDITY1.dta")) %>%
  mutate(mom_morbidity1_2012_observed = 1)

# mom morbidity2 (2012)
mom_morbidity2_2012 <- as_tibble(read_dta("../Data/2012/Survey/MORBIDITY2.dta")) %>%
  mutate(mom_morbidity2_2012_observed = 1)

# person (2012)
mom_person_2012 <- as_tibble(read_dta("../Data/2012/Survey/PERSON.dta")) %>%
  filter(relmom == 0) %>%
  mutate(mom_person_2012_observed = 1)

# CRP (2012)
mom_crp_2012 <- as_tibble(read_dta('../Data/Mothers aging data for Elijah/AGING GRANT ROUND 1 DATA FILES/biomarker/From NU lab/NIA_CRP.dta')) %>%
  mutate(mom_crp_2012_observed = 1)

#hh wealth pcs
wealth_pc_1983 <- as_tibble(read_dta("../Data/cebupc.dta")) %>% filter(year==1983) %>% select(-pc2, -pc3, -site)
wealth_pc_1991 <- as_tibble(read_dta("../Data/cebupc.dta")) %>% filter(year==1991) %>% select(-pc2, -pc3, -site)
wealth_pc_1994 <- as_tibble(read_dta("../Data/cebupc.dta")) %>% filter(year==1994) %>% select(-pc2, -pc3, -site)
wealth_pc_1998 <- as_tibble(read_dta("../Data/cebupc.dta")) %>% filter(year==1998) %>% select(-pc2, -pc3, -site)
wealth_pc_2002 <- as_tibble(read_dta("../Data/cebupc.dta")) %>% filter(year==2002) %>% select(-pc2, -pc3, -site)
wealth_pc_2005 <- as_tibble(read_dta("../Data/cebupc.dta")) %>% filter(year==2005) %>% select(-pc2, -pc3, -site)
