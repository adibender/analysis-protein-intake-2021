library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())
library(pammtools)
library(mgcv)
source("scripts/helpers.R")

patient <- readRDS("data/patient.Rds")
daily <- readRDS("data/daily.Rds")

dir.create("models/icu")

ll_fun <- function(t, tz) { t >= tz + 4 & t <= tz * 3 + 10 }
ll <- get_laglead(
  0:60,
  tz = 1:11,
  ll_fun = ll_fun)

####################### Model I ###############
formula = as.formula('ped_status ~ s(int_mid, bs = "ps") +
  Year +
  ApacheIIScore + ApacheIIScore:int_mid +
  s(Age, by = int_mid, bs = "ps") +
  s(BMI, by = int_mid, bs = "ps") + DiagID2 + AdmCatID + Gender +
  inMV2_4 + Propofol2_4 + OralIntake2_4 + PN2_4 +
  s(CombinedicuID, bs = "re", by = icuByDummy) +
  te(t, tz, by = I(LL * proteinCat3),
  bs = "ps", m = list(c(2, 1), c(2, 1)), id = "cal") +
  te(t, tz, by = I(LL * proteinCat2),
  bs = "ps", m = list(c(2, 1), c(2, 1)), id = "cal")
  ')

############## Model 1 ###############
### 60 days
### Subdist hazard 60 days: Admin censoring after 60 days, Protein
### This model sets ALL survivors (more prec. dischargees) to the maximal event
### time of 60.
### Model A: Death
patient_1_A <- patient
patient_1_A$time <- patient_1_A$surv_icu0to60
patient_1_A$time[patient_1_A$surv_icu_status == 1] <- 61
patient_1_A$status <- 1L * (patient_1_A$surv_icu_status == 2)
ped_1_A <- as_ped(
  data    = patient_1_A,
  formula = Surv(time, status) ~ Year + DiagID2 + AdmCatID + Gender +
    ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_1_A <- ped_1_A %>%
  add_cumulative_eff_vec(daily, "proteinCat2", "proteinCat3", LL = ll)
ped_1_A$int_mid <- 0.5 * (ped_1_A$tstart + ped_1_A$tend)
ped_1_A <- ped_1_A[ped_1_A$tend >= 5, ]
model_1_A <- bam(
  formula = formula,
  data = ped_1_A,
  family = poisson(),
  offset = offset)
saveRDS(model_1_A, "models/icu/m_1_A.Rds")

############## Model 2 ###############
### Cause spec hazard 60 days: Admin censoring after 60 days, Protein
### time of 60.
### Model A: Death
patient_2_A <- patient
patient_2_A$time <- as.numeric(patient_2_A$surv_icu0to60)
patient_2_A$status <- 1L * (patient_2_A$surv_icu_status == 2)
ped_2_A <- as_ped(
  data    = patient_2_A,
  formula = Surv(time, status) ~ Year + DiagID2 + AdmCatID + Gender +
    ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_A <- ped_2_A %>%
  add_cumulative_eff_vec(daily, "proteinCat2", "proteinCat3", LL = ll)
ped_2_A$int_mid <- 0.5 * (ped_2_A$tstart + ped_2_A$tend)
saveRDS(ped_2_A, "data/ped-data-death.Rds")
# one of the exclusion criteria was that subjects had to be on ICU for at least
# 4 (calendar) days and survived more than 4 days. Thus we also subset the ped data
ped_2_A <- ped_2_A[ped_2_A$tend >= 5, ]
saveRDS(ped_2_A, "data/ped-data-death-sub5.Rds")

model_2_A <- bam(
  formula = formula,
  data    = ped_2_A,
  family  = "poisson",
  offset  = offset)
saveRDS(model_2_A, "models/icu/m_2_A.Rds")

### Model B: Discharge
patient_2_B <- patient
patient_2_B$time <- patient_2_B$surv_icu0to60
patient_2_B$status <- 1L * (patient_2_B$surv_icu_status == 1)
ped_2_B <- as_ped(
  data    = patient_2_B,
  formula = Surv(time, status) ~ Year + DiagID2 + AdmCatID +
    Gender + ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_B <- ped_2_B %>%
  add_cumulative_eff_vec(daily, "proteinCat2", "proteinCat3", LL = ll)
ped_2_B$int_mid <- 0.5 * (ped_2_B$tstart + ped_2_B$tend)
ped_2_B <- ped_2_B[ped_2_B$tend >= 5, ]
saveRDS(ped_2_B, "data/ped-data-discharge-sub5.Rds")

model_2_B <- bam(
  formula = formula,
  data    = ped_2_B,
  family  = "poisson",
  offset  = offset)

saveRDS(model_2_B, "models/icu/m_2_B.Rds")
