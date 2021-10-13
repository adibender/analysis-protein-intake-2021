library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())
library(pammtools)
library(mgcv)

# create folder
dir.create("models/bmi")

source("scripts/helpers.R")

patient <- readRDS("data/patient.Rds")
daily <- readRDS("data/daily.Rds")

ll_fun <- function(t, tz) { t >= tz + 4 & t <= tz * 3 + 10 }
ll <- get_laglead(
  0:60,
  tz = 1:11,
  ll_fun = ll_fun)
ll_dynamic <- gg_laglead(ll) + scale_x_discrete(guide = guide_axis(angle = 90))

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


############## Model 2 ###############
### Cause spec hazard 60 days: Admin censoring after 60 days, Protein
### This model sets ALL survivors (more prec. dischargees) to the maximal event
### time of 60.
### Model A: Death
ped_2_A <- as_ped(
  data    = filter(patient, BMI >30),
  formula = Surv(event, PatientDied) ~ Year + DiagID2 + AdmCatID + Gender +
    ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_A <- ped_2_A %>%
  add_cumulative_eff_vec(daily, "proteinCat2", "proteinCat3", LL = ll)
ped_2_A$int_mid <- 0.5 * (ped_2_A$tstart + ped_2_A$tend)
ped_2_A <- ped_2_A[ped_2_A$tend >= 5, ]
saveRDS(ped_2_A, "data/ped-data-death-hosp-bmi-subgroup.Rds")

model_2_A <- bam(
  formula = formula,
  data    = ped_2_A,
  family  = "poisson",
  offset  = offset)
saveRDS(model_2_A, "models/bmi/m_2_A-bmi-subgroup.Rds")

### Model B: Discharge
ped_2_B <- as_ped(
  data    = filter(patient, BMI > 30),
  formula = Surv(event, PatientDischarged) ~ Year + DiagID2 + AdmCatID +
    Gender + ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_B <- ped_2_B %>%
  add_cumulative_eff_vec(daily, "proteinCat2", "proteinCat3", LL = ll)
ped_2_B$int_mid <- 0.5 * (ped_2_B$tstart + ped_2_B$tend)
ped_2_B <- ped_2_B[ped_2_B$tend >= 5, ]
saveRDS(ped_2_B, "data/ped-data-discharge-hosp-bmi-subgroup.Rds")

model_2_B <- bam(
  formula = formula,
  data    = ped_2_B,
  family  = "poisson",
  offset  = offset)
saveRDS(model_2_B, "models/bmi/m_2_B-bmi-subgroup.Rds")
