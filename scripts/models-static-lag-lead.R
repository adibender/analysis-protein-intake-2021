library(dplyr)
library(tidyr)
# library(devtools)
#devtools::install_github("https://github.com/adibender/pammtools")
library(pammtools)
library(mgcv)

dir.create("models/static")

source("scripts/helpers.R")

patient <- readRDS("data/patient.Rds")
daily <- readRDS("data/daily.Rds")
ll_fun_static <- function(t, tz) { t >= tz + 4 }
ll_static <- get_laglead(
  0:60,
  tz = 1:11,
  ll_fun = ll_fun_static)
ll_static <- gg_laglead(ll_static) + scale_x_discrete(guide = guide_axis(angle = 90))
ggsave("results/figures/ll-static.png")
saveRDS(ll_fun_static, "data/ll-static.Rds")
# gg_laglead(ll) + scale_x_discrete(guide = guide_axis(angle = 90))
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

ped_static_A <- as_ped(
  data    = patient,
  formula = Surv(event, PatientDied) ~ Year + DiagID2 + AdmCatID + Gender +
    ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut = 0:60,
  id = "CombinedID")
ped_static_A <- ped_static_A %>%
  add_cumulative_eff_vec(daily, "proteinCat2", "proteinCat3", LL = ll_static)
ped_static_A$int_mid <- 0.5 * (ped_static_A$tstart + ped_static_A$tend)
ped_static_A <- ped_static_A[ped_static_A$tend >= 5, ]

model_static_A <- bam(
  formula = formula,
  data    = ped_static_A,
  family  = "poisson",
  offset  = offset)
saveRDS(model_static_A , "models/static/m_static_A.Rds")



### Model B: Discharge
ped_static_B <- as_ped(
  data    = patient,
  formula = Surv(event, PatientDischarged) ~ Year + DiagID2 + AdmCatID +
    Gender + ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_static_B <- ped_static_B %>%
  add_cumulative_eff_vec(daily, "proteinCat2", "proteinCat3", LL = ll_static)
ped_static_B$int_mid <- 0.5 * (ped_static_B$tstart + ped_static_B$tend)
ped_static_B <- ped_static_B[ped_static_B$tend >= 5, ]

model_static_B <- bam(
  formula = formula,
  data    = ped_static_B,
  family  = "poisson",
  offset  = offset)
saveRDS(model_static_B, "models/static/m_static_B.Rds")
