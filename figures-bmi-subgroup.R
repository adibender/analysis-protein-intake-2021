library(patchwork)
library(ggplot2)
theme_set(theme_bw())
library(pammtools)
library(purrr)

# set maximal follow-up
max.follow <- 60
# for preprocessing, set survival time of patients with survival > max.follow to max.time
max.time <- 60.1
# interval break points (for piece-wise exponential model format)
brks <- c(0:floor(max.time))

## for pattern_label function
# define number of days time-dependent covariate was provided
maxdays.nutri <- 11
# define protocols for comparisons
source("scripts/helpers.R")
m_2_A   <- readRDS("models/bmi/m_2_A-bmi-subgroup.Rds")
m_2_B   <- readRDS("models/bmi/m_2_B-bmi-subgroup.Rds")
daily   <- readRDS("data/daily.Rds")
patient <- readRDS("data/patient.Rds")
ll      <- readRDS("data/ll.Rds")
comparisons_df <- readRDS("results/comparisons_df.Rds")

f6 <- make_six_frames(m_2_A, patient, daily, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)
f6_discharge <- make_six_frames(m_2_B, patient, daily, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)

protocol1 <- rep("< 0.8 g/kg", maxdays.nutri)
protocol2 <- c(rep("< 0.8 g/kg", 4), rep("0.8 - 1.2 g/kg", 7))
protocol3 <- rep("0.8 - 1.2 g/kg", maxdays.nutri)
protocol4 <- c(rep("0.8 - 1.2 g/kg", 4), rep("> 1.2 g/kg", 7))
protocol5 <- rep("> 1.2 g/kg", maxdays.nutri)


# 6-comparisons
list_p_death <- purrr::map2(
  .x = seq_along(f6),
  .y = seq_len(nrow(comparisons_df)),
  .f = ~single_plot(f6[[.x]], comparisons_df$comparison_text[.y]))

list_p_discharge <- purrr::map2(
  .x = seq_along(f6_discharge),
  .y = seq_len(nrow(comparisons_df)),
  .f = ~single_plot(f6_discharge[[.x]], comparisons_df$comparison_text[.y]))

p6_death <- Reduce("+", list_p_death)
ggsave(p6_death, file = "results/figures/supplement/bmi-subgroup-6panels-death.png")
p6_discharge <- Reduce("+", list_p_discharge)
ggsave(p6_discharge, file = "results/figures/supplement/bmi-subgroup-6panels-discharge.png")


## number of subjects in subgroup
ped_death_bmi_subgroup <- readRDS("data/ped-data-death-hosp-bmi-subgroup.Rds")
length(unique(ped_death_bmi_subgroup$CombinedID))
