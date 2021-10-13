library(dplyr)
library(purrr)

ped_death <- readRDS("data/ped-data-death-hosp-sub5.Rds")
ped_discharge <- readRDS("data/ped-data-discharge-hosp-sub5.Rds")

dir.create("results/main")

# Table of patient at risk/dying/discharged per interval
smry_death <- ped_death %>%
  group_by(interval) %>%
  summarize(at_risk = n(), dying = sum(ped_status))
smry_discharge <- ped_discharge %>%
  group_by(interval) %>%
  summarize(at_risk = n(), discharge = sum(ped_status)) %>%
  select(-at_risk)

smry_tab <- smry_death %>% left_join(smry_discharge, by = c("interval"))

readr::write_excel_csv(smry_tab, "results/supplement-table-death-discharge-per-interval.csv")


# HR per interval (death and discharge) for all comparisons
frames2A <- readRDS("models/m_2_A_6frames.Rds")
frames2B <- readRDS("models/m_2_B_6frames.Rds")

frames2A <- map(
  .x = frames2A,
  .f = ~{
    .x %>%
    mutate_at(
      c("fit", "lo", "hi"),
      ~round(.x, 2)) %>%
    mutate(CI = paste0("[", lo, ", ", hi, ",]"))
  }
)

frames2B <- map(
  .x = frames2B,
  .f = ~{
    .x %>%
    mutate_at(
      c("fit", "lo", "hi"),
      ~round(.x, 2)) %>%
    mutate(CI = paste0("[", lo, ", ", hi, ",]"))
  }
)

iwalk(frames2A, ~ readr::write_excel_csv(.x, paste0("results/main/tab-death-", .y, ".csv")))
iwalk(frames2B, ~ readr::write_excel_csv(.x, paste0("results/main/tab-discharge-", .y, ".csv")))
