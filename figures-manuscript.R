library(patchwork)
library(ggplot2)
theme_set(theme_bw())
library(pammtools)
library(purrr)

dir.create("results/figures/", recursive = TRUE)

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
m_1_A   <- readRDS("models/m_1_A.Rds")
m_2_A   <- readRDS("models/m_2_A.Rds")
m_2_B   <- readRDS("models/m_2_B.Rds")
daily   <- readRDS("data/daily.Rds")
patient <- readRDS("data/patient.Rds")
ll      <- readRDS("data/ll.Rds")

f6 <- make_six_frames(m_2_A, patient, daily, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)
f6_discharge <- make_six_frames(m_2_B, patient, daily, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)
f6_subdistribution <- make_six_frames(m_1_A, patient, daily, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)


protocol1 <- rep("< 0.8 g/kg", maxdays.nutri)
protocol2 <- c(rep("< 0.8 g/kg", 4), rep("0.8 - 1.2 g/kg", 7))
protocol3 <- rep("0.8 - 1.2 g/kg", maxdays.nutri)
protocol4 <- c(rep("0.8 - 1.2 g/kg", 4), rep("> 1.2 g/kg", 7))
protocol5 <- rep("> 1.2 g/kg", maxdays.nutri)

list_p_effect <- purrr::map(seq_along(f6), ~single_plot(f6[[.x]], NULL))
list_p_discharge <- purrr::map(
  seq_along(f6_discharge),
  ~single_plot(f6_discharge[[.x]], NULL))

comparisons_df <-
  data.frame(
    comparison = LETTERS[1:6],
    p1 = c("protocol1", "protocol2", "protocol1", "protocol3", "protocol4", "protocol3"),
    p2 = c("protocol2", "protocol3", "protocol3", "protocol4", "protocol5", "protocol5"),
    p1.name = c("exclusively low", "late standard", "exclusively low", "early standard", "late high", "early standard"),
    p2.name = c("late standard", "early standard", "early standard", "late high", "early high", "early high"),
    p1.count = c(1, 2, 1, 3, 4, 3),
    p2.count = c(2, 3, 3, 4, 5, 5)
  )
saveRDS(comparisons_df, "results/comparisons_df.Rds")
names(f6) <- names(f6_discharge) <- names(f6_subdistribution) <-
  paste0(comparisons_df$p1.name, " vs. ", comparisons_df$p2.name)
saveRDS(f6, "models/m_2_A_6frames.Rds")
saveRDS(f6_discharge, "models/m_2_B_6frames.Rds")
saveRDS(f6_subdistribution, "models/m_1_A_6frames.Rds")

list_p_scheme <- purrr::map(
  .x = seq_len(nrow(comparisons_df)),
  .f = ~ggcomparison_scheme(
    protocol1 = get(comparisons_df[.x,"p1"]),
    protocol2 = get(comparisons_df[.x, "p2"]),
    protocol1.name = comparisons_df[.x, "p1.name"],
    protocol2.name = comparisons_df[.x, "p2.name"],
    legend = TRUE
  )
)
saveRDS(list_p_scheme, "data/list_p_scheme.Rds")

p_1 <- list_p_scheme[[1]] + ggtitle("Comparison") +
  list_p_effect[[1]] + ggtitle("Death") +
  list_p_discharge[[1]] + ggtitle("Discharge") +
  list_p_scheme[[2]] + list_p_effect[[2]] + list_p_discharge[[2]] +
  list_p_scheme[[3]] + list_p_effect[[3]] + list_p_discharge[[3]] +
  plot_layout(nrow = 3, byrow = TRUE)
p_1
ggsave(p_1, file = "results/figures/protein_effect_comparisons_1.png", height = 9, width = 9)

p_2 <- list_p_scheme[[4]] + ggtitle("Comparison") +
  list_p_effect[[4]] + ggtitle("Death") +
  list_p_discharge[[4]] + ggtitle("Discharge") +
  list_p_scheme[[5]] + list_p_effect[[5]] + list_p_discharge[[5]] +
  list_p_scheme[[6]] + list_p_effect[[6]] + list_p_discharge[[6]] +
  plot_layout(nrow = 3, byrow = TRUE)
p_2

ggsave(p_2, file = "results/figures/protein_effect_comparisons_2.png", height = 9, width = 9)
