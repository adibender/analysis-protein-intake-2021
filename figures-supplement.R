library(pammtools)
library(ggplot2)
theme_set(theme_bw())
library(patchwork)
library(magrittr)

# preproc
dir.create("results/figures/supplement", recursive = TRUE)

# code
source("plot-helpers.R")

# global variables

surv_max <- 60

## objects
# data frames that contain hazard ratios for 6 comparisons of different nutrition
# protocols for the different models
f6                 <- readRDS("models/m_2_A_6frames.Rds")
f6_discharge       <- readRDS("models/m_2_B_6frames.Rds")
f6_subdistribution <- readRDS("models/m_1_A_6frames.Rds")
# data frame that summarizes the 6 comparisons
comparisons_df <- readRDS("results/comparisons_df.Rds") %>%
  mutate(comparison_text = paste0(comparison, ": ", p1.name, " vs. ", p2.name))

# model objects
mod_2_A <- readRDS("models/m_2_A.Rds") # cause specific: death
mod_2_B <- readRDS("models/m_2_B.Rds") # cause specific: discharge

########################## main model (cause specific ) ########################

# smooth effects baseline, age, bmi
p_4_panel_2_A <- gg_4_panel(mod_2_A)
p_4_panel_2_B <- gg_4_panel(mod_2_B)

ggsave(p_4_panel_2_A, file = "results/figures/supplement/p_4_panel_2_A.png")
ggsave(p_4_panel_2_B, file = "results/figures/supplement/p_4_panel_2_B.png")

# fixed coefficients
p_coef_tab_2_A_B <- coefPlotGAM(
  models = list(mod_2_A, mod_2_B),
  modelnames = c("Risk of death", "Risk of discharge"))

ggsave(p_coef_tab_2_A_B, file="results/figures/supplement/fixed-coefs-main.png")


# # 6-comparisons
# list_p_death <- purrr::map2(
#   .x = seq_along(f6),
#   .y = seq_len(nrow(comparisons_df)),
#   .f = ~single_plot(f6[[.x]], comparisons_df$comparison_text[.y]))
# list_p_death <- purrr::map2(
#   .x = seq_along(f6_discharge),
#   .y = seq_len(nrow(comparisons_df)),
#   .f = ~single_plot(f6_discharge[[.x]], comparisons_df$comparison_text[.y]))

# p6_death <- Reduce("+", list_p_death)
# ggsave(p6_death, file = "results/figures/supplement/main-6panels-death.png")
# p6_discharge <- Reduce("+", list_p_discharge)
# ggsave(p6_discharge, file = "results/figures/supplement/main-6panels-discharge.png")


############################# main model (subdistribution) #####################
# 6-comparisons
list_p_subdistribution <- purrr::map2(
  .x = seq_along(f6_subdistribution),
  .y = seq_len(nrow(comparisons_df)),
  .f = ~single_plot(f6_subdistribution[[.x]], comparisons_df$comparison_text[.y]))
p6_subdistribution <- Reduce("+", list_p_subdistribution)
ggsave(p6_subdistribution,
  file = "results/figures/supplement/main-subdistribution-6panels-death.png",
  width = 12, height = 9, dpi = 600)



## Cumulative Incidence Function
# ==============================================================================

only_low_ped <- prepare_predict_frame(
  patient,
  daily,
  ll_fun   = ll,
  scheme   = list(C2 = 0, C3 = 0),
  surv_max = surv_max,
  type     = 1,
  var1     = "proteinCat2",
  var2     = "proteinCat3")

only_med_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = 1, C3 = 0),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")

low_to_med_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = c(rep(0, 4), rep(1, 7)), C3 = 0),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")

med_to_high_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = c(rep(1, 4), rep(0, 7)), C3 = c(rep(0, 4), rep(1, 7))),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")

only_high_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = 0, C3 = 1),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")



# create list of comparisons
comp1 <- list(only_low_ped, low_to_med_ped)
names(comp1) <- c("low", "low_to_med")
comp2 <- list(low_to_med_ped, only_med_ped)
names(comp2)<- c("low_to_med", "med")
comp3 <-  list(only_low_ped, only_med_ped)
names(comp3)<- c("low", "med")
comp4 <-  list(only_med_ped, med_to_high_ped)
names(comp4)<- c("med", "med_to_high")
comp5 <-  list(med_to_high_ped, only_high_ped)
names(comp5)<- c("med_to_high", "high")
comp6 <-  list(only_med_ped, only_high_ped)
names(comp6)<- c("med", "high")

comparisons_list <- list(comp1, comp2, comp3, comp4, comp5, comp6)
names(comparisons_list) <- purrr::map_chr(comparisons_list, ~paste0(names(.x), collapse = " vs. "))

res_df_cif <- purrr::imap_dfr(
  .x = comparisons_list,
  .f = ~{
    get_cif_df(.x, m_2_A, m_2_B, comparison_name = .y)
  })

res_df_cif <- res_df_cif %>%
  mutate(comparison = factor(comparison,
    levels = c("low vs. low_to_med", "low_to_med vs. med", "low vs. med",
      "med vs. med_to_high", "med_to_high vs. high", "med vs. high"),
    labels = comparisons_df$comparison_text))
res_df_cif <- res_df_cif %>%
  mutate(
    protocol = factor(
      protocol,
      levels = c("low", "low_to_med", "med", "med_to_high", "high"),
      labels = c("exclusively low", "late standard", "early standard", "late high", "early high"))
  )

res_df_cif_sd <- purrr::imap_dfr(
  .x = comparisons_list,
  .f = ~{
    get_cif_df_sd(.x, m_1_A, comparison_name = .y)
  })

res_df_cif_sd <- res_df_cif_sd %>%
  mutate(comparison = factor(comparison,
    levels = c("low vs. low_to_med", "low_to_med vs. med", "low vs. med",
      "med vs. med_to_high", "med_to_high vs. high", "med vs. high"),
    labels = paste0("Comparison ", LETTERS[1:6])))

res_df_cif_sd <- res_df_cif_sd %>%
  mutate(
    protocol = factor(protocol,
      levels = c("low", "low_to_med", "med", "med_to_high", "high"),
      labels = c("exclusively low", "late standard", "early standard", "late high", "early high"))
  )

p_cs <- ggplot(
  data = res_df_cif,
  aes(x = int_mid, y = cif, col = protocol)) +
  geom_line() +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = protocol), alpha = .3) +
  facet_grid(comparison ~ cause) +
  ylab("CIF") + xlab("Days after ICU admission") +
  ggtitle("Cause specific")
p_cs_death <- ggplot(
  data = res_df_cif %>% filter(cause == "death"),
  aes(x = int_mid, y = cif, col = protocol)) +
  geom_line() +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = protocol), alpha = .3) +
  facet_wrap(~comparison) +
  ylab("CIF") + xlab("Days after ICU admission") +
  ggtitle("Cause specific (in-hospital death)") +
  labs(color = "hypothetical diet", fill = "hypothetical diet")

p_cs_discharge <- ggplot(
  data = res_df_cif %>% filter(cause == "discharge"),
  aes(x = int_mid, y = cif, col = protocol)) +
  geom_line() +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = protocol), alpha = .3) +
  facet_wrap(~comparison) +
  ylab("CIF") + xlab("Days after ICU admission") +
  ggtitle("Cause specific (life hospital discharge)") +
  labs(color = "hypothetical diet", fill = "hypothetical diet")

p_sd <- ggplot(
  data = res_df_cif_sd,
  aes(x = int_mid, y = cif, col = protocol)) +
  geom_line() +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = protocol), alpha = .3) +
  facet_grid(comparison ~ cause) +
  ylab("CIF") + xlab("Days after ICU admission") +
  ggtitle("Subdistribution") +
  labs(color = "hypothetical diet", fill = "hypothetical diet")
# only death needed for subdistribution hazards
# (cif is modeled directly in subdistribution models)
p_sd_death <- ggplot(
  data = res_df_cif_sd %>% filter(cause == "death"),
  aes(x = int_mid, y = cif, col = protocol)) +
  geom_line() +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = protocol), alpha = .3) +
  facet_wrap(~comparison) +
  ylab("CIF") + xlab("Days after ICU admission") +
  ggtitle("Subdistribution (in-hospital death)") +
  labs(color = "hypothetical diet", fill = "hypothetical diet")

ggsave(p_cs_death, filename = "results/figures/supplement/p-cs-cif-death.png",
  height = 6, width = 9)
ggsave(p_cs_discharge, filename = "results/figures/supplement/p-cs-cif-discharge.png",
  height = 6, width = 9)
ggsave(p_sd_death, filename = "results/figures/supplement/p-sd-cif-death.png",
  height = 6, width = 9)


########################### sensitivity static #################################
# 6 panel plot static
m_2_A_static <- readRDS("models/static/m_static_A.Rds")
m_2_B_static <- readRDS("models/static/m_static_B.Rds")
ll_static <- readRDS("data/ll-static.Rds")
f6_static <- make_six_frames(m_2_A_static, patient, daily, ll_static,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = 60)
f6_static_discharge <- make_six_frames(m_2_B_static, patient, daily, ll_static,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = 60)

# 6-comparisons
# death
list_p_static <- purrr::map2(
  .x = seq_along(f6_static),
  .y = seq_len(nrow(comparisons_df)),
  .f = ~single_plot(f6_static[[.x]], comparisons_df$comparison_text[.y]))
p6_static <- Reduce("+", list_p_static)
ggsave(p6_static,
  file = "results/figures/supplement/main-static-6panels-death.png",
  width = 12, height = 9, dpi = 600)
# discharge
list_p_static_discharge <- purrr::map2(
  .x = seq_along(f6_static_discharge),
  .y = seq_len(nrow(comparisons_df)),
  .f = ~single_plot(f6_static_discharge[[.x]], comparisons_df$comparison_text[.y]))
p6_static_discharge <- Reduce("+", list_p_static_discharge)
ggsave(p6_static_discharge,
  file = "results/figures/supplement/main-static-6panels-discharge.png",
  width = 12, height = 9, dpi = 600)



####################### sensitivity in-icu death/discharge #####################
m_2_A_icu <- readRDS("models/icu/m_2_A.Rds")
m_2_B_icu <- readRDS("models/icu/m_2_B.Rds")
list_p_scheme <- readRDS("data/list_p_scheme.Rds")

f6_icu <- make_six_frames(m_2_A_icu, patient, daily, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)
f6_icu_discharge <- make_six_frames(m_2_B_icu, patient, daily, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)

list_p_effect <- purrr::map(seq_along(f6_icu), ~single_plot(f6_icu[[.x]], NULL))
list_p_discharge <- purrr::map(
  seq_along(f6_icu_discharge),
  ~single_plot(f6_icu_discharge[[.x]], NULL))

p_1_icu <- list_p_scheme[[1]] + ggtitle("Comparison") +
  list_p_effect[[1]] + ggtitle("In-ICU death") +
  list_p_discharge[[1]] + ggtitle("Live ICU discharge") +
  list_p_scheme[[2]] + list_p_effect[[2]] + list_p_discharge[[2]] +
  list_p_scheme[[3]] + list_p_effect[[3]] + list_p_discharge[[3]] +
  plot_layout(nrow = 3, byrow = TRUE)
p_1_icu
ggsave(p_1_icu, file = "results/figures/supplement/icu-protein_effect_comparisons_1.png")

p_2 <- list_p_scheme[[4]] + ggtitle("Comparison") +
  list_p_effect[[4]] + ggtitle("In-ICU death") +
  list_p_discharge[[4]] + ggtitle("Live ICU discharge") +
  list_p_scheme[[5]] + list_p_effect[[5]] + list_p_discharge[[5]] +
  list_p_scheme[[6]] + list_p_effect[[6]] + list_p_discharge[[6]] +
  plot_layout(nrow = 3, byrow = TRUE)
p_2

ggsave(p_2, file = "results/figures/supplement/icu-protein_effect_comparisons_2.png")


# correlation Protein intake vs. Calorie intake
daily_all <- readRDS("data/daily_all.Rds")
ped_death <- readRDS("data/ped-data-death-hosp-sub5.Rds")
daily_all <- filter(daily_all, CombinedID %in% unique(ped_death$CombinedID))
daily_all <- daily_all %>%
  mutate(
    cal_tot  = EN_Calories + PN_Calories,
    prot_tot = EN_Protein + PN_Protein)

p_corr <- ggplot(daily_all, aes(x = cal_tot, y = prot_tot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Protein") +
  xlab("Calories") +
  ylim(c(0, 500)) + xlim(c(0, 6000))
ggsave(p_corr, file = "results/figures/supplement/cor-protein-calories.png")
cor(daily_all$cal_tot, daily_all$prot_tot, use="pairwise.complete.obs")
