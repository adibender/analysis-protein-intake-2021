# Improvised function adding cumulative effects to the PED.
add_cumulative_eff_vec <- function(
  ped,
  daily_set,
  varname1,
  varname2,
  LL) {

  ped           <- ped[ped$CombinedID %in% daily_set$CombinedID, ]
  daily_set     <- daily_set[daily_set$CombinedID %in% ped$CombinedID, ]
  grouped_ped   <- ped %>% group_by(CombinedID)
  nested_ped    <- grouped_ped %>% tidyr::nest()
  grouped_daily <- daily_set %>% group_by(CombinedID)
  nested_daily  <- grouped_daily %>% tidyr::nest()
  t_matrix      <- matrix(rep(ped$tend), ncol = 11, nrow = nrow(ped)) - 0.5
  tz_matrix     <- matrix(1:11, ncol = 11, nrow = nrow(ped), byrow = TRUE) + 1
  var1_matrix   <- matrix(1, nrow = nrow(ped), ncol = 11)
  var2_matrix   <- matrix(1, nrow = nrow(ped), ncol = 11)
  LL_matrix     <- matrix(0, nrow = nrow(ped), ncol = 11)
  max_LL_matrix <- matrix(c(LL$LL), nrow = length(unique(LL$t)), byrow = TRUE)
  max_LL_matrix <- max_LL_matrix[2:(max(ped$tend) + 1), ]
  ids <- unique(ped$CombinedID)
  for (i in 1:length(ids)) {
    number <- which(attr(grouped_ped, "groups")$CombinedID == ids[i])
    cur_rows <- attr(grouped_ped, "groups")$.rows[[number]]
    matching <- which(nested_daily$CombinedID == ids[i])
    cur_data <- nested_daily$data[[matching]]
    var1_vec <- rep(NA, 11)
    var2_vec <- rep(NA, 11)
    var1_vec[1:nrow(cur_data)] <- dplyr::pull(cur_data, varname1)
    var2_vec[1:nrow(cur_data)] <- dplyr::pull(cur_data, varname2)
    var1_matrix[cur_rows, ] <- matrix(rep(var1_vec, length(cur_rows)),
                                      nrow = length(cur_rows), byrow = TRUE)
    var2_matrix[cur_rows, ] <- matrix(rep(var2_vec, length(cur_rows)),
                                      nrow = length(cur_rows), byrow = TRUE)
    LL_matrix[cur_rows, ] <- max_LL_matrix[1:length(cur_rows), ]
  }
  t_matrix[is.na(var1_matrix)] <- 0
  tz_matrix[is.na(var1_matrix)] <- 0
  var1_matrix[is.na(var1_matrix)] <- 0
  var2_matrix[is.na(var2_matrix)] <- 0
  ped$t <- t_matrix
  ped$tz <- tz_matrix
  ped$var1 <- var1_matrix
  ped$var2 <- var2_matrix
  ped$LL <- LL_matrix
  colnames(ped)[colnames(ped) == "var1"] <- varname1
  colnames(ped)[colnames(ped) == "var2"] <- varname2
  ped
}

# Compute effect differences as in Bender et al. 2018
# Taken from ELRA paper
getEffectDiffs <- function(
  m,
  newdata1,
  newdata2,
  effectname = "calCat") {

  ## deal with fits created with old mgcv before ti() was implemented
  ## when loaded with newer versions that have ti(): add the missing "mc"-slot
  if(!is.null(mgcv::ti)){
    m$smooth <- lapply(m$smooth, function(trm){
      if("margin" %in% names(trm) & is.null(trm$mc)){
        trm$mc <- rep(TRUE, length(trm$margin))
        #message("added mc slot to ", trm$label)
      }
      trm
    } )
  }

  useColumns <- grep(effectname, names(m$coefficients))
  covCoefs <- m$Vp[useColumns, useColumns]
  X1 <- predict(m, newdata = newdata1, type = "lpmatrix")[,useColumns]
  X2 <- predict(m, newdata = newdata2, type = "lpmatrix")[,useColumns]
  X <- X2 - X1
  fit <- drop(X %*% m$coefficients[useColumns])
  se <- sqrt(rowSums((X %*% covCoefs) * X))
  hi <- fit + 2 * se
  lo <- fit - 2 * se
  cbind(intmid = newdata1$int_mid, fit = exp(fit), se = se, hi = exp(hi),
        lo = exp(lo))

}

# Prepare data frame for prediction
#
# Arguments:
# p: the patient data set (one row per patient)
# d: the daily data set (one row per day in ICU)
# scheme: the nutrition protocol
# surv_max: the maximally considered follow up
prepare_predict_frame <- function(
  p,
  d,
  ll_fun,
  scheme = list(C2 = NULL, C3 = NULL),
  surv_max = 60,
  type,
  var1,
  var2) {
  p <- p[1, ]
  p$CombinedID = 1
  p$PatientDied <- 1
  p$PatientDischarged <- 1
  p$event <- surv_max
  d <- d[1:11, ]
  d$CombinedID <- 1
  d$Study_Day <- 1:11
  d$OralIntake <- 0
  if (var1 %in% c("calCat2", "calCat3")) {
    d$calCat2 <- scheme[[1]]
    d$calCat3 <- scheme[[2]]
  } else {
    d$proteinCat2 <- scheme[[1]]
    d$proteinCat3 <- scheme[[2]]
  }
  ll <- get_laglead(0:surv_max,
                    tz = 1:11,
                    ll_fun = ll_fun)
  if (type == "1") {
    final_ped <- as_ped(
      data    = p,
      formula = Surv(event, PatientDied) ~ Year + DiagID2 + AdmCatID + Gender +
        ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
          CombinedicuID + icuByDummy,
      cut     = 0:surv_max, id = "CombinedID") %>%
      add_cumulative_eff_vec(d, var1, var2, LL = ll)
  } else {
    final_ped <- as_ped(data    = p,
                        formula = Surv(event, PatientDischarged) ~ Year + DiagID2 + AdmCatID + Gender +
                          ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
                          CombinedicuID + icuByDummy,
                        cut     = 0:surv_max, id = "CombinedID") %>%
      add_cumulative_eff_vec(daily_set = d, var1, var2, LL = ll)
  }
  final_ped$icuByDummy <- 0
  final_ped$int_mid <- 0.5 * (final_ped$tstart + final_ped$tend)
  final_ped$CombinedicuID <- 11325L
  final_ped <- final_ped[final_ped$tend > 4, ]
  final_ped
}

# Make a single comparison plot
single_plot <- function(df, number) {
  p <- ggplot(df, aes(x = intmid, y = fit)) +
    geom_hline(yintercept = 1, col = "red") +
    geom_step(aes(x = intmid, y = fit)) +
    geom_step(aes(x = intmid, y = lo), linetype = "dotted") +
    geom_step(aes(x = intmid, y = hi), linetype = "dotted") +
    ylab("Hazard ratio") + xlab("Days after ICU admission") +
    ggtitle(paste(number)) +
    scale_y_continuous(limits = c(0.15, 4), trans = "log1p",
                       breaks = c(0.15, 0.5, 1, 1.5, 2, 4))
  p
}

#
six_plots <- function(df1, df2, df3, df4, df5, df6) {
  p1 <- single_plot(df1, "A")
  p2 <- single_plot(df2, "B")
  p3 <- single_plot(df3, "C")
  p4 <- single_plot(df4, "D")
  p5 <- single_plot(df5, "E")
  p6 <- single_plot(df6, "F")
  grid.arrange(
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    nrow = 2,
    top = "Comparison of different diets"
  )
}

make_six_frames <- function(
  m, patient, daily, ll_fun, type,
  var1 = "calCat2", var2 = "calCat3",
  effect = "calcat", surv_max = 60, min_surv = 5) {
  only_high_ped <- prepare_predict_frame(
    patient, daily, ll_fun,
    scheme = list(C2 = 0, C3 = 1),
    surv_max = surv_max, type, var1, var2)

  only_med_ped <- prepare_predict_frame(
    patient, daily, ll_fun,
    scheme = list(C2 = 1, C3 = 0),
    surv_max = surv_max, type, var1, var2)

  only_low_ped <- prepare_predict_frame(
    patient, daily, ll_fun,
    scheme = list(C2 = 0, C3 = 0),
    surv_max = surv_max, type, var1, var2)

  low_to_med_ped <- prepare_predict_frame(
    patient, daily, ll_fun,
    scheme = list(C2 = c(rep(0, 4), rep(1, 7)), C3 = 0),
    surv_max = surv_max, type, var1, var2)

  med_to_high_ped <- prepare_predict_frame(
    patient, daily, ll_fun,
    scheme = list(C2 = c(rep(1, 4), rep(0, 7)), C3 = c(rep(0, 4), rep(1, 7))),
    surv_max = surv_max, type, var1, var2)

  comp1 <- as.data.frame(
    getEffectDiffs(m, only_low_ped, low_to_med_ped, effectname = effect))
  comp2 <- as.data.frame(
    getEffectDiffs(m, low_to_med_ped, only_med_ped, effectname = effect))
  comp3 <- as.data.frame(
    getEffectDiffs(m, only_low_ped, only_med_ped, effectname = effect))
  comp4 <-
    as.data.frame(
      getEffectDiffs(m, only_med_ped, med_to_high_ped, effectname = effect))
  comp5 <- as.data.frame(
    getEffectDiffs(m, med_to_high_ped, only_high_ped, effectname = effect))
  comp6 <- as.data.frame(
    getEffectDiffs(m, only_med_ped, only_high_ped, effectname = effect))
  res <- list(comp1, comp2, comp3, comp4, comp5, comp6)
  attr(res, "exampledf") <- only_high_ped
  res
}

three_plots <- function(df1, df2, df3, names = c("A", "B", "C"),
                        main = "Comparison of different diets") {
  p1 <- single_plot(df1, names[1])
  p2 <- single_plot(df2, names[2])
  p3 <- single_plot(df3, names[3])
  grid.arrange(
    p1,
    p2,
    p3,
    nrow = 1,
    top = main
  )
}


ggcomparison_scheme <- function(
    protocol1,
    protocol2,
    protocol1.name,
    protocol2.name,
    legend=TRUE,
    coding=data.frame(
        value=c("< 0.8 g/kg", "0.8 - 1.2 g/kg", "> 1.2 g/kg"),
        code=c(.49, .99, 1.41))) {
    library(dplyr)

    df.scheme  <- cbind.data.frame(
        p1 = as.character(protocol1),
        p2 = as.character(protocol2),
        te = seq_along(protocol1))
    mdf.scheme <- tidyr::gather(df.scheme, variable, value, -te)
    mdf.scheme <- left_join(mdf.scheme, coding)
    mdf.scheme$variable <- factor(
        mdf.scheme$variable,
        levels=unique(mdf.scheme$variable),
        labels=paste0(c(protocol1.name, protocol2.name)))

   gg.scheme <- ggplot(mdf.scheme, aes(x=te, y=code, fill=variable)) +
        geom_bar(stat="identity", position="dodge") +
        # geom_hline(yintercept=3, col="grey12") +
        ylab("g protein/kg day") +
        xlab("Nutrition day") +
        scale_fill_manual(name=NULL, values=c("grey70", "grey30")) +
        scale_y_continuous(breaks=seq(0, 2, by = .2), limits = c(0, 2)) +
        scale_x_continuous(breaks=1:11)
    if(legend) {
        gg.scheme + theme(
          legend.position=c(0, 1),
          legend.justification=c(0, 1),
          legend.direction="horizontal",
          legend.background=element_rect(colour="black"),
          legend.key.size=unit(0.3, "cm"))
    } else {
        gg.scheme + theme(legend.position="none")
    }

}
