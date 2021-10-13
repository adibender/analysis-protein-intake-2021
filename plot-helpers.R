## 4er panel
# - baseline hazard, age, BMI, ICU Random effects
gg_4_panel <- function(mod) {
  smooths <- tidy_smooth(mod) %>%
    dplyr::filter(xlab %in% c("int_mid", "Age", "BMI"))
  p_re <- gg_re(mod)

  df_baseline <- dplyr::filter(smooths, xlab == "int_mid")
  p_baseline <- ggplot(
      dplyr::filter(smooths, xlab == "int_mid"),
      aes(x = x, y = fit)) +
    geom_line() +
    geom_line(aes(y = ci_lower), lty = 2) +
    geom_line(aes(y = ci_upper), lty = 2) + xlab("t") + ylab(df_baseline$ylab[1])
  df_age <- dplyr::filter(smooths, xlab == "Age")
  p_age <- ggplot(
      dplyr::filter(smooths, xlab == "Age"),
      aes(x = x, y = fit)) +
    geom_line() +
    geom_line(aes(y = ci_lower), lty = 2) +
    geom_line(aes(y = ci_upper), lty = 2) + xlab("Age") + ylab(df_age$ylab[1])
  df_bmi <- dplyr::filter(smooths, xlab == "BMI")
  p_bmi <- ggplot(
      dplyr::filter(smooths, xlab == "BMI"),
      aes(x = x, y = fit)) +
    geom_line() +
    geom_line(aes(y = ci_lower), lty = 2) +
    geom_line(aes(y = ci_upper), lty = 2) + xlab("BMI") + ylab(df_bmi$ylab[1])

    p_baseline + p_age + p_bmi + p_re

}


## adapted for mgcv gams from https://gist.github.com/dsparks/818976

coefPlotGAM <- function(
  models = list(),
  data,
  alpha       = 0.05,
  modelnames  = NULL,
  col0        = "black",
  intercept   = FALSE,
  nrow.legend = 2,
  legend.position = "bottom",
  mod.cols    = c("#E41A1C", "#377EB8"),
  exclude     = NULL,
  prettify    = TRUE) {
    # models must be a list()

  library(ggplot2)

  models <- lapply(models, function(z) {
    if( is.character(z) ) {
      summary(readRDS(z))
    }
    else {
      if( !(class(z)[1] == "summary.gam") ) {
        summary(z)
      }
      else z
    }
  })

  Multiplier <- qnorm(1 - alpha / 2)

    ## prep coefficient table, exclude not needed vars
  CoefficientTables <- lapply(models, "[[", i = "p.table")
  if( !intercept )
    CoefficientTables <- lapply(CoefficientTables, function(z) z[-1, ])
  if( !is.null(exclude) ) {
    CoefficientTables <- lapply(CoefficientTables, function(z) {
      z[-grep(exclude, row.names(z)), ]
    })
  }

  for(i in seq_along(CoefficientTables)) {
    CoefficientTables[[i]] <- cbind.data.frame(CoefficientTables[[i]],
      vars=rownames(CoefficientTables[[i]]))
    if( prettify ) {
      df.i <- prettify_labels(CoefficientTables[[i]])
      CoefficientTables[[i]] <- merge(CoefficientTables[[i]], df.i, by="vars")
    } else{
      colnames(CoefficientTables[[i]])["vars"] <- "IV"
    }
    CoefficientTables[[i]]$Model <- factor(modelnames[i], rev(modelnames))
  }

  matrix.of.models <- do.call(rbind, CoefficientTables)
  matrix.of.models <- matrix.of.models[, c("Estimate", "Std. Error", "IV", "Model")]
  colnames(matrix.of.models) <- c("Estimate", "StandardError", "IV", "Model")
  matrix.of.models[["Multiplier"]] <- Multiplier
    # matrix.of.models$IV <- factor(matrix.of.models$IV, levels = matrix.of.models$IV)
    # matrix.of.models[, -c(1, 6)] <- apply(matrix.of.models[, -c(1, 6)], 2,
    #         function(x){as.numeric(as.character(x))})

    # for(i in levels(matrix.of.models$IV))
    #   matrix.of.models$IV <- relevel(matrix.of.models$IV, i)
  p <- ggplot(data = matrix.of.models) +
  geom_hline(yintercept = 0, alpha = 0.8, lty = 1, col = col0 ) +
  theme_bw() +
  xlab("") +
  coord_flip()

  if(length(models) > 1) {
    p + geom_pointrange(aes(y = Estimate, x = IV,
      color = Model,
      ymax = Estimate + Multiplier * StandardError,
      ymin = Estimate - Multiplier * StandardError),
    position=position_dodge(width=.6)) +
    scale_color_manual("Model", values = rev(mod.cols)) +
    guides(color = guide_legend(nrow = nrow.legend, reverse=TRUE)) +
    theme(legend.position=legend.position)


  } else {
    p + geom_pointrange(aes(
      y = Estimate, x = IV,
      ymax = Estimate + Multiplier * StandardError,
      ymin = Estimate - Multiplier * StandardError)) +
    theme(legend.position="none") +
    scale_color_manual(values=mod.col)
  }


}


prettify_labels <- function(coeftab) {

  labels <- c(
    "Intercept",
    "Year of therapy: 2008",
    "Year of therapy: 2009",
    "Year of therapy: 2011",
    "Year of therapy: 2013",
    "Year of therapy: 2014",
    "Apache II Score",
    "Admission diagnosis: Cardio-Vascular",
    "Admission diagnosis: Respiratory",
    "Admission diagnosis: Gastrointestinal",
    "Admission diagnosis: Neurologic",
    "Admission diagnosis: Sepsis",
    "Admission diagnosis: Orthopedic/Trauma",
    "Admission diagnosis: Metabolic",
    "Admission diagnosis: Renal",
    "Admission: Surgical/Elective",
    # "Admission category: Medical",
    "Admission category: Surgical/Emergency",
    "Gender: Male",
    "#days with MV\n (up to day 4)",
    "#days with Propofol\n (up to day 4)",
    "#days with Oral Intake\n (up to day 4)",
    "#days Parenteral Nutrition\n (up to day 4)",
    "Apache II Score:t")

  vars <- c(
    "(Intercept)",
    "Year2008",
    "Year2009",
    "Year2011",
    "Year2013",
    "Year2014",
    "ApacheIIScore",
    "DiagID2Cardio-Vascular",
    "DiagID2Respiratory",
    "DiagID2Gastrointestinal",
    "DiagID2Neurologic",
    "DiagID2Sepsis",
    "DiagID2Orthopedic/Trauma",
    "DiagID2Metabolic",
    "DiagID2Renal",
    "AdmCatIDSurgical/Elective",
    # "AdmCatIDMedical",
    "AdmCatIDSurgical/Emeregency",
    "GenderMale",
    "inMV2_4",
    "Propofol2_4",
    "OralIntake2_4",
    "PN2_4",
    "ApacheIIScore:int_mid")

  df.labs <- cbind.data.frame(vars = as.factor(vars), IV = as.factor(labels))
  for( i in labels ) {
    df.labs$IV <- relevel(df.labs$IV, ref = i)
  }

  df.labs[match(rownames(coeftab), df.labs$vars), ]

}

# Differences between nutriton protocols
getEffectDiffs <- function(
  m,
  newdata1,
  newdata2,
  effectname = "ProteinCat") {

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
  cbind(int_mid = newdata1$int_mid, fit = fit, se = se, hi = hi, lo = lo)

}

ggcomparison_scheme <- function(
    protocol1,
    protocol2,
    protocol1.count,
    protocol2.count,
    legend=TRUE,
    coding=data.frame(
        value=c("low", "mid", "full"),
        code=c(1,2,3))) {

    library(reshape2)
    library(dplyr)

    df.scheme  <- cbind.data.frame(
        p1 = as.character(protocol1),
        p2 = as.character(protocol2),
        te = seq_along(protocol1))
    mdf.scheme <- melt(df.scheme, id.vars="te")
    mdf.scheme <- left_join(mdf.scheme, coding)
    mdf.scheme$variable <- factor(
        mdf.scheme$variable,
        levels=unique(mdf.scheme$variable),
        labels=paste0("plan #", c(protocol1.count, protocol2.count)))

   gg.scheme <- ggplot(mdf.scheme, aes(x=te, y=code, fill=variable)) +
        geom_bar(stat="identity", position="dodge") +
        # geom_hline(yintercept=3, col="grey12") +
        ylab("Category of caloric intake") +
        xlab("Nutrition day") +
        scale_fill_manual(name=NULL, values=c("grey70", "grey30")) +
        scale_y_continuous(breaks=1:3, labels=c("I", "II", "III"),
            limits=c(0, 3.5)) +
        scale_x_continuous(breaks=1:11)
    if(legend) {
        gg.scheme + theme(
          legend.position=c(0, 1.05),
          legend.justification=c(0, 1),
          legend.direction="horizontal",
          legend.background=element_rect(colour="black"),
          legend.key.size=unit(0.3, "cm"))
    } else {
        gg.scheme + theme(legend.position="none")
    }

}

#' Comparison plots
#'
#' Given 2 nutrition protocols, creates comparison plot.
#'
#' @rdname ggcomparison
#' @export
ggcomparison <- function(
  comp.data,
  comp.effect=TRUE,
  plot.legend=TRUE,
  y.lim=c(-2.25, 2.25),
  brks=0:30) {

    gg.diff <- ggplot(comp.data,
      aes(x = int_mid, y = fit)) +
    geom_hline(yintercept = 0, lty=1, col="grey50") +
    geom_line(aes(linetype=Model)) +
    geom_line(data=filter(comp.data, Model==levels(comp.data$Model)[1]),
        aes(y = lo), lty=2) +
    geom_line(data=filter(comp.data, Model==levels(comp.data$Model)[1]),
        aes(y = hi), lty=2) +
    ylab("Hazard ratio") +
    xlab("Days after ICU admission") +
    scale_linetype_manual(values=c(1, 3))+
    scale_fill_manual(values=c(NA, ribbon.col)) +
    scale_x_continuous(minor_breaks = c(4, int_info(brks=brks)$tend)) +
    guides(color = guide_legend(order = 0)) +
    geom_text(data=data.frame(x=45, y=1, label=unique(comp.data$comparison)),
        aes(x=x, y=y, label=label)) +
    scale_y_continuous(
        breaks = log(sy <- c(0.25, 0.5, 0.75, 1, 1.25, 2, 4)),
        labels = sy,
        limits = c(-2, 2))

    if(comp.effect) {
        gg.diff <- gg.diff +
            geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model)) +
             geom_line(data=subset(comp.data, Model==levels(comp.data$Model)[2]),
                colour="grey20", lty=3)
    }
    if(plot.legend) {
        gg.diff + theme(legend.position="top")
    }
    else {
        gg.diff + theme(legend.position="none")
    }

}

#' @rdname ggcomparison
#' @export
ggcomparison_single <- function(
    comp.data,
    comp.effect=NULL,
    plot.legend=TRUE,
    y.lim=c(-1.5, 1.5),
    legend.labels=c("Hazard ratio", "95% CI"),
    brks=0:30)  {

    comp.data <- filter(comp.data, Model==levels(Model)[1])
    comp.data <- melt(comp.data[, c("fit", "int_mid", "hi", "lo", "label",
        "comparison")], id.vars=c("int_mid", "label", "comparison"))
    comp.data$type <- factor(comp.data$variable=="fit", levels=c(TRUE, FALSE),
        labels=legend.labels)

    gg.diff <- ggplot(comp.data, aes(x = int_mid, y = value)) +
        geom_hline(yintercept = 0, lty=1, col="grey50") +
        # geom_step(aes(group=variable, lty=type)) +
        geom_line(aes(group=variable, lty=type)) +
        # geom_point(aes(group=variable, shape=type)) +
        ylab("Hazard ratio") +
        xlab("Days after ICU admission") +
        scale_linetype_manual(name=NULL, values=c(1, 2))+
        scale_x_continuous(minor_breaks = c(4, int_info(brks=brks)$tend)) +
        # scale_shape_manual(name=NULL, values=c(19, NA)) +
        # geom_text(data=data.frame(x=45, y=1, label=unique(comp.data$comparison)),
        #     aes(x=x, y=y, label=label)) +
        scale_y_continuous(
            breaks = log(sy <- c(0.25, 0.5, 0.75, 1, 1.25, 2, 4)),
            labels = sy) +
        coord_cartesian(ylim=y.lim)
    if(plot.legend) {
        gg.diff + theme(
          legend.position      = c(0, 1),
          legend.justification = c(0, 1),
          legend.direction     = "horizontal",
          legend.background    = element_rect(colour = "grey30"),
          legend.key.size      = unit(0.5, "cm"),
          legend.title         = element_blank())
    }
    else {
        gg.diff + theme(legend.position="none")
    }

}

#' @rdname ggcomparison
#' @export
ggplot_comparisons3 <- function(
  diffs.data,
  pc,
  pdf,
  pdf.scheme  = data.frame(
    protocol    = c("low", "lowmid", "mid", "midfull", "full"),
    label       = paste0("plan #", 1:5),
    title = c("complete, severly hypocaloric", "delayed, mildly hypocaloric",
      "early, mildly hypocaloric", "delayed, near target",
      "early, near target")),
  comp.effect  = TRUE,
  plot.legend  = FALSE,
  debug        = FALSE,
  brks         = 0:30,
  title.prefix = "",
  y.lim        = c(-1.5, 1.5)) {

  comparisons <- levels(as.factor(droplevels(diffs.data)$comparison))
  pc <- as.matrix(pc)

  comp.list <- list()
  if(comp.effect) {
    gg_comp <- ggcomparison
  } else {
    gg_comp <- ggcomparison_single
  }
  for(i in seq_along(comparisons)) {

    pc.i <- pc[pc[, "comparison"]==comparisons[i] , drop=TRUE]
    scheme <- ggcomparison_scheme(
      pdf[, pc.i[1]],
      pdf[, pc.i[2]],
      pdf.scheme[pdf.scheme[,1]==pc.i[1], 2],
      pdf.scheme[pdf.scheme[,1]==pc.i[2], 2],
      legend=plot.legend)
    comp   <- gg_comp(
      filter(diffs.data, comparison==comparisons[i]),
      plot.legend = plot.legend,
      comp.effect = comp.effect,
      brks=brks,
      y.lim=y.lim)
    comp.list[[i]] <- arrangeGrob(scheme, comp,
      bottom=textGrob("\n\n"),
      top=textGrob(
        paste0(
          title.prefix, letters[i], ": ", paste(pdf.scheme[pdf.scheme[,1]==pc.i[1], 3],
          pdf.scheme[pdf.scheme[,1]==pc.i[2], 3],
          sep = " vs. ")),
        gp=gpar(fontface="bold")),
      nrow=1,
      heights = c(6))

  }

  if(debug) return(comp.list)
  grid.arrange(
    grobs   = comp.list,
    ncol    = 1,
    heights = c(7, 7, 7))

}

#' Label patterns
#'
#' @importFrom utils head
#' @keywords internal
#' @export
pattern_label <- function(
    pattern,
    maxdays.tdc,
    effectname="ProteinCat") {

    patt.string <- {
        rles      <- rle(as.character(pattern))
        start     <- ifelse(maxdays.tdc==12, 0, 1)
        startdays <- head(c(start, start + cumsum(rles$lengths)), -1)
        stopdays  <- start + cumsum(rles$lengths) - 1
        tmp       <- paste0(startdays, "-", stopdays, ": ", rles$values)
        tmp[1]    <- paste("days", tmp[1])
        lvls      <- if(grepl("ProteinCat", effectname)) {
            c("C I", "C II", "C III")
        } else {
            c("<0.6 g/kg", "0.6-1.2 g/kg", ">1.2 g/kg")
        }
        paste(
            gsub(pattern="low",  replacement=lvls[1],
                gsub(pattern="mid",  replacement=lvls[2],
                    gsub(pattern="full", replacement=lvls[3], x=tmp))),
            collapse="; ")
    }

    # return
    patt.string
}



#' get differences of effects on log hazard
#'
#' @importFrom reshape2 melt
#' @keywords internal
#' @export
get_comp_diffs <- function(
  patients.list,
  protocols.df,
  protocols.to.compare,
  model,
  effectname = "proteinCat") {

  maxdays.tdc <- nrow(protocols.df)

  diff.labs <- apply(protocols.to.compare, 1, function(z) {
    paste(pattern_label(protocols.df[, z[1]], maxdays.tdc=maxdays.tdc), "vs. \n",
      pattern_label(protocols.df[, z[2]], maxdays.tdc=maxdays.tdc))
  })

  effect.diffs <- apply(protocols.to.compare, 1, function(z) {
    diff.i <- as.data.frame(
      getEffectDiffs(
        m=model,
        patients.list[[z[1]]],
        patients.list[[z[2]]],
        effectname = effectname))
  })

  names(effect.diffs) <- diff.labs
  m.effect.diffs <- melt(effect.diffs, id.vars=c("fit", "int_mid", "se", "hi", "lo"))
  m.effect.diffs$label <- factor(m.effect.diffs$L1, levels = diff.labs,
    labels=diff.labs)
  m.effect.diffs$comparison <- factor(m.effect.diffs$L1,
    levels=diff.labs, labels=levels(protocols.to.compare[, "comparison"]))

  ## return
  m.effect.diffs

}

#' @keywords internal
#' @export
pattern_pat <- function(
  pattern,
  ped,
  m,
  median.patient,
  protoPat = NULL,
  effectname    = "ProteinCat") {

  gc <- grep(effectname, names(coef(m)), value = TRUE)

  effectname <- sub(".*I\\(", "", gc)
  effectname <- unique(sub("Tot.*", "", effectname))

  low.var    <- paste0(effectname, "Tot0to30")
  mid.var    <- paste0(effectname, "Tot30To70")
  full.var   <- paste0(effectname, "TotAbove70")

  # pick a patient that remained under risk all the way to overwrite their data:
  if(is.null(protoPat)) {
    protoPat <- {
      ind <- which(ped$int_mid == max(ped$int_mid))[1]
      id <- ped$CombinedID[ind]
      subset(ped, CombinedID == id)
    }
  }
  # overwrite with median data for confounders:
  for(var in colnames(median.patient)){
    protoPat[, var] <- median.patient[1, var]
  }

  maxdays.tdc <- ncol(ped$DaysMat)
  max_per_int <- rowSums(protoPat$LHartlDynf)

  stopifnot(length(pattern) == maxdays.tdc,
    all(pattern %in% c("low", "mid", "full")))

  # switch
  protoPat[[low.var]]  <- outer(rep(1, nrow(protoPat)), pattern == "low")
  protoPat[[mid.var]]  <- outer(rep(1, nrow(protoPat)), pattern == "mid")
  protoPat[[full.var]] <- outer(rep(1, nrow(protoPat)), pattern == "full")

  protoPat

}


#' The modus of a categorical variable
#'
#' @keywords internal
#' @export
modus <- function(f){
  freqs <- table(f)
  names(freqs)[which.max(freqs)]
}

#' Extract median covariate information from columns
#'
#' @importFrom stats median
#' @importFrom stats coefficients
#' @keywords internal
#' @export
median_patient <- function(
    data,
    patient.vars  = c("ApacheIIScore", "Age", "BMI", "Year", "DiagID2",
        "AdmCatID", "Gender", "inMV2_4", "Propofol2_4", "OralIntake2_4",
        "PN2_4"),
    random.effect = "CombinedicuID",
    re.byvar      = "icuByDummy") {

    median.x <- data.frame(
        lapply(
            subset(data, select=patient.vars),
            function(x) {
                if(is.factor(x)){
                    modus(x)
                } else {
                    median(x, na.rm=TRUE)
                }
            }))
    if(!is.null(random.effect)) {
    # pic any RE factor level and set according by variable to 0
        median.x[[random.effect]] <- factor(
            x = levels(data[[random.effect]])[1])
        median.x[[re.byvar]] <- 0
    }

    # return
    median.x

}


# Cumulative incidence function

get_cif_cs <- function(
  newdata,
  object1,
  object2,
  protocol_name = "",
  n_sim = 500L,
  time_var = "int_mid",
  alpha = 0.05,
  cause_names = c("death", "discharge")
  ) {

  newdata <- newdata %>%
    mutate(intlen = tend - tstart)

  coefs_1        <- coef(object1)
  V_1            <- object1$Vp
  sim_coef_mat_1 <- mvtnorm::rmvnorm(n_sim, mean = coefs_1, sigma = V_1)

  coefs_2        <- coef(object2)
  V_2            <- object2$Vp
  sim_coef_mat_2 <- mvtnorm::rmvnorm(n_sim, mean = coefs_2, sigma = V_2)

  hazards <- purrr:::map2(
    .x = list(object1, object2),
    .y = list(sim_coef_mat_1, sim_coef_mat_2),
    .f = ~{
      .df <- newdata %>% arrange(.data[[time_var]], .by_group = TRUE)
      X <- predict(.x, .df, type = "lpmatrix")
      apply(.y, 1, function(z) exp(X %*% z))
    }
  )

  overall_survivals <- apply(
    Reduce("+", hazards),
    2,
    function(z) exp(-cumsum(z * newdata[["intlen"]])))

  hps <- purrr::map(hazards, ~ .x * (overall_survivals - 1e-10))
  cifs <- purrr::map(hps, ~ apply(.x, 2, function(z) cumsum(z * newdata[["intlen"]])))
  names(cifs) <- cause_names
  cifs_df <- purrr::imap_dfr(
    .x = cifs,
    .f = ~{
      newdata[["cif"]]       <- rowMeans(.x)
      newdata[["cif_lower"]] <- apply(.x, 1, quantile, alpha/2)
      newdata[["cif_upper"]] <- apply(.x, 1, quantile, 1-alpha/2)
      newdata[["cause"]] <- .y
      newdata[["protocol"]] <- protocol_name
      newdata
    }
  )


}
get_cif_sd <- function(
  newdata,
  object1,
  protocol_name = "",
  n_sim = 500L,
  time_var = "int_mid",
  alpha = 0.05,
  cause_names = c("death")
  ) {

  newdata <- newdata %>%
    mutate(intlen = tend - tstart)

  coefs_1        <- coef(object1)
  V_1            <- object1$Vp
  sim_coef_mat_1 <- mvtnorm::rmvnorm(n_sim, mean = coefs_1, sigma = V_1)

  hazards <- purrr:::map2(
    .x = list(object1),
    .y = list(sim_coef_mat_1),
    .f = ~{
      .df <- newdata %>% arrange(.data[[time_var]], .by_group = TRUE)
      X <- predict(.x, .df, type = "lpmatrix")
      apply(.y, 1, function(z) exp(X %*% z))
    }
  )

  overall_survivals <- purrr::map(
    .x = hazards,
    .f = ~ apply(.x, 2, function(z) exp(-cumsum(z * newdata[["intlen"]])))
  )

  hps <- purrr::map2(hazards, overall_survivals, ~.x * (.y - 1e-10))
  cifs <- purrr::map(hps, ~ apply(.x, 2, function(z) cumsum(z * newdata[["intlen"]])))
  names(cifs) <- cause_names
  cifs_df <- purrr::imap_dfr(
    .x = cifs,
    .f = ~{
      newdata[["cif"]]       <- rowMeans(.x)
      newdata[["cif_lower"]] <- apply(.x, 1, quantile, alpha/2)
      newdata[["cif_upper"]] <- apply(.x, 1, quantile, 1-alpha/2)
      newdata[["cause"]] <- .y
      newdata[["protocol"]] <- protocol_name
      newdata
    }
  )
}

get_cif_df <- function(
  comparison_list,
  object1,
  object2,
  comparison_name = "",
  n_sim = 500L,
  time_var = "int_mid",
  alpha = 0.05,
  cause_names = c("death", "discharge")) {

  .out <- purrr::imap_dfr(
    .x = comparison_list,
    .f = ~{
      get_cif_cs(.x, object1 = object1, object2 = object2, protocol_name = .y,
        n_sim = n_sim, time_var = time_var, alpha = alpha, cause_names = cause_names)
  })
  .out[["comparison"]] <- comparison_name
  .out

}

get_cif_df_sd <- function(
  comparison_list,
  object1,
  comparison_name = "",
  n_sim = 500L,
  time_var = "int_mid",
  alpha = 0.05,
  cause_names = c("death")) {

  .out <- purrr::imap_dfr(
    .x = comparison_list,
    .f = ~{
      get_cif_sd(.x, object1 = object1, protocol_name = .y,
        n_sim = n_sim, time_var = time_var, alpha = alpha, cause_names = cause_names)
  })
  .out[["comparison"]] <- comparison_name
  .out

}
