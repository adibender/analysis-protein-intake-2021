# @Author: kopperphilipp@gmail.com
# @Date:   2020-03-09 17:19:23
# @Description:
# This script modifies an earlier script of @andreas.bender@stat.uni-muenchen.de
# to fit the current data situation.

## load packages with additional functionality
library(lubridate) # handling date formats/variables
library(checkmate) # checking inputs
library(dplyr)     # data manipulation library
library(haven)

# create folders
dir.create("data")

source("preproc/preproc-helpers.R")

icu_sas <-
  read_sas("preproc/data/icuinfo.sas7bdat", encoding = "latin1")
patient_sas <-
  read_sas("preproc/data/patientinfo.sas7bdat", encoding = "latin1")
daily_sas <-
  read_sas("preproc/data/obsinfo.sas7bdat", encoding = "latin1")

ICU <- icu_sas
patient <- patient_sas
daily <- daily_sas

###################################### ICU #####################################
icu.vars.old <- c(
  "year"          , "combinedicuid"       , "countryid"     ,
  "hosp_type"     , "hosp_beds"           , "multi_icu"     ,
  "icu_type"      , "medical_director_yn" , "icu_beds"      ,
  "dietitianwork" , "feeding_protocol_yn" , "grv"           ,
  "hosp_name"     , "icu_name"   )


icu.vars.new <- c(
  "Year"          , "CombinedicuID"       , "CountryID"     ,
  "InstTypeID"    , "hospsize"            , "multipleICUs"  ,
  "ICUTypeID"     , "MedDirector"         , "NumBeds"       ,
  "dietitianwork" , "HasProtocols"        , "GastricVolume" ,
  "hosp_name"     , "icu_name"  )


ICU <- ICU[, icu.vars.old]
for( i in seq_along(icu.vars.new) ) {
  colnames(ICU)[colnames(ICU) == icu.vars.old[i]] <- icu.vars.new[i]
}

################################### PATIENT ####################################
#### rename variables
## define variables (the ones which appear in previous evaluations)
patient.vars.old <- c(
  "year"             , "combinedicuid" ,
  "study_id"         , "sex"           ,
  "age"              ,
  "mv_start_yn"      , "admission_type",
  "DiagID"           ,
  "apache_score"     , "height"        ,
  "weight"           , "bmi"           ,
  "calsperkg"        , "protperkg",
  "patientdied"      ,
  "pt_die_icu_yn"    ,
  "mv_discontinued_icu_yn"
  )

patient.vars.new <- c(
  "Year"             , "CombinedicuID"     ,
  "CombinedID"       , "Gender"            ,
  "Age"              ,
  "MechVentPrior"    , "AdmCatID"          ,
  "DiagID"           ,
  "ApacheIIScore"    , "Height"            ,
  "Weight"           , "BMI"               ,
  "Calories"         , "Protein"           ,
  "PatientDied"      ,
  "icuDischargeDeath",
  "VentDisICU"
  )

#### time variables
td.vars.old <- c(
  "mvstartdt"         , "endt"              , "pndt"             ,
  "hospadmindt"       , "icuadmindt"        , "deathdt"          ,
  "icudischargedt"    , "hospdischargedt"   , "mvdiscontinuedt")

# patient.vars.new old dataset
td_vars_new <- c(
  "MechVent"          , "enTime"            , "pnTime"           ,
  "HospAdmissionDate" , "AdmissionDate"     , "PatientDeathTime" ,
  "icuDischargeDate"  , "HospDischargeDate" , "ventDisTime")

vars.to.rename <- c(patient.vars.old, td.vars.old)
patient.td_vars_new <- c(patient.vars.new, td_vars_new)

## first subset patient data
patient <- patient[, vars.to.rename]

# rename to fit old script, via loop since ordering is crucial
for(i in seq_along(vars.to.rename) ) {
  ind <- which(colnames(patient) == vars.to.rename[i])
  colnames(patient)[ind] <- patient.td_vars_new[i]
}

# change time variables to POSIX format
lct <- Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "C")

# check integrity of variables containing date information
sapply(td_vars_new, function(x) {
  tmp  <- patient[[x]]
  tmp2 <- as.POSIXct(tmp, tz="CET")
  if(any(is.na(tmp2)!= is.na(tmp))) browser()
  patient[[x]] <<- tmp2
  return(NULL)
})

Sys.setlocale("LC_TIME", lct)

patient$CombinedID <- as.integer(gsub("_", "", patient$CombinedID))

####################################### DAILY ##################################
## fix variable names
daily.old.vars <- c(
  "year"             , "combinedicuid"  ,
  "study_id"         , "study_day"       ,
  "en_given_yn"      , "pn_nutrition_yn" ,
  "oral_nutrition_yn", "morning_glucose" ,
  "propofol_yn"      , "propofolcalories",
  "encalories"       , "enprotein"       ,
  "pn_calories"      , "pn_protein"      ,
  "AdequacyCalsEN"   , "AdequacyProtEN"  ,
  "AdequacyCalsTot"  , "AdequacyProtTot")


daily.new.vars <- c(
  "Year"          ,   "CombinedicuID",    "CombinedID"   ,    "Study_Day",
  "EN"              , "PN"             ,  "OralIntake"     ,  "Morning_BS",
  "Propofol"        , "PropofolCal"    ,  "EN_Calories"    ,  "EN_Protein",
  "PN_Calories"    ,  "PN_Protein"     ,  "AdequacyCalsEN",
  "AdequacyProtEN"  , "AdequacyCalsTot",  "AdequacyProtTot")

## subset daily dataset
daily <- select(daily, one_of(daily.old.vars))

## change names
for( i in seq_along(daily.old.vars) ) {
  colnames(daily)[colnames(daily) == daily.old.vars[i]] <- daily.new.vars[i]
}

daily$CombinedID <- as.integer(gsub("_", "", daily$CombinedID))
daily_all <- daily
saveRDS(daily_all, "data/daily_all.Rds")
#### check td_vars_new for inconsistencies
## mech vent discontinuation before mech vent start
# we start with: 21100
nrow(patient)
ind.start.before.end <- patient$MechVent > patient$ventDisTime &
  !(is.na(patient$MechVent > patient$ventDisTime))
patient[ind.start.before.end, c("MechVent", "ventDisTime", "HospAdmissionDate",
                                "AdmissionDate")]

removed_patients <- data.frame(id = patient$CombinedID[ind.start.before.end],
                               cause = "inconsistency",
                               cause2 = "mech vent disc. before start")

patient <- patient[!ind.start.before.end, ]
## mech vent discontinuation equal to mech vent start?
ind.start.equal.end <- patient$MechVent == patient$ventDisTime &
  !(is.na(patient$MechVent == patient$ventDisTime))
pEq <- patient[ind.start.equal.end, c("CombinedID", td_vars_new[c(1, 5, 7, 9)])]
dim(pEq)
## 10 patients with equal mech vent start and discontinuation time;
# we assume that this is ok and leave the patients in the data set

## en/pn start before Admission?
en.before.adm <- patient$enTime - patient$AdmissionDate
sum(en.before.adm < 0, na.rm = TRUE)
hist(as.numeric(
  en.before.adm[en.before.adm < 0 & (!is.na(en.before.adm))]) / 360 ,
     xlab = "hours of discrp.", main = "EN")
# 188 patients with EN before Admission!
# Partly within 24 hours or within 48 hours, for these we assumed this is OK.
# The rest is removed (129 obs)
sum(en.before.adm < - (48 * 360), na.rm = TRUE)

removed_patients <-
  rbind(
    removed_patients,
    data.frame(id = patient$CombinedID[!is.na(en.before.adm)][
      (en.before.adm[!is.na(en.before.adm)] < - (48 * 360))],
      cause = "inconsistency",
      cause2 = "EN before Admission"))

patient <- patient[(!en.before.adm < - (48 * 360)) | is.na(en.before.adm), ]

pn.before.adm <- patient$pnTime - patient$AdmissionDate
sum(pn.before.adm < 0, na.rm = TRUE)
hist(as.numeric(
  pn.before.adm[pn.before.adm < 0 & (!is.na(pn.before.adm))]) / 360 ,
  xlab = "hours of discrp.", main = "PN")
# 47 patients with PN before Admission! (in total there were 52 before EN)
# Some within 24 hours or within 48 hours.
# The rest is removed (22 obs)
sum(pn.before.adm < - (48 * 360), na.rm = TRUE)

removed_patients <- rbind(
  removed_patients,
  data.frame(id = patient$CombinedID[!is.na(pn.before.adm)][
    (pn.before.adm[!is.na(pn.before.adm)] < - (48 * 360))],
    cause = "inconsistency",
    cause2 = "EN before Admission")
)

patient <- patient[(!pn.before.adm < - (48 * 360)) | is.na(pn.before.adm), ]

## death before en/pn?
assert_tdiff(patient, "PatientDeathTime", "enTime", threshold=0)# no!
assert_tdiff(patient, "PatientDeathTime", "pnTime", threshold=0)# no!

## en/pn before hospital admission?
check_time_consistency(patient, "enTime", "HospAdmissionDate", threshold=24)
# several cases where EN before Hosp. admission, but mostly the same day or
# less than 24h difference, thus we assumed this is OK
check_time_consistency(patient, "pnTime", "HospAdmissionDate", threshold=24)
# same as enTime before HospAdmissionDate

## icu discharge before admission?
assert_tdiff(patient, "icuDischargeDate", "AdmissionDate", threshold=0)#no!

## icu admission before hospital admission?
#assert_tdiff(patient, "AdmissionDate", "HospAdmissionDate", threshold=0)
sum(patient$AdmissionDate < patient$HospAdmissionDate, na.rm = TRUE)
## only 4 obs: we can drop them

removed_patients <- rbind(
  removed_patients,
  data.frame(
    id = patient$CombinedID[((patient$AdmissionDate < patient$HospAdmissionDate)) &
                              !(is.na(patient$AdmissionDate < patient$HospAdmissionDate))],
    cause = "inconsistency",
    cause2 = "icu admission befor hosp admission")
)

patient <- patient[(!(patient$AdmissionDate < patient$HospAdmissionDate)) |
                     is.na(patient$AdmissionDate < patient$HospAdmissionDate), ]

## hospital discharge before icu discharge?
check_time_consistency(patient, "HospDischargeDate", "icuDischargeDate")
# 5 patients discharged from hospital before being discharged from icu.
# 2 seem to be just inaccurateley timed in the system. Keep them.
# The rest is removed.

removed_patients <- rbind(
  removed_patients,
  data.frame(
    id = c(20101705, 20101927, 20102315),
    cause = "inconsistency",
    cause2 = "icu discharge before hosp discharge")
)

patient <- patient[!(patient$CombinedID %in% c(20101705, 20101927, 20102315)), ]

## posthumous intubation?
mechvent.after.death <- check_time_consistency(patient, "PatientDeathTime", "MechVent")
## No.

## remove patients without prescribed calories or proteins (72)
sum((is.na(patient$Calories) | is.na(patient$Protein)))

removed_patients <- rbind(
  removed_patients,
  data.frame(
    id = patient$CombinedID[(is.na(patient$Calories) | is.na(patient$Protein))],
    cause = "inconsistency",
    cause2 = "No prescribed calories / protein")
)


patient <- filter(patient, !(is.na(Calories) | is.na(Protein)))

## patient death time before admission?
assert_tdiff(patient, "PatientDeathTime", "AdmissionDate", threshold=0)#no!

## mech vent (long) after admission?
sum(patient$MechVent - patient$AdmissionDate > 96 * 60 * 60, na.rm = TRUE)
diff.mv.adm <- difftime(patient$MechVent, patient$AdmissionDate, unit = "hours")

## mech vent discontinuation before admission?
id.mvend.before.adm <- patient %>% filter(ventDisTime < AdmissionDate) %>%
  select(CombinedID)
id.mvend.before.adm
## exclude 1 patient with mv end before admission

removed_patients <- rbind(
  removed_patients,
  data.frame(
    id = as.numeric(id.mvend.before.adm),
    cause = "inconsistency",
    cause2 = "Mech vent end before admission")
)

patient <- filter(patient, !(CombinedID %in% id.mvend.before.adm$CombinedID))

#### No patient with HospDischargeDate < Admission Date
ind.disch.before.adm <- patient$HospDischargeDate < patient$AdmissionDate
nrow(patient[ind.disch.before.adm & !is.na(ind.disch.before.adm), td_vars_new])

#### remove (10) patients with missing indication of death/survival
nrow(filter(patient, is.na(PatientDied)))
patient <- filter(patient, !(is.na(PatientDied)))

#### 103 without BMI
nrow(filter(patient, is.na(BMI)))

removed_patients <- rbind(
  removed_patients,
  data.frame(
    id = patient$CombinedID[is.na(patient$BMI)],
    cause = "inconsistency",
    cause2 = "No BMI")
)

patient <- filter(patient, !(is.na(BMI)))

#### 354 without AII Score
nrow(filter(patient, is.na(ApacheIIScore)))

removed_patients <- rbind(
  removed_patients,
  data.frame(
    id = patient$CombinedID[is.na(patient$ApacheIIScore)],
    cause = "inconsistency",
    cause2 = "No AIIScore")
)

patient <- filter(patient, !(is.na(ApacheIIScore)))

#### 3 without sex information (only binary coding applied here).
nrow(filter(patient, is.na(Gender)))

removed_patients <- rbind(
  removed_patients,
  data.frame(
    id = patient$CombinedID[is.na(patient$Gender)],
    cause = "inconsistency",
    cause2 = "No Sex")
)

patient <- filter(patient, !(is.na(Gender)))

### extremely high / low EN Protein / PN Protein?                               ##### REMOVE!?
summary(daily$EN_Protein)
summary(daily$PN_Protein)
# investigate
daily[daily$PN_Protein > 300, c("CombinedID", "Study_Day", "PN_Protein")]
ids_PN_too_high <- daily[daily$PN_Protein > 300, c("CombinedID")] %>%
  unlist()
daily[daily$CombinedID %in% ids_PN_too_high,
      c("CombinedID", "Study_Day", "PN_Protein")]
# more likely 70.8 (71) and 46 respectively
daily[daily$PN_Protein > 300, "PN_Protein"] <- c(71, 46)

## check PN / EN
summary(daily$EN)
summary(daily$PN)
# one case with PN == 2, make it 1
daily$PN[daily$PN == 2] <- 1L

## exclude patients with missing admission date and weight (none!)
nrow(filter(patient, is.na(AdmissionDate)))
nrow(filter(patient, is.na(Weight)))
patient <- filter(patient, !(is.na(AdmissionDate) | is.na(Weight)))

# check variables
assert_numeric(as.numeric(patient$AdmissionDate), any.missing=FALSE,
               finite=TRUE, min.len=1)
assert_numeric(patient$Weight, any.missing=FALSE, finite=TRUE, lower=30,
               upper=350, min.len=1)

print(paste("Excluded", nrow(removed_patients), "obs due to data quality."))
### %%%%%%%%%%%%%%%%%%%%%% END SANITY CHECKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ###


#### redefine Diagnosis variable (see report docu)
patient$DiagID <- as.factor(patient$DiagID)
(levs.diag <- levels(patient$DiagID))
addmargins(table(patient$DiagID))
new.diag <- numeric(nrow(patient))
diag_list <- list(cv     = c(1:9, 50:54, 52.1, 52.2),
                  resp   = c(10:19, 57:58),
                  gastro = c(21:27, 61:68),
                  neuro  = c(29:35, 70:74),
                  sepsis = c(37, 38, 82:85),
                  trauma = c(39, 40, 76, 77, 81),
                  meta   = c(41:43),
                  renal  = c(47, 78))
new.diag[patient$DiagID %in% diag_list[[1]]] <- 1 # CV
new.diag[patient$DiagID %in% diag_list[[2]]] <- 2 # Resp
new.diag[patient$DiagID %in% diag_list[[3]]] <- 3 # Gastro
new.diag[patient$DiagID %in% diag_list[[4]]] <- 4 # neuro
new.diag[patient$DiagID %in% diag_list[[5]]] <- 5 # sepsis
new.diag[patient$DiagID %in% diag_list[[6]]] <- 6 # trauma
new.diag[patient$DiagID %in% diag_list[[7]]] <- 7 # metabolic
new.diag[patient$DiagID %in% diag_list[[8]]] <- 8 # renal
new.diag[!(new.diag %in% 1:8)] <- 9 #other

patient$DiagID2 <- factor(new.diag, labels = c(
  "Cardio-Vascular", "Respiratory",  "Gastrointestinal",
  "Neurologic", "Sepsis", "Orthopedic/Trauma",
  "Metabolic", "Renal", "Other"))

# remove old definition
patient$DiagID <- NULL

#### calculate survival time
#### We define discharge from hospital as the LATEST measured discharge event
#### in case there were discrepancies. However discharge from hospital must
#### have been tracked.
s0to60 <- patient$PatientDeathTime - patient$AdmissionDate

discharge_hosp <- difftime(
  patient$HospDischargeDate,
  patient$AdmissionDate,
  unit = "days")
discharge_icu <- difftime(
  patient$icuDischargeDate,
  patient$AdmissionDate,
  unit = "days")
ind_na_death <- is.na(s0to60)
ind_na_icu <- is.na(discharge_icu)
ind_na_hosp <- is.na(discharge_hosp)
table(ind_na_icu, ind_na_hosp)
# patients with death time but PatientDied == 0
pat_sub  <- patient[which(!ind_na_death & patient$PatientDied == 0),
  c("CombinedID", "AdmissionDate", "PatientDied", "PatientDeathTime", "HospDischargeDate", "icuDischargeDate") ]
pat_sub  %>% print(n = 100)
pat_sub$PatientDeathTime - pat_sub$AdmissionDate  # all more than 60 days, so would be censored at 60

discharge_hosp[is.na(discharge_hosp)] <- discharge_icu[is.na(discharge_hosp)]
#discharge_icu[is.na(discharge_icu)] <- discharge_hosp[is.na(discharge_icu)]
#discharge_time <- pmax(discharge_hosp, discharge_icu)
#discharge_time[discharge_time == 0L] <- NA
discharge_time <- discharge_hosp
discharge_time[patient$PatientDied == 1] <- NA
d0to60 <- discharge_time
patient$Surv0To60 <- as.numeric(s0to60)
patient$Disc0To60 <- as.numeric(d0to60)

# in the following we consider time until death in ICU. subjects that die
# after discharge from ICU are therefore considered censored w.r.t. to death on ICU
# release from ICU will be coded as cause =1
# death in ICU as cause = 2
# if subjects remain in ICU for t > 60 they will be considered as censored
surv_icu0to60 <- s0to60
surv_icu_status <- NA_integer_
# if death time is given and < 60 days, these are deaths, cause will be set to 2
# death time > 60 days will be administratively censored
surv_icu_status[!is.na(surv_icu0to60)] <- 2L
# patient discharged from ICU before death are considered a competing risk for death
# the time until event is set to discharge time and cause to 1
ind_disc_icu_before_death <- which(discharge_icu < surv_icu0to60)
surv_icu0to60[ind_disc_icu_before_death] <- discharge_icu[ind_disc_icu_before_death]
surv_icu_status[ind_disc_icu_before_death] <- 1L
# subjects where discharge time from icu is given but no information on death time
# time until event is set to discharge time, cause is set to 1
ind_disc_icu <- which(!is.na(discharge_icu) & is.na(surv_icu0to60))
surv_icu0to60[ind_disc_icu] <- discharge_icu[ind_disc_icu]
surv_icu_status[ind_disc_icu] <- 1L
# remaining NAs are assumed to not have been discharged from ICU during the observation period
# time until discharge will be set to 61, cause will be set to 0, i.e. censoring
ind_na_remaining <- which(is.na(surv_icu0to60))
length(ind_na_remaining) # # these are 1117 subjects that were not discharged from icu until end of follow up
surv_icu0to60[ind_na_remaining] <- 61
surv_icu_status[ind_na_remaining] <- 0L
# in case there are further observations with time until event > 60
# set time until event to 61, cause to 0, i.e. censored
surv_icu0to60[surv_icu0to60 > 60] <- 61
surv_icu_status[surv_icu0to60 > 60] <- 0L
patient$surv_icu0to60 <- surv_icu0to60
patient$surv_icu_status <- surv_icu_status


### check how many missing events exist:
sum(is.na(s0to60) & is.na(d0to60))
### There are more than 1000 "missing" events (1117).
missing_events <- patient[is.na(s0to60) & is.na(d0to60), ]
summary(missing_events)
### For all we find a missing discharge date, however, we know that ALL survived.
### We use last_date_icu for these missing ones (all receive a 61L).
patient$Disc0To60[is.na(s0to60) & is.na(d0to60)] <- 61L
patient$AdCens <- ifelse(is.na(s0to60) & is.na(d0to60), 1, 0)
### if we measured longer survival /discharge, we censor it to 61L.
### We leave an indicator that we already censored these values.
### (assumed to have survived, Note that we do not make further use of the
### indicator.)
patient$Disc0To60[patient$Surv0To60 > 61 & (!is.na(patient$Surv0To60))] <- 61L
patient$Surv0To60[patient$Surv0To60 > 61 & (!is.na(patient$Surv0To60))] <- 61L
patient$PatientDied[patient$Surv0To60 > 61 & (!is.na(patient$Surv0To60))] <- 0L
patient$Disc0To60[patient$Disc0To60 > 61 & (!is.na(patient$Disc0To60))] <- 61L
patient$PatientDischarged <- ifelse(patient$PatientDied == 0 & patient$Disc0To60 < 61L & !is.na(patient$Disc0To60), 1, 0)

#check
table(patient$PatientDied, patient$surv_icu_status)
# 1248 patients died but were released from icu prior to death -> status = 1
# 3771 patients died while at icu -> status 2
# 5 patients who died but after 60 days -> admistrative censoring
patient %>%
  filter(PatientDied == 1, surv_icu_status == 0) %>%
  select(PatientDied, AdmissionDate, PatientDeathTime, icuDischargeDate, surv_icu0to60)


### Compute days in ICU
# First create a variable that give last date time in ICU
last_date_icu <- patient$icuDischargeDate
last_date_icu[is.na(last_date_icu) & patient$PatientDied == 1] <-
  patient$PatientDeathTime[is.na(last_date_icu) & patient$PatientDied == 1]

last_date_icu[is.na(last_date_icu)] <-
  patient$AdmissionDate[is.na(last_date_icu)] + days(61L)
patient$last_date_icu <- last_date_icu
patient$DaysInICU <-
  as.numeric(with(patient,
                  difftime(last_date_icu, AdmissionDate, units = "days")))
# check variable
assert_numeric(patient$DaysInICU, lower = 0, min.len = 1, any.missing = FALSE,
  finite = TRUE)

# create event variable
patient$event <- patient$Surv0To60
patient$event[is.na(patient$event)] <- patient$Disc0To60[is.na(patient$event)]

# Either > 60 or tiny difference, negligible
filter(patient, event < surv_icu0to60) %>%
  select(PatientDied, surv_icu0to60, event)

# calendar days in ICU, will be needed for descriptive analysis and application
# of exclusion criteria
patient$calendarDaysInICU <-
  round(as.numeric(with(patient,
                        difftime(
                          ceiling_date(last_date_icu, unit = "day"),
                          floor_date(AdmissionDate, unit = "day"),
                          units = "days"))))
# check variable
assert_numeric(patient$calendarDaysInICU, lower = 0, min.len = 1,
               any.missing = FALSE, finite = TRUE)

## 2445 patients died before day twelve
sum(patient$Surv0To60 < 12, na.rm = TRUE)
## 356 out of these were only though "early dischargess"
## 318 if the first 8 hours are not considered.
## This seems reasonable
died_early <- filter(patient, Surv0To60 < 12)
sum(difftime(died_early$PatientDeathTime,
             died_early$icuDischargeDate, units = "hours") > 0, na.rm = TRUE)
sum(difftime(died_early$PatientDeathTime,
             died_early$icuDischargeDate, units = "hours") > 8, na.rm = TRUE)

# 11554 patients died or realeased from hospital before 12 days
sum(patient$surv_icu0to60 < 12, na.rm = TRUE)
# 2089 of which experienced an event before 12 days
sum(patient$surv_icu0to60 < 12 & patient$surv_icu_status == 2, na.rm = TRUE)

## Compute calender days to event
patient$calendarEvent <- 0
patient$calendarEvent[is.na(patient$Disc0To60)] <-
  with(patient[is.na(patient$Disc0To60), ],
       round(difftime(ceiling_date(PatientDeathTime, unit = "day"),
                      floor_date(AdmissionDate, unit = "day"),
                      units = "days")))
patient$calendarEvent[is.na(patient$Surv0To60)] <-
  with(patient[is.na(patient$Surv0To60), ],
       round(difftime(ceiling_date(HospDischargeDate, unit = "day"),
                      floor_date(AdmissionDate, unit = "day"),
                      units = "days")))
sum(is.na(patient$calendarEvent))
patient$calendarEvent[is.na(patient$calendarEvent)] <- 62L
# set NAs to maximal follow-up for patients who survived
patient$calendarEvent[patient$calendarEvent >= 62L] <- 62L
assert_numeric(patient$calendarEvent, lower = 0, upper = 62,
               any.missing = FALSE, min.len = 1)

# if mech vent occured prior to admission, we set starting point to
# time of admission.
# We set beginMechVent to 99 (after max follow up), when there was no prior mech
# vent and mech vent is missing
patient$beginMechVent <- with(patient, {
  tmp <- as.numeric(difftime(MechVent, AdmissionDate, units = "days"))
  tmp[!is.na(MechVentPrior)] <- 0
  tmp[is.na(MechVentPrior) & is.na(MechVent)] <- 99
  tmp
})
assert_numeric(patient$beginMechVent, lower = 0, upper = 99, min.len = 1,
               any.missing = FALSE, finite = TRUE)

patient$calendarBeginMechVent <- round(with(patient, {
  tmp <- as.numeric(
    difftime(
      ceiling_date(MechVent, unit = "day"),
      floor_date(AdmissionDate, unit = "day"),
      units = "days"))
  tmp[!is.na(MechVentPrior)] <- 1
  tmp[is.na(MechVentPrior) & is.na(MechVent)] <- 99
  tmp
}))
assert_numeric(patient$calendarBeginMechVent, lower = 0, upper = 99,
               min.len = 1, any.missing = FALSE, finite = TRUE)

# same as begin mech vent
patient$endMechVent <- with(patient, {
  tmp <- as.numeric(difftime(ventDisTime, AdmissionDate, units = "days"))
  tmp[tmp > 60] <- 61
  #if no ventDisTime given MechVent presumed to go on until death or censoring
  tmp[is.na(ventDisTime)] <- event[is.na(ventDisTime)]
  tmp
})
assert_numeric(patient$endMechVent, lower = 0, upper = 61, min.len = 1,
               any.missing = FALSE, finite = TRUE)

patient$calendarEndMechVent <- round(with(patient, {
  tmp <- as.numeric(
    difftime(
      ceiling_date(ventDisTime, unit = "day"),
      floor_date(AdmissionDate, unit = "day"),
      units = "days"))
  tmp[tmp > 60] <- 61
  #if no ventDisTime given MechVent presumed to go on until death or censoring
  tmp[is.na(ventDisTime)] <- calendarEvent[is.na(ventDisTime)]
  tmp
}))
# upper can be 62, because often calendar days are one more than 24h periods
assert_numeric(patient$calendarEndMechVent, lower = 0, upper = 62, min.len = 1,
               any.missing = FALSE, finite = TRUE)

# calculate actual time with mechanical ventilation
patient$DaysMechVent <- with(patient, endMechVent - beginMechVent)
patient$DaysMechVent[patient$beginMechVent == 99] <- 0
patient$DaysMechVent[patient$DaysMechVent < 0] <- 0

#### ICU ids correct
summary(ICU)
sum(table(ICU$CombinedicuID, ICU$Year) > 1)
nrow(unique((ICU[, !(colnames(ICU) %in% "CombinedicuID")])))
nrow(ICU)
# some ICUs are doubled, however there is only one true obvious error:
duplic_ICU <- ICU[duplicated(ICU[, - c(2, 14)], fromLast = FALSE), ]
duplic_ICU_2 <- ICU[duplicated(ICU[, - c(2, 14)], fromLast = TRUE), ]
duplic_ICU <- rbind(duplic_ICU, duplic_ICU_2)
data.frame(duplic_ICU[order(duplic_ICU$hospsize), c(1, 2, 5, 9, 13, 14)])
ICU$CombinedicuID[ICU$CombinedicuID == 11306] <- 2306
patient$CombinedicuID[patient$CombinedicuID == 11306] <- 2306
ICU <- ICU[, 1:12]

######################## daily data   ##########################################
# remove all observations with unspecified patient ID
# none here
daily <- filter(daily, !is.na(CombinedID))
# remove all observations from patients that are not in patient data
sum(!(daily$CombinedID %in% patient$CombinedID))
## 6363 patient days -- most due to the earlier deletion of associated patients.
daily <- filter(daily, CombinedID %in% patient$CombinedID)
# check if set of CombinedID is equal in patient and daily data set
assert_set_equal(unique(daily$CombinedID), patient$CombinedID)
## sort to have proper study_day ordering
## later on we need to convert the types back
daily <- arrange(daily, CombinedID, Study_Day, Year)
daily$Year <- as.integer(daily$Year)
daily$CombinedicuID <- as.integer(daily$CombinedicuID)
patient$Year <- as.integer(patient$Year)
patient$CombinedicuID <- as.integer(patient$CombinedicuID)

# merge patient and daily data
daily <- inner_join(daily, patient,
                    by = c("CombinedID" = "CombinedID",
                           "CombinedicuID" = "CombinedicuID",
                           "Year" = "Year"))
# create variable that indicates maximal number of days with nutrition protocol
# per patient (maxday)
daily <- daily %>%
  group_by(CombinedID) %>%
  mutate(maxday = max(Study_Day)) %>% ungroup

## new variable for total caloric intake (EN + PN + Propofol)
daily$EN_Calories[is.na(daily$EN_Calories)] <- 0
daily$PN_Calories[is.na(daily$PN_Calories)] <- 0
daily$PropofolCal[is.na(daily$PropofolCal)] <- 0
daily$caloriesIntake <- select(daily, EN_Calories, PN_Calories, PropofolCal) %>%
  rowSums(na.rm = TRUE)
# check variable
assert_numeric(daily$caloriesIntake, lower = 0, finite = TRUE,
               any.missing = FALSE, min.len = 1)

## in general: we need to make EN_/PN_Protein and EN/PN 0 where NA
daily$EN[is.na(daily$EN)] <- 0
daily$PN[is.na(daily$PN)] <- 0
daily$EN_Protein[is.na(daily$EN_Protein)] <- 0
daily$PN_Protein[is.na(daily$PN_Protein)] <- 0

# investigate mutual inconsistencies:
ids_PN_wrong <- unique(daily$CombinedID[daily$PN > 0 &
  daily$PN_Protein == 0 & daily$PN_Calories == 0])
ids_EN_wrong <- unique(daily$CombinedID[daily$EN > 0 &
  daily$EN_Protein == 0 & daily$EN_Calories == 0])
daily[daily$CombinedID %in% ids_EN_wrong, ][19:28, c("CombinedID", "Study_Day", "EN", "EN_Protein", "EN_Calories", "DaysInICU")]
# make last day intentionally incomplete if inconsistent OR event time indicates
# that incomplete
# Most extreme case: patient is admitted on Jan 1 23:59
# and discharged on Jan 8 00:01
# -> Jan 8 is her 8th Study Day while she only spent six complete days in ICU
# -> two runs of this seems justified, all further IDs seem to be intrinsically
# problematic.
daily <- daily[daily$CombinedID %in% patient$CombinedID, ]
daily$uni <- 1:nrow(daily)
helper_daily <- daily %>%
  group_by(CombinedID) %>% filter(row_number() %in% (n() - 1):n()) %>%
  select(c(CombinedID, Study_Day, EN, PN, EN_Protein, PN_Protein,
           EN_Calories, PN_Calories, calendarDaysInICU, uni))
inconsistent_EN <- (helper_daily$EN > 0) & (helper_daily$EN_Calories == 0) &
  (helper_daily$EN_Protein == 0)
inconsistent_PN <- (helper_daily$PN > 0) & (helper_daily$PN_Calories == 0) &
  (helper_daily$PN_Protein == 0)
inconsistent_last_day_icu <- helper_daily$Study_Day >= helper_daily$calendarDaysInICU
existing_intake <- (helper_daily$EN_Protein + helper_daily$PN_Protein) > 0
inconsistent_last_day <- inconsistent_EN | inconsistent_PN | (inconsistent_last_day_icu & !existing_intake)
helper_daily <- helper_daily[inconsistent_last_day, ]
daily <- daily[!(daily$uni %in% helper_daily$uni), ]
daily$uni <- NULL
#removed
removed_patients <-
  rbind(
    removed_patients,
    data.frame(id = patient$CombinedID[!(patient$CombinedID %in% daily$CombinedID)],
               cause = "inconsistency",
               cause2 = "EN/PN values inherently inconsistent"))
patient <- patient[patient$CombinedID %in% daily$CombinedID, ]
# still have to remove obs of inconclusive last days that are inconsistent
ids_PN_wrong <- unique(daily$CombinedID[daily$PN > 0 & daily$PN_Calories == 0])
ids_EN_wrong <- unique(daily$CombinedID[daily$EN > 0 & daily$EN_Calories == 0])
daily[daily$CombinedID %in% ids_EN_wrong, c("CombinedID", "Study_Day", "EN", "EN_Protein", "EN_Calories", "calendarDaysInICU", "DaysInICU", "OralIntake")]
# if no EN/PN calories have been reported --> assoc. EN/PN are zero
daily$EN <- ifelse(daily$EN_Calories == 0, 0, daily$EN)
daily$PN <- ifelse(daily$PN_Calories == 0, 0, daily$PN)
daily$EN_Protein <- ifelse(daily$EN == 0, 0, daily$EN_Protein)
daily$PN_Protein <- ifelse(daily$PN == 0, 0, daily$PN_Protein)

## new variable for total protein intake (EN + PN)
# no protein through propofol
daily$proteinIntake <- select(daily, EN_Protein, PN_Protein) %>%
  rowSums(na.rm = TRUE)
assert_numeric(daily$proteinIntake, lower = 0, finite = TRUE,
               any.missing = FALSE, min.len = 1)
# adjusted protein intake (see Email WH_2012_12_11_Protein)
daily <- mutate(daily,
                proteinIntakeAdjusted =
                  rowSums(cbind(EN_Protein, PN_Protein * 0.83), na.rm = TRUE))
# check variable
assert_numeric(daily$proteinIntakeAdjusted, lower=0, finite=TRUE,
               any.missing = FALSE, min.len = 1)

## assert that all unique IDs match:
assert_set_equal(patient$CombinedID, unique(daily$CombinedID))
# check icu id
ICU <- filter(ICU, CombinedicuID %in% unique(patient$CombinedicuID))
assert_set_equal(unique(patient$CombinedicuID), unique(ICU$CombinedicuID))
assert_set_equal(unique(daily$CombinedicuID), unique(ICU$CombinedicuID))

#
daily$PropofolCal[is.na(daily$Propofol) & (daily$PropofolCal != 0L) &
                    !is.na(daily$PropofolCal)] <- 1L
daily$Propofol[is.na(daily$Propofol)] <- 0L
daily$Propofol <- factor(daily$Propofol, exclude = NULL)
## oral intake is na --> "NO"
summary(daily$OralIntake)
daily$OralIntake[is.na(daily$OralIntake)] <- 0L
# check variable
assert_factor(daily$Propofol, levels = c("0", "1"), ordered=FALSE,
              any.missing=FALSE)

# CountryID in daily and patient not the same as in ICU data
# remove as it's not needed anyway
patient$CountryID <- daily$CountryID <- NULL
patientdaily <- daily

# Year variable to factor
patient$Year      <- as.factor(patient$Year)
ICU$Year          <- as.factor(ICU$Year)
daily$Year        <- as.factor(patientdaily$Year)
patientdaily$Year <- as.factor(patientdaily$Year)
ICU$CombinedicuID <- as.integer(ICU$CombinedicuID)
patientdaily$CombinedicuID <- as.integer(patientdaily$CombinedicuID)

# merge patientdaily with ICU information
data <- inner_join(patientdaily, ICU,
                  by = c("CombinedicuID"="CombinedicuID", "Year"="Year"))
nrow(data) == nrow(patientdaily)#TRUE
levels(data$Year)

## add indicator for mechanical ventilation (for each study day)
maxdays.nutri <- 12L
data <-
  mutate(data,
         inMV =
           as.integer((Study_Day >= pmax(1, calendarBeginMechVent) &
                         Study_Day <= pmin(maxdays.nutri, calendarEndMechVent))))
assert_integer(data$inMV, lower = 0, upper = 1, any.missing = FALSE, min.len = 1)
data %>%
  select(one_of(c("AdmissionDate", "icuDischargeDate",
                  "CombinedID"            , "MechVent"            , "ventDisTime"       ,
                  "calendarBeginMechVent" , "calendarEndMechVent" , "calendarDaysInICU" ,
                  "Study_Day"             , "maxday"              , "inMV"))) %>%
  print(n = 40, width = Inf)

data <- mutate(data,
               inProtocol=as.integer(Study_Day <= pmin(maxdays.nutri, calendarDaysInICU)))

# compute relative intake adequacy
# Calories/Protein is the prescribed amount (at admission)
assert_numeric(data$Calories, any.missing=FALSE, finite=TRUE, min.len=1, lower=0)
assert_numeric(data$Protein,  any.missing=FALSE, finite=TRUE, min.len=1, lower=0)
data$caloriesPercentage        <- (data$caloriesIntake/ (data$Calories * data$Weight)) * 100
data$proteinPercentage         <- (data$proteinIntake / (data$Protein * data$Weight)) * 100
data$proteinAdjustedPercentage <- data$proteinIntakeAdjusted / data$Protein
data$calproKg                  <- data$caloriesIntake / data$Weight

# needed in transformToPEM.R
daily$caloriesPercentage        <- data$caloriesPercentage
daily$proteinPercentage         <- data$proteinPercentage
daily$proteinAdjustedPercentage <- data$proteinAdjustedPercentage
daily$calproKg                  <- data$calproKg


# check consistency of factor variables in data
for(i in colnames(data)[sapply(data, is.factor)]) {
  print(i)
  assert_factor(data[[i]], empty.levels.ok=FALSE, min.levels=2)
}

# check consistency of factor variables in patient
for(i in colnames(patient)[sapply(patient, is.factor)]) {
  print(i)
  assert_factor(patient[[i]], empty.levels.ok=FALSE, min.levels=2)
}

## New definition of protein variable
data$proteinGproKG <- daily$proteinGproKG <-
  data$proteinIntakeAdjusted / data$Weight

# check variable
sum(data$proteinGproKG < 0)
sum(data$proteinGproKG > 7)
assert_numeric(data$proteinGproKG, lower=0, upper=7, any.missing=FALSE,
               finite=TRUE)

# categorical version of proteinGproKG
data$proteinCat <- factor(rep("mid", nrow(data)),
                          levels = c("lower", "mid", "upper"))
data$proteinCat[data$proteinGproKG < 0.8] <- "lower"
data$proteinCat[data$proteinGproKG >= 1.2] <- "upper"
data$proteinCat[is.na(data$proteinGproKG)] <- NA
data$proteinCat2 <- ifelse(data$proteinCat == "mid", 1, 0)
data$proteinCat3 <- ifelse(data$proteinCat == "upper", 1, 0)
data$protCat <- ifelse(data$proteinCat == "lower" & data$OralIntake == 1,
                       "mid", as.character(data$proteinCat))
data$protCat <- ifelse(data$proteinCat == "mid" & data$OralIntake == 1,
                       "upper", data$protCat)
data$protCat <- as.factor(data$protCat)
daily$proteinCat <- data$proteinCat
daily$protCat <- data$protCat
daily$proteinCat2 <- data$proteinCat2
daily$proteinCat3 <- data$proteinCat3
data$protCat2 <- ifelse(data$protCat == "mid", 1, 0)
data$protCat3 <- ifelse(data$protCat == "upper", 1, 0)
daily$protCat2 <- data$protCat2
daily$protCat3 <- data$protCat3

## Variables for descriptive analyes
# calories categorical before OI
data$caloriesCat <- factor("30 <= CA < 70",
                           levels = c("CA < 30", "30 <= CA < 70", "CA >= 70"))
data$caloriesCat[data$caloriesPercentage < 30] <- "CA < 30"
data$caloriesCat[data$caloriesPercentage >= 70] <- "CA >= 70"
data$caloriesCat[is.na(data$caloriesPercentage)] <- NA
daily$caloriesCat <- data$caloriesCat

### move one cat up if oral intake
daily$calCat <- NULL
daily$calCat <-
  ifelse(daily$OralIntake == 1 & daily$caloriesCat == "30 <= CA < 70",
         "CA >= 70", as.character(daily$caloriesCat))
daily$calCat <-
  ifelse(daily$OralIntake == 1 & daily$caloriesCat == "CA < 30",
         "30 <= CA < 70", as.character(daily$calCat))
daily$calCat <- as.factor(daily$calCat)
data$calCat <- daily$calCat

daily$calCat2 <- ifelse(daily$calCat == "30 <= CA < 70", 1, 0)
daily$calCat3 <- ifelse(daily$calCat == "CA >= 70", 1, 0)
daily$caloriesCat2 <- ifelse(daily$caloriesCat == "30 <= CA < 70", 1, 0)
daily$caloriesCat3 <- ifelse(daily$caloriesCat == "CA >= 70", 1, 0)

# add variables to daily data set / data
daily$inMV   <- data$inMV
data$calCat2 <- daily$calCat2
data$calCat3 <- daily$calCat3

data$DiedBeforeOrAtMax <- data$calendarEvent <= maxdays.nutri
data <- data %>%
  group_by(CombinedID) %>%
  mutate(maxday = max(Study_Day)) %>% ungroup
data$incomplete <- with(data,
                        (DiedBeforeOrAtMax & (maxday < calendarEvent) & (calendarDaysInICU <= maxday)) |
                          (!DiedBeforeOrAtMax & (maxday < maxdays.nutri) & (calendarDaysInICU <= maxday)))
daily$incomplete <- data$incomplete
daily$maxday <- data$maxday

##### final cleaning or trafos
removed_patients <-
  rbind(
    removed_patients,
    data.frame(id = patient$CombinedID[patient$ApacheIIScore > 71],
               cause = "inconsistency",
               cause2 = "AII Score > 71"))

patient <- patient[patient$ApacheIIScore <= 71, ] # last cleaning step (9 patients)
data <- data[data$CombinedID %in% patient$CombinedID, ]
daily <- daily[daily$CombinedID %in% patient$CombinedID, ]
nrow(patient)
length(unique(data$CombinedID))

reduced <-
  data[, c("CombinedID", "Study_Day", "Propofol",
           "inMV", "OralIntake", "PN", "EN")]
reduced$Propofol <- as.integer(as.character(reduced$Propofol))
reduced_14 <- reduced[reduced$Study_Day <= 4, ]
patient_14 <-
  data.frame(aggregate(reduced_14[c("Propofol", "inMV", "OralIntake", "PN", "EN")],
                       by = list(CombinedID = reduced_14$CombinedID), FUN = sum))

reduced_24 <- reduced[reduced$Study_Day >= 2 & reduced$Study_Day <= 4, ]
patient_24 <-
  data.frame(aggregate(reduced_24[c("Propofol", "inMV", "OralIntake", "PN", "EN")],
                       by = list(reduced_24$CombinedID), FUN = sum))
patient$Propofol1_4 <- patient_14$Propofol
patient$Propofol2_4 <- patient_24$Propofol
patient$inMV1_4 <- patient_14$inMV
patient$inMV2_4 <- patient_24$inMV
patient$OralIntake1_4 <- patient_14$OralIntake
patient$OralIntake2_4 <- patient_24$OralIntake
patient$PN1_4 <- patient_14$PN
patient$PN2_4 <- patient_24$PN
patient$EN1_4 <- patient_14$EN
patient$EN2_4 <- patient_24$EN
patient$icuByDummy <- 1
patient$AdmCatID <- as.factor(patient$AdmCatID)
patient$Gender <- as.factor(patient$Gender)
patient$AdCens <- as.factor(patient$AdCens)
patient$AdmCatID <- as.character(patient$AdmCatID)
patient$AdmCatID[patient$AdmCatID == "1"] <- "Medical"
patient$AdmCatID[patient$AdmCatID == "2"] <- "Surgical/Elective"
patient$AdmCatID[patient$AdmCatID == "3"] <- "Surgical/Emeregency"
patient$AdmCatID <- as.factor(patient$AdmCatID)
patient$DiagID2 <- relevel(patient$DiagID2, ref = "Other")
patient$Gender <- as.character(patient$Gender)
patient$Gender[patient$Gender == "1"] <- "Male"
patient$Gender[patient$Gender == "2"] <- "Female"
patient$Gender <- as.factor(patient$Gender)
patient$Year <- factor(as.character(patient$Year))
patient$Year <- relevel(patient$Year, ref = "2007")
patient$CombinedicuID <- as.factor(patient$CombinedicuID)

nrow(patient)
ids_no_PN2_4 <- patient$CombinedID[(is.na(patient$PN2_4))] # 0
patient <- filter(patient, !(is.na(PN2_4)))
ids_no_proteinCat <- unique(daily$CombinedID[is.na(daily$proteinCat)]) # 0
daily <- filter(daily, !(is.na(proteinCat)))
data <- data[data$event + 1 > data$Study_Day, ]
daily <- daily[daily$CombinedID %in% data$CombinedID, ]
patient <- patient[patient$CombinedID %in% data$CombinedID, ]
daily <- daily[daily$CombinedID %in% patient$CombinedID, ]
data <- data[data$CombinedID %in% patient$CombinedID, ]

# if no PN and EN but OI --> cat II for prot and cals,
# approx. 5 % of patient days!
daily$caloriesCat2[daily$PN == 0 & daily$EN == 0 & daily$OralIntake == 1] <- 1
daily$caloriesCat3[daily$PN == 0 & daily$EN == 0 & daily$OralIntake == 1] <- 0
daily$calCat2[daily$PN == 0 & daily$EN == 0 & daily$OralIntake == 1] <- 1
daily$calCat3[daily$PN == 0 & daily$EN == 0 & daily$OralIntake == 1] <- 0
daily$proteinCat2[daily$PN == 0 & daily$EN == 0 & daily$OralIntake == 1] <- 0
daily$proteinCat3[daily$PN == 0 & daily$EN == 0 & daily$OralIntake == 1] <- 0
daily$protCat2[daily$PN == 0 & daily$EN == 0 & daily$OralIntake == 1] <- 1
daily$protCat3[daily$PN == 0 & daily$EN == 0 & daily$OralIntake == 1] <- 0

daily_pre <- daily

# incomplete prototcols
# DOCU MISSING
daily$incomplete_day <- 0
library(tidyverse)
daily_l <- daily %>% group_by(CombinedID) %>% nest()
data_list <- daily_l$data
for (i in 1:length(data_list)) {
  nrow_d <- nrow(data_list[[i]])
  if (nrow_d != 12) {
    add_frame <- data_list[[i]][rep(nrow_d, 12 - nrow_d), ]
    add_frame[, c("calCat2", "calCat3", "protCat2", "protCat3",
                   "proteinCat2", "proteinCat3", "incomplete_day")] <-
      matrix(rep(c(0, 1, 1, 0, 1, 0, 1), 12 - nrow_d), nrow = 12 - nrow_d, byrow = TRUE)
    add_frame$Study_Day <- (nrow_d + 1):12
    data_list[[i]] <- rbind(data_list[[i]], add_frame)
  }
  if ((i %% 2500) == 0) print(paste(round((i / length(data_list)) * 100, 0), "%", " completed.", sep = ""))
}
daily_l$data <- data_list
daily <- daily_l %>% unnest(cols = c(data)) %>% ungroup()
table(daily$Study_Day)
# in very few cases: order disturbed
for (i in 2:nrow(daily)) {
  if ((daily$Study_Day[i] != (daily$Study_Day[i - 1] + 1)) & (daily$Study_Day[i] != 1)) {
    daily$Study_Day[i] <- daily$Study_Day[i] - 1
  }
}
table(daily$Study_Day)
daily$incomplete_day <- as.logical(daily$incomplete_day)

## we identified one special case which we manually change:
d <- daily %>% group_by(CombinedID) %>% summarise(m = mean(Study_Day))
problematic_ids <- unlist(d[d$m != 6.5, "CombinedID"])
daily[daily$CombinedID == problematic_ids[1], "Study_Day"] <- 1:12
table(daily$Study_Day)

## set parameters for exclusion
# parameters regarding time scale of survival (days after ICU admission)
mindays.surv  = 5 # minimal number of days after ICU admission survived
mindays.icu   = 4 # minimal number of days spent on ICU after ICU admission
# parameters regarding time scale of nutrition (calendar days/logged protocol days)
mindays.nutri = 4L # minimal number of calendar/protocol days with certain nutrition
# calendar/protocol days
## apply following exclusion criteria
# - survival time shorter than 4 days or discharged alive from ICU within 4 days
# - end of nutrition protocol before the end of calendar day 12 while patient
# still on ICU (we consider this to be missing data, the last nutrition days
# are hard to examine, though: Mostly missingness is due to discharge on the
# respective calendar day --> buffer.)
# - neither enteral nor parenteral nutrition during the first 4 days of protocol,
# but additional oral intake during first 4 days
# - no mechenical ventilation during first 4 days
# - patients that have not reached age of at least 18

ids_no_mechvent <- patient$CombinedID[patient$inMV2_4 == 0]
patient <- patient[!(patient$CombinedID %in% ids_no_mechvent), ]
ids_only_oi <- patient$CombinedID[(patient$EN2_4 == 0 & patient$PN2_4 == 0 &
                                     patient$OralIntake2_4 > 0)]
patient <- patient[!(patient$CombinedID %in% ids_only_oi), ]
ids_died_too_soon <- patient$CombinedID[patient$calendarEvent < 5]
ids_died_too_soon <- unique(c(ids_died_too_soon,
                              patient$CombinedID[patient$event <= 4]))
patient <- patient[!(patient$CombinedID %in% ids_died_too_soon), ]
ids_too_few_icudays <- patient$CombinedID[patient$DaysInICU < 4]
patient <- patient[!(patient$CombinedID %in% ids_too_few_icudays), ]
ids_too_young <- patient$CombinedID[patient$Age < 18]
patient <- patient[!(patient$CombinedID %in% ids_too_young), ]
ids_BMI_too_low <- patient$CombinedID[patient$BMI <= 13]
patient <- patient[!(patient$CombinedID %in% ids_BMI_too_low), ]
ids_2010 <- patient$CombinedID[patient$Year == 2010]
patient <- patient[!(patient$CombinedID %in% ids_2010), ]

daily <- daily[daily$CombinedID %in% patient$CombinedID, ]
ids_end_protocol_before_end_icu <-
  unique(daily$CombinedID[daily$calendarDaysInICU > 12 & daily$incomplete_day & daily$Study_Day < 12])
ids_end_protocol_before_end_icu <-
  unique(
    c(ids_end_protocol_before_end_icu,
      unique(daily$CombinedID[daily$incomplete_day & daily$Study_Day < 5])
    ))
patient <- patient[!(patient$CombinedID %in% ids_end_protocol_before_end_icu), ]

removed_patients <-
  rbind(
    removed_patients,
    data.frame(id = ids_end_protocol_before_end_icu,
               cause = "Consistency",
               cause2 = "Protocol too short."),
    data.frame(id = ids_no_mechvent,
               cause = "Exclusion",
               cause2 = "No mechanical ventilation in first days."),
    data.frame(id = ids_only_oi,
               cause = "Exclusion",
               cause2 = "Only oral intake in the first days (no art. nutri.)."),
    data.frame(id = ids_died_too_soon,
               cause = "Exclusion",
               cause2 = "Patient died before day 5"),
    data.frame(id = ids_too_few_icudays,
               cause = "Exclusion",
               cause2 = "Patient remained too briefly in ICU."),
    #data.frame(id = ids_too_young, ## is empty
    #           cause = "Exclusion",
    #           cause2 = "Patient too young"),
    data.frame(id = ids_BMI_too_low,
               cause = "Exclusion",
               cause2 = "BMI too low."),
    data.frame(id = ids_2010,
               cause = "Exclusion",
               cause2 = "Year (2010) excluded.")
    )

nrow(patient)
data <- data[data$CombinedID %in% patient$CombinedID, ]
daily <- daily[daily$CombinedID %in% patient$CombinedID, ]

# exclude mutually

patient <- patient[, c("CombinedID", "Year", "DiagID2", "AdmCatID", "Gender",
                       "ApacheIIScore", "BMI", "Propofol2_4", "inMV2_4", "PatientDied",
                       "PatientDischarged", "OralIntake2_4", "PN2_4", "Age", "EN2_4",
                       "Surv0To60", "Disc0To60", "surv_icu0to60", "surv_icu_status",
                       "event", "CombinedicuID",
                       "icuByDummy", "AdCens", "Calories", "Protein",
                       "DaysInICU", "DaysMechVent", "Weight")]
data <- merge(patient,
              daily[, c("CombinedID", "EN", "PN", "OralIntake", "calproKg",
                        "caloriesIntake", "caloriesPercentage", "proteinIntake",
                        "proteinGproKG", "Study_Day", "calCat2", "calCat3",
                        "caloriesPercentage", "proteinGproKG")],
              by = "CombinedID")

#saveRDS(data, file = "data/int/mergedAndCleanedData.Rds")

daily <- daily[, c("CombinedID", "Study_Day", "OralIntake", "incomplete",
                   "incomplete_day", "calproKg", "EN_Protein", "PN_Protein",
                   "calCat2", "calCat3", "caloriesCat2", "caloriesCat3",
                   "protCat2", "protCat3", "proteinCat2", "proteinCat3",
                   "caloriesPercentage", "proteinGproKG",
                   "proteinAdjustedPercentage",
                   "Propofol", "PN", "EN", "caloriesIntake")]

# Adapt to Study definition
daily$Study_Day <- daily$Study_Day - 1
daily <- daily[daily$Study_Day != 0, ]

data <- merge(patient, daily, by = "CombinedID")
data <- arrange(data, CombinedID, Year, Study_Day)
nrow(patient)

saveRDS(daily,   file = "data/daily.Rds")
saveRDS(patient, file = "data/patient.Rds")
saveRDS(ICU,     file = "data/ICU.Rds")
saveRDS(data,    file = "data/mergedAndCleanedData.Rds")

saveRDS(removed_patients, file = "data/excluded_ids.Rds")
write_csv(removed_patients, file = "data/excluded.csv")
