library(dplyr)

patient           <- readRDS("data/patient.Rds")
daily             <- readRDS("data/daily.Rds")
ped_icu_death     <- readRDS("data/ped-data-death-sub5.Rds")
ped_icu_discharge <- readRDS("data/ped-data-discharge-sub5.Rds")
ped_death         <- readRDS("data/ped-data-death-hosp-sub5.Rds")
ped_discharge     <- readRDS("data/ped-data-discharge-hosp-sub5.Rds")

# uid <- unique(ped$CombinedID)
uid <- unique(ped_death$CombinedID)
# number of subjects in final model
length(unique(uid))
# number discharged before day 60
ped_discharge %>% pull(ped_status) %>% sum()
ped_death %>% pull(ped_status) %>% sum()


# Year
table(patient$Year)
round((table(patient$Year) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$Year)))
names(res) <- unique(patient$Year)
i <- 0
for (y in unique(patient$Year)) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$Year == y])
}
round(res * 100, 1)

# Gender
table(patient$Gender)
round((table(patient$Gender) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$Gender)))
names(res) <- unique(patient$Gender)
i <- 0
for (y in unique(patient$Gender)) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$Gender == y])
}
round(res * 100, 1)

# MV
MV <- patient$DaysMechVent * 24
MV <- cut(MV , c(0, 24, 48, max(MV)), right = TRUE, include.lowest = TRUE)
table(MV)
round((table(MV) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(MV)))
names(res) <- unique(MV)
i <- 0
for (y in unique(MV)) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[MV == y])
}
round(res * 100, 1)

#OI
table(patient$OralIntake2_4)
round((table(patient$OralIntake2_4) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$OralIntake2_4)))
names(res) <- sort(unique(patient$OralIntake2_4))
i <- 0
for (y in sort(unique(patient$OralIntake2_4))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$OralIntake2_4 == y])
}
round(res * 100, 1)

#PF
table(patient$Propofol2_4)
round((table(patient$Propofol2_4) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$Propofol2_4)))
names(res) <- sort(unique(patient$Propofol2_4))
i <- 0
for (y in sort(unique(patient$Propofol2_4))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$Propofol2_4 == y])
}
round(res * 100, 1)

#PN
table(patient$PN2_4)
round((table(patient$PN2_4) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$PN2_4)))
names(res) <- sort(unique(patient$PN2_4))
i <- 0
for (y in sort(unique(patient$PN2_4))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$PN2_4 == y])
}
round(res * 100, 1)

#EN # empty
table(patient$EN2_4)
round((table(patient$EN2_4) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$EN2_4)))
names(res) <- sort(unique(patient$EN2_4))
i <- 0
for (y in sort(unique(patient$EN2_4))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$EN2_4 == y])
}
round(res * 100, 1)

#Admin cat
table(patient$AdmCatID)
round((table(patient$AdmCatID) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$AdmCatID)))
names(res) <- names(table(patient$AdmCatID))
i <- 0
for (y in unique(patient$AdmCatID)) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$AdmCatID == y])
}
round(res * 100, 1)

# Diag
table(patient$DiagID2)
round((table(patient$DiagID2) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$DiagID2)))
names(res) <- names(table(patient$DiagID2))
i <- 0
for (y in names(table(patient$DiagID2))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$DiagID2 == y])
}
round(res * 100, 1)

# summaries

summary(patient[, c("Age", "ApacheIIScore", "BMI")])
