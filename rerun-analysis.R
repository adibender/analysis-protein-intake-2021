# Preprocessing
source("preproc/preproc-data.R", echo = TRUE)

# Fit the models
source("scripts/models-main.R", echo = TRUE)
source("scripts/models-icu.R", echo = TRUE)
source("scripts/models-static-lag-lead.R", echo = TRUE)
source("scripts/models-bmi-subgroup.R", echo = TRUE)

# create publication figures
source("figures-manuscript.R", echo = TRUE)

# create supplement figures
source("figures-supplement.R", echo = TRUE)

# numbers, tables and summaries
source("numbers-manuscript.R", echo = TRUE)
source("numbers-supplement.R", echo = TRUE)
