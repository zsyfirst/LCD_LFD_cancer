library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(msm)
library(purrr)
library(CMAverse)

setwd("lowcarb_cancer/mediation")
pool <- fread("lowcarb_cancer/data/pool_allcancer.csv")

exp = c("lcssvq","alcssvq","vlcssvq","hlcssvq", "ulcssvq",
        "lfssvq","alfssvq","vlfssvq","hlfssvq", "ulfssvq")
#define mediators (time-varying BMI)
med <- c("bmicon") 
#define covariates
covar <- c("agecon","cohort","intervalp","calorsvq", 
           "smkstatus", "alc", "actc",  "bmib3cat", 
           "white", "mv", "asp", "cafh", 
           "pmhever", "mam", "endo2c") 

outcome <- "tobeca"
event <- "can_obese"

## keep only the variables we need and drop rows with missing
vars_needed <- c(outcome, event, exp, med, covar)

pool <- pool %>%
  dplyr::select(all_of(vars_needed)) %>%
  dplyr::filter(complete.cases(.))

# make sure event is 0/1 numeric
pool[[event]] <- as.integer(pool[[event]])

## Controlled mediator values (median)
metab_dataset <- pool %>% dplyr::select(all_of(med))

controlled_med <- metab_dataset %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable) %>%
  summarise(median = median(value, na.rm = TRUE), .groups = "drop")

mval_bmi <- controlled_med$median[controlled_med$variable == "bmicon"]
# ------------------------------------------------------------
#    Run loops for weight-change mediator
#    (diet score -> bmi -> can_obese)
# ------------------------------------------------------------
run_cmest <- function(data, exposure, mediator, mval,
                      astar = 1,  # reference level (e.g., Q1)
                      a     = 5)  # exposed level (e.g., Q5)
{
  cmest(data      = data,
        model     = "rb",
        outcome   = outcome,
        event     = event,
        exposure  = exposure, mediator  = mediator, basec     = covar,
        EMint     = FALSE,
        mreg      = list("linear"),   # linear model for continuous mediator
        yreg      = "coxph",
        astar     = astar,
        a         = a,
        mval      = list(mval),
        estimation = "imputation",
        inference  = "bootstrap",
        # nboot      = 1000, # we use bootstrap n = 1000 for the results
        nboot      = 3,
        boot.ci.type="bca")  # notes: 1 boot needs almost 1 hour
}

cmest_wt_list <- lapply(exp, function(e) {
  fit <- run_cmest(data = pool, exposure = e,
                   mediator = "bmicon", mval = mval_bmi,
                   astar = 1,   # Q1
                   a     = 5)   # Q5
  list(exposure = e, fit = fit)})

# extract summary tables and combine
res_bmi <- lapply(cmest_wt_list, function(x) {
  s <- as.data.frame(as.matrix(summary(x$fit)$summarydf))
  s$exposure <- x$exposure
  s$mediator <- "BMI"
  s
}) %>% bind_rows()

res_bmi$Outcome <- "Obesity-related cancer"

write.csv(res_bmi,
          file = "lowcarb_cancer/mediation/BMI_Med_allscores.csv",
          row.names = FALSE)

