
# load libraries
rm( list = ls ())
if (! require (readxl)) install . packages ("readxl"); require (readxl)
if (! require ( tidyverse )) install . packages (’tidyverse ’); require (tidyverse)
if (! require ( reshape2 )) install . packages ("reshape2"); require (reshape2)
if (! require (rstan)) install . packages ("rstan"); require (rstan)
if (! require ( ggmcmc )) install . packages (ggmcmc ’); require (ggmcmc)
if (! require ( ggfortify )) install . packages ("ggfortify"); require (ggfortify)
if (! require (JM)) install . packages ("JM"); require (JM)
if (! require ( JMbayes )) install . packages (JMbayes ’); require (JMbayes)
if (! require ( lubridate )) install . packages ("lubridate"); require (lubridate)
if (! require ( survival )) install . packages ("survival"); require (survival)
setwd("/Users/eric/Library/Mobile Documents/com~apple~CloudDocs/MSc Statisitcs with Data Science - University of Edinburgh/Consultancy Style projects /Project 1 - Liver cancer /Data")

df.JM <- read_excel("Screening_cohort_clean.xlsx",sheet = 1,
                 col_types = c("text",rep("numeric",13)))
df.HCC.JM <-read_excel("Long_HCC_clean_final.xlsx",sheet = 1,
                                 col_types = c("text",rep("numeric",12)))
full.data.long <- read_excel("Full_data_final.xlsx",sheet = 1,
                             col_types = c("text",rep("numeric",11)))

#------------------------------ preparing data for Joint Model package -----------------------------

# set time at sample
df.JM$Obstime <- df.JM$Sample_time-df.JM$Time

# Add censoring status
df.JM$Status <- rep(0,length(df.JM$ID))

# Transform AFP values to log(AFP+0.01)
df.JM$AFP <- log(df.JM$AFP+0.01,base = 10)
df.HCC.JM$AFP <- log(df.HCC.JM$AFP+0.01,base = 10)

# Select important covariates
data.JM.cohort <- df.JM %>% dplyr::select(ID:Age,Time,ALD:Other,AFP,Obstime,Status)
data.JM.HCC <- df.HCC.JM %>% dplyr::select(ID:Age,TimeEvent,ALD:Other,AFP,Obstime,Status)

# Combine the two datasets
colnames(data.JM.cohort)[colnames(data.JM.cohort) == 'Time'] <- 'TimeEvent'
data.JM.cohort$TimeEvent <- as.numeric(-data.JM.cohort$TimeEvent)
full.data.long.test <- bind_rows(data.JM.HCC,data.JM.cohort)

# Remove inconsitencies in data reporting (patient missing, wrond date of entry etc...)
full.data.long.test <- full.data.long.test %>% dplyr::select(ID:Age,TimeEvent,ALD:Other,AFP,Obstime,Status)%>%
  filter_at(vars(Obstime), any_vars(.>=0))
#--------------------------------------------------------------------------------------------------
#------------------------------ preparing data for Survival Model ---------------------------------
# convert screening cohort to wide format
full.data.wide.test <- full.data.long.test %>% group_by(ID) %>% filter(row_number(ID) == 1)

data.long.test <- full.data.long.test
data.long.test$TimeEvent <- data.long.test$TimeEvent+1

data.wide.test <-data.long.test %>% group_by(ID) %>% filter(row_number(ID) == 1)

stopifnot(data.long.test$TimeEvent>data.long.test$Obstime)
stopifnot(data.long.test$TimeEvent>0)
stopifnot(data.long.test$Obstime>=0)

### sorting the data for the Joint Model
sortnames <- c("ID","TimeEvent")
DATA.long <- full.data.long[do.call("order", full.data.long[sortnames]),]
DATA.wide <- DATA.long[!duplicated(DATA.long$ID),]

#------------------------------ Linear Mixed Model ------------------------------------------------

ctrl <- lmeControl(opt='optim')
lmeFit <- lme(AFP ~ Obstime,
                   random = ~ Obstime | ID,data = data.long.test)
summary(lmeFit)

# or using spline

lmeFit2 <- lme(AFP ~ ns(Obstime, 2), data = full.data.long,
              random = ~ ns(Obstime, 2) | ID)
#------------------------------ Cox model ---------------------------------------------------------
coxFit <- coxph(Surv(TimeEvent, Status) ~ Sex+Age,
                     data = data.wide.test, x = TRUE)
summary(coxFit)

#or
coxFit2 <- coxph(Surv(TimeEvent, Status) ~ Sex +Age,
                data = full.data.wide, x = TRUE)
#------------------------------ Joint Model -------------------------------------------------------
jointFit <- jointModel(lmeFit, coxFit,
                            timeVar = "Obstime", method = "piecewise-PH-GH")

summary(jointFit)
#### --------------------------- Bayesian Approach to Joint Model ------------------------

jointFitBayes <- jointModelBayes(lmeFit, coxFit, timeVar = "Obstime",baseHaz ="P-splines")

## test on smaller data set
lmeFit.test <- lme(AFP ~ Obstime,
              random = ~ Obstime | ID,data = data.JM.HCC)
coxfit.test <- coxph(Surv(TimeEvent, Status) ~ Sex+Age+ALD,
                     data = aaa, x = TRUE)
jointFitBayes.test <- jointModelBayes(lmeFit.test, coxfit.test, timeVar = "Obstime")
summary(jointFitBayes.test)

lmeFit.test2 <- lme(AFP ~ Obstime,
                    random = ~ Obstime | ID,data = DATA.long)

coxfit.test2 <- coxph(Surv(TimeEvent, Status) ~ Sex + Age,
                     control=ctrl.test,data = DATA.wide, x = TRUE)
jointFitBayes.test2 <- jointModelBayes(lmeFit.test2, coxfit.test2, timeVar = "Obstime")
summary(jointFitBayes.test2)

inference.length <- 80000
M1.cov.results <- coda.samples(jointFitBayes.test$mcmc)

## GGplot graphics
plot.JMbayes(jointFitBayes.test2)
plot.survfit.JMbayes(jointFitBayes.test2)

M1.cov.graphics <- ggs(jointFitBayes.test)
ggs_density(jointFitBayes.test$priors$priorTau.betas)
ggs_traceplot(jointFitBayes.test)

#plots of conditional probabilities of survival.
ND <- DATA.long[DATA.wide$ID == "A13", ] # the data of Subject 2
survPreds <- vector("list", nrow(ND))
for (i in 1:nrow(ND)) {
  survPreds[[i]] <- survfitJM(jointFitBayes.test2, newdata = ND[1:i, ],idVar = "ID")
}
par(mfrow = c(2, 2), oma = c(0, 2, 0, 2))
for (i in c(1,3,5,7)) {
  plot(survPreds[[i]], estimator = "median", conf.int = TRUE,
       include.y = TRUE, main = paste("Follow-up time:",
                                      round(survPreds[[i]]$last.time, 1)), ylab = "", ylab2 = "")
}
