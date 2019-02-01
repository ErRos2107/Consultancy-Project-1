#library("xlsx")
library("dplyr")
library("readxl")
library("reshape2")
library("ggplot2")
library("mice")
library("Amelia")
library("lattice")
library("gdata")
library("lubridate")
library("pec")
library("rms")
library("randomForestSRC")
library("ggRandomForests")
library("party")
library("survival")
library("survminer")
library("prodlim")
library("stargazer")
library("caret")
#
setwd("/Users/eric/Library/Mobile Documents/com~apple~CloudDocs/MSc Statisitcs with Data Science - University of Edinburgh/Consultancy Style projects /Project 1 - Liver cancer /Data")

data(GBSG2)


GBSG2$status <- GBSG2$cens
GBSG2 <- GBSG2[order(GBSG2$time,-GBSG2$status),]
GBSG2$grade.bin <- as.factor(as.numeric(GBSG2$tgrade!="I"))

levels(GBSG2$grade.bin) <- c("I","II/III")

fitformGBSG <- Surv(time,status)~age+tsize+pnodes+progrec+estrec+grade.bin
fitcox <- cph(fitformGBSG,data=GBSG2,surv=TRUE,se.fit=FALSE)
set.seed(17)
fitrsf.GBSG2 <- rfsrc(fitformGBSG,data=GBSG2,forest=TRUE,ntree=100)

set.seed(17)
fitcforest <- pecCforest(fitformGBSG, data=GBSG2,
                            controls=cforest_classical(ntree=100))

extends <- function(...)TRUE
set.seed(2006)

fitpec <- pec(list("Cox"=fitcox,"rsf"=fitrsf,"cforest"=fitcforest),
              formula=Surv(time,cens)~age+tsize+grade.bin+pnodes+progrec+estrec, data=GBSG2,
              cens.model="cox", splitMethod="Boot632plus", maxtime=2000, B=5, keep.index=TRUE,
              keep.matrix=TRUE)


plot.variable(fitrsf, plots.per.page = 3)

plot.variable(fitrsf, plots.per.page = 2, xvar.names = c("tsize", "time", "age"))
plot.variable(fitrsf, surv.type = "surv", nvar = 1, time = 200)
plot.variable(fitrsf, surv.type = "surv", partial = TRUE, smooth.lines = TRUE)
plot.variable(fitrsf, surv.type = "rel.freq", partial = TRUE, nvar = 2)
## example of plot.variable calling a pre-processed plot.variable object
p.v <- plot.variable(fitrsf, surv.type = "surv", partial = TRUE, smooth.lines = TRUE)
plot.variable(p.v)
p.v$plots.per.page <- 1
p.v$smooth.lines <- FALSE
plot.variable(p.v)
## competing risks
data(follic, package = "randomForestSRC")
follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
plot.variable(follic.obj, target = 2)

plot.survival(fitrsf.GBSG2, subset = c(1,3,6)) # patient 1,3 and 6
plot.survival(fitrsf.GBSG2, subset = 3, k = 100)




cox2 <- coxph(Surv(time,cens)~tgrade+age+tsize+pnodes,data=GBSG2)
newdata <- data.frame(tgrade=c("I","II","III"),age=50,tsize=30,pnodes=8)
plotPredictSurvProb(cox2,
                    sort(unique(GBSG2$time)),
                    newdata=newdata,
                    col=c("darkgreen","darkorange","red"),
                    legend.title="Tumor grade",
                    legend.legend=c("I","II","III"))
mtext("Individualized survival curves from multiple Cox regression
      age=50, tumor size= 30, no. positive lymph nodes=8",line=1.5,cex=1.3)







##### Random forest survival model

df <- read_excel("Screening_cohort_clean.xlsx",sheet = 1,
                 col_types = c("text",rep("numeric",14)))
# convert screening cohort to wide format
df.cohort <- df %>% group_by(ID) %>% filter(row_number(ID) == 1) 
df.HCC <-read_excel("HCC_in_screening_clean.xlsx",sheet = "Long_format",
                    col_types = c("text","numeric","date",rep("numeric",10),"date",rep("numeric",3)))
df.HCC$`Age at sample` <- as.integer(df.HCC$`Age at sample`)
# adding censoring/event status column for both data sets
df.cohort$status <- rep(0,length(df.cohort$ID))
df.HCC$status <- rep(1,length(df.HCC$ID))

#### prepare data for randomForestSRC package
df.temp <- df.HCC%>% select(ID,Sex,`Time relative`,`Age at sample`,ALD:Other,status) %>%
  filter_all(all_vars(.!='NA'))
# remove left censored data
df.HCC.rf <- df.temp%>% select(ID,Sex,`Time relative`,`Age at sample`,ALD:Other,status)%>%
  filter_at(vars("Time relative"), any_vars(.<0))

df.temp.screening <- df.cohort %>% select(ID,Sex,"Time relative","Age  at sample",ALD:Other,status) %>%
  filter_all(all_vars(.!='NA'))
df.cohort.rf <- df.temp.screening%>% select(ID,Sex,"Time relative","Age  at sample",ALD:Other,status)%>%
  filter_at(vars("Time relative"), any_vars(.<0))

colnames(df.HCC.rf)[colnames(df.HCC.rf) == 'Age at sample'] <- 'Age'
colnames(df.cohort.rf)[colnames(df.cohort.rf) == "Age  at sample"] <- 'Age'

full.data <- bind_rows(df.HCC.rf,df.cohort.rf)
colnames(full.data)[colnames(full.data) == "Time relative"] <- 'time'

full.data <- full.data[order(full.data$time,-full.data$status),]

full.data$time <- as.numeric(-(full.data$time))

#### Estimating the median follow-up time
quantile(prodlim(Hist(time,status)~1,data=full.data,reverse=TRUE))
#The median potential follow-up time of the HCC study was 1645 days 
#(IQR: [1100 days;1714 days]). This means that 50% of the patients 
#would have been observed for at least 1645 days had there been no events.


#### Kaplan-Meier graph
km0 <- prodlim(Hist(time,status)~1,data=full.data)
par(mar=c(7,7,5,5), # margin of figure
    mgp=c(4,1,0))
plot(km0,
     xlab="Years",  # label for x-axis
     axis1.at=c(0,1461,2922,4383,5844), # time grid for x-axis
     axis1.labels=c(0,4,8,12,14), # time labels for x-axis
     axis2.las=2, # rotate labels of y-axis
     atrisk.dist=2) # adjust numbers below the figur) 
#Published table
publish(km0,times=seq(0,2900,365.25),org=TRUE)


#### Fitting the survival function
fitformHCC <- Surv(time,event=status)~Age+Sex+ALD+NAFLD+HepB+HepC+Other
fit <- surv_fit(fitformHCC,data = full.data)
summary(fit)


#### Fitting Cox regression
fitcox <- cph(fitformHCC,data=full.data,surv=TRUE,se.fit=FALSE)
dd <- datadist(full.data)
options(datadist="dd")
factor.effect.summary <- summary(fitcox)
stargazer(factor.effect.summary)

# # Fit a stratified model, clustered on patients 
# full.data.age <- full.data[full.data$Age<60,]
# fitcox.age <- coxph(Surv(time, status) ~ (Sex+ALD+NAFLD+HepB+HepC+Other) * strata(Age) + 
#         cluster(ID), full.data.age)
# summary(fitcox.age)

# predict the 10 year survival probabilities for Cox model
pcox <- predictSurvProb(fitcox,newdata=newData,times=10*365.25)


############################ Fit randomForest #############################
#Data Partitionning 
train.idx <- createDataPartition(full.data$status, p=0.8)$Resample1
train.data <- full.data[train.idx,]
test.data <- full.data[-train.idx,]

set.seed(17)
fitrsf <- rfsrc(fitformHCC,data =train.data ,forest=TRUE,ntree=1000,mtry = 3,
                tree.err = TRUE,importance=TRUE)
summary(fitrsf)
print.rfsrc(fitrsf)
# predict the 10 year survival probabilities for RF model on test set
prsf <-predict.rfsrc(fitrsf,newdata=test.data,times=10*365.25,importance = TRUE)

#generalization error estimate
plot(gg_error(fitrsf))+ggtitle("Generalization error estimate")

#predicted survival from our RSF model for training data
ggRFsrc <- plot(gg_rfsrc(fitrsf), alpha = 0.2) +
   #theme(legend.position = "right") +
  labs(y = "Survival Probability", x = "Time (days)") +
  coord_cartesian(ylim = c(-0.01, 1.01))+
  scale_shape_discrete(name="Status",breaks=c(0, 1),
  labels=c("Dead", "Censored"))+
  ggtitle("Survival curve for every patient")
  
show(ggRFsrc)

# plot survival curves for first 10 individuals: direct way
matplot(fitrsf$time.interest, 100 * t(fitrsf$survival[1:10, ]),
        xlab = "Time", ylab = "Survival", type = "l", lty = 1)

# predicted survival from our RSF model for test data
plot(gg_rfsrc(prsf), alpha=.2) +
     theme(legend.position = "right") 
     labs(y = "Survival Probability", x = "Time (days)") +
     coord_cartesian(ylim = c(-0.01, 1.01))

#compares the predicted survival between male and female
plot(gg_rfsrc(fitrsf, by = "Sex")) +
   theme(legend.position = c(0.2, 0.2)) +
  labs(y = "Survival Probability", x = "Time (days)") +
  coord_cartesian(ylim = c(-0.01, 1.01))


##### Variable importance (VIMP)


# Extracts variable importance (Mean Decrease in Gini Index)
# Sorts by variable importance and relevels factors to match ordering
var_importance <- data_frame(variable=setdiff(colnames(full.data), c("status","time","ID")),
                             importance=as.vector(vimp(fitrsf)$importance))
var_importance <- arrange(var_importance, desc(importance))
var_importance$variable <- factor(var_importance$variable, levels=var_importance$variable)

p <- ggplot(var_importance, aes(x=variable, weight=importance))
p <- p + geom_bar() + ggtitle("Variable Importance from Random Forest Fit")
p <- p + xlab("Covariates") + ylab("Variable Importance (Mean Decrease in Gini Index)")
p + theme(axis.text.x=element_blank(),
          axis.text.y=element_text(size=12),
          axis.title=element_text(size=12),
          plot.title=element_text(size=12),
          legend.title=element_text(size=16),
          legend.text=element_text(size=12)) +
          coord_flip()
  
### Minimum depth
varselect.fitrst <- var.select(fitrsf)
gg_md <- gg_minimal_depth(varselect.fitrst, lbls = st.labs)
print(gg_md)

plot(gg_minimal_vimp(gg_md)) +
   theme(legend.position=c(0.8, 0.2)) +
   ggtitle("Variable selection comparison")


# Variable dependence
gg_v <- gg_variable(fitrsf, time = c(1, 5),
                     time.labels = c("1 Year", "5 Years"))
plot(gg_v, xvar = "Age", alpha = 0.4) + 
   theme(legend.position = "right") +
   coord_cartesian(ylim = c(-0.01, 1.01))

xvar <- c("Sex", "Age", "ALD", "NAFLD", "HepB","HepC")

   plot(gg_v, xvar = xvar[-1], panel = TRUE, alpha = 0.4) +
   labs(y = "Survival") +
   theme(legend.position = "right") +
   coord_cartesian(ylim = c(-0.05, 1.05))
#Loess smooth curve with
#shaded 95% confidence band indicates decreasing survival with increasing age.

##calculating all pair- wise minimal depth interactions
ggint <- gg_interaction(fitrsf)
plot(ggint)
 
#Variable dependence coplot of survival

# Get variable dependence at 1 year
ggvar <- gg_variable(fitrsf, time = 1)

# For labeling coplot membership
ggvar$ALD <- paste("ALD = ", ggvar$ALD, sep = "")

# Plot with linear smooth (method argument)
var_dep <- plot(ggvar, xvar = "Age", alpha = 0.5) +
     #  geom_smooth(method = "glm",se = FALSE) +
     labs(y = "Survival", x = "Age") +
     theme(legend.position = "right") +
     #scale_color_manual(values = strCol, labels = event.labels) +
     #scale_shape_manual(values = event.marks, labels = event.labels) +
     coord_cartesian(y = c(-.01,1.01))+

var_dep + facet_grid(~ALD)
######## Fit conditional inference forest 
set.seed(17)
fitcforest <- pecCforest(fitformHCC, data=full.data,
                         controls=cforest_classical(ntree=1000))
summary(fitcforest)

#predict the 10 year survival probabilities for conditional RF model
extends <- function(...)TRUE
pcf <- predictSurvProb(fitcforest,newdata=test.set,times=10*365.25)
#set.seed(2006)

# prediction error curve using boot362 method to compare pec between RF and Cox model
fitpec <- pec(list("Cox"=fitcox,"rsf"=fitrsf,"cforest"=fitcforest),
              formula=fitformHCC, data=full.data,
              cens.model="cox", splitMethod="Boot632plus", maxtime=2000, B=5, keep.index=TRUE,
              keep.matrix=TRUE)
summary(fitpec)
crps(fitpec,times=1000) # between min time and 1000 days
plot(fitpec, predErr="Boot632plusErr",xlim=c(0,10*365.25), 
     axis1.at=seq(0,10*365.25, 2*365.25), axis1.label=seq(0,10,2))
# Survival curves

ggsurvplot(fit, data = full.data,
           title = "Survival Curves",    # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines  # Change legend label               # Use JCO journal color palette
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE               # Hide tables y axis text
)

# Distribution of events'time
surv <- Surv(full.data$time, full.data$status)
ggsurvevents(surv)

# 
fit.test1 <- survfit(Surv(time,status)~Sex,data = full.data)
fit.test2 <- survfit(Surv(time,status)~1,data = full.data)

ggsurvplot(fit.test1, data = full.data,
           title = "Survival Curves",
           pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines              
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE,
           ggtheme = theme_grey()              # Hide tables y axis text
)
# Plot cumulative events
ggsurvplot(fit.test1,
           risk.table = TRUE, risk.table.col = "strata",
           fun = "event")
# Adjusted survival curves for cox proportional hazards
ggcoxadjustedcurves(fitcox, data = full.data,
                   individual.curves = TRUE)

#  
my_variables <- c("Sex", "ALD", "NAFLD","HepB","HepC","Other")
my_formulas <- list()
for(variable in my_variables){
  my_formulas[[variable]] <- paste0("Surv(time, status) ~ ", variable) %>%
    as.formula()
}
my_formulas
my_fits <- surv_fit(my_formulas, data = full.data)
ggsurvs <- ggsurvplot(my_fits, risk.table = TRUE, pval = TRUE)
names(ggsurvs)

# Print survival curves for all variables at once
ggsurvs
#
ggsurvplot_group_by(fit, data=full.data, group.by=c("HepB","ALD"))

# Print survival curves for the variable Sex
ggsurvs$`full.data::Other`
#
ggforest(fitcox)


