library("survival")
library("survminer")
library("readxl")
library("car")
library("stargazer")
library("ggplot2")
library("grid")
require(gridExtra)

setwd("/Users/eric/Library/Mobile Documents/com~apple~CloudDocs/MSc Statisitcs with Data Science - University of Edinburgh/Consultancy Style projects /Project 1 - Liver cancer /Data")
df.long <- read_excel("Full_data_final.xlsx",sheet = 1,
                             col_types = c("text",rep("numeric",11)))
df.wide <- df.long[!duplicated(df.long$ID),]
## -------------------- Unadjusted Cox regression models -----------------------

# Univariate Cox regression

#categorize age
df.wide$age.cat <- recode(df.wide$Age, " lo:53=0; 53:hi=1 ")


covariates <- c("Sex","age.cat","ALD", "NAFLD", "HepB","HepC","Other")
univ.formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(TimeEvent, Status)~', x)))

univ.models <- lapply( univ.formulas, function(x){coxph(x, data = df.wide)})

# Extract summary for each univariate model
univ.results <- lapply(univ.models,
                       function(x){ c
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ.results, check.names = FALSE))
univ.res <- as.data.frame(res)
stargazer(univ.res,summary = FALSE)

## ---------------------- Multivariate Cox regression model ----------------------------

res.cox <- coxph(Surv(TimeEvent, Status) ~ Sex+age.cat+Other+HepC+ALD+NAFLD+HepB, data =df.wide)
sum.res.cox <- summary(res.cox)
sum.res.cox$coefficients
test.ph <- cox.zph(res.cox, transform="km", global=TRUE)
ggcoxzph(test.ph)

fit1 <- surv_fit(Surv(TimeEvent, Status) ~ Sex, data = df.wide)
surv_median(fit1)

# Diagnostic plots

par(mfrow=c(2, 2))
plot(cox.zph(res.cox))
#Systematic departures from a horizontal line are indicative of
#non-proportional hazards

# Testing influential observations

ggcoxdiagnostics(res.cox, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

# Forest plot
ggforest(res.cox)
# --------------------------- Optional ------------------------------------------------

#categorize age (0-50 and 50+ )
test.cat <- df.wide
test.cat$age.cat <- recode(df.wide$Age, " lo:50=1; 50:hi=2 ")

res.cox2 <- coxph(Surv(TimeEvent, Status) ~ Sex+age.cat+ALD+NAFLD+Other+HepB, data =test.cat)
summary(res.cox2)
cox.zph(res.cox2)

#need to stratify by age
res.cox3 <- coxph(Surv(TimeEvent, Status) ~ Sex+strata(age.cat)+ALD+NAFLD+Other+HepB, data =test.cat)
summary(res.cox3)
cox.zph(res.cox3)
#--------------------------------------------------------------------------------------


#Estimated Survival Function Distribution
ggsurvplot(survfit(res.cox), color = "#2E9FDF",
           ggtheme = theme_minimal(),data = df.wide)

#assess the impact of the sex on the estimated survival probability

  # Create the new data
aetiologies.df <- with(df.wide, data.frame(Sex = c(1, 1,1,1,1),age.cat=c(1,1,1,1,1),
                        ALD = c(1,0,0,0,0),NAFLD = c(0,1,0,0,0),HepB = c(0,0,1,0,0),
                        HepC=c(0,0,0,1,0),Other = c(0,0,0,0,1)))

# predicted survival curves
fit <- survfit(res.cox, newdata = aetiologies.df)
ggsurv <- ggsurvplot(fit,data = df.wide, conf.int = FALSE,legend.labs=c("ALD=1", "NAFLD=2",
                                                             "HepB=3","Hepc=4","Other=5"),
           ggtheme = theme_minimal()) +ggtitle("Pre-conditional disease impact")

ggsurv$plot+facet_wrap(~ Sex)

test.fit <- survfit(res.cox, newdata = df.wide)
ggsurvplot(test.fit, data = df.wide,
           palette = "jco", pval = TRUE,conf.int = FALSE,legend.labs=c("ALD=1", "NAFLD=2",
                                                                       "HepB=3","Hepc=4","Other=5"))



test.fit <- survfit(res.cox, newdata = aetiologies.df)
ggsurv <- ggsurvplot(test.fit, data = df.wide,
                     palette = "jco",conf.int = FALSE,legend.labs=c("ALD=1", "NAFLD=2",
                      "HepB=3","Hepc=4","Other=5"),
                     ggtheme = theme_minimal())

a <- ggsurv$plot+ggtitle("Aetiology effect on HCC diagnosis")

aetiologies.df2 <- with(df.wide, data.frame(Sex = c(0,0,0,0,0),age.cat=c(1,1,1,1,1),
                                            ALD = c(1,0,0,0,0),NAFLD = c(0,1,0,0,0),HepB = c(0,0,1,0,0),
                                            HepC=c(0,0,0,1,0),Other = c(0,0,0,0,1)))



testfit.female <- survfit(res.cox, newdata = aetiologies.df2)
ggsurv.female <- ggsurvplot(testfit.female, data = df.wide,
                     palette = "jco",conf.int = FALSE,legend.labs=c("ALD=1", "NAFLD=2",
                                                                    "HepB=3","Hepc=4","Other=5"),
                     ggtheme = theme_minimal())
b <- ggsurv.female$plot+ggtitle("Aetiology effect on HCC diagnosis-Male")

# Plot survival prob for males vs females
grid.arrange(a, b, nrow = 2)

# Same plot with sharing legend
grid.arrange.shared.legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}
grid_arrange_shared_legend(a, b)
