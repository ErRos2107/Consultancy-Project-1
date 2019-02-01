library('PEB')
library("dplyr")
library("readxl")
library("reshape2")
library("ggplot2")
library("broom")
beta <- multilevel::ICC1(aov(df$`Result AFP`~df$ID)) #0.6673252
beta.hcc <- multilevel::ICC1(aov(df.HCC$`AFP value`~df.HCC$ID)) #0.672548
setwd("/Users/eric/Library/Mobile Documents/com~apple~CloudDocs/MSc Statisitcs with Data Science - University of Edinburgh/Consultancy Style projects /Project 1 - Liver cancer /Data")

df <- read_excel("Screening_cohort_clean.xlsx",sheet = 1,
                 col_types = c("text",rep("numeric",14)))

# checking condition for PEB
# 1) Normality of between patient AFP measurements
# 2) Normality within patient 
# 3) constant variance
#1)
Check1 <- df %>% group_by(ID) %>% 
  summarise(mean=mean(`Result AFP`),variance=var(`Result AFP`))

check2 <- as.data.frame(df %>% group_by(ID) %>% filter(row_number(`Result AFP`)==1)) 
colnames(check2)[colnames(check2) == "Result AFP"] <- 'logAFP'
ggplot(check2,aes(logAFP)) + geom_density(stat = "density",adjust=3) +
  ggtitle("Distribution of first AFP measure per patient")
# or 
check3 <- as.data.frame(df)
colnames(check3)[colnames(check3) == "Result AFP"] <- 'logAFP'
ggplot(check3,aes(logAFP)) + geom_density(stat = "density",adjust=3) +
  ggtitle("Distribution of AFP levels")
#2)
check4 <- df %>% group_by(ID) %>% filter(row_number(`Result AFP`)>=45)
unique(check4$ID)

df1<-data.frame(logAFP=df$`Result AFP` [df$ID == 16])
df2<-data.frame(logAFP=df$`Result AFP` [df$ID == 237])
df3<-data.frame(logAFP=df$`Result AFP` [df$ID == 1266])
df4<-data.frame(logAFP=df$`Result AFP` [df$ID == 1363])

ggplot(df1,aes(logAFP))+geom_density(aes(color="ID 16"),adjust=3)+
  geom_density(data=df2,aes(color="ID 237"),adjust=3)+
  geom_density(data=df3,aes(color="ID 1266"),adjust=3)+
  geom_density(data=df4,aes(color="ID 1363"),adjust=3)+
  labs(color="Density curve")+
  ggtitle("Within patient AFP levels")



df$ID <- as.factor(df$ID)
df$`Result AFP` <- log(df$`Result AFP`+0.01,base = 10)
AFP.per.patients <- df %>% group_by(ID,"Result AFP") %>% summarize(count=n())

params <- PEBparams(x=df$`Result AFP`,
                    id=df$ID)
# Accuracy of specificity level
params_iterative <- PEBparams(x=df$`Result AFP`,
            id=df$ID,
            method='iterative',iterations = 2)


(y_past = df$`Result AFP` [df$ID == 9]) # let's use person '9' as an example
qpeb(p = 0.95,
     params_iterative,
     n=length(y_past),
     ybar=mean(y_past),
     conf.level=0.9)

#  and a personalized threshold with 95% specificity can be set via:
qpeb(p=0.95, # specificity level
     n=length(y_past),
     ybar=mean(y_past),
     params) #confidence interval for the false positive rate is fairly narrow (good)
#> .0.7948526
threshold <- NULL
#personalised threshold for every controls:
for (i in unique(df$ID)){
  y_past=df$`Result AFP`[df$ID==i]
  threshold[i] <- qpeb(p=0.5,n=length(y_past),
       ybar=mean(y_past),
       params)
}
threshold <- NULL
for (i in 1:AFP.per.patients$count[3]){
  y_past <- df$`Result AFP`[df$ID==100]
  threshold[i] <- qpeb(p=0.95,n=length(y_past[1:i]),
                       ybar=mean(y_past[1:i]),
                       params)
}

# check sensitivity of PEB for last AFP measurement for each patient 
Last.AFP.measurement <- df %>%
  group_by(ID) %>%
  slice(n()) %>%
  ungroup() %>%
  select("ID","Result AFP")
sum((Last.AFP.measurement$`Result AFP`>threshold)+0)

####### Using the HCC data set ####
df.HCC <-read_excel("HCC_in_screening_clean.xlsx",sheet = "Sheet1",
                    col_types = c("text","numeric","date",rep("numeric",8),"date",rep("numeric",2)))

df.HCC$ID <- as.factor(df.HCC$ID)
df.HCC$`AFP value` <- log(df.HCC$`AFP value`+0.01,base = 10)
df.HCC <-df.HCC %>% filter_at(vars("AFP value"), any_vars(.!='NA'))
AFP.per.patients.HCC <- df.HCC %>% group_by(ID,"AFP value") %>% summarise(count=n())

# test for one patient in HCC
threshold <- NULL
for (i in 1:AFP.per.patients.HCC$count[8]){
  y_past <- df.HCC$`AFP value`[df.HCC$ID=="A16"]
  threshold[i] <- qpeb(p=0.95,n=length(y_past[1:i]),
                       ybar=mean(y_past[1:i]),
                       params)
}

#test for subset of patients (linear trend in AFP over time: R^2>=0.7)
AFP.temp <- NULL
AFP.thresholds <- list()
test.1 <- df.HCC[1:60,]
AFP.subset <- AFP.r.squared %>% filter_at(vars("r.squared"), any_vars(.>=0.8))
for (i in unique(AFP.subset$ID)){
  for (j in 1:AFP.per.patients.HCC$count[AFP.per.patients.HCC$ID==i]){
    y_past <- df.HCC$`AFP value`[df.HCC$ID==i]
    AFP.temp[j] <- qpeb(p=0.95,n=length(y_past[1:j]),
                        ybar=mean(y_past[1:j]),
                        params)
  }
  AFP.thresholds[[i]] <- AFP.temp #store threshold values for patient i
  AFP.temp <- NULL # reset the vector
}

# test for all patients
AFP.thresholds <- list()
AFP.temp <- NULL
for (i in unique(df.HCC$ID)){
  for (j in 1:AFP.per.patients.HCC$count[AFP.per.patients.HCC$ID==i]){
    y_past <- df.HCC$`AFP value`[df.HCC$ID==i]
   AFP.temp[j] <- qpeb(p=0.95,n=length(y_past[1:j]),
                                ybar=mean(y_past[1:j]),
                                params)
  }
  AFP.thresholds[[i]] <- AFP.temp
  AFP.temp <- NULL
}
  
A <- unlist(AFP.thresholds)
df.HCC$threshold <- unlist(AFP.thresholds)

df.HCC$diff <- df.HCC$threshold-df.HCC$`AFP value`
A <- df.HCC %>% group_by(ID) %>% summarise(count=n(),flag.number=sum(diff<0))

A %>% filter(flag.number!=0) # number of flagged patient 179 out of 348 so TPR=51.4%

B <
# flagging AFP > Threshold




# Checking linearity of AFP vs Time 
AFP.linear <- df.HCC %>% group_by(ID) %>% do(model = lm(`AFP value` ~ `Age at sample/years`,
                                               data = .,na.action =na.exclude ))
tidy(AFP.linear, model)
AFP.r.squared <- glance(AFP.linear,model)
AFP.subset <- AFP.r.squared %>% filter_at(vars("r.squared"), any_vars(.>=0.8))

# remove factor
droplevels()

