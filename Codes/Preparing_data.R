library("dplyr")
library("readxl")
library("xlsx")
library("reshape2")

setwd("/Users/eric/Library/Mobile Documents/com~apple~CloudDocs/MSc Statisitcs with Data Science - University of Edinburgh/Consultancy Style projects /Project 1 - Liver cancer /Data")

######################################################################################################
#------------------------------------- HCC in screening dataset --------------------------------------
######################################################################################################
######################################################################################################

HCC.data <-read_excel("HCC_in_screening_clean.xlsx",sheet = "Sheet1",
                      col_types = c("text","numeric","date",rep("numeric",8),"date","numeric","numeric"))

# reorder dataframe
HCC.data <- HCC.data[order(HCC.data$ID),]
## cleaning the data 
count.HCC <- HCC.data %>% group_by(ID) %>% summarize(count=n())
#check ordering
stopifnot(count.HCC$ID==unique(HCC.data$ID))

#Impute values into long format
sex.HCC <- HCC.data %>%  dplyr::select(Sex)%>% filter(Sex!="NA") 
aetiology.HCC <- HCC.data %>%dplyr:: select(ALD:Other) %>%
  filter_all(all_vars(.!='NA'))
HCC.data$Sex <- rep(sex.HCC$Sex,times=count.HCC$count)
for (i in 1:dim(aetiology.HCC)[2]){
  HCC.data[,6+i] <- rep(pull(aetiology.HCC[,i]),times=count.HCC$count)
}

# Add status (event)
HCC.data$Status <- rep(1,length(HCC.data$ID))

# Add observation time (time point for each measurements) per patient
Obstime <- list()
obstime.temp <- NULL
for (i in unique(HCC.data$ID)){
  for (j in 1:count.HCC$count[count.HCC$ID==i]){
    patient.dates <- HCC.data$"Sample_date"[HCC.data$ID==i]
    init.date <- patient.dates[1]
    obstime.temp[j] <- difftime(patient.dates[j],init.date,units = "days")
  }
  Obstime[[i]] <- obstime.temp
  obstime.temp <- NULL
}

HCC.data$Obstime <- unlist(Obstime)

# Recalculate time to event 
data.sample.temp <- HCC.data %>% group_by(ID) %>% filter(row_number(ID) == 1) 
date.sample <- ymd(as.character(data.sample.temp$Sample_date))
start.date <- date.sample[!is.na(date.sample)]
date.diagnosis <- ymd(as.character(HCC.data$Diagnosis))
end <- date.diagnosis[!is.na(date.diagnosis)]

time.to.event <- difftime(end,start.date,units="days")
HCC.data$TimeEvent <- rep(time.to.event,times=count.HCC$count)

# Remove NAN 
HCC.temp1 <- HCC.data%>% dplyr::select(ID,Sex,Age_sample,Time,AFP,ALD:Other,
                                       Status,Obstime,TimeEvent)%>%
  filter_all(all_vars(.!='NA'))

# rename column
colnames(HCC.temp1)[colnames(HCC.temp1) == 'Age_sample'] <- 'Age'
colnames(HCC.temp1)[colnames(HCC.temp1) == 'HEPC'] <- 'HepC'

#remove left censored event
HCC.temp2 <- HCC.temp1%>%dplyr:: select(ID,Sex,Age,AFP,Time,ALD:Other
                                       ,Status,Obstime,TimeEvent)%>%
  filter_at(vars(TimeEvent), any_vars(.>0))

#Remove longitudinal measurements after event/diagnosis
HCC.temp3 <- HCC.temp2%>%dplyr:: select(ID,Sex,Age,Time,AFP,ALD:Other
                                      ,Status,Obstime,TimeEvent)%>%
filter_at(vars(Time), any_vars(.<0))

# check that TimeEvent > Obstime 
stopifnot(HCC.temp3$TimeEvent > HCC.temp3$Obstime)

indices <- which(HCC.temp3$TimeEvent < HCC.temp3$Obstime)
# remove 
HCC.temp3 <- HCC.temp3[-indices,]

HCC.data <- HCC.temp3 %>% dplyr:: select(ID,Sex,Age,AFP,ALD:Other
                                        ,Status,Obstime,TimeEvent)

######################################################################################################
#------------------------------------- Cohort in screening dataset --------------------------------------
######################################################################################################
######################################################################################################

df.JM <- read_excel("Screening_cohort_clean.xlsx",sheet = 1,
                    col_types = c("text",rep("numeric",13)))

# set time at sample
df.JM$Obstime <- df.JM$Sample_time-df.JM$Time

# Add censoring status
df.JM$Status <- rep(0,length(df.JM$ID))

# Select important covariates
data.JM.cohort <- df.JM %>% dplyr::select(ID:Age,Time,ALD:Other,AFP,Obstime,Status)

# Filter patients with at least two measurements 
df.JM.temp <- data.JM.cohort %>% group_by(ID) %>% filter(row_number(ID) >= 2)

#Change sign of Time to event
df.JM.temp$Time <- as.numeric(-df.JM.temp$Time)

# rename column
colnames(df.JM.temp)[colnames(df.JM.temp) == 'Time'] <- 'TimeEvent'


#Filter NA values
df.JM.temp <- df.JM.temp %>% dplyr::select(ID:Age,TimeEvent,ALD:Other,AFP,Obstime,Status) %>%
  filter_all(all_vars(.!='NA'))

#remove left censored event
df.JM.temp2 <- df.JM.temp %>% dplyr:: select(ID:Age,TimeEvent,ALD:Other,AFP,Obstime,Status)%>%
  filter_at(vars(TimeEvent), any_vars(.>0))

# check that Time >= Obstime 
stopifnot(df.JM.temp2$TimeEvent >= df.JM.temp2$Obstime)

# Remove inconsitencies in data reporting (patient missing, wrond date of entry etc...)
JM.cohort <- df.JM.temp2 %>% dplyr::select(ID:Age,TimeEvent,ALD:Other,AFP,Obstime,Status)%>%
  filter_at(vars(Obstime), any_vars(.>=0))

#------------------------------------- Joining the two datasets --------------------------------------

full.data.long <- bind_rows(HCC.data,JM.cohort)

#Transform AFP levels to lo(AFP+0.01)
full.data.long$AFP <- log(full.data.long$AFP+0.01,base = 10)
# check that Time >= Obstime 
stopifnot(full.data.long$TimeEvent >= full.data.long$Obstime)


# export cleaned data as spreasheet
write.xlsx(full.data.long, "Full_data_final.xlsx")

