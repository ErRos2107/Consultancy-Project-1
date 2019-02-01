
# load libraries
rm( list = ls ())
if (! require (readxl)) install . packages ("readxl"); require (readxl)
if (! require ( tidyverse )) install . packages (’tidyverse ’); require (tidyverse)
if (! require ( reshape2 )) install . packages ("reshape2"); require (reshape2)

setwd("/Users/eric/Library/Mobile Documents/com~apple~CloudDocs/MSc Statisitcs with Data Science - University of Edinburgh/Consultancy Style projects /Project 1 - Liver cancer /Data")


HCC.data <-read_excel("HCC_in_screening_clean.xlsx",sheet = "Sheet1",
                      col_types = c("text","numeric","date",rep("numeric",8),"date","numeric","numeric"))

# reorder dataframe
HCC.data <- HCC.data[order(HCC.data$ID),]

## ----------------- cleaning the data --------------------

# Number of measurement per patient/ID
count.HCC <- HCC.data %>% group_by(ID) %>% summarize(count=n())

#check vector length
stopifnot(count.HCC$ID==unique(HCC.data$ID))

# Fill missing longitudonal measurements for variable 'Sex'
sex.HCC <- HCC.data %>%  dplyr::select(Sex)%>% filter(Sex!="NA") # remove NA
aetiology.HCC <- HCC.data %>%dplyr:: select(ALD:Other) %>%
  filter_all(all_vars(.!='NA')) #remove NA
HCC.data$Sex <- rep(sex.HCC$Sex,times=count.HCC$count)
for (i in 1:dim(aetiology.HCC)[2]){
  HCC.data[,6+i] <- rep(pull(aetiology.HCC[,i]),times=count.HCC$count)
}

# Add status (event) , 1 for HCC 0 for no HCC
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
HCC.temp1 <- HCC.data%>% dplyr::select(ID,Sex,Age_sample,AFP,Time,ALD:Other,
                                       Status,Obstime,TimeEvent)%>%
  filter_all(all_vars(.!='NA'))

#remove left censored event
HCC.temp2 <- HCC.data%>%dplyr:: select(ID,Sex,Age_sample,AFP,Time,ALD:Other
                                       ,Status,Obstime,TimeEvent)%>%
  filter_at(vars(TimeEvent), any_vars(.>0))

HCC.temp3 <- HCC.data%>%dplyr:: select(ID,Sex,Age_sample,AFP,Time,ALD:Other
                                       ,Status,Obstime,TimeEvent)%>%
  filter_at(vars(Time), any_vars(.<0))

# check that TimeEvent > Obstime
stopifnot(HCC.temp3$TimeEvent > HCC.temp3$Obstime)

indices <- which(HCC.temp3$TimeEvent < HCC.temp3$Obstime)

# remove observations where TimeEvent < Obstime
HCC.temp3 <- HCC.temp3[-indices,]

# rename columns
colnames(HCC.temp3)[colnames(HCC.temp2) == 'Age_sample'] <- 'Age'

length(unique(HCC.data$ID))

# export cleaned data as spreasheet
write.xlsx(HCC.temp3, "Long_HCC_clean_final.xlsx")
