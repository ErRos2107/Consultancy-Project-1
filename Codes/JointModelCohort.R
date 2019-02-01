# load libraries
rm( list = ls ())
if (! require (readxl)) install . packages ("readxl"); require (readxl)
if (! require ( tidyverse )) install . packages (’tidyverse ’); require (tidyverse)
if (! require ( reshape2 )) install . packages ("reshape2"); require (reshape2)
if (! require ( ggplot2 )) install . packages ("ggplot2"); require (ggplot2)
if (! require ( JM )) install . packages ("JM"); require (JM)


setwd("/Users/eric/Library/Mobile Documents/com~apple~CloudDocs/MSc Statisitcs with Data Science - University of Edinburgh/Consultancy Style projects /Project 1 - Liver cancer /Data")


df.JM <- read_excel("Screening_cohort_clean.xlsx",sheet = 1,
                    col_types = c("text",rep("numeric",13)))

# set time at sample
df.JM$Obstime <- df.JM$Sample_time-df.JM$Time

# Add censoring status
df.JM$Status <- rep(0,length(df.JM$ID))

# Select important covariates
data.JM.cohort <- df.JM %>% dplyr::select(ID:Age,Time,ALD:Other,AFP,Obstime,Status)

# Filter patients with at least two measurements
df.JM.temp <- data.JM.cohort %>% group_by(ID) %>% filter(n() >= 2)

#Filter NA values
df.JM.temp <- df.JM.temp %>% dplyr::select(ID:Age,Time,ALD:Other,AFP,Obstime,Status) %>%
  filter_all(all_vars(.!='NA'))

#Change sign of Time to event
df.JM.temp$Time <- as.numeric(-df.JM.temp$Time)

#remove left censored event
df.JM.temp2 <- df.JM.temp%>%dplyr:: select(ID:Age,Time,ALD:Other,AFP,Obstime,Status)%>%
  filter_at(vars(Time), any_vars(.>0))

# check that Time >= Obstime
stopifnot(df.JM.temp2$Time >= df.JM.temp2$Obstime)

# Remove inconsitencies in data reporting (patient missing, wrond date of entry etc...)
JM.cohort <- df.JM.temp2 %>% dplyr::select(ID:Age,Time,ALD:Other,AFP,Obstime,Status)%>%
  filter_at(vars(Obstime), any_vars(.>=0))
# rename column
colnames(JM.cohort)[colnames(JM.cohort) == 'Time'] <- 'TimeEvent'

# export cleaned data as spreasheet
write.xlsx(JM.cohort, "Long_Cohort_clean_final.xlsx")
