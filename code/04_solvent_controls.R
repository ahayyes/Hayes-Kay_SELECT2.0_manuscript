####### Solvent Control for Paper#########

# Purpose of this script is to analyse the solvent controls generated as
# part of the this study


####################### Load Libraries ###########################

library(tidyverse)
library(ggplot2)
library(reshape2)
library(MetBrewer)
library(plotrix)
library(gcplyr)


##################### Confidence Interval ############################

conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}


################### NaOH Control ########################

# read in file
df <- read.delim("data/naoh_control.txt", header=TRUE)


### in this plate we had A - NaOH a1-A6
### we had positive control E1 - E6 
#### We had a negative control which was fine

# remove blanks 

df <- subset(df, select = c(A1:A6, E1:E6))

# rename columns
names(df)[names(df)=="A1"] <- "NaOH R1"
names(df)[names(df)=="A2"] <- "NaOH R2"
names(df)[names(df)=="A3"] <- "NaOH R3"
names(df)[names(df)=="A4"] <- "NaOH R4"
names(df)[names(df)=="A5"] <- "NaOH R5"
names(df)[names(df)=="A6"] <- "NaOH R6"

names(df)[names(df)=="E1"] <- "PosControl R1"
names(df)[names(df)=="E2"] <- "PosControl R2"
names(df)[names(df)=="E3"] <- "PosControl R3"
names(df)[names(df)=="E4"] <- "PosControl R4"
names(df)[names(df)=="E5"] <- "PosControl R5"
names(df)[names(df)=="E6"] <- "PosControl R6"


#Time column
names(df)[names(df)=="Kinetic.read"] <- "Reading"

df$Reading <- seq(0, 1440, by=10)


#Melt dataframe
df_melt <- melt(df, id=c("Reading"))

#Separate
df_melt <- separate(df_melt, variable, c("Treatment", "Replicate"), sep=" ")


# summarise dataframe 

summ <- df_melt %>%
  group_by(Reading, Treatment) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


# plot growth in exponential phase (hour 12 = 720 minutes), 

ggplot(subset(summ, Reading<800), aes(x=Reading, y=avg, color=Treatment)) +
  geom_ribbon(aes(ymin=avg-CI95, ymax=avg+CI95, fill=Treatment),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment")+
  labs(fill = "Treatment")+
  theme(text=element_text(size=20))+
  theme(legend.position="bottom")+
  scale_fill_manual(values=met.brewer("Hiroshige", 3), name = "Treatment")+
  scale_color_manual(values=met.brewer("Hiroshige", 3), name = "Treatment")+
  labs(x="Time (minutes)")+
  theme(strip.background = element_rect(fill = "white", colour = "black"))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))

