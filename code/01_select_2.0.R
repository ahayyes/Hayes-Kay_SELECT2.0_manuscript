################ SELECT 1.0, 2.0 Analysis - Antibiotics ###########################

# Purpose of this script is to use the EC1 for all antibiotics tested in this study

## our data in these are in mg/L. To calculate ug/L multiply by 1000


####################### Load Libraries ###########################

library(tidyverse)
library(ggplot2)
library(reshape2)
library(MetBrewer)
library(plotrix)
library(gcplyr)
library(drc)
library(dunn.test)


##################### Confidence Interval ############################

conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}


###################### Amoxicillin  ##############################

# read in file
df <- read.csv("data/amoxicillin_221223.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<54), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) #  8.284e-13

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  3.255e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # 1.483e-15

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) #5.321e-10

dunn.test(df$X60, df$conc)

# SELECT 1.0 LOEC = 0.0625    mg/L or 62.5ug/L



########## now do EC1

# exponential phase df
xpo <- subset(df_melt, variable<54)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="confidence", add = TRUE)


#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED1 and ED 10 
# also displaying the 95% confidence intervals
# we also want to see if the CIs get smaller with a smaller EC

ED <- ED(drc, c(1, 10), interval="delta")

# Estimated effective doses
# 
# Estimate Std. Error    Lower    Upper
# e:1:1  0.144553   0.054278 0.036651 0.252454
# e:1:10 0.503950   0.086973 0.331053 0.676847

# ED1 =  0.144553  mg/L  = 144.553 ug/L

write.csv(ED, "readouts/ED_Amoxicillin.csv")




################### Ampicillin ########################

# read in file
df <- read.csv("data/ampicillin_050923.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)


# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))

## graph growth curve
ggplot(subset(summ, variable>2 & variable<55), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=26))+
  theme(axis.title = element_text(size = 20))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
#want to find the point in exponential phase with the largest dose response. this is hour five for this data
# for your data you will want to test all time points within exponential phase (so normally something between hours 4-12 and see which has the largest correlation coefficient)

# pearsons = normal
# spearmans = non normal


## list of hour time points when taking 10 minute readings
# 0, #60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 660
# in our numbers they are 
# 0, 6, 12, 18 ,24, 30, 36, 42, 48, 54, 60, 66

# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X60) # 2.307e-16
cor.test(df$conc, df$X60, method = c("spearman")) #1.556e-12


shapiro.test(df$X42) #9.279e-14
cor.test(df$conc, df$X42, method = c("spearman")) #< 4.33e-16

shapiro.test(df$X66) #0.0001144
cor.test(df$conc, df$X66, method = c("spearman")) #1.376e-09

shapiro.test(df$X72) #0.0006799
cor.test(df$conc, df$X72, method = c("spearman")) #0.01114

shapiro.test(df$X48) # 4.499e-13
cor.test(df$conc, df$X48, method = c("spearman")) #< 2.2e-16 ## this one 

dunn.test(df$X48, df$conc)

# SELECT 1.0 LOEC is 1000ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<55)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="confidence", add = TRUE)


# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED1 and ED 10 
# also displaying the 95% confidence intervals
# we also want to see if the CIs get smaller with a smaller EC

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Ampicillin.csv")




############## Azithromycin #######################

# read in file
df <- read.csv("data/azithromycin_160823.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<58), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=26))+
  theme(axis.title = element_text(size = 30))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
#want to find the point in exponential phase with the largest dose response. this is hour five for this data
# for your data you will want to test all time points within exponential phase (so normally something between hours 4-12 and see which has the largest correlation coefficient)

# pearsons = normal
# spearmans = non normal


## list of hour time points when taking 10 minute readings
# 0, #60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 660
# in our numbers they are 
# 0, 6, 12, 18 ,24, 30, 36, 42, 48, 54, 60, 66

# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X60) # 1.958e-10 N Normal (NN)
cor.test(df$conc, df$X60, method = c("spearman")) #< 2.2e-16


shapiro.test(df$X54) # (NN)
cor.test(df$conc, df$X54, method = c("spearman")) #< 2.2e-16

shapiro.test(df$X42) #1.022e-10 (NN)
cor.test(df$conc, df$X42, method = c("spearman")) #< 2.2e-16

shapiro.test(df$X66) #0.0001144 NN 
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) #0.0006799
cor.test(df$conc, df$X72, method = c("spearman")) #< 2.2e-16

shapiro.test(df$X48) # 4.499e-13
cor.test(df$conc, df$X48, method = c("spearman")) #< 2.2e-16 ## this one 


# if all same do hour 11 
dunn.test(df$X66, df$conc)

#LOEC is 0.78125 mg/L which is 781.25 ug/L

############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<58)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)


plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

# looks like an alright fit actually

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

# v good fit there

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Azithromycin.csv")



######################## Cefotaxime #################################

# read in file
df <- read.csv("data/cefotaxime_090823.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 
summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<54), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X60) # 2.307e-16
cor.test(df$conc, df$X60, method = c("spearman")) #1.556e-12

shapiro.test(df$X42) #9.279e-14
cor.test(df$conc, df$X42, method = c("spearman")) #< 4.33e-16

shapiro.test(df$X66) #0.0001144
cor.test(df$conc, df$X66, method = c("spearman")) #1.376e-09

shapiro.test(df$X72) #0.0006799
cor.test(df$conc, df$X72, method = c("spearman")) #0.01114

shapiro.test(df$X48) # 4.499e-13
cor.test(df$conc, df$X48, method = c("spearman")) #< 2.2e-16 ## this one 


# do hour 66, latest but with greatest difference
dunn.test(df$X48, df$conc)


# 1.0 LOEC is 0.03125 mg/L or 31.25ug/L




############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<54)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

# lets look at the data first
plot(maxDens ~ log(conc+.1), data = select, main="Logarithmic Dose Scale")

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

# looks like an alright fit actually

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

# v good fit there

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ec1 and ec10
# also dispolying the 95% confidence intervals


ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Cefotaxime.csv")



####################### Ceftriofur ###########################

# read in file
df <- read.csv("data/ceftiofur_110823.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 
summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve

ggplot(subset(summ, variable>2 & variable<55), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X60) # 2.307e-16
cor.test(df$conc, df$X60, method = c("spearman")) #1.556e-12

shapiro.test(df$X42) #9.279e-14
cor.test(df$conc, df$X42, method = c("spearman")) #< 4.33e-16

shapiro.test(df$X66) #0.0001144
cor.test(df$conc, df$X66, method = c("spearman")) #1.376e-09

shapiro.test(df$X72) #0.0006799
cor.test(df$conc, df$X72, method = c("spearman")) #0.01114

shapiro.test(df$X48) # 4.499e-13
cor.test(df$conc, df$X48, method = c("spearman")) #< 2.2e-16 ## this one 


# do hour 66, latest but with greatest difference
dunn.test(df$X60, df$conc)

# 1.0 LOEC is 62.5



############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))



select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

# lets look at the data first
plot(maxDens ~ log(conc+.1), data = select, main="Logarithmic Dose Scale")

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

# looks like an alright fit actually

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

# v good fit there

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also dispolying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Ceftiofur.csv")




################# Cetriaxone ######################

# read in file
df <- read.csv("data/ceftriaxone_060223.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<66), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))




############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) #9.789e-08
cor.test(df$conc, df$X42, method = c("spearman")) #< 2.2e-16

shapiro.test(df$X48) # 8.888e-07
cor.test(df$conc, df$X48, method = c("spearman")) # < 2.2e-16e 

shapiro.test(df$X54) #4.464e-05
cor.test(df$conc, df$X54, method = c("spearman")) #< 2.2e-16

shapiro.test(df$X60) #0.0002633
cor.test(df$conc, df$X60, method = c("spearman"))# < 2.2e-16

shapiro.test(df$X66) #0.03838
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) #0.02255
cor.test(df$conc, df$X72, method = c("spearman")) # < 2.2e-16

# could use any

# do hour 12, latest but with greatest difference
dunn.test(df$X72, df$conc)


# 1.0 LOEC is 0.03125 mg/L or 31.25ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<66)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

# lets look at the data first
plot(maxDens ~ log(conc+.1), data = select, main="Logarithmic Dose Scale")

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

# v good fit there

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also dispolying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Cefriaxone.csv")



###################### Chloramphenicol ##############################

# read in file
df <- read.csv("data/chloramphenicol_100823.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<66), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) #< 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) # < 2.2e-16e 

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # 2.669e-11

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman"))# < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("spearman")) # 3.113e-08


# do hour 11, latest but with greatest difference
dunn.test(df$X66, df$conc)


# 1.0 LOEC is 1 mg/L or 1000 ug/L




############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<66)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

# lets look at the data first
plot(maxDens ~ log(conc+.1), data = select, main="Logarithmic Dose Scale")

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

# v good fit there

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also dispolying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Chloramphenicol.csv")


###################### Ciprofloxacin ##############################


# read in file
df <- read.csv("data/ciprofloxacin_020623.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<72), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) #< 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) # < 2.2e-16e 

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # 2.669e-11

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman"))# < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("spearman")) # 3.113e-08


# do hour 11, latest but with greatest difference
dunn.test(df$X72, df$conc)


# 1.0 LOEC is 1 mg/L or 1000 ug/L




############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))



select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

# lets look at the data first
plot(maxDens ~ log(conc+.1), data = select, main="Logarithmic Dose Scale")

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

# v good fit there

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Ciprofloxacin.csv")





###################### Clarithromycin ##############################

# read in file
df <- read.csv("data/clarithromycin_280623.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve

ggplot(subset(summ, variable>2 & variable<72), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # N
cor.test(df$conc, df$X42, method = c("pearson")) # 1.073e-14

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) # < 2.2e-16e 

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # 2.669e-11

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman"))# < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("spearman")) #  7.051e-10


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)


# 1.0 LOEC is 4 mg/L or 4000 ug/L



############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

# lets look at the data first
plot(maxDens ~ log(conc+.1), data = select, main="Logarithmic Dose Scale")

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

# v good fit there

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)

# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Clarithromycin.csv")



###################### Colistin Sulfate ##############################

# read in file
df <- read.csv("data/colistinsulfate_180624.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")

# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 
summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # N
cor.test(df$conc, df$X42, method = c("pearson")) # 1.162e-10

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) # < 2.2e-16e 

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # 9.587e-13

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # 5.066e-08 

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # 2.951e-10

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("spearman")) #   1.806e-10


# do 11 latest but largest correlation
dunn.test(df$X48, df$conc)


# 1.0 LOEC is 0.0625mg/L or 62.5 ug/L

############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

# lets look at the data first
plot(maxDens ~ log(conc+.1), data = select, main="Logarithmic Dose Scale")

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

# v good fit there

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Colistin.csv")



###################### Doripenem  ##############################

# read in file
df <- read.csv("data/doripenem_270624.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<66), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  3.459e-10

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16 

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("spearman")) #   < 2.2e-16


# do 11 latest but largest correlation
dunn.test(df$X72, df$conc)


# 1.0 LOEC is 0.015625   mg/L or 15.625 ug/L



############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

# lets look at the data first
plot(maxDens ~ log(conc+.1), data = select, main="Logarithmic Dose Scale")

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Doripenem.csv")



###################### Doxycycline  ##############################

# read in file
df <- read.csv("data/doxycycline_150124.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve

ggplot(subset(summ, variable>2 & variable<50), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # N
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16e

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) # < 2.2e-16e 

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16e

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16e

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16e

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("spearman")) #   < 2.2e-16e


# do 11 latest but largest correlation
dunn.test(df$X72, df$conc)


# 1.0 LOEC is 0.0625mg/L or 62.5 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<50)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

# lets look at the data first
plot(maxDens ~ log(conc+.1), data = select, main="Logarithmic Dose Scale")

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)



# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Doxycycline.csv")



###################### Enrofloxacin  ##############################


# read in file
df <- read.csv("data/enrofloxacin_17112025.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<48), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # N
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16e

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) # < 2.2e-16e 

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16e

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16e

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) #  2.642e-08

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("spearman")) #   0.0345


# do 11 latest but largest correlation
dunn.test(df$X60, df$conc)


# 1.0 LOEC is 1.25 mg/L or 125 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<48)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))

## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Enrofloxacin.csv")



###################### Erythromycin  ##############################

# read in file
df <- read.csv("data/erythromycin_270623.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # N
cor.test(df$conc, df$X42, method = c("spearman")) # 3.898e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) # < 2.2e-16e 

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16e

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16e

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16e

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) # 9.514e-10


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)


# 1.0 LOEC is 0.0625mg/L or 62.5 ug/L

############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<66)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)



# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Erythromycin.csv")



###################### Florfenicol  ##############################

# read in file
df <- read.csv("data/florfenicol_281122.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<30), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # N
cor.test(df$conc, df$X42, method = c("spearman")) # 5.503e-08

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) # 1.625e-07

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) #  5.87e-07

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) #  5.87e-07

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # 1.84e-09

# shapiro.test(df$X72) # NN
# cor.test(df$conc, df$X72, method = c("pearson")) # 9.514e-10


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)


# 1.0 LOEC is 0.0625mg/L or 62.5 ug/L

############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<30)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Florfenicol.csv")



###################### Gentamicin  ##############################

# read in file
df <- read.csv("data/gentamicin_260723.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<72), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) #3.077e-05

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  2.094e-11

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) #  1.879e-13

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) #  7.135e-14

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # 1.863e-10

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) # < 2.2e-16


# do 11 latest but largest correlation
dunn.test(df$X72, df$conc)


# 1.0 LOEC is 0.5 mg/L or 500 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<66)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve


## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Gentamicin.csv")


###################### Imipenem  ##############################

# read in file
df <- read.csv("data/imipenem_250723.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve

ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) # < 2.2e-16


# do 11 latest but largest correlation
dunn.test(df$X72, df$conc)


# 1.0 LOEC is 0.0625 mg/L or 62.5ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))




########## now do ED10

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Imipenem.csv")



###################### Kanamycin  ##############################

# read in file
df <- read.csv("data/kanamycin_070524.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) # 1.328e-10


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)


# 1.0 LOEC is 0.5 mg/L or 500 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Kanamycin.csv")



###################### Meropenem  ##############################

# read in file
df <- read.csv("data/meropenem_070223.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve


ggplot(subset(summ, variable>2 & variable<54), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) # 6.128e-09


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)


###### if we did the itme point we cut by for growth then we get a similar LOEC 

# 1.0 LOEC is 0.125   mg/L which is 125 ug/L 



############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<54)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
# 
# drc <- drm(maxDens ~ conc, data = select, 
#            fct = LL.4())

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)



# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Meropenem.csv")



###################### Nitrofurantoin  ##############################

# read in file
df <- read.csv("data/nitrofurantoin_200624.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve

ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))




############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

# shapiro.test(df$X72) # NN
# cor.test(df$conc, df$X72, method = c("pearson")) # 6.128e-09


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)

# SELECT 1.0 LOEC = 500 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<66)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)



# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Nitrofurantoin.csv")



###################### Norfloxacin  ##############################

# read in file
df <- read.csv("data/norfloxacin_310523.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve

ggplot(subset(summ, variable>2 & variable<66), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # 2.506e-12

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) #  4.65e-13

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # 4.943e-08

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) # 3.975e-09


# do 11 latest but largest correlation
dunn.test(df$X54, df$conc)

# SELECT 1.0 LOEC = 0.03125 mg/L or 31.25 ug/L



############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Norfloxacin.csv")



###################### Ofloxacin  ##############################


# read in file
df <- read.csv("data/ofloxacin_010623.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve

ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) #  4.286e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) # 5.103e-05


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)

# SELECT 1.0 LOEC = 0.007812  mg/L or 7.812 ug/L



############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)

# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Ofloxacin.csv")



###################### Oxytetracycline  ##############################

# read in file
df <- read.csv("data/oxytetracycline_241123.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<54), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) # 0.003051


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)

# SELECT 1.0 LOEC = 0.03125   mg/L or 31.25 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<54)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")



#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)



# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Oxytetracycline.csv")


###################### Penicillin ##############################

# read in file
df <- read.csv("data/penicillin_221123.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve

ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # 1.72e-13

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) #  7.021e-08


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)

# SELECT 1.0 LOEC = 0.5   mg/L or 500 ug/L



############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<54)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Penicillin.csv")


###################### Streptomycin ##############################

# read in file
df <- read.csv("data/streptomycin_290323.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve

ggplot(subset(summ, variable>2 & variable<66), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # 1.72e-13

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # 9.257e-14

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # 9.257e-14

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) #  0.0001008


# do 11 latest but largest correlation
dunn.test(df$X48, df$conc)

# SELECT 1.0 LOEC = 0.03125   mg/L or 31.25 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<66)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")



#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Streptomycin.csv")


###################### Sulfamethoxazole ##############################

# read in file
df <- read.csv("data/sulfamethoxazole_190224.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve


ggplot(subset(summ, variable>2 & variable<54), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) #  0.0009939


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)

# SELECT 1.0 LOEC = 0.25   mg/L or 250 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<54)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Sulfamethoxazole.csv")




###################### Sulfadiazine ##############################

# read in file
df <- read.csv("data/To Deposit/sulfadiazine_090124.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve


ggplot(subset(summ, variable>2 & variable<54), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("spearman")) #  < 2.2e-16


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)

# SELECT 1.0 LOEC = 0.234375     mg/L or 234.375 ug/L



############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<54)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

## now do ED10

select$conc <- as.numeric(as.character(select$conc))

## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Sulfadiazine.csv")


###################### Sulfapyridine ##############################

# read in file
df <- read.csv("data/sulfapyridine_220124.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # NN
cor.test(df$conc, df$X72, method = c("pearson")) #  3.923e-07


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)

# SELECT 1.0 LOEC = 4   mg/L or 4000 ug/L



############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Sulfapyridine.csv")



###################### Tetracycline ##############################

# read in file
df <- read.csv("data/tetracycline_141023.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X72) # N
cor.test(df$conc, df$X72, method = c("pearson")) #   3.415e-07


# do 11 latest but largest correlation
dunn.test(df$X66, df$conc)

# SELECT 1.0 LOEC = 4   mg/L or 4000 ug/L

############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Tetracycline.csv")


###################### Thiamphenicol ##############################

# read in file
df <- read.csv("data/thiamphenicol_240723.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # 3.202e-13

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

# shapiro.test(df$X72) # N
# cor.test(df$conc, df$X72, method = c("spearman")) #   < 2.2e-16


# did time point 10 because time point 11 had difficult to manage trend of significance 
dunn.test(df$X60, df$conc)

# SELECT 1.0 LOEC = 1   mg/L or 1000 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Thiamphenicol.csv")



###################### Trimethoprim ##############################

# read in file
df <- read.csv("data/trimethoprim_100323.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve

ggplot(subset(summ, variable>2 & variable<54), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))


############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) # < 2.2e-16

# shapiro.test(df$X72) # N
# cor.test(df$conc, df$X72, method = c("spearman")) #   < 2.2e-16


# did time point 10 because time point 11 had difficult to manage trend of significance 
dunn.test(df$X60, df$conc)

# SELECT 1.0 LOEC = 0.03125    mg/L or 31.25 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<54)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")


#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)


# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Trimethoprim.csv")


###################### Vancomycin ##############################

# read in file
df <- read.csv("data/vancomycin_130224.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

#### Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")


# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe 

summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))


summ$variable <- as.numeric(as.character(summ$variable))


## graph growth curve
ggplot(subset(summ, variable>2 & variable<60), aes(x=variable, y=avg, color=conc)) + #can subset to whichever hour you prefer
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=conc),
              color=NA, alpha=0.1) +
  geom_line(linewidth = 1) +
  theme_bw()+
  labs(y='OD600', x="Time") +
  labs(colour = "Treatment (mg/L)")+
  labs(fill = "Treatment (mg/L)")+
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black", fill="white"))+
  scale_color_manual(values=met.brewer("Hiroshige", 16))+
  scale_fill_manual(values=met.brewer("Hiroshige", 16))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size = 22))+
  ylab(bquote('Absorbance '~'('~OD[600]~')'))



############### SELECT 1.0 

#shapiro test - normal if the p value is above 0.05
# pearsons = normal
# spearmans = non normal


# make df numeric for conc and remove ntc
df <- subset(df, conc !="ntc")
df$conc <- as.numeric(as.character(df$conc)) 

shapiro.test(df$X42) # NN
cor.test(df$conc, df$X42, method = c("spearman")) # 3.202e-13

shapiro.test(df$X48) # NN
cor.test(df$conc, df$X48, method = c("spearman")) #  < 2.2e-16

shapiro.test(df$X54) # NN
cor.test(df$conc, df$X54, method = c("spearman")) # < 2.2e-16

shapiro.test(df$X60) # NN 
cor.test(df$conc, df$X60, method = c("spearman")) # 8.75e-14

shapiro.test(df$X66) # NN
cor.test(df$conc, df$X66, method = c("spearman")) #5.913e-10

# shapiro.test(df$X72) # N
# cor.test(df$conc, df$X72, method = c("spearman")) #   < 2.2e-16


# did time point 10 because time point 11 had difficult to manage trend of significance 
dunn.test(df$X54, df$conc)

# SELECT 1.0 LOEC = 8   mg/L or 8000 ug/L


############## SELECT 2.0 #

# exponential phase df
xpo <- subset(df_melt, variable<60)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

## get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))


select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

## make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this 

plot(drc, broken=TRUE, bty="l",
     xlab="Concentration of Antibiotic", ylab="Maximum OD", type="all")

#Test for lack of fit
modelFit(drc)

# plot residuals
qqnorm(residuals(drc))
qqline(residuals(drc))

plot(residuals(drc) ~ fitted(drc), main="Residuals vs Fitted")
abline(h=0)

# get ED 10 and ED 5
# also displaying the 95% confidence intervals

ED <- ED(drc, c(1, 10), interval="delta")
write.csv(ED, "readouts/ED_Vancomycin.csv")

