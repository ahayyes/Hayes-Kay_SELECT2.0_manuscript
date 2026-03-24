####################### SELECT 2.0 Walkthrough #############################

## This script will allow you to estimate effect concentrations for antimicrobials 
## and visualise growth of a microbial community 

####################### Load Libraries ###########################

library(tidyverse)
library(ggplot2)
library(reshape2)
library(MetBrewer)
library(plotrix)
library(drc)

##################### Confidence Interval ############################

conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}

############################ Step 1: read in data ########################

# read in file
df <- read.csv("data/Antibiotic_OD_Files/amoxicillin.csv", header=TRUE)

# paste conc and rep
df2 <- df %>% unite(conc_rep, "conc", "replicate", remove = TRUE)

# Melt the dataframe 
df_melt <- melt(df2, id=c("conc_rep"))

df_melt <- df_melt %>%
  separate(conc_rep, c("conc", "rep"), sep = "_")

# remove the X in the variable column 
df_melt$variable <- gsub("X", "", df_melt$variable)

# summarise dataframe across replicates per each treatment
summ <- df_melt %>%
  group_by(conc, variable) %>%
  summarise(N=length(value), 
            avg=mean(value), 
            CI95 = conf_int95(value),
            sd=sd(value))

# make sure this is numeric 
summ$variable <- as.numeric(as.character(summ$variable))



########## Step 2: Visualise the growth curves

# graph this from time point 2 and within exponential phase 
# the first time points can often be noisy due to condensation and machine settling

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



############### Step 3: Find the EC1 #############################

# make a dataframe of exponential phase data only 
xpo <- subset(df_melt, variable<54)

# remove ntc (these are media blanks)
xpo <- subset(xpo, conc != "ntc")

# get maximum growth rates
select <- xpo %>%
  group_by(conc, rep) %>%
  summarise(maxDens = max(value))

# make sure this is numeric 
select$conc <- as.numeric(as.character(select$conc))


## fit a four parameter log logistic curve

# make the model
drc <- drm(maxDens ~ conc, data = select, 
           fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(drc)

# plot this in base R
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

