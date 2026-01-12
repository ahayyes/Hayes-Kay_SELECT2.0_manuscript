##################### Script 3 - Environmental Risk Assessments ###############


# The purpose of this script is to generate environmental risk assessments (ERAs)
# using data from the UBA database for global risk, and the UKWIR Chemicals Investigations
# Programme for the England and Wales specific ERA

####################### Load Libraries ###########################

library(tidyverse)
library(ggplot2)
library(reshape2)
library(MetBrewer)
library(patchwork)
library(plotrix)
library(stringr)


##################### Confidence Interval ############################


conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}



################ Environmental Risk Assessments (ERAs) Global ######################

# These ERAs will use data collated from the UBA database
# this is freely available 

uba <- read.csv("data/uba_mecs_transposed.csv", header=TRUE)
#exported data was exported on the 08042024 at 2:34pm

# make R friendly
uba_melt <- melt(uba, id=c("water_type", "antibiotic"))

# remove variable column because unnecessary
uba_melt <- uba_melt %>% select(-variable)

# uba_melt <- merge(uba_melt, ec10s, by = "antibiotic")

#Add EC10 numbers 
uba_melt$EC10 <- if_else(uba_melt$antibiotic=="Amoxicillin", "1000",
                         if_else(uba_melt$antibiotic=="Ampicillin", "2000",
                                 if_else(uba_melt$antibiotic=="Azithromycin", "781", 
                                         if_else(uba_melt$antibiotic=="Cefotaxime", "15.6",
                                                 if_else(uba_melt$antibiotic=="Ceftiofur", "125",
                                                         if_else(uba_melt$antibiotic=="Ceftriaxone", "7.8",   
                                                                 if_else(uba_melt$antibiotic=="Chloramphenicol", "500", 
                                                                         if_else(uba_melt$antibiotic=="Ciprofloxacin","3.9",
                                                                                 if_else(uba_melt$antibiotic=="Clarithromycin","2000",
                                                                                         if_else(uba_melt$antibiotic=="Doxycycline","250", 
                                                                                                 if_else(uba_melt$antibiotic=="Enrofloxacin","4.88",     
                                                                                                         if_else(uba_melt$antibiotic=="Erythromycin", "2000",
                                                                                                                 if_else(uba_melt$antibiotic=="Florfenicol","937.5",
                                                                                                                         if_else(uba_melt$antibiotic=="Norfloxacin","31.25",  
                                                                                                                                 if_else(uba_melt$antibiotic=="Orfloxacin","7.81",
                                                                                                                                         if_else(uba_melt$antibiotic=="Oxytetracycline","125",     
                                                                                                                                                 if_else(uba_melt$antibiotic=="Streptomycin","250",
                                                                                                                                                         if_else(uba_melt$antibiotic=="Sulfamethoxazole","925",
                                                                                                                                                                 if_else(uba_melt$antibiotic=="Sulfapyridine","1250",
                                                                                                                                                                         if_else(uba_melt$antibiotic=="Tetracycline","250",
                                                                                                                                                                                 if_else(uba_melt$antibiotic=="Thiamphenicol","1000", "31.3")))))))))))))))))))))



# Make EC10 numeric
uba_melt$EC10 <- as.numeric(uba_melt$EC10)

# Make PNEC
uba_melt$pnec <- uba_melt$EC10/10

uba_melt$pnec <- as.numeric(uba_melt$pnec)

# Make RQ
uba_melt$RQ <- uba_melt$value / uba_melt$pnec

# Add level of risk
uba_melt$risk <- if_else(uba_melt$RQ>=1, "high risk", if_else(uba_melt$RQ>=0.1, "medium risk", "low risk"))

# remove rows with NAs in the value column
uba_melt <- uba_melt %>% drop_na(value)

# Reorder water type column for easy plotting
uba_melt$water_type <- factor(uba_melt$water_type, levels=c("Influent", "Effluent"))

# Generate summary table
rqsummary = uba_melt %>% 
  group_by(antibiotic, water_type) %>%
  summarise(MedianRQ = median(RQ), MaxRQ = max(RQ)) %>%
  gather(Risk, RQ, -antibiotic, -water_type)

write.csv(rqsummary, "data/RQ_Summary_UBA.csv")

# make RQ summary with only high risk things
rqsummary_high <- subset(rqsummary, RQ>=1)

write.csv(rqsummary_high, "data/RQ_Summary_UBA_highrisk.csv")

############ Make UBA Global RQ Plot

# specify log values for different RQ cut offs
lowrq = log10(0.1)+1
highrq = log10(1)+1

# boxplots of all data points log(RQs), facet wrapped by antibiotic for both sample types (colour). 
# Text labels and lines for different risk (RQ) groupings


# add antibiotic class to uba data
unique(uba_melt$antibiotic)

uba_melt <- merge(uba_melt, antibiotic_class, by = "antibiotic")

# make plot
uba_plot <- ggplot(uba_melt, aes(x=water_type, y=RQ+1)) +
  geom_hline(yintercept= 1.1, linetype= "longdash", linewidth = 1.2, colour = "darkred") + # the low line
  geom_hline(yintercept= 2, colour = "darkred", linewidth = 1.2) + # the high line
  geom_point(aes(colour = water_type), position = position_jitter(), alpha=0.3, size = 6, show.legend = FALSE)+
  facet_wrap(~class, scales = "free") +
  ylab(bquote(log[10](RQ+1)))+
  labs(x="Wastewater Type", 
       title = "A) Global Environmental Risk") +
  theme_bw()+
  scale_color_manual(values=met.brewer("Hiroshige", 2))+
  scale_fill_manual(values=met.brewer("Hiroshige", 2))+
  theme(legend.position="none",
        strip.background=element_rect(colour="black",
                                      fill="white"))+
  theme(text=element_text(size=26))+
  theme(axis.title = element_text(size = 30))+
  scale_y_log10(labels = scales::comma)

uba_plot

ggsave("figures/Figure 4A.tiff", height = 8, width = 12, dpi = 300)


############################ ERAs for England and Wales data #######################

cip3 <- read.csv("data/cip3_mecs.csv", header=TRUE)
#exported data was exported on the 08042024 at 2:34pm

#Filter to include Treatment Influent and Treatment Effluent only 
cip3 <- subset(cip3, SampleLocationName=="Treatment Effluent" | SampleLocationName=="Treatment Influent")

# remove Treatment from these cells as just need Influent and Effluent
cip3 <- cip3 %>%
  mutate(SampleLocationName = recode(SampleLocationName,
                                     `Treatment Effluent` = "Effluent",
                                     `Treatment Influent` = 'Influent'))

#Filter to only include above min detection
cip3 <- subset(cip3, BelowMinReading!="Yes")

# rename antibiotic in cip3 
names(cip3)[names(cip3)=="NameDeterminandName"] <- "antibiotic"

#Add PNEC numbers 
cip3$EC10 <- if_else(cip3$antibiotic=="erythromycin", "2000",
                     if_else(cip3$antibiotic=="azithromycin", "390", 
                             if_else(cip3$antibiotic=="clarithromycin","2000",
                                     if_else(cip3$antibiotic=="ciprofloxacin","3.9",
                                             if_else(cip3$antibiotic=="Sulfamethoxazole","925",
                                                     if_else(cip3$antibiotic=="Trimethoprim","31.3",
                                                             if_else(cip3$antibiotic=="Sulfadiazine (Silvadene)","468.8",
                                                                     if_else(cip3$antibiotic=="Florfenicol","937.5", "31.25"))))))))


#Make EC10 numerica
cip3$EC10 <- as.numeric(cip3$EC10)

#Make PNEC
cip3$pnecr <- cip3$EC10/10

cip3$pnecr <- as.numeric(cip3$pnecr)

#Make RQ
cip3$RQ <- cip3$SampleValue / cip3$pnecr

#Add level of risk
cip3$risk <- if_else(cip3$RQ>=1, "high risk", if_else(cip3$RQ>=0.1, "medium risk", "low risk"))


#Generate summary table
cip3_rqsummary = cip3 %>% 
  group_by(antibiotic, SampleLocationName) %>%
  summarise(MedianRQ = median(RQ), MaxRQ = max(RQ)) %>%
  gather(Risk, RQ, -antibiotic, -SampleLocationName)

write.csv(cip3_rqsummary, "data/RQ_Summary_cip3.csv")

#Generate summary table high risk 
cip3_rqsummary_high <- subset(cip3_rqsummary, RQ>=1)

write.csv(cip3_rqsummary, "data/RQ_Summary_cip3_high.csv")

#order the Sample type to Inflow then Effluent, for plotting
cip3$SampleLocationName <- factor(cip3$SampleLocationName, levels=c("Influent", "Effluent"))


# add antibiotic clas
# cip3 <- merge(cip3, antibiotic_class, by = "antibiotic")
cip_class <- merge(cip3, antibiotic_class, by = "antibiotic", all.x = T)

## add name for Facet wrap labels
antibiotic.label_all <- c("Azithromycin", "Ciprofloxacin", "Clarithromycin", "Erythromycin", "Florfenicol", 
                          "Sulfadiazine (Silvadene)", "Sulfamethoxazole", "Trimethoprim")
names(antibiotic.label_all) <- c("azithromycin", "ciprofloxacin", "clarithromycin", "erythromycin",
                                 "Florfenicol", 
                                 "Sulfadiazine (Silvadene)", "Sulfamethoxazole", "Trimethoprim")


## make plot
cip3_plot <- ggplot(cip_class, aes(x=SampleLocationName, y=RQ)) +
  geom_hline(yintercept= 0.1, linetype= "longdash", linewidth = 1.2, colour = "darkred") + # the low line
  geom_hline(yintercept= 1, colour = "darkred", linewidth = 1.2) + # the high line
  geom_point(aes(colour = SampleLocationName), position = position_jitter(), alpha=0.3, size = 6, show.legend = FALSE)+
  facet_wrap(~antibiotic, scales = "free", labeller = labeller(antibiotic = antibiotic.label_all)) +
  ylab(bquote(log[10](RQ)))+
  labs(x="Wastewater Type", 
       title = "B) England and Wales Environmental Risk") +
  theme_bw()+
  scale_color_manual(values=met.brewer("Hiroshige", 2))+
  scale_fill_manual(values=met.brewer("Hiroshige", 2))+
  theme(legend.position="none",
        strip.background=element_rect(colour="black",
                                      fill="white"))+
  theme(text=element_text(size=26))+
  theme(axis.title = element_text(size = 30))+
  scale_y_log10(labels = scales::comma)

cip3_plot

######################### Patchwork ERA Plots Figure 4 ############################

uba_plot / cip3_plot

ggsave("figures/Figure 4.tiff", dpi = 300, height = 16, width = 12)

