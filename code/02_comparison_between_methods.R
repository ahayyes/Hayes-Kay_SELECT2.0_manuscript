######################### SELECT 2.0 Analysis ###########################

## this script will compare the PNECRs generated using the SELECT 1.0 and 
## SELECT 2.0 methods

####################### Load Libraries ###########################

library(tidyverse)
library(ggplot2)
library(reshape2)
library(MetBrewer)
library(patchwork)
library(plotrix)

##################### Confidence Interval ############################

conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}

################### Read in data and organise ########################


# read in file
df <- read.csv("data/results.csv", header=TRUE)

# make a dataframe with those antibiotics that include classes that we have data added to

# make PNECRs
df$select1.0_pnecr <- (df$select1.0/2)/10 

df$select2.0_pnecr <- df$select2.0/10 

## do fold hcnage
df$foldchange <- df$select1.0_pnecr / df$select2.0_pnecr

## read this out for inclusion in manuscript
write.csv(df, "data/Table 2.csv")

######## Make plot Figure 1 SELECT 2.0 EC10s boxplots per class #################

ggplot(df, aes(x=fct_reorder(class, select2.0_pnecr), y=select2.0_pnecr)) +
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, colour = "grey20")+
  geom_point(size = 6, alpha = 0.7, position = position_jitter(width = 0.1), colour = "grey20") +
  theme_bw() +
  labs(x = "Antibiotic Class",
       y = "log10(EC1 PNECR)") +
  ylab(bquote(log[10](`EC1 PNECR`)))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme(text=element_text(size=16))+
  theme(axis.title = element_text(size=20))+
  scale_y_log10() 

ggsave("figures/Figure 1.tiff", height = 8, width = 12, dpi = 300)

####################### Plot SELECT 1.0 vs SELECT 2.0 ###########################

## melt dataframe to make for ease of comparison

df_small <- df %>%
  select(c(antibiotic, select1.0_pnecr, select2.0_pnecr)) %>%
  pivot_longer(!antibiotic, names_to = "method", values_to = "pnecr")

# take class data from df
antibiotic_class <- df %>%
  select(c(antibiotic, class))

# merge them together
df_comp <- merge(df_small, antibiotic_class, by = "antibiotic")

fig2a <- ggplot(df_comp, aes(x=fct_reorder(class, pnecr), y=pnecr, colour = method)) +
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(size = 6, alpha = 0.9, position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0)) +
  theme_bw() +
  labs(x = "Antibiotic Class",
       y = "log10(PNECR)",
       title = "A) Comparison between SELECT 1.0 and SELECT 2.0 methods") +
  ylab(bquote(log[10](PNECR)))+
  scale_colour_manual(values=met.brewer("Hiroshige", 2), name = "Method", 
                      labels = c("SELECT 1.0", "SELECT 2.0"))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme(text=element_text(size=16))+
  theme(axis.title = element_text(size=20))+
  scale_y_log10(labels = scales::comma) 

fig2a

ggsave("figures/Figure 2a.tiff", height = 8, width = 12, dpi = 300)

############################# Bland Altman Analysis ######################

# can't have NAs in the dataframe 
# remove excess columns 
# want to use pnecRs not just the LOECs/EC10s

df_BA <- df %>%
  na.omit(df) %>%
  select(c(antibiotic, class, select1.0_pnecr, select2.0_pnecr))

## Manually make means and BA

df_BA$means <- (df_BA$select1.0_pnecr + df_BA$select2.0_pnecr) / 2

df_BA$diff <- df_BA$select1.0_pnecr - df_BA$select2.0_pnecr

# Average difference (aka the bias)
bias <- mean(df_BA$diff)

# Sample standard deviation
sd <- sd(df_BA$diff)

# Limits of agreement
upper_loa <- bias + 2 * sd
lower_loa <- bias - 2 * sd

#bias line (average difference) = grey 
#higher and lower limits of agreement are in dashed black
#zero line = no dash black

## Add confidence intervals
n <- nrow(df_BA)

# We want 95% confidence intervals
conf_interval <- 0.95

# Variance
var <- sd**2

# Standard error of the bias
se_bias <- sqrt(var / n)

# Standard error of the limits of agreement
se_loas <- sqrt(3 * var / n)

# Confidence intervals
confidence_intervals <- conf_int95(df_BA$means)
confidence_intervals_diff <- conf_int95(df_BA$diff)

fig2b <- ggplot(df_BA, aes(x = means, y = diff, colour = class))+
  geom_hline(yintercept = 0, col = "black")+ 
  geom_hline(yintercept = upper_loa, col = "black", lty = "dashed") +
  geom_hline(yintercept = bias, col = "azure4", lty = "dashed") +
  geom_hline(yintercept = lower_loa, col = "black", lty = "dashed") +
  geom_ribbon(aes(ymin=upper_loa+confidence_intervals_diff, ymax=upper_loa-confidence_intervals_diff, fill="azure2"),
              color=NA, alpha=0.1) +
  geom_ribbon(aes(ymin=lower_loa+confidence_intervals_diff, ymax=lower_loa-confidence_intervals_diff, fill="azure2"),
              color=NA, alpha=0.1) +
  geom_point(shape = 19, alpha = 0.9, size = 6) +
  scale_colour_manual(values=met.brewer("Hiroshige", 12), name = "Antibiotic\nClass")+
  labs(title = "B) Bland-Altman Analysis", 
       y = "Difference", 
       x = "Mean of Measurements")+
  theme_classic()  +
  theme(text=element_text(size=16),
        axis.title = element_text(size=20))+
  guides(fill = "none") 

fig2b

ggsave("figures/Bland Altman.tiff", dpi = 300, height = 12, width = 12)

####################### Patchwork the plots for figure 2######################

fig2a / fig2b

# write out dataframe to get PNECRs
ggsave("figures/Figure 2.tiff", dpi = 300, height = 12, width = 12)

####################### Comparison with meta-analysis data ######################

# need to load in data from the Murray et al., 2024 analysis paper
df_murray24 <- read.csv("data/results_murray2024.csv", header=TRUE)

# add column for 2024 study
df_murray24$study <- "2024"

unique(df_murray24$antibiotic)
# 16 antibiotics in this study

# add our data to this using select 2.0 data

# filter df_comp to have select 2.0 pnects
df_select2.0_pnec <- subset(df_comp, method == "select2.0_pnecr")

# add column for our study 
df_select2.0_pnec$study <- "2025"

# remove method column
df_select2.0_pnec <- df_select2.0_pnec %>% select(-c(method))

unique(df_select2.0_pnec$antibiotic)
# 32 antibiotics in this study

# filter our dataframe to only contain those antibiotics we can compare to
df_murray24_small <- filter(df_murray24, antibiotic %in% df_select2.0_pnec$antibiotic)

df_select2.0_pnec_small <- filter(df_select2.0_pnec, antibiotic %in% df_murray24$antibiotic)

# merge the two dataframes
df_2.0_meta_comp <- rbind(df_select2.0_pnec_small, df_murray24_small)

## graph these comparisons

# graph by antibiotic individually
ggplot(df_2.0_meta_comp, aes(x=fct_reorder(antibiotic, pnecr), y=pnecr, colour = study)) +
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(size = 6, alpha = 0.9, position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0)) +
  theme_bw() +
  labs(x = "Antibiotic",
       y = "log10(PNECR)",
       colour = "Study") +
  ylab(bquote(log[10](PNECR)))+
  scale_colour_manual(values=met.brewer("Hiroshige", 2), name = "Study",
                      labels = c("Murray et al.,\n2024", "This study\n2025"))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme(text=element_text(size=16))+
  theme(axis.title = element_text(size=20))+
  scale_y_log10(labels = scales::comma) 

ggsave("figures/Figure 3_sep.tiff", dpi = 300, height = 8, width = 12)




