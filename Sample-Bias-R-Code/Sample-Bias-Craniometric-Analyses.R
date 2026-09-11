######################################################
# Title: Sample Bias - Craniometric Analyses 
# Author: Rhian Dunn & Joe Hefner
# Last Update: 27 April 2026
# Updates:
# 1. 
# 2. 
# 3. 
# 4. 
# notes: 
#
#######################################################

# START WITH A CLEAN WORKING DIRECTORY 
rm(list=ls()) #will erase everything 

#######################################################

#LOAD LIBRARIES
library(dplyr) # For data partitioning 
library(tidyr) # For %>% 
library(psych) # For describeBy function 
library(mice)  # For data imputation
library(VIM) # For missing data analysis
library(corrplot) # For correlation plot 
library(candisc) # For Canonical Variate Analysis
library(ggplot2) # For plotting 
library(ggrepel) # For plotting 
library(reshape2) # For melting 
library(geosphere) # For GPS coordinate analysis
library(vegan) # For procrustes analysis
library(plotly) # For 3D plot of Mahalanobis Distances 
library(caret) # For data partitioning 
library(MASS) # Model fitting and StepAIC
library(klaR) # for stepwise variable selection
library(factoextra) # for scree plot and clustering
library(pROC) # for ROC curve plot

#######################################################

##IMPORTING DATA
data <- read.csv("data", header = TRUE) #replace "data" with name of data file
# View the first few rows of the data
head(data)

#######################################################
set.seed(1234)
#######################################################

#######################################################
#Analysis of Data Within and Between Collections
#######################################################

#Table of Contents
# 1- Data Preparation  (line 56) 
# 2- MANOVAs and CVA (line 248)
# 3- Mahalanobis Distances (line 510)
# 4- Classification (line 908)

#######################################################
# Filtering data
#######################################################
# Remove measurements not needed for analysis
reduced_mms <- data %>% dplyr::select(GOL, BNL, BBH, XCB, XFB, WFB, ZYB, AUB, ASB, 
                                      BPL, NPH, NLH, JUB, NLB, MAL, OBH, OBB, DKB, 
                                      WNB, ZMB, EKB, FRC, PAC, OCC, FOL, FOB, MOW, 
                                      UFBR, UFHT)
# Pull out just demographic data (first 37 columns)
demographic_data <- data[,c(1:37)]
# Combine back into same dataset
data <- cbind(demographic_data,reduced_mms)
#######################################################
# Plots to find outliers
#######################################################
# Boxplots for all measurement columns (38:66)
data %>%
  dplyr::select(all_of(38:66)) %>%
  pivot_longer(cols = everything(), names_to = "Measurement", values_to = "Value") %>%
  ggplot(aes(x = Measurement, y = Value)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots for Identifying Outliers") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot variables in paired scatterplots
pairs(~ GOL + BNL + BBH + XCB + XFB + WFB + ZYB + AUB + ASB + BPL + NPH + NLH + JUB + NLB + MAL + 
        OBH + OBB + DKB + WNB + ZMB + EKB + FRC + PAC + OCC + FOL + FOB + MOW + UFBR + UFHT, 
      data = data)
# Identify which MOW measurements are outliers
outlier_values<-boxplot.stats(data$MOW)$out # Replace with variable of interest 
outlier_values
#######################################################
# Data imputation 
#######################################################
# Perform multiple imputation
imputed_data <- mice(data[, 38:66], m = 5, method = 'pmm', seed = 123)
# Plot the convergence traceplots
plot(imputed_data)
# Get the completed dataset
complete_data <- complete(imputed_data)
#compile new dataset with demographic data and full measurements
clean_data <- cbind(demographic_data, complete_data)
write.csv(clean_data,file="Dunn_Dissertation_Imputated Dataset.csv")
# Check missing data patterns
md.pattern(data[, 38:66])
# Missing data pattern visualization with the VIM package
aggr(data[, 38:66], col = c("#18453B", "#CFD5D2"), numbers = TRUE, sortVars = TRUE,
     labels = names(data[, 38:66]), cex.axis = .7, gap = 3, 
     ylab = c("Missing data pattern", "Count"))
#######################################################
# Cooks Distance to find outliers
######################################################## 
# calculate GLM for cooks distance
mod <- glm(GOL + BNL + BBH + XCB + XFB + WFB + ZYB + AUB + ASB + BPL + NPH + NLH + JUB + NLB + MAL + 
             OBH + OBB + DKB + WNB + ZMB + EKB + FRC + PAC + OCC + FOL + FOB + MOW + UFBR + UFHT~Sex, data=clean_data)
cooksd <- cooks.distance(mod)
# plot cooks distance
plot(cooksd, pch=19, cex=1, main="Outlier detection using Cooks distance")  # plot cook's distance
abline(h = 10*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>10*mean(cooksd, na.rm=T),names(cooksd),""), col="black", adj=1,pos=3,offset=.5)  # add labels
# Find the observation(s) 10*mean
mean_cooksd <- mean(cooksd)
outliers_10x_mean <- which(cooksd > 10 * mean_cooksd)
outlier_data_10x <- clean_data[outliers_10x_mean, ]
print(outlier_data_10x) 
# Find the observation(s) with the largest Cook's D
outlier_index <- which.max(cooksd)
print(data[outlier_index, ])
# Identify top outlier 
car::outlierTest(mod) 
clean_data <- clean_data %>% filter(ID != "278")
#######################################################
# Descriptive stats
#######################################################
# Select measurement columns and calculate min, 25% quantile, mean, 75% quantile, and max for each
summary_stats <- clean_data %>%
  summarise(across(all_of(38:66), 
                   list(min = ~min(.x, na.rm = TRUE), 
                        q25 = ~quantile(.x, 0.25, na.rm = TRUE),  # 25th percentile
                        mean = ~mean(.x, na.rm = TRUE), 
                        q75 = ~quantile(.x, 0.75, na.rm = TRUE),  # 75th percentile
                        max = ~max(.x, na.rm = TRUE))))
# Pivot the result to long format for readability
summary_stats_long <- summary_stats %>%
  pivot_longer(cols = everything(), names_to = c("Measurement", "Statistic"), names_sep = "_") 
# View the long-format summary statistics
print(summary_stats_long)

# Descriptive statistics (by Collection)
desc.data<-describeBy(clean_data,clean_data$Collection)
write.csv(desc.data$JAWDHSC,file="desc.data.jawdhsc.csv")
write.csv(desc.data$MLOC,file="desc.data.mloc.csv")
write.csv(desc.data$MMDSC,file="desc.data.mmdsc.csv")
write.csv(desc.data$`MSUFAL DSC`,file="desc.data.msufal.csv")
write.csv(desc.data$`STAFS DSC`,file="desc.data.stafs.csv")
write.csv(desc.data$TXSTDSC,file="desc.data.txstdsc.csv")
write.csv(desc.data$`UTK DSC`,file="desc.data.utk.csv")
write.csv(desc.data$`WMed STARS`,file="desc.data.wmed.csv")
#######################################################
# T-tests of measurements by Collection with Bonferroni Correction
#######################################################
# List of variables left:
# GOL, BNL, BBH, XCB, XFB, WFB, ZYB, AUB, ASB,BPL, NPH, NLH, JUB, NLB, MAL, 
# OBH, OBB, DKB, WNB, ZMB, EKB, FRC, PAC, OCC, FOL, FOB, MOW, UFBR, UFHT

#Subset of Males for test 
m_clean_data <- clean_data[clean_data$Sex == "M", ]
pairwise.t.test(m_clean_data$GOL, m_clean_data$Collection, p.adjust.method="bonferroni", print =T)
pairwise.t.test(m_clean_data$BNL, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$BBH, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$XCB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$XFB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$WFB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$ZYB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$AUB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$ASB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$BPL, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$NPH, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$NLH, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$JUB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$NLB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$MAL, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$OBH, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$OBB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$DKB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$WNB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$ZMB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$EKB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$FRC, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$OCC, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$FOL, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$FOB, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$MOW, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$UFBR, m_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(m_clean_data$UFHT, m_clean_data$Collection, p.adjust.method="bonferroni")

#Subset of Females for test
f_clean_data <- clean_data[clean_data$Sex == "F", ]
pairwise.t.test(f_clean_data$GOL, f_clean_data$Collection, p.adjust.method="bonferroni", print =T)
pairwise.t.test(f_clean_data$BNL, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$BBH, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$XCB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$XFB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$WFB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$ZYB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$AUB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$ASB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$BPL, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$NPH, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$NLH, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$JUB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$NLB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$MAL, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$OBH, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$OBB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$DKB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$WNB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$ZMB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$EKB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$FRC, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$OCC, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$FOL, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$FOB, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$MOW, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$UFBR, f_clean_data$Collection, p.adjust.method="bonferroni")
pairwise.t.test(f_clean_data$UFHT, f_clean_data$Collection, p.adjust.method="bonferroni")
#######################################################
# Center and scale data 
#######################################################
#Pull out measurements 
measurements <- clean_data[,38:66]
#center and scale data
center_scale<-function(x){
  scale(x,scale = FALSE)
}
#center means
measurements<-center_scale(measurements)
measurements<-as.matrix(measurements)
#Remake demo data as outlier was removed
demographic_data <- clean_data[,1:37]
scaled_data <- cbind(demographic_data,measurements)
#######################################################



#######################################################
# Prior MANOVA and CVA - All data by Collection
#######################################################
# One way MANOVA to set the scene
prior_manova_model <- manova(as.matrix(scaled_data[, 38:66]) ~ Collection, data = scaled_data)
# Summary of MANOVA results
summary(prior_manova_model, test = "Pillai")
# Univariate ANOVA for each dependent variable after significant MANOVA
summary.aov(prior_manova_model)
# For post-hoc pairwise comparisons (Tukey HSD) on a specific dependent variable, for example, the first one
prior_aov_model <- aov(DKB ~ Collection, data = scaled_data) 
TukeyHSD(prior_aov_model)
####
# Perform canonical variates analysis on the MANOVA model for 'Collection'
prior_cva_collection <- candisc(prior_manova_model, term = "Collection")
summary(prior_cva_collection)
# Extract the canonical scores
prior_cva_collection_scores <- as.data.frame(prior_cva_collection$scores)
# Calculate the means for Can1 and Can2 for each Collection group
prior_cva_collection_group_means <- prior_cva_collection_scores %>%
  group_by(Collection) %>%
  summarise(Mean_Can1 = mean(Can1, na.rm = TRUE),
            Mean_Can2 = mean(Can2, na.rm = TRUE),
            Mean_Can3 = mean(Can3, na.rm = TRUE),
            Mean_Can4 = mean(Can4, na.rm = TRUE),
            Mean_Can5 = mean(Can5, na.rm = TRUE),
            Mean_Can6 = mean(Can6, na.rm = TRUE))
# Plot the Canonical Variates and overlay group means
ggplot(prior_cva_collection_group_means, aes(x = Mean_Can1, y = Mean_Can2)) +
  geom_point(size = 4, shape = 21, fill = "white") +
  geom_text(aes(label = Collection), vjust = -1, color = "black", size = 5) +  # Label the group means
  theme_classic() +
  labs(title = "Canonical Variates Plot for Collection",
       x = "Canonical Variate 1 (34.4%)",
       y = "Canonical Variate 2 (22.0%)") 
#######################################################
# MANOVA1 and CVA - All data by Collection and Sex 
#######################################################
# Perform MANOVA with selected variables
manova_model <- manova(as.matrix(scaled_data[, 38:66]) ~ Collection+Sex, data = scaled_data)
# Summary of MANOVA results
summary(manova_model, test = "Pillai")
# Univariate ANOVA for each dependent variable after significant MANOVA
summary.aov(manova_model)
# For post-hoc pairwise comparisons (Tukey HSD) on a specific dependent variable, for example, the first one
aov_model <- aov(DKB ~ Collection, data = scaled_data) # Test each measurement (listed below)
TukeyHSD(aov_model)
####
# Perform canonical variates analysis on the MANOVA model for 'Collection'
cva_collection <- candisc(manova_model, term = "Collection")
summary(cva_collection)
# Extract the canonical scores
cva_collection_scores <- as.data.frame(cva_collection$scores)
# Calculate the means for Can1 and Can2 for each Collection group
cva_collection_group_means <- cva_collection_scores %>%
  group_by(Collection) %>%
  summarise(Mean_Can1 = mean(Can1, na.rm = TRUE),
            Mean_Can2 = mean(Can2, na.rm = TRUE),
            Mean_Can3 = mean(Can3, na.rm = TRUE),
            Mean_Can4 = mean(Can4, na.rm = TRUE),
            Mean_Can5 = mean(Can5, na.rm = TRUE),
            Mean_Can6 = mean(Can6, na.rm = TRUE))
# Plot the Canonical Variates and overlay group means
ggplot(cva_collection_group_means, aes(x = Mean_Can1, y = Mean_Can2)) +
  geom_point(size = 4, shape = 21, fill = "white") +
  geom_text(aes(label = Collection), vjust = -1, color = "black", size = 5) +  # Label the group means
  theme_classic() +
  labs(title = "Canonical Variates Plot for Collection (Controlling for Sex)",
       x = "Canonical Variate 1 (34.4%)",
       y = "Canonical Variate 2 (22.1%)") 
#######################################################
# MANOVA2 and CVA - All data by Collection, Population, and Sex
#######################################################
# Perform MANOVA with selected variables
manova_model2 <- manova(as.matrix(scaled_data[, 38:66]) ~ Collection+Sex+Population.Short, data = scaled_data)
# Summary of MANOVA results
summary(manova_model2, test = "Pillai")
# Univariate ANOVA for each dependent variable after significant MANOVA
summary.aov(manova_model2)
# For post-hoc pairwise comparisons (Tukey HSD) on a specific dependent variable, for example, the first one
aov_model2 <- aov(DKB ~ Collection, data = scaled_data) 
TukeyHSD(aov_model2)
####
# Perform canonical variates analysis on the MANOVA model for 'Collection'
cva_collection2 <- candisc(manova_model2, term = "Collection")
summary(cva_collection2)
# Extract the canonical scores
cva_collection_scores2 <- as.data.frame(cva_collection2$scores)
# Calculate the means for Can1 and Can2 for each Collection group
cva_collection_group_means2 <- cva_collection_scores2 %>%
  group_by(Collection) %>%
  summarise(Mean_Can1 = mean(Can1, na.rm = TRUE),
            Mean_Can2 = mean(Can2, na.rm = TRUE),
            Mean_Can3 = mean(Can3, na.rm = TRUE),
            Mean_Can4 = mean(Can4, na.rm = TRUE),
            Mean_Can5 = mean(Can5, na.rm = TRUE),
            Mean_Can6 = mean(Can6, na.rm = TRUE))
# Plot the Canonical Variates and overlay group means
ggplot(cva_collection_group_means2, aes(x = Mean_Can1, y = Mean_Can2)) +
  geom_point(size = 4, shape = 21, fill = "white") +
  geom_text(aes(label = Collection), vjust = -1, color = "black", size = 5) +  # Label the group means
  theme_classic() +
  labs(title = "Canonical Variates Plot for Collection (Controlling for Sex and Population)",
       x = "Canonical Variate 1 (37.7%)",
       y = "Canonical Variate 2 (20.5%)") 
#######################################################
# MANCOVA1 and CVA - All data by Collection, Population, Sex, and Age
#######################################################
# Perform MANCOVA with selected variables
mancova_model <- manova(as.matrix(scaled_data[, 38:66]) ~ Collection+Sex+Population.Short+Age, data = scaled_data)
# Summary of MANCOVA results
summary(mancova_model, test = "Pillai")
# Univariate ANOVA for each dependent variable after significant MANCOVA
summary.aov(mancova_model)
# For post-hoc pairwise comparisons (Tukey HSD) on a specific dependent variable, for example, the first one
aov_model3 <- aov(DKB ~ Collection, data = scaled_data) 
TukeyHSD(aov_model3)
####
# Perform canonical variates analysis on the MANCOVA model for 'Collection'
cva_collection3 <- candisc(mancova_model, term = "Collection")
summary(cva_collection3)
# Extract the canonical scores
cva_collection_scores3 <- as.data.frame(cva_collection3$scores)
# Calculate the means for Can1 and Can2 for each Collection group
cva_collection_group_means3 <- cva_collection_scores3 %>%
  group_by(Collection) %>%
  summarise(Mean_Can1 = mean(Can1, na.rm = TRUE),
            Mean_Can2 = mean(Can2, na.rm = TRUE),
            Mean_Can3 = mean(Can3, na.rm = TRUE),
            Mean_Can4 = mean(Can4, na.rm = TRUE),
            Mean_Can5 = mean(Can5, na.rm = TRUE),
            Mean_Can6 = mean(Can6, na.rm = TRUE))
# Plot the Canonical Variates and overlay group means
ggplot(cva_collection_group_means3, aes(x = Mean_Can1, y = Mean_Can2)) +
  geom_point(size = 4, shape = 21, fill = "white") +
  geom_text(aes(label = Collection), vjust = -1, color = "black", size = 5) +  # Label the group means
  theme_classic() +
  labs(title = "Canonical Variates Plot for Collection (Controlling for Sex, Population, and Age)",
       x = "Canonical Variate 1 (39.4%)",
       y = "Canonical Variate 2 (20.2%)") 
#######################################################
# MANCOVA2 and CVA - All data by Collection, Population, Sex, Age, and Birth Year
#######################################################
# Perform MANCOVA with selected variables
mancova_model2 <- manova(as.matrix(scaled_data[, 38:66]) ~ Collection+Sex+Population.Short+Age+BirthYear, data = scaled_data)
# Summary of MANCOVA results
summary(mancova_model2, test = "Pillai")
# Univariate ANOVA for each dependent variable after significant MANCOVA
summary.aov(mancova_model2)
# For post-hoc pairwise comparisons (Tukey HSD) on a specific dependent variable, for example, the first one
aov_model4 <- aov(DKB ~ Collection, data = scaled_data) 
TukeyHSD(aov_model4)
####
# Perform canonical variates analysis on the MANCOVA model for 'Collection'
cva_collection4 <- candisc(mancova_model2, term = "Collection")
summary(cva_collection4)
# Extract the canonical scores
cva_collection_scores4 <- as.data.frame(cva_collection4$scores)
# Calculate the means for Can1 and Can2 for each Collection group
cva_collection_group_means4 <- cva_collection_scores4 %>%
  group_by(Collection) %>%
  summarise(Mean_Can1 = mean(Can1, na.rm = TRUE),
            Mean_Can2 = mean(Can2, na.rm = TRUE),
            Mean_Can3 = mean(Can3, na.rm = TRUE),
            Mean_Can4 = mean(Can4, na.rm = TRUE),
            Mean_Can5 = mean(Can5, na.rm = TRUE),
            Mean_Can6 = mean(Can6, na.rm = TRUE))
# Plot the Canonical Variates and overlay group means
ggplot(cva_collection_group_means4, aes(x = Mean_Can1, y = Mean_Can2)) +
  geom_point(size = 4, shape = 21, fill = "white") +
  geom_text(aes(label = Collection), vjust = -1, color = "black", size = 5) +  # Label the group means
  theme_classic() +
  labs(title = "Canonical Variates Plot for Collection (Controlling for Sex, Population, Age, and Birth Year)",
       x = "Canonical Variate 1 (29.9%)",
       y = "Canonical Variate 2 (23.7%)") 
#######################################################
# MANCOVA3 and CVA - All data by Collection, Population, Sex, Age, Birth Year, and State of Birth
#######################################################
# Perform MANCOVA with selected variables
mancova_model3 <- manova(as.matrix(scaled_data[, 38:66]) ~ Collection+Sex+Population.Short+Age+BirthYear+Birthplace_State, data = scaled_data)
# Summary of MANCOVA results
summary(mancova_model3, test = "Pillai")
# Univariate ANOVA for each dependent variable after significant MANCOVA
summary.aov(mancova_model3)
# For post-hoc pairwise comparisons (Tukey HSD) on a specific dependent variable, for example, the first one
aov_model5 <- aov(DKB ~ Collection, data = scaled_data) # Test each measurement (listed below)
TukeyHSD(aov_model5)
####
# Perform canonical variates analysis on the MANCOVA model for 'Collection'
cva_collection5 <- candisc(mancova_model3, term = "Collection")
summary(cva_collection5)
# Extract the canonical scores
cva_collection_scores5 <- as.data.frame(cva_collection5$scores)
# Calculate the means for Can1 and Can2 for each Collection group
cva_collection_group_means5 <- cva_collection_scores5 %>%
  group_by(Collection) %>%
  summarise(Mean_Can1 = mean(Can1, na.rm = TRUE),
            Mean_Can2 = mean(Can2, na.rm = TRUE),
            Mean_Can3 = mean(Can3, na.rm = TRUE),
            Mean_Can4 = mean(Can4, na.rm = TRUE),
            Mean_Can5 = mean(Can5, na.rm = TRUE),
            Mean_Can6 = mean(Can6, na.rm = TRUE))
# Plot the Canonical Variates and overlay group means
ggplot(cva_collection_group_means5, aes(x = Mean_Can1, y = Mean_Can2)) +
  geom_point(size = 4, shape = 21, fill = "white") +
  geom_text(aes(label = Collection), vjust = -1, color = "black", size = 5) +  # Label the group means
  theme_classic() +
  labs(title = "Canonical Variates Plot for Collection (Controlling for Sex, Population, Age, Birth Year, and State of Birth)",
       x = "Canonical Variate 1 (35.6%)",
       y = "Canonical Variate 2 (21.8%)") 
#######################################################
# MANCOVA4 and CVA - Reduced MMDSC data by Collection, Population, Sex, Age, Birth Year, and State of Birth
#######################################################
new_scaled_data <- scaled_data
#Remove transfers of unidentified from SOM
new_scaled_data <- new_scaled_data %>% filter(DonationType != "School of Medicine")
#Remove doc doe
new_scaled_data <- new_scaled_data %>% 
  filter(!Donor_ID %in% c("18", "74", "120", "124", "126", "129", "153"))
#Remove donations without donation forms
new_scaled_data <- new_scaled_data %>% 
  filter(!Donor_ID %in% c("19", "94", "103", "109", "111", "112", "151", 
                          "162", "164", "167", "168", "169", "171", "177", 
                          "183", "188", "202", "211", "215", "218", "220", 
                          "222", "223", "224", "226", "227", "228", "229", 
                          "230", "231", "238"))
# Perform MANCOVA with selected variables
new_mancova_model <- manova(as.matrix(new_scaled_data[, 38:66]) ~ Collection+Sex+Population.Short+Age+BirthYear+Birthplace_State, data = new_scaled_data)
# Summary of MANOVA results
summary(new_mancova_model, test = "Pillai")
# Univariate ANOVA for each dependent variable after significant MANOVA
summary.aov(new_mancova_model)
# For post-hoc pairwise comparisons (Tukey HSD) on a specific dependent variable, for example, the first one
new_aov_model <- aov(XCB ~ Collection, data = new_scaled_data) # Test each measurement 
TukeyHSD(new_aov_model)
####
# Perform canonical variates analysis on the MANOVA model for 'Collection'
new_cva_collection <- candisc(new_mancova_model, term = "Collection")
summary(new_cva_collection)
# Extract the canonical scores
new_cva_collection_scores <- as.data.frame(new_cva_collection$scores)
# Calculate the means for Can1 and Can2 for each Collection group
new_cva_collection_group_means <- new_cva_collection_scores %>%
  group_by(Collection) %>%
  summarise(Mean_Can1 = mean(Can1, na.rm = TRUE),
            Mean_Can2 = mean(Can2, na.rm = TRUE),
            Mean_Can3 = mean(Can3, na.rm = TRUE),
            Mean_Can4 = mean(Can4, na.rm = TRUE),
            Mean_Can5 = mean(Can5, na.rm = TRUE),
            Mean_Can6 = mean(Can6, na.rm = TRUE))
# Plot the Canonical Variates and overlay group means
ggplot(new_cva_collection_group_means, aes(x = Mean_Can1, y = Mean_Can2)) +
  geom_point(size = 4, shape = 21, fill = "white") +
  geom_text(aes(label = Collection), vjust = -1, color = "black", size = 5) +  # Label the group means
  theme_classic() +
  labs(title = "Canonical Variates Plot with Modified MMDSC Data for Collection (Controlling for Sex, Population, Age, Birth Year, and State of Birth)",
       x = "Canonical Variate 1 (35.6%)",
       y = "Canonical Variate 2 (22.5%)") 
#######################################################



#######################################################
# Mahalanobis Distances - All data by Collection 
#######################################################
# specific measurement columns
measurement_columns <- colnames(scaled_data)[38:66]
# Calculate the mean vector for each collection
group_means <- scaled_data %>%
  group_by(Collection) %>%
  summarise(across(all_of(measurement_columns), \(x) mean(x, na.rm = TRUE)))
# Calculate the pooled covariance matrix (based on group means, not raw data)
cov_matrix <- cov(scaled_data %>% dplyr::select(all_of(measurement_columns)), use = "complete.obs")
# Create an empty matrix to store Mahalanobis distances
mahalanobis_distances <- matrix(NA, nrow = nrow(group_means), ncol = nrow(group_means))
# Add labels to the distance matrix (collection names as row/column labels)
rownames(mahalanobis_distances) <- group_means$Collection
colnames(mahalanobis_distances) <- group_means$Collection
# Calculate Mahalanobis distances between group means
for (i in 1:nrow(group_means)) {
  for (j in 1:nrow(group_means)) {
    # Convert to numeric vectors and calculate Mahalanobis distance
    mahalanobis_distances[i, j] <- mahalanobis(
      as.numeric(group_means[i, -1]),  # Group i mean vector (excluding 'Collection' column)
      as.numeric(group_means[j, -1]),  # Group j mean vector (excluding 'Collection' column)
      cov_matrix  # Covariance matrix
    )
  }
}
# Display the Mahalanobis distance matrix
mahalanobis_distances
# Copy for plot 
mahalanobis_distances_copy <- mahalanobis_distances
# Set the lower triangular and diagonal values to NA
mahalanobis_distances_copy[lower.tri(mahalanobis_distances_copy, diag = TRUE)] <- NA
# Convert the distance matrix into a data frame suitable for ggplot2
melted_mahalanobis_distances <- melt(as.matrix(mahalanobis_distances_copy), na.rm = TRUE)  
# Create a heatmap with the distance values displayed in each cell with only the upper triangular part
ggplot(data = melted_mahalanobis_distances, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +                                                
  geom_text(aes(label = round(value, 2)), color = "black") +  
  scale_fill_gradient(low = "white", high = "#18453B") +            
  theme_minimal() +                                            
  labs(title = "Heatmap of Mahalanobis Distances for Collection Groups", 
       x = "Groups", y = "Groups", fill = "Mahalanobis Distance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))     
####
# Perform Multidimensional Scaling (MDS)
mds_coords <- cmdscale(mahalanobis_distances, k = 2) 
# Plot the MDS results
plot(mds_coords, 
     type = "n", 
     xlab = "MDS Dimension 1", 
     ylab = "MDS Dimension 2", 
     main = "MDS Plot of Mahalanobis Distances"); text(mds_coords, labels = rownames(mahalanobis_distances), cex = 0.8)
# Perform hierarchical clustering
hclust_res <- hclust(as.dist(mahalanobis_distances))
# Plot the dendrogram
plot(hclust_res, main = "Dendrogram of Mahalanobis Distances by Collection")
# Perform MDS in 3D
dist_matrix <- dist(mahalanobis_distances)
mds_result <- cmdscale(dist_matrix, k = 3)  
# Create a dataframe for MDS coordinates
mds_coords_3d <- data.frame(mds_result)
# Rename columns for clarity
colnames(mds_coords_3d) <- c("x", "y", "z")  
mds_coords_3d$Collection <- group_means$Collection  
# Create a 3D scatter plot with labels
plot_ly(mds_coords_3d, 
        x = ~x, y = ~y, z = ~z, 
        type = "scatter3d", 
        mode = "markers",  
        marker = list(size = 5),  
        color = ~Collection,  
        text = ~Collection,  
        hoverinfo = 'text') %>%  
  layout(title = "3D MDS Plot of Mahalanobis Distances by Collection",
         scene = list(xaxis = list(title = "MDS Dimension 1"),
                      yaxis = list(title = "MDS Dimension 2"),
                      zaxis = list(title = "MDS Dimension 3")))
#######################################################
# Mahalanobis Distances - All data by Collection GPS Location
#######################################################
# Create a lookup table for GPS coordinates based on Collection names
gps_coordinates <- data.frame(
  Collection = c("JAWDHSC", "MLOC","MMDSC","MSUFAL DSC","STAFS DSC","TXSTDSC","UTK DSC","WMed STARS"),
  Latitude = c(35.3090,21.2964,35.0843,42.7251,30.7136,29.8884,35.9544,42.2896),
  Longitude = c(-83.1864,-157.8637,-106.6198,-84.4791,-95.5473,-97.9384,-83.9295,-85.5795)
)
gps_coordinates$Latitude <- as.numeric(gps_coordinates$Latitude)
gps_coordinates$Longitude <- as.numeric(gps_coordinates$Longitude)
#  Find centroids for analysis 
collection_centroids <- gps_coordinates %>%
  group_by(Collection) %>%
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude))
# Calculate the pairwise geodesic distances (Haversine)
geo_distance_matrix <- distm(collection_centroids[, c("Longitude", "Latitude")], fun = distHaversine)
# Assign row and column names as collection names
rownames(geo_distance_matrix) <- collection_centroids$Collection
colnames(geo_distance_matrix) <- collection_centroids$Collection
# View the geodesic distance matrix
print(geo_distance_matrix)
####
geo_distance_matrix_copy <- geo_distance_matrix
# Set the diagonal and lower triangular part to NA (to remove them from the plot)
geo_distance_matrix_copy[lower.tri(geo_distance_matrix_copy, diag = TRUE)] <- NA
# Convert the distance matrix into a data frame suitable for ggplot2
melted_geo_distance_matrix <- melt(as.matrix(geo_distance_matrix_copy), na.rm = TRUE)  
# Create a heatmap with the distance values displayed in each cell with only the upper triangular part
ggplot(data = melted_geo_distance_matrix, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +                                                
  geom_text(aes(label = round(value, 2)), color = "black") +   
  scale_fill_gradient(low = "white", high = "#18453B") +            
  theme_minimal() +                                            
  labs(title = "Heatmap of Geographic Distances for Physical Collection Locations", 
       x = "Groups", y = "Groups", fill = "Distance (Meters)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))     
####
# Perform Multidimensional Scaling (MDS)
geo_mds_coords <- cmdscale(geo_distance_matrix, k = 2)  
# Plot the MDS results
plot(geo_mds_coords, 
     type = "n", 
     xlab = "MDS Dimension 1", 
     ylab = "MDS Dimension 2", 
     main = "MDS Plot of Distances by GPS Coordinates"); text(geo_mds_coords, labels = rownames(geo_distance_matrix), cex = 0.8)
# Perform hierarchical clustering
geo_hclust_res <- hclust(as.dist(geo_distance_matrix))
# Plot the dendrogram
plot(geo_hclust_res, main = "Dendrogram of Mahalanobis Distances Based on GPS Coordinates")
# Perform MDS in 3D
geo_dist_matrix <- dist(geo_distance_matrix)
geo_mds_result <- cmdscale(geo_dist_matrix, k = 3) 
# Create a dataframe for MDS coordinates
geo_mds_coords_3d <- data.frame(geo_mds_result)
# Rename columns for clarity
colnames(geo_mds_coords_3d) <- c("x", "y", "z")  
geo_mds_coords_3d$Collection <- collection_centroids$Collection  
# Create a 3D scatter plot with labels
plot_ly(geo_mds_coords_3d, 
        x = ~x, y = ~y, z = ~z, 
        type = "scatter3d", 
        mode = "markers",  
        marker = list(size = 5),  
        color = ~Collection,  
        text = ~Collection,  
        hoverinfo = 'text') %>%  
  layout(title = "3D MDS Plot of Mahalanobis Distances by Collection",
         scene = list(xaxis = list(title = "MDS Dimension 1"),
                      yaxis = list(title = "MDS Dimension 2"),
                      zaxis = list(title = "MDS Dimension 3")))
#######################################################
# Plot Collection and GPS Collection Together
#######################################################
# need mds_coords and geo_mds_coords from above
# Create a vector of collection labels
collection_labels <- gps_coordinates$Collection
# Align geographic coordinates to craniometric coordinates
procrustes_result <- procrustes(mds_coords, geo_mds_coords, scale = TRUE)
# Extract transformed points
transformed_geo <- procrustes_result$Yrot  # The transformed geographic coordinates
original_mah <- procrustes_result$X        # The original craniometric coordinates
# Plot both sets of points together
plot(original_mah, col = "blue", pch = 16, xlim = range(c(original_mah[,1], transformed_geo[,1])), 
     ylim = range(c(original_mah[,2], transformed_geo[,2])), xlab = "Dimension 1", ylab = "Dimension 2",
     main = "Procrustes-Transformed Coordinates")
points(transformed_geo, col = "red", pch = 17)
# Add collection labels to the points
text(original_mah, labels = collection_labels, pos = 3, cex = 0.8, col = "blue")
text(transformed_geo, labels = collection_labels, pos = 3, cex = 0.8, col = "red")
# Add a legend to differentiate between craniometric and geographic data
legend("topleft", legend = c("Craniometric Distances", "Geographic Distances"), 
       col = c("blue", "red"), pch = c(16, 17))
# Create a data frame with collection labels and point type
plot_data <- data.frame(
  LD1 = c(original_mah[, 1], transformed_geo[, 1]),
  LD2 = c(original_mah[, 2], transformed_geo[, 2]),
  LD1_end = c(transformed_geo[, 1], transformed_geo[, 1]),  
  LD2_end = c(transformed_geo[, 2], transformed_geo[, 2]),
  Collection = rep(collection_labels, 2),
  Type = rep(c("Craniometric Distance", "Geographic Distance"), each = length(collection_labels))
)
# Plot the data with lines
ggplot(plot_data) +
  geom_point(aes(x = LD1, y = LD2, color = Type, shape = Type), size = 3) +
  geom_segment(aes(x = LD1, y = LD2, xend = LD1_end, yend = LD2_end),
               linetype = "dotted", show.legend = FALSE) +
  geom_text_repel(aes(x = LD1, y = LD2, label = Collection, color = Type),
                  size = 3) +
  scale_color_manual(values = c("#18453B", "#7F988E")) +  
  scale_shape_manual(values = c(16, 17), labels = c("Craniometric Distance", "Geographic Distance")) +
  labs(title = "Procrustes-Transformed Plot (Craniometric vs. Geographic Distance)",
       x = "Dimension 1", y = "Dimension 2", shape = "Data Type", color = "Data Type") +
  theme_classic() +
  theme(legend.position = "right")
#######################################################
# Mahalanobis Distances - All data by Collection, Sex, and Pop
#######################################################
#Reduce to only groups with >4 individuals per 
collpopsex_scaled_data <- scaled_data %>% filter(!is.na(Collection.Pop.Sex.Red2) & Collection.Pop.Sex.Red2 != "")
# specific measurement columns
collpopsex_measurement_columns <- colnames(collpopsex_scaled_data)[38:66]
# Calculate the mean vector for each collection
collpopsex_group_means <- collpopsex_scaled_data %>%
  group_by(Collection.Pop.Sex.Red2) %>% 
  summarise(across(all_of(collpopsex_measurement_columns), \(x) mean(x, na.rm = TRUE)))
# Calculate the pooled covariance matrix (based on group means, not raw data)
collpopsex_cov_matrix <- cov(collpopsex_scaled_data %>% dplyr::select(all_of(collpopsex_measurement_columns)), use = "complete.obs")
# Create an empty matrix to store Mahalanobis distances
collpopsex_mahalanobis_distances <- matrix(NA, nrow = nrow(collpopsex_group_means), ncol = nrow(collpopsex_group_means))
# Add labels to the distance matrix (collection names as row/column labels)
rownames(collpopsex_mahalanobis_distances) <- collpopsex_group_means$Collection.Pop.Sex.Red2
colnames(collpopsex_mahalanobis_distances) <- collpopsex_group_means$Collection.Pop.Sex.Red2
# Calculate Mahalanobis distances between group means
for (i in 1:nrow(collpopsex_group_means)) {
  for (j in 1:nrow(collpopsex_group_means)) {
    # Convert to numeric vectors and calculate Mahalanobis distance
    collpopsex_mahalanobis_distances[i, j] <- mahalanobis(
      as.numeric(collpopsex_group_means[i, -1]),  # Group i mean vector (excluding 'Collection' column)
      as.numeric(collpopsex_group_means[j, -1]),  # Group j mean vector (excluding 'Collection' column)
      collpopsex_cov_matrix  # Covariance matrix
    )
  }
}
# Display the Mahalanobis distance matrix
collpopsex_mahalanobis_distances
write.csv(collpopsex_mahalanobis_distances, file = "Collection Pop Sex MD Matrix.csv")
####
collpopsex_mahalanobis_distances_copy <- collpopsex_mahalanobis_distances
# Set the diagonal and lower triangular part to NA (to remove them from the plot)
collpopsex_mahalanobis_distances_copy[lower.tri(collpopsex_mahalanobis_distances_copy, diag = TRUE)] <- NA
# Convert the distance matrix into a data frame suitable for ggplot2
melted_collpopsex_mahalanobis_distances <- melt(as.matrix(collpopsex_mahalanobis_distances_copy), na.rm = TRUE) 
# Create a heatmap with the distance values displayed in each cell with only the upper triangular part
ggplot(data = melted_collpopsex_mahalanobis_distances, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +                                                
  geom_text(aes(label = round(value, 2)), color = "black") +   
  scale_fill_gradient(low = "white", high = "#18453B") +           
  theme_minimal() +                                            
  labs(title = "Heatmap of Mahalanobis Distances for Collection, Population, and Sex Specific Groups", 
       x = "Groups", y = "Groups", fill = "Mahalanobis Distance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))     
####
# Perform Multidimensional Scaling (MDS)
collpopsex_mds_coords <- cmdscale(collpopsex_mahalanobis_distances, k = 2)  
# Plot the MDS results
plot(collpopsex_mds_coords, 
     type = "n", 
     xlab = "MDS Dimension 1", 
     ylab = "MDS Dimension 2", 
     main = "MDS Plot of Mahalanobis Distances"); text(collpopsex_mds_coords, labels = rownames(collpopsex_mahalanobis_distances), cex = 0.8)
# Perform hierarchical clustering
collpopsex_hclust_res <- hclust(as.dist(collpopsex_mahalanobis_distances))
# Plot the dendrogram
plot(collpopsex_hclust_res, main = "Dendrogram of Mahalanobis Distances")
# Perform MDS in 3D
collpopsex_dist_matrix <- dist(collpopsex_mahalanobis_distances)
collpopsex_mds_result <- cmdscale(collpopsex_dist_matrix, k = 3)
# Create a dataframe for MDS coordinates
collpopsex_mds_coords_3d <- data.frame(collpopsex_mds_result)
# Rename columns for clarity
colnames(collpopsex_mds_coords_3d) <- c("x", "y", "z")  
collpopsex_mds_coords_3d$Collection.Pop.Sex.Red2 <- collpopsex_group_means$Collection.Pop.Sex.Red2  
# Create a 3D scatter plot with labels
plot_ly(collpopsex_mds_coords_3d, 
        x = ~x, y = ~y, z = ~z, 
        type = "scatter3d", 
        mode = "markers",  
        marker = list(size = 5),  
        color = ~Collection.Pop.Sex.Red2,  
        text = ~Collection.Pop.Sex.Red2,  
        hoverinfo = 'text') %>%  
  layout(title = "3D MDS Plot of Mahalanobis Distances by Collection",
         scene = list(xaxis = list(title = "MDS Dimension 1"),
                      yaxis = list(title = "MDS Dimension 2"),
                      zaxis = list(title = "MDS Dimension 3")))
#######################################################
# Mahalanobis Distances - All European Individuals by Collection
#######################################################
# Pull out only European individuals
eur_scaled_data <- scaled_data[scaled_data$Population.Short == "W", ]
# specific measurement columns
eur_measurement_columns <- colnames(eur_scaled_data)[38:66]
# Calculate the mean vector for each collection
eur_group_means <- eur_scaled_data %>%
  group_by(Collection) %>%
  summarise(across(all_of(eur_measurement_columns), \(x) mean(x, na.rm = TRUE)))
# Calculate the pooled covariance matrix (based on group means, not raw data)
eur_cov_matrix <- cov(eur_scaled_data %>% dplyr::select(all_of(eur_measurement_columns)), use = "complete.obs")
# Create an empty matrix to store Mahalanobis distances
eur_mahalanobis_distances <- matrix(NA, nrow = nrow(eur_group_means), ncol = nrow(eur_group_means))
# Add labels to the distance matrix (collection names as row/column labels)
rownames(eur_mahalanobis_distances) <- eur_group_means$Collection
colnames(eur_mahalanobis_distances) <- eur_group_means$Collection
# Calculate Mahalanobis distances between group means
for (i in 1:nrow(eur_group_means)) {
  for (j in 1:nrow(eur_group_means)) {
    # Convert to numeric vectors and calculate Mahalanobis distance
    eur_mahalanobis_distances[i, j] <- mahalanobis(
      as.numeric(eur_group_means[i, -1]),  # Group i mean vector (excluding 'Collection' column)
      as.numeric(eur_group_means[j, -1]),  # Group j mean vector (excluding 'Collection' column)
      eur_cov_matrix  # Covariance matrix
    )
  }
}
# Display the Mahalanobis distance matrix
eur_mahalanobis_distances
####
eur_mahalanobis_distances_copy <- eur_mahalanobis_distances
# Set the diagonal and lower triangular part to NA (to remove them from the plot)
eur_mahalanobis_distances_copy[lower.tri(eur_mahalanobis_distances_copy, diag = TRUE)] <- NA
# Convert the distance matrix into a data frame suitable for ggplot2
melted_eur_mahalanobis_distances <- melt(as.matrix(eur_mahalanobis_distances_copy), na.rm = TRUE)  
# Create a heatmap with the distance values displayed in each cell with only the upper triangular part
ggplot(data = melted_eur_mahalanobis_distances, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +                                                
  geom_text(aes(label = round(value, 2)), color = "black") +   
  scale_fill_gradient(low = "white", high = "#18453B") +            
  theme_minimal() +                                            
  labs(title = "Heatmap of Mahalanobis Distances for Collection Groups (European American Samples)", 
       x = "Groups", y = "Groups", fill = "Mahalanobis Distance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))     
####
# Perform Multidimensional Scaling (MDS)
eur_mds_coords <- cmdscale(eur_mahalanobis_distances, k = 2)  
# Plot the MDS results
plot(eur_mds_coords, 
     type = "n", 
     xlab = "MDS Dimension 1", 
     ylab = "MDS Dimension 2", 
     main = "MDS Plot of Mahalanobis Distances"); text(eur_mds_coords, labels = rownames(eur_mahalanobis_distances), cex = 0.8)
# Perform hierarchical clustering
eur_hclust_res <- hclust(as.dist(eur_mahalanobis_distances))
# Plot the dendrogram
plot(eur_hclust_res, main = "Dendrogram of Mahalanobis Distances")
# Perform MDS in 3D
eur_dist_matrix <- dist(eur_mahalanobis_distances)
eur_mds_result <- cmdscale(eur_dist_matrix, k = 3)
# Create a dataframe for MDS coordinates
eur_mds_coords_3d <- data.frame(eur_mds_result)
# Rename columns for clarity
colnames(eur_mds_coords_3d) <- c("x", "y", "z")  
eur_mds_coords_3d$Collection <- eur_group_means$Collection  
# Create a 3D scatter plot with labels
plot_ly(eur_mds_coords_3d, 
        x = ~x, y = ~y, z = ~z, 
        type = "scatter3d", 
        mode = "markers",  
        marker = list(size = 5),  
        color = ~Collection,  
        text = ~Collection,  
        hoverinfo = 'text') %>%  
  layout(title = "3D MDS Plot of Mahalanobis Distances by Collection",
         scene = list(xaxis = list(title = "MDS Dimension 1"),
                      yaxis = list(title = "MDS Dimension 2"),
                      zaxis = list(title = "MDS Dimension 3")))
#######################################################
# Plot European Collection and GPS Collection Together
#######################################################
# Align geographic coordinates to craniometric coordinates
eur_procrustes_result <- procrustes(eur_mds_coords, geo_mds_coords, scale = TRUE)
# Extract transformed points
eur_original_mah <- eur_procrustes_result$X        # The original craniometric coordinates
# Plot both sets of points together
plot(eur_original_mah, col = "blue", pch = 16, xlim = range(c(eur_original_mah[,1], transformed_geo[,1])), 
     ylim = range(c(eur_original_mah[,2], transformed_geo[,2])), xlab = "Dimension 1", ylab = "Dimension 2",
     main = "Procrustes-Transformed Coordinates with European Only Samples")
points(transformed_geo, col = "red", pch = 17)
# Add collection labels to the points
text(eur_original_mah, labels = collection_labels, pos = 3, cex = 0.8, col = "blue")
text(transformed_geo, labels = collection_labels, pos = 3, cex = 0.8, col = "red")
# Add a legend to differentiate between craniometric and geographic data
legend("topright", legend = c("Craniometric Distances", "Geographic Distances"), 
       col = c("blue", "red"), pch = c(16, 17))
# Create a data frame with collection labels and point type
eur_plot_data <- data.frame(
  LD1 = c(eur_original_mah[, 1], transformed_geo[, 1]),
  LD2 = c(eur_original_mah[, 2], transformed_geo[, 2]),
  LD1_end = c(transformed_geo[, 1], transformed_geo[, 1]),  # Coordinates for segment end points
  LD2_end = c(transformed_geo[, 2], transformed_geo[, 2]),
  Collection = rep(collection_labels, 2),
  Type = rep(c("Craniometric Distance", "Geographic Distance"), each = length(collection_labels))
)
# Plot the data with lines
ggplot(eur_plot_data) +
  geom_point(aes(x = LD1, y = LD2, color = Type, shape = Type), size = 3) +
  geom_segment(aes(x = LD1, y = LD2, xend = LD1_end, yend = LD2_end),
               linetype = "dotted", show.legend = FALSE) +
  geom_text_repel(aes(x = LD1, y = LD2, label = Collection, color = Type),
                  size = 3) +
  scale_color_manual(values = c("#18453B", "#7F988E")) +  
  scale_shape_manual(values = c(16, 17), labels = c("Craniometric Distance", "Geographic Distance")) +
  labs(title = "Procrustes-Transformed Plot (Craniometric vs. Geographic Distance of European American Samples)",
       x = "Dimension 1", y = "Dimension 2", shape = "Data Type", color = "Data Type") +
  theme_classic() +
  theme(legend.position = "right")
#######################################################



#######################################################
# Classification - Reduced data by Collection
#######################################################
# Using reduced sample with just collections that have at least 100 individuals 
reduced_scaled_data <- scaled_data %>%
  filter(!(Collection %in% c("MSUFAL DSC", "WMed STARS")))
reduced_scaled_data$Collection <- droplevels(reduced_scaled_data$Collection)
####
# Fit the LDA model
model <- lda(Collection ~ ., data = reduced_scaled_data[, c(2, 38:66)], CV = TRUE) # Selecting Collection and measurements 
model
predicted_classes <- model$class
actual_classes <- reduced_scaled_data$Collection
confusion_matrix <- confusionMatrix(predicted_classes, actual_classes)
print(confusion_matrix)
#######################################################
# Classification - All data by Population/Sex instead of Collection
#######################################################
# Create df of balanced pop.sex groups
popsex_scaled_data <- scaled_data
popsex_scaled_data$Pop.Sex <- as.factor(popsex_scaled_data$Pop.Sex)
# Count the individuals in each pop.sex group
popsex_counts <- popsex_scaled_data %>%
  group_by(Pop.Sex) %>%
  summarise(count = n())
# Filter out groups with less than 7
popsex_scaled_data <- popsex_scaled_data %>%
  inner_join(popsex_counts %>% filter(count > 6), by = "Pop.Sex")
#Identify min size of Pop.Sex groups
popsex_min_size <- popsex_scaled_data %>%
  group_by(Pop.Sex) %>%
  summarise(count = n()) %>%
  pull(count) %>%
  min()
# Downsample each category to have the same number of observations
popsex_scaled_data <- popsex_scaled_data %>%
  group_by(Pop.Sex) %>%
  sample_n(popsex_min_size) %>%
  ungroup()
# View the updated data
head(popsex_scaled_data)
popsex_scaled_data$Pop.Sex <- droplevels(popsex_scaled_data$Pop.Sex)
####
popsex_model <- lda(Pop.Sex ~ ., data = popsex_scaled_data[, c(23, 38:66)], CV = TRUE) # Selecting Pop.Sex and measurements 
popsex_model
popsex_predicted_classes <- popsex_model$class
popsex_actual_classes <- popsex_scaled_data$Pop.Sex
popsex_confusion_matrix <- confusionMatrix(popsex_predicted_classes, popsex_actual_classes)
print(popsex_confusion_matrix) 
#######################################################
# LDFA Results Plot
#######################################################
# Create a dataframe for both models' performance metrics
model_performance <- data.frame(
  Metric = c("Accuracy", "Kappa", "Sensitivity", "Specificity"),
  Model1 = c(0.3201, 0.1502, mean(c(0.076, 0.32, 0.455, 0.14, 0.436, 0.290)), 
             mean(c(0.960, 0.948, 0.780, 0.919, 0.720, 0.819))),
  Model2 = c(0.4896, 0.4167, mean(c(0.416, 0.333, 0.5, 0.583, 0.583, 0.5, 0.333, 0.666)), 
             mean(c(0.916, 0.928, 0.928, 0.916, 0.964, 0.928, 0.916, 0.916)))
)

# Reshape data for ggplot
model_performance_long <- reshape2::melt(model_performance, id = "Metric")

# Plot Bar Chart
ggplot(model_performance_long, aes(x = Metric, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "LDFA Performance Comparison", x = "Metric", y = "Value") +
  scale_fill_manual(values = c("#A7B7B0", "#18453B"), labels = c("Collection", "Population and Sex")) +
  theme_classic()
#######################################################