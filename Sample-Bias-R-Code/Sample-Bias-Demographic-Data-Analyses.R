######################################################
# Title: Sample Bias - Demographic Analyses
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

library(tidyverse)     
library(sf)
library(usmap)
library(chisq.posthoc.test)

#######################################################
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
# 1- Sex (line 53)
# 2- Birth Year (line 74)
# 3- Year of Death (line 110)
# 4- Age (line 142)
# 5- Population (line 198)
# 6- Height (line 224)
# 7- SES (line 296)
# 8- Birthplace (line 349) 
# 9- Place at Time of Death (line 393)

#######################################################
#SEX OF INDIVIDUALS IN COLLECTIONS
#######################################################
# Create sex contingency table
sex_contingency_table <- data %>%
  count(Collection, Sex) %>%
  xtabs(n ~ Collection + Sex, data = .)
####
# Perform the Chi-squared test for sex
chisq.test(sex_contingency_table)
# Post hoc analysis
chisq.posthoc.test(sex_contingency_table, method = "bonferroni")
####
# Stacked bar plot for Sex by Collection
ggplot(data, aes(x = Collection, fill = Sex)) +
  scale_fill_manual(values = c("#18453B", "#7F988E")) + 
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "Proportion of Males and Females by Collection",
       x = "Collection",
       y = "Proportion")
#######################################################
#BIRTH YEAR OF INDIVIDUALS 
#######################################################
# Compute summary statistics by Collection
data %>%
  filter(!is.na(BirthYear)) %>%
  group_by(Collection) %>%
  summarise(
    Mean = mean(BirthYear),
    Median = median(BirthYear),
    SD = sd(BirthYear),
    Count = n(),
    .groups = "drop"
  )
####
# Perform ANOVA for Birth Year by Collection
birthyear_anova <- aov(BirthYear ~ Collection, data = data)
summary(birthyear_anova)
# Post hoc test (Tukey's HSD)
TukeyHSD(birthyear_anova)
####
# Boxplot for birth year by Collection
ggplot(data, aes(x = Collection, y = BirthYear)) +
  geom_boxplot(fill = "#A7B7B0") +
  labs(title = "Birth Year by Collection",
       x = "Collection",
       y = "Birth Year") +
  theme_classic() +
  theme(legend.position = "none")         
# Voilin plot for birth year by Collection
ggplot(data, aes(x = Collection, y = BirthYear, fill = Collection)) +
  geom_violin(fill = "#A7B7B0") +
  labs(title = "Violin Plot of Birth Year by Collection",
       x = "Collection",
       y = "Birth Year") +
  theme_minimal()
#######################################################
#YEAR OF DEATH OF INDIVIDUALS 
#######################################################
# Compute summary statistics by Collection
data %>%
  filter(!is.na(DeathYear)) %>%
  group_by(Collection) %>%
  summarise(
    Mean = mean(DeathYear),
    Median = median(DeathYear),
    SD = sd(DeathYear),
    Count = n(),
    .groups = "drop"
  )
####
# Perform ANOVA for Year of Death by Collection
deathyear_anova <- aov(DeathYear ~ Collection, data = data)
summary(deathyear_anova)
# Post hoc test (Tukey's HSD)
TukeyHSD(deathyear_anova)
####
# Boxplot for Year of Death by Collection
ggplot(data, aes(x = Collection, y = DeathYear)) +
  geom_boxplot(fill = "#A7B7B0") +
  theme_minimal()
# Voilin plot for Year of Death by Collection
ggplot(data, aes(x = Collection, y = DeathYear, fill = Collection)) +
  geom_violin(fill = "#A7B7B0") +
  labs(title = "Violin Plot of Year of Death by Collection",
       x = "Collection",
       y = "Year of Death") +
  theme_minimal()
#######################################################
#AGE OF INDIVIDUALS ACROSS COLLECTIONS
#######################################################
# Compute summary statistics by Collection
data %>%
  filter(!is.na(Age)) %>%
  group_by(Collection) %>%
  summarise(
    Mean = mean(Age),
    Median = median(Age),
    SD = sd(Age),
    Count = n(),
    .groups = "drop"
  )
####
# Perform ANOVA for Age by Collection
age_anova <- aov(Age ~ Collection, data = data)
summary(age_anova)
# Post hoc test (Tukey's HSD)
TukeyHSD(age_anova)
####
# Perform ANOVA for Age by Collection and Sex
agesex_anova <- aov(Age ~ Sex*Collection, data = data)
summary(agesex_anova)
# Post hoc test (Tukey's HSD)
TukeyHSD(agesex_anova)
####
# Boxplot for Age by Collection
ggplot(data, aes(x = Collection, y = Age)) +
  geom_boxplot(fill = "#A7B7B0") +
  theme_minimal()
# Voilin plot for Age by Collection
ggplot(data, aes(x = Collection, y = Age, fill = Collection)) +
  geom_violin(fill = "#A7B7B0") +
  labs(title = "Violin Plot of Age by Collection",
       x = "Collection",
       y = "Age") +
  theme_minimal()
# Boxplot for Age by Collection and Sex
ggplot(data, aes(x = Collection, y = Age, fill = Sex)) +
  geom_boxplot() +
  scale_fill_manual(values = c("F" = "#18453B", "M" = "#7F988E")) +
  labs(
    title = "Interaction Between Sex and Collection on Age",
    x = "Collection",
    y = "Age"
  ) +
  theme_minimal()
# Voilin plot for Age by Collection and Sex
ggplot(data, aes(x = Collection, y = Age, fill = Sex)) +
  geom_violin() +
  scale_fill_manual(values = c("F" = "#18453B", "M" = "#7F988E")) +
  labs(title = "Violin Plot of Age by Collection and Sex",
       x = "Collection",
       y = "Age") +
  theme_minimal()
######################################################
#POPULATION OF INDIVIDUALS IN COLLECTIONS
#######################################################
# Filter out individuals with two or more population groups for chi-squared analysis
pop_data <- data %>%
  filter(Population.Short %in% c("W", "B", "H", "I", "A"))
####
# Create a contingency table for population data
pop_table <- pop_data %>%
  filter(!is.na(Population.Short), Population.Short != "") %>%
  count(Collection, Population.Short) %>%
  xtabs(n ~ Collection + Population.Short, data = .)
pop_table
####
# Perform the Chi-squared test for population
chisq.test(pop_table, simulate.p.value = TRUE)
# Post hoc analysis
chisq.posthoc.test(pop_table, method = "bonferroni")
####
# Stacked bar plot for Populations in each Collection
ggplot(pop_data, aes(x = Collection, fill = Population.Short)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "Proportion of Population Groups by Collection",
       x = "Collection",
       y = "Proportion")
######################################################
#HEIGHT OF INDIVIDUALS IN THE COLLECTION 
#######################################################
# Compute summary statistics by Collection
data %>%
  mutate(Height_CM = as.numeric(Height_CM)) %>%
  filter(!is.na(Height_CM), Height_CM >= 122, !is.na(Sex)) %>%
  group_by(Collection, Sex) %>%
  summarise(
    Mean = mean(Height_CM),
    Median = median(Height_CM),
    SD = sd(Height_CM),
    Count = n(),
    .groups = "drop"
  )
####
# Perform ANOVA for Birth Year by Collection
height_anova <- aov(Height_CM ~ Collection, data = data)
summary(height_anova)
# Post hoc test (Tukey's HSD)
TukeyHSD(height_anova)
####
# Perform ANOVA for Birth Year by Collection + Sex + Age + Birth Year
height_anova2 <- aov(Height_CM ~ Collection + Sex + BirthYear, data = data)
summary(height_anova2)
####
# Boxplot for height by Collection
ggplot(data, aes(x = Collection, y = Height_CM, fill = Collection)) +
  geom_boxplot(fill = "#A7B7B0") +
  theme_minimal()
# Boxplot of height by Collection and Sex
ggplot(data, aes(x = Collection, y = Height_CM, color = Sex)) +
  geom_boxplot() +
  labs(title = "Interaction Between Sex and Collection on Height",
       x = "Collection",
       y = "Height (cm)") +
  theme_minimal()
# Voilin plot for height by Collection
ggplot(data, aes(x = Collection, y = Height_CM, fill = Collection)) +
  geom_violin() +
  labs(title = "Violin Plot of Height by Collection",
       x = "Collection",
       y = "Height (cm)") +
  theme_minimal()
# Voilin plot for height by Collection and sex
ggplot(data, aes(x = Collection, y = Height_CM, fill = Sex)) +
  scale_fill_manual(values = c("F" = "#18453B", "M" = "#7F988E")) +
  geom_violin() +
  labs(title = "Violin Plot of Height by Collection and Sex",
       x = "Collection",
       y = "Height (cm)") +
  theme_minimal()
# Plot of height by birth year 
ggplot(data, aes(x = BirthYear, y = Height_CM)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    title = "Relationship Between Height and Birth Year",
    x = "Birth Year",
    y = "Height (cm)"
  ) +
  theme_minimal()
# Plot of height by birth year by sex
ggplot(data, aes(x = BirthYear, y = Height_CM, color = Sex)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  labs(
    title = "Height by Birth Year and Sex",
    x = "Birth Year",
    y = "Height (cm)"
  ) +
  theme_minimal()
#######################################################
#SOCIOECONOMIC STATUS (Childhood and Adulthood)
#######################################################
####
# CHILDHOOD
####
# Get childhood SES data
child_ses_data <- data %>%
  filter(!is.na(ChildSES), ChildSES != "")
# Make table 
child_ses_table <- child_ses_data %>%
  count(Collection, ChildSES, name = "Count") %>%
  xtabs(Count ~ Collection + ChildSES, data = .)
child_ses_table
# chi-squared analysis 
chisq.test(child_ses_table, simulate.p.value = TRUE)
####
# Stacked bar plot for Childhood SES Category by Collection
ggplot(child_ses_data, aes(x = Collection, fill = ChildSES)) +
  geom_bar(position = "fill") +
  labs(
    title = "Proportion of Childhood SES Categories by Collection",
    x = "Collection",
    y = "Proportion"
  ) +
  theme_minimal()
####
# ADULT
####
# Get adulthood SES data
adult_ses_data <- data %>%
  filter(!is.na(AdultSES), AdultSES != "")
# Make table 
adult_ses_table <- adult_ses_data %>%
  count(Collection, AdultSES, name = "Count") %>%
  xtabs(Count ~ Collection + AdultSES, data = .)
adult_ses_table
####
#Chi-square analysis
adult_ses_chisq <- chisq.test(adult_ses_table, simulate.p.value = TRUE)
adult_ses_chisq
# Post hoc analysis
chisq.posthoc.test(adult_ses_table, method = "bonferroni")
####
# Stacked bar plot for Adulthood SES Category by Collection
ggplot(adult_ses_data, aes(x = Collection, fill = AdultSES)) +
  geom_bar(position = "fill") +
  labs(
    title = "Proportion of Adult SES Categories by Collection",
    x = "Collection",
    y = "Proportion"
  ) +
  theme_minimal()
#######################################################
#BIRTHPLACE OF INDIVIDUALS
#######################################################
# Filter the dataset to keep only rows with state codes
birthplace_data <- data %>%
  filter(!is.na(Birthplace_State),
         Birthplace_State %in% state.abb)
# Counts by collection and birthplace state
birthplace_counts <- birthplace_data %>%
  count(Collection, Birthplace_State, name = "Count")
birthplace_counts
####
#Get within-collection birth state proportions 
birthplace_props <- birthplace_counts %>%
  group_by(Collection) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  rename(state = Birthplace_State)
birthplace_props
####
# Plot birthplace for each collection
plot_usmap(
  data = birthplace_props,
  values = "Proportion",
  regions = "states"
) +
  scale_fill_gradient(
    low = "#edf2f4",
    high = "#18453B",
    name = "Proportion"
  ) +
  facet_wrap(~ Collection) +
  labs(title = "State of Birth") +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )
####
# Chi-squared test for birthplace 
birthplace_table <- birthplace_counts %>%
  xtabs(Count ~ Collection + Birthplace_State, data = .)
chisq.test(birthplace_table, simulate.p.value = TRUE)
# Posthoc analysis 
chisq.posthoc.test(birthplace_table, method = "bonferroni")
#######################################################
#PLACE OF DEATH OF INDIVIDUALS
#######################################################
# Filter the dataset to keep only rows with state codes
deathplace_data <- data %>%
  filter(!is.na(Death_State),
         Death_State %in% state.abb)
# Counts by collection and state at time of death
deathplace_counts <- deathplace_data %>%
  count(Collection, Death_State, name = "Count")
deathplace_counts
####
# Get within-collection death state proportions 
deathplace_props <- deathplace_counts %>%
  group_by(Collection) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  rename(state = Death_State)
deathplace_props
####
# Plot place of death for each collection
plot_usmap(
  data = deathplace_props,
  values = "Proportion",
  regions = "states"
) +
  scale_fill_gradient(
    low = "#edf2f4",
    high = "#18453B",
    name = "Proportion"
  ) +
  facet_wrap(~ Collection) +
  labs(title = "State at Time of Death") +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )
####
# Chi-squared test
deathplace_table <- deathplace_counts %>%
  xtabs(Count ~ Collection + Death_State, data = .)
chisq.test(deathplace_table, simulate.p.value = TRUE)
# Post hoc analysis 
chisq.posthoc.test(deathplace_table, method = "bonferroni")
#######################################################