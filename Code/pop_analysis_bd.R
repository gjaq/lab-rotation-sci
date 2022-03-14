# Load packages ----
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)

#---
# Load the data ----
# Dataset used for the analysis of blood draws
n_test_tot_df <- read.csv("~/Data/n_test.csv")

# Dataset with the demographic information
demog_df <- read.csv("~/Data/demographic_info.csv")


########################################################################################################## #
#### Construct the dataset with all necessary information for the analysis ####  
########################################################################################################## #
demog_df_red <- demog_df %>% select(c("patientennummer", "AgeAtDOI", "Sex"))
n_test_tot_df_red <- n_test_tot_df %>% select(c("patientennummer", "Grade", "Plegia"))
murnau_demog <- merge(x = demog_df_red, y = n_test_tot_df_red, by = "patientennummer")


########################################################################################################## #
#### Check general characteristics ####  
########################################################################################################## #
# Sex distribution
sex_counts <- table(murnau_demog$Sex)

# Age distribution
age_mean <- mean(murnau_demog$AgeAtDOI, na.rm = T)
age_sd <- sd(murnau_demog$AgeAtDOI, na.rm = T)
age_min <- min(murnau_demog$AgeAtDOI, na.rm = T)
age_max <- max(murnau_demog$AgeAtDOI, na.rm = T)

# Injury level
para <- sum(murnau_demog$Plegia == "para")
tetra <- sum(murnau_demog$Plegia == "tetra")

# ASIA grade
grade_counts <- table(murnau_demog$Grade)


########################################################################################################## #
#### Analysis for paraplegic patients ####  
########################################################################################################## #
# Make the dataset
murnau_demog_para <- murnau_demog %>% filter(Plegia == "para")
mean(murnau_demog_para$AgeAtDOI, na.rm = T)
sd(murnau_demog_para$AgeAtDOI, na.rm = T)

table(murnau_demog_para$Sex)
table(murnau_demog_para$Sex)/147

########################################################################################################## #
#### Analysis for tetraplegic patients ####  
########################################################################################################## #
# Make the dataset
murnau_demog_tetra <- murnau_demog %>% filter(Plegia == "tetra")
mean(murnau_demog_tetra$AgeAtDOI, na.rm = T)
sd(murnau_demog_tetra$AgeAtDOI, na.rm = T)

table(murnau_demog_tetra$Sex)
table(murnau_demog_tetra$Sex)/152
