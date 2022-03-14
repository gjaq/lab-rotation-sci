library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(UpSetR)
library(gtools)
library(data.table)
library(ggpubr)
library(lattice)
library(FSA)

#---
# Load the data ----
df_murnau <- read.csv("~/Data/HematologicalBiomark_DATA_2021-01-07_0729.csv",
                      na.strings = c(""))
path <- "~/Data/HematologicalBiomarkersForSCI_DataDictionary_2019-10-01.csv"
dictionary_murnau <- read.csv(path)

# Change columns names of the dictionary (remove dots)
names(dictionary_murnau) <- gsub("\\.", "", names(dictionary_murnau))


########################################################################################################## #
#### Useful functions (the plot functions are in the different sections of this script) ####  
########################################################################################################## #

# Function to extract the days from the dataset
extract_days = function(df, name_df, rm_time_in_day = T) {
  test_points <- names(df)
  
  # Remove the name of the marker, so that we only get the day
  to_remove <- paste0(name_df, "_")
  test_days <- gsub(to_remove, "", test_points)
  
  if(rm_time_in_day) {
    # Remove the last character of each name, indicating whether this is the first, second, ...
    # test done on one specific day
    test_days <- gsub(".{1}$", "", test_days)
  }
  
  return(test_days)
}

# Function to make the sum of the NAs in vector x
sum_na = function(x) {
  return(sum(is.na(x)))
}

# Function that adds a columns "grade" (for the very acute and the acute stages) to the dataframe
# By default, the grade determined in the acute stage is taken.
# If however the grade at the acute stage is NA or ND, then the grade at the very acute stage is used
add_grade_to_df = function(df, df_murnau_ = df_murnau) {
  df$Grade_va <- df_murnau_$va_ais
  df$Grade_ai <- df_murnau_$ai_ais
  
  df <- df %>% mutate(Grade = ifelse(Grade_va == "ND" | is.na(Grade_va), Grade_ai, Grade_va)) %>%
    select(-c("Grade_va", "Grade_ai"))
  return(df)
}
#mutate(Grade = ifelse(Grade_ai == "ND" | is.na(Grade_ai), Grade_va, Grade_ai))

########################################################################################################## #
#### Create the 28 datasets with the 28 routinely assessed blood markers ####
########################################################################################################## #

df_amylase <- df_murnau %>% select(contains("amylase")) %>% select(!contains("neg"))
df_ap <- df_murnau %>% select(starts_with("ap")) %>% select(!contains("neg"))
df_calcium <- df_murnau %>% select(contains("calcium")) %>% select(!contains("neg"))
df_che <- df_murnau %>% select(starts_with("che")) %>% select(!contains("neg"))
df_creatinin <- df_murnau %>% select(contains("creatinin")) %>% select(!contains("neg"))
df_crp <- df_murnau %>% select(contains("crp")) %>% select(!contains("neg"))
df_ery <- df_murnau %>% select(starts_with("ery")) %>% select(!contains("neg")) %>% 
  select(!contains("sed"))
df_gamma_gt <- df_murnau %>% select(contains("gamma_gt")) %>% select(!contains("neg"))
df_bilirubin <- df_murnau %>% select(contains("gesamt_bilirubin")) %>% select(!contains("neg"))
df_gesamteiweiss <- df_murnau %>% select(contains("gesamteiweiss")) %>% select(!contains("neg"))
df_glucose <- df_murnau %>% select(starts_with("glucose")) %>% select(!contains("neg"))
df_harnstoff <- df_murnau %>% select(contains("harnstoff")) %>% select(!contains("neg"))
df_hb <- df_murnau %>% select(starts_with("hb_")) %>% select(!contains("neg"))
df_hbe <- df_murnau %>% select(contains("hbe")) %>% select(!contains("neg")) %>% select(!contains("mch"))
df_hk <- df_murnau %>% select(contains("hk")) %>% select(!contains("neg"))
df_inr <- df_murnau %>% select(contains("inr")) %>% select(!contains("neg"))
df_kalium <- df_murnau %>% select(contains("kalium")) %>% select(!contains("neg")) %>%
  select(!contains("urin"))
df_ldh <- df_murnau %>% select(contains("ldh")) %>% select(!contains("neg"))
df_leuco_nl <- df_murnau %>% select(contains("leuco_nl")) %>% select(!contains("neg"))
df_lipase <- df_murnau %>% select(contains("lipase")) %>% select(!contains("neg"))
df_mchc <- df_murnau %>% select(contains("mchc")) %>% select(!contains("neg"))
df_mcv <- df_murnau %>% select(contains("mcv")) %>% select(!contains("neg"))
df_natrium <- df_murnau %>% select(contains("natrium")) %>% select(!contains("neg")) %>%
  select(!contains("urin"))
df_ptt <- df_murnau %>% select(contains("ptt")) %>% select(!contains("neg"))
df_quick <- df_murnau %>% select(contains("quick")) %>% select(!contains("neg"))
df_thrombo <- df_murnau %>% select(starts_with("thrombo_")) %>% select(!contains("neg"))
df_got <- df_murnau %>% select(contains("got")) %>% select(!contains("neg"))
df_gpt <- df_murnau %>% select(contains("gpt")) %>% select(!contains("neg"))

list_names_murnau = c("Amylase", "Alkaline phosphatase", "Calcium", "Cholinesterase", "Creatinin",
                      "CRP", "Erythrocytes", "Gamma GT", "Total Bilirubin", "Proteins", "Glucose", 
                      "Blood urea", "Hemoglobin", "Hemoglobin per erythrocyte", "Hematocrit", "INR", 
                      "Potassium", "Lactate dehydrogenase", "Leucocytes", "Lipase","MCHC", "MCV", 
                      "Sodium", "Prothrombin time", "Quick test", "Thrombocytes", "ASAT", "ALAT")

list_murnau = list(df_amylase, df_ap, df_calcium, df_che, df_creatinin, df_crp, df_ery, df_gamma_gt, 
                   df_bilirubin, df_gesamteiweiss, df_glucose, df_harnstoff, df_hb, df_hbe, df_hk, df_inr, 
                   df_kalium, df_ldh, df_leuco_nl, df_lipase, df_mchc, df_mcv, df_natrium, df_ptt, 
                   df_quick, df_thrombo, df_got, df_gpt)

list_df_names_murnau = c("amylase", "ap", "calcium", "che", "creatinin", "crp", "ery", "gamma_gt", 
                         "gesamt_bilirubin", "gesamteiweiss", "glucose", "harnstoff", "hb", "hbe",
                         "hk", "inr", "kalium", "ldh", "leuco_nl", "lipase","mchc", "mcv", "natrium", 
                         "ptt", "quick", "thrombo", "got", "gpt")

# Convert all values to numeric
#list_murnau_numeric <- sapply(list_murnau, as.numeric)

# Remove columns that only contain NAs
remove_NA_columns = function(df) {
  return(Filter(function(x) !all(is.na(x)), df))
} 
list_murnau_no_na_col <- sapply(list_murnau, remove_NA_columns)

# Sort the datasets according to the lab tests orders
sort_dataframe = function(df) {
  df <- df[ , str_sort(names(df), numeric = TRUE)]
  return(df)
}

list_murnau_sorted <- sapply(list_murnau_no_na_col, sort_dataframe)


########################################################################################################## #
#### Count the number of AIS grades at each stage ####
########################################################################################################## #

grades <- c("A", "B", "C", "D", "E", "ND", "NA")
grades_df <- data.frame(Very = df_murnau$va_ais, AcuteI = df_murnau$ai_ais, 
                        AcuteII = df_murnau$aii_aiis, AcuteIII = df_murnau$aiii_aiiis,
                        Chronic = df_murnau$c_cs)
grades_df <- add_grade_to_df(grades_df)
df_summary_grades = data.frame("AIS_grade" = grades, 
                               "Very_acute" = count(grades_df, Very)$n, 
                               "AcuteI" = count(grades_df, AcuteI)$n,
                               "AcuteII" = count(grades_df, AcuteII)$n,
                               "AcuteIII" = count(grades_df, AcuteIII)$n,
                               "Chronic" = count(grades_df, Chronic)$n,
                               "Grade_given" = count(grades_df, Grade)$n)
rownames(df_summary_grades) <- grades


########################################################################################################## #
#### Create the list of dataframes for the first two weeks ####
########################################################################################################## #

# From now on: focus on the first two weeks
# For each marker, create a dataset for the first two weeks
subset_df_by_time = function(df, name_df, n_days) {
  # Extract the days from the dataset
  days <- extract_days(df, name_df)
  
  # Select the columns whose day is less than n_days
  # Selects eg the first 14 days, first month, ...
  inds <- as.numeric(days) <= n_days
  df_month <- df[, inds]
  
  return(df_month)
}

## Apply the function to each dataset to select data for the first two weeks
n_days <- rep(14, 28)
list_murnau_2_weeks <- mapply(subset_df_by_time, list_murnau_sorted, list_df_names_murnau, n_days)


########################################################################################################## #
#### Melt all 28 dataset ####
########################################################################################################## #

# Construct a function that melt all values of a dataframe
# So that we obtain a two-columns dataframe, with the first column recording the day and the patient
# (i.e. 0a_1) and the first column the test value
bin_melt_df_day_patient_value = function(df, name_df) {
  # Make a binary dataset: 0 NA, 1 value not NA
  df <- df %>% mutate_all(function(x) ifelse(is.na(x), 0, 1))
  
  # Change the column names so that we only have the day
  days <- extract_days(df, name_df, F)
  names(df) <- days
  
  # Add columns "Grade" and "Patient" to the dataframe
  df$Patient <- 1:363
  df <- add_grade_to_df(df)
  
  # Melt the dataframe: keep columns grade and patient
  # Successively, remove column DayTime
  df_melt <- pivot_longer(df, cols = -c("Grade", "Patient"), names_to = "DayTime", values_to = name_df) %>% 
    select(-DayTime)
  
  # Change the name of the columns, so that all dataframes will have four columns, the first three named
  # "Patient", "Grade", "Cases" and "DayTime" respectively, while the fourth will have the name of the marker
  names(df_melt) <- c("Patient", "Grade", name_df)
  
  return(df_melt)
}

list_df_melted <- mapply(bin_melt_df_day_patient_value, list_murnau_2_weeks, list_df_names_murnau,
                         SIMPLIFY = F)


########################################################################################################## #
#### Make the dataset for the analysis: dataset of 2-weeks counts ####
########################################################################################################## #
count_tests_per_patient = function(df, name_df) {
  df <- df %>% 
    group_by(Patient) %>% 
    mutate(N_tests = sum(!!sym(name_df))) %>% 
    distinct(Patient, Grade, .keep_all = TRUE) %>% 
    select(-all_of(name_df)) %>% 
    ungroup
  names(df) <- c("Patient", "Grade", name_df)
  return(df)
}

list_df_counts <- mapply(count_tests_per_patient, list_df_melted, list_df_names_murnau,
                         SIMPLIFY = F)

# Combine all dataset together
df_counts <- join_all(list_df_counts, by = c("Grade", "Patient"), type = "inner")

# Count the total number of tests done on one patient, considering all markers
df_counts_tot <- df_counts %>% mutate(Total = select(., amylase:gpt) %>% rowSums())

# Visualise on boxplot
df_counts_tot %>% filter(Grade %in% c("A", "B", "C", "D")) %>% 
  ggplot(aes(x = Grade, y = Total, fill = Grade)) +
  geom_boxplot() +
  geom_point(aes(fill = Grade), size = 1, shape = 21, position = position_jitterdodge()) +
  labs(title = "Distribution of the n° of tests per patient \n according to their AIS grade", 
       x = "", y = "N° tests in the first 2 weeks")
ggsave(filename = "boxplot_28_tot.png", path = "~/Plots/Plots_28_LM_tests")

df_counts_tot %>% filter(Grade %in% c("A", "B", "C", "D")) %>% 
  ggplot(aes(x = Grade, y = Total, fill = Grade)) +
  geom_violin() +
  geom_point(aes(fill = Grade), size = 1, shape = 21, position = position_jitterdodge()) +
  labs(title = "Distribution of the n° per patient according \n to their AIS grade", 
       x = "", y = "N° tests in the first 2 weeks")
ggsave(filename = "violin_28_tot.png", path = "~/Plots/Plots_28_LM_tests")

# Plot a boxplot for each marker
# Do this to see whether there is one marker that is more appropriate to use in order to carry out ANOVA
# for (i in 1:28) {
#   marker <- list_df_names_murnau[[i]]
#   df <- df_counts_tot %>% select(c("Grade", marker))
#   
#   df %>% filter(Grade %in% c("A", "B", "C", "D")) %>% 
#     ggplot(aes(x = Grade, y = !!sym(marker), fill = Grade)) +
#     geom_boxplot() +
#     geom_point(aes(fill = Grade), size = 1, shape = 21, position = position_jitterdodge()) +
#     labs(title = "Distribution of the n° per patient according \n to their AIS grade", 
#          x = "", y = "N° tests in the first 2 weeks")
#   
#   f1 <- paste0(marker, "_boxplot.png")
#   ggsave(filename = f1, path = "~/Code/Plots_28_LM_tests")
# }


########################################################################################################## #
#### Check if the counts are normally distributed and with equal variances ####
########################################################################################################## #
df_counts_a <- df_counts_tot %>% filter(Grade == "A") %>% select(c("Patient", "Grade", "Total"))
df_counts_b <- df_counts_tot %>% filter(Grade == "B") %>% select(c("Patient", "Grade", "Total"))
df_counts_c <- df_counts_tot %>% filter(Grade == "C") %>% select(c("Patient", "Grade", "Total"))
df_counts_d <- df_counts_tot %>% filter(Grade == "D") %>% select(c("Patient", "Grade", "Total"))

# Rule of thumb in order to check whether variances are equal
# Check if the highest variance is less than 4 times the smallest variance
var_a <- var(df_counts_a$Total)
var_b <- var(df_counts_b$Total)
var_c <- var(df_counts_c$Total)
var_d <- var(df_counts_d$Total)
variances <- c(var_a, var_b, var_c, var_d)
max(variances) < 4 * min(variances)
# Using this rule of thumb, we can say that equality of variances assumption is met.

# Now check whether the 4 datasets are normally distributed
# First, make visual checking
# Density plot in order to judge whether the distribution is bell shaped
a <- ggdensity(df_counts_a$Total, main = "Density plot of number of tests for AIS grade A", xlab = "N° tests")
b <- ggdensity(df_counts_b$Total, main = "Density plot of number of tests for AIS grade B", xlab = "N° tests")
c <- ggdensity(df_counts_c$Total, main = "Density plot of number of tests for AIS grade C", xlab = "N° tests")
d <- ggdensity(df_counts_d$Total, main = "Density plot of number of tests for AIS grade D", xlab = "N° tests")
p1 <- plot_grid(plotlist = list(a, b, c, d), nrow = 2, ncol = 2)
save_plot(filename = "counts_density.png", plot = p1, path = "~/Plots/Plots_28_LM_tests",
          nrow = 2, ncol = 2)
# From the plots, none of the groups looks normally distributed.

# Since I'm not sure that plotting a continuous distribution is wise in this case: plot histogram
df_counts_tot_ <- filter(df_counts_tot, Grade %in% c("A", "B", "C", "D"))
histogram(~ Total | Grade, data = df_counts_tot_, layout = c(1,4))

# QQplot: plots one quantile against the other
# If the distribution is normal, the obtained plot will be a straight line, meaning that the quantiles are equal
# (since the density curve is bell-shaped)
a <- ggqqplot(df_counts_a$Total, title = "Q-Q plot for AIS grade A")
b <- ggqqplot(df_counts_b$Total, title = "Q-Q plot for AIS grade B")
c <- ggqqplot(df_counts_c$Total, title = "Q-Q plot for AIS grade C")
d <- ggqqplot(df_counts_d$Total, title = "Q-Q plot for AIS grade D")
p1 <- plot_grid(plotlist = list(a, b, c, d), nrow = 2, ncol = 2)
save_plot(filename = "counts_qq_plot.png", plot = p1, path = "~/Plots/Plots_28_LM_tests",
          nrow = 2, ncol = 2)

# Statistical test: Shapiro-Wilk's method
# Tests if the four distributions are normally distributed
shapiro.test(df_counts_a$Total)
shapiro.test(df_counts_b$Total)
shapiro.test(df_counts_c$Total)
shapiro.test(df_counts_d$Total)

# From the p-values obtained, we can tell that D population (and maybe also A) is not normally distributed.
# Shapiro test tends to say that populations with small size are normally distributed
# but from the plots, we should probably conclude that neither B or C are normally distributed.

########################################################################################################## #
#### Perform ANOVA: all 28 ####
########################################################################################################## #
df_counts_tot_an <- df_counts_tot %>% filter(Grade %in% c("A", "B", "C", "D"))
one.way <- aov(che ~ Grade, data = df_counts_tot_an)
summary_an <- summary(one.way)

tukey.one.way <- TukeyHSD(one.way)
tukey.one.way

# The difference between the group means is significant and is given by grade D.
# However, normality assumption is not met, thus probably ANOVA is not a good choice here.


########################################################################################################## #
#### Another implementation of ANOVA ####
########################################################################################################## #
# This implementation of ANOVA doesn't assume that variances are equal
df_counts_tot_an <- df_counts_tot %>% filter(Grade %in% c("A", "B", "C", "D"))
oneway <- oneway.test(Total ~ Grade, data = df_counts_tot_an)


########################################################################################################## #
#### Kruskal-Wallis test ####
########################################################################################################## #
# Non parametric test, doesn't assume that data are normally distributed
# Null hypothesis: the mean ranks of the groups are the same (groups are sampled from the same population)
# Assumptions: observations in each group come from populations with the same shape of distribution 
# (which is not my case: D is right skewed, the other are more left skewed)
# Designed for ranked variables
kruskal.test(Total ~ Grade, data = df_counts_tot_an) 

PT = dunnTest(Total ~ Grade, data = df_counts_tot_an, method = "bh")



