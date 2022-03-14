# Load packages ----
library(dplyr)
library(plyr)
library(sjmisc)
library(ggplot2)
library(gtools)
library(reshape2)
library(UpSetR)
library(tidyr)
library(stringr)
library(ggpubr)
library(lattice)
library(FSA)


# Load the data ----
df_murnau <- read.csv("~/Data/HematologicalBiomark_DATA_2021-01-07_0729.csv",
                      na.strings = c(""))
path <- "~/Data/HematologicalBiomarkersForSCI_DataDictionary_2019-10-01.csv"
dictionary_murnau <- read.csv(path)

# Change columns names of the dictionary (remove dots)
names(dictionary_murnau) <- gsub("\\.", "", names(dictionary_murnau))


########################################################################################################## #
#### Useful functions ####                                                                                   
########################################################################################################## #

extract_days = function(test_names, rm_time_in_day = T) {
  # Removes everything before and including the last _
  test_days <- gsub(".*_", "", test_names)
  
  if(rm_time_in_day) {
    # Remove the last character of each name, indicating whether this is the first, second, ...
    # test done on one specific day
    test_days <- gsub(".{1}$", "", test_days)
  }
  
  return(test_days)
}


# Function that adds a columns "grade" (for the very acute and the acute stages) to the dataframe
# By default, the grade determined in the very acute stage is taken.
# If however the grade at the acute stage is NA or ND, then the grade at the acute stage is used
return_grade = function(df_murnau_ = df_murnau) {
  df <- data.frame(Grade_va = df_murnau_$va_ais, Grade_ai = df_murnau_$ai_ais)
  
  df <- df %>% mutate(Grade = ifelse(Grade_va == "ND" | is.na(Grade_va), Grade_ai, Grade_va))
  return(df$Grade)
}


########################################################################################################## #
#### Count the number of AIS grades at each stage ####
########################################################################################################## #
grades <- c("A", "B", "C", "D", "E", "ND", "NA")
grades_df <- data.frame(Very = df_murnau$va_ais, AcuteI = df_murnau$ai_ais, 
                        AcuteII = df_murnau$aii_aiis, AcuteIII = df_murnau$aiii_aiiis,
                        Chronic = df_murnau$c_cs)
grade_given <- return_grade()
grades_df$Grade_given <- grade_given
df_summary_grades = data.frame("AIS_grade" = grades, 
                               "Very_acute" = count(grades_df, Very)$n, 
                               "AcuteI" = count(grades_df, AcuteI)$n,
                               "AcuteII" = count(grades_df, AcuteII)$n,
                               "AcuteIII" = count(grades_df, AcuteIII)$n,
                               "Chronic" = count(grades_df, Chronic)$n,
                               "Grade_given" = count(grades_df, Grade_given)$n)
rownames(df_summary_grades) <- grades


########################################################################################################## #
#### Extract all the unique lab tests made from the Murnau dataset ####       
########################################################################################################## #

# Use the dictionary to do it
lab_test_df <- filter(dictionary_murnau, grepl("lab", FormName))
lab_test <- lab_test_df$VariableFieldName
lab_test_df <- select(df_murnau, any_of(lab_test))

# Remove the columns that only contain NAs
lab_test_df_no_na_col <- Filter(function(x) !all(is.na(x)), lab_test_df)

# Replace all values that are empty strings with NAs
lab_test_df_no_na_col <- mutate_all(lab_test_df_no_na_col, list(~na_if(.,"")))

# Only select the columns that contains tests done after the injury
lab_test_df_no_na_col <- lab_test_df_no_na_col[, str_contains("neg", 
                                                              names(lab_test_df_no_na_col), switch = TRUE) == F]

# Find all unique tests 
test_names <- sub("[0-9]+.*", "", names(lab_test_df_no_na_col))
test_names_unique <- unique(test_names)
# Remove the last character: _
test_names_unique <- substr(test_names_unique, 1, nchar(test_names_unique)-1)


########################################################################################################## #
#### Prepare the dataset for ANOVA: CK-MB ####       
########################################################################################################## #
# Select all tests for creatinin kinase MB
df_ck <- lab_test_df_no_na_col %>% select(starts_with("ck_mb") & !contains("prozent"))

# Select all tests done in the first 2 weeks
days <- extract_days(names(df_ck), T)
# Select the columns whose day is less than 7
inds <- as.numeric(days) <= 7
df_ck_1_week <- df_ck[, inds]
# Order the columns by mixed order (numbers treat as numbers)
df_ck_1_week <- df_ck_1_week[, mixedorder(names(df_ck_1_week))]

# Count the total number of tests done for each patient
df_ck_1_week <- df_ck_1_week %>% 
  mutate_all(function(x) ifelse(is.na(x), 0, 1)) %>%
  mutate(Grade = return_grade()) %>% 
  mutate(N_tests = select(., ck_mb_0a:ck_mb_7a) %>% rowSums()) %>% 
  mutate(Patient = 1:363)


########################################################################################################## #
#### Test if CK-MB fulfills ANOVA assumptions ####       
########################################################################################################## #
df_counts_a <- df_ck_1_week %>% filter(Grade == "A") %>% select(c("Patient", "Grade", "N_tests"))
df_counts_b <- df_ck_1_week %>% filter(Grade == "B") %>% select(c("Patient", "Grade", "N_tests"))
df_counts_c <- df_ck_1_week %>% filter(Grade == "C") %>% select(c("Patient", "Grade", "N_tests"))
df_counts_d <- df_ck_1_week %>% filter(Grade == "D") %>% select(c("Patient", "Grade", "N_tests"))

# Check if the highest variance is less than 4 times the smallest variance
var_a <- var(df_counts_a$N_tests)
var_b <- var(df_counts_b$N_tests)
var_c <- var(df_counts_c$N_tests)
var_d <- var(df_counts_d$N_tests)
variances <- c(var_a, var_b, var_c, var_d)
max(variances) < 4 * min(variances)
# TRUE

# Plot the boxplots of the four distributions
df_ck_1_week %>% filter(Grade %in% c("A", "B", "C", "D")) %>% 
  ggplot(aes(x = Grade, y = N_tests, fill = Grade)) +
  geom_boxplot() +
  geom_point(aes(fill = Grade), size = 1, shape = 21, position = position_jitterdodge()) +
  labs(title = "Distribution of the n° per patient according \n to their AIS grade: CK-MB", 
       x = "", y = "N° tests in the first 2 weeks")
ggsave(filename = "boxplot_ck_mb.png", path = "~/Plots/Plots_all_LM_stat_tests")

# Density plot to check whether the distributions are bell-shaped
a <- ggdensity(df_counts_a$N_tests, main = "Density plot of number of tests for AIS grade A", xlab = "N° tests")
b <- ggdensity(df_counts_b$N_tests, main = "Density plot of number of tests for AIS grade B", xlab = "N° tests")
c <- ggdensity(df_counts_c$N_tests, main = "Density plot of number of tests for AIS grade C", xlab = "N° tests")
d <- ggdensity(df_counts_d$N_tests, main = "Density plot of number of tests for AIS grade D", xlab = "N° tests")
p1 <- plot_grid(plotlist = list(a, b, c, d), nrow = 2, ncol = 2)
save_plot(filename = "counts_density.png", plot = p1, path = "~/Plots/Plots_all_LM_stat_tests",
          nrow = 2, ncol = 2)
# From the plots, none of the groups looks normally distributed.

# Since I'm not sure that plotting a continuous distribution is wise in this case: plot histogram
df_ck_1_week_ <- filter(df_ck_1_week, Grade %in% c("A", "B", "C", "D"))
histogram(~ N_tests | Grade, data = df_ck_1_week_, layout = c(1, 4))

# Q-Q plot: plots one quantile against the other
# If the distribution is normal, the obtained plot will be a straight line, meaning that the quantiles are equal
# (since the density curve is bell-shaped)
a <- ggqqplot(df_counts_a$N_tests, title = "Q-Q plot for AIS grade A")
b <- ggqqplot(df_counts_b$N_tests, title = "Q-Q plot for AIS grade B")
c <- ggqqplot(df_counts_c$N_tests, title = "Q-Q plot for AIS grade C")
d <- ggqqplot(df_counts_d$N_tests, title = "Q-Q plot for AIS grade D")
p1 <- plot_grid(plotlist = list(a, b, c, d), nrow = 2, ncol = 2)
save_plot(filename = "counts_qq_plot.png", plot = p1, path = "~/Plots/Plots_all_LM_stat_tests",
          nrow = 2, ncol = 2)
# Also from the Q-Q plot the distributions don't look normal

# Statistical test: Shapiro-Wilk's method
shapiro.test(df_counts_a$N_tests)
shapiro.test(df_counts_b$N_tests)
shapiro.test(df_counts_c$N_tests)
shapiro.test(df_counts_d$N_tests)
# From the Shapiro-Wilk test of normality, only the hypothesis that the distribution of B is normal
# is not rejected.
# However, the size of group B is fairly small and the plots suggest that B is not normal.


########################################################################################################## #
#### ANOVA: creatinin kinase MB ####       
########################################################################################################## #
df_ck_an <- df_ck_1_week %>% filter(Grade %in% c("A", "B", "C", "D"))
one.way <- aov(N_tests ~ Grade, data = df_ck_an)
summary(one.way)

tukey.one.way <- TukeyHSD(one.way)
tukey.one.way


########################################################################################################## #
#### Kolmogorov-Smirnov test: creatinin kinase MB ####       
########################################################################################################## #
# Non-paramteric pairwise test
# Null hypothesis: the groups were sampled from populations with identical distributions
# Less power than Mann-Whitney test in detecting shifts in the median between the two groups
# More power in detecting shifts in the shape of the distribution
# Doesn't really work for discrete distributions...
ks.test(df_counts_a$N_tests, df_counts_d$N_tests)

# Mann-Whitney U test: requires that the shape of the distribution are similar, if we want to compare medians
# insted of comparing mean ranks


