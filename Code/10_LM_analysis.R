# Load packages ----
library(dplyr)
library(gtools)
library(stringr)
library(ggplot2)
library(tibble)
library(tidyr)
library(corrplot)
library(data.table)

# Load the data ----
# Murnau dataset
df_murnau <- read.csv("~/Data/HematologicalBiomark_DATA_2021-01-07_0729.csv",
                      na.strings = c(""))

# Murnau dictionary
path <- "~/Data/HematologicalBiomarkersForSCI_DataDictionary_2019-10-01.csv"
dictionary_murnau <- read.csv(path)
# Change columns names of the dictionary (remove dots)
names(dictionary_murnau) <- gsub("\\.", "", names(dictionary_murnau))

# Murnau patients information
df_demog <- read.csv("~/Data/demographic_info.csv")

# Load the 2 weeks dataset
df_2_weeks <- read.csv("~/Lab_tests_two_weeks.csv")


########################################################################################################## #
#### Extract the 10 middle markers ####       
########################################################################################################## #
list_names_middle_markers <- c("a_luc", "ckmb_prozent", "ck_mb", "a_basophile", "a_eosinophile", "a_monocyten",
                               "a_lymphocyten", "atrophile", "nrbc", "gesamt_ck")
list_middle_markers <- c("LUC", "Creatin kinase MB", "Creatin kinase MB %", "Basophile", "Eosinophile",
                         "Monocytes", "Lymphocytes", "Atrophile", "NRBC", "Total CK")

df_middle_markers <- dplyr::select(df_2_weeks, contains(list_names_middle_markers))


########################################################################################################## #
#### Check if white blood cells are always tested together ####       
########################################################################################################## #
# Extract the tests for white blood cells from the dataset
df_wbcs <- df_middle_markers %>% 
  dplyr::select(contains(c("a_basophile", "a_eosinophile", "a_monocyten", "a_lymphocyten"))) %>% 
  add_column(Patient_number = df_demog$patientennummer, .before = "a_basophile_0a")

# Remove patients that were never tested for WBCs
# There are 237 patients that were tested at least once for one WBC
df_wbcs_tested <- df_wbcs[rowSums(is.na(df_wbcs)) != (ncol(df_wbcs)-1), ]

# Melt the dataset and add columns day and marker, in order to facilitate the following pivot_wider
df_wbcs_tested_molten <- df_wbcs_tested %>% pivot_longer(cols = 2:101) %>%
  mutate(DayTime = sub(".*_", "", name)) %>% 
  mutate(Marker = sub("_[^_]+$", "", name)) %>% 
  dplyr::select(-name)

df_wbcs_tested_final <- pivot_wider(df_wbcs_tested_molten, id_cols = c("Patient_number", "DayTime"),
                                    names_from = Marker, values_from = value)

# Dataframe recording the times in which all 4 white blood cells were tested
# Done by removing all rows that contain at least one NA (because if they contain a NA value, is for sure in
# either a_basophile, a_eosinophile, a_monocyten or a_lymphocyten columns)
df_wbcs_tested_all <- na.omit(df_wbcs_tested_final)

# Dataframe made by discarding the rows that contain only NAs in columns a_basophile, a_eosinophile,
# a_monocyten and a_lymphocyten.
# Check if the size of this dataset is the same as the size of the previous dataset.
df_wbcs_tested_some <- df_wbcs_tested_final[rowSums(is.na(df_wbcs_tested_final[c("a_basophile", 
                                                                                 "a_eosinophile", 
                                                                                 "a_monocyten", 
                                                                                 "a_lymphocyten")])) != 4, ]

# There is one patient that was once tested for only some WBCs (but not all). Check which one:
df_wbcs_tested_one <- df_wbcs_tested_some[rowSums(is.na(df_wbcs_tested_some[c("a_basophile", 
                                                                                 "a_eosinophile", 
                                                                                 "a_monocyten", 
                                                                                 "a_lymphocyten")])) > 0, ]

## Conclusion: white blood cells are always tested together, except for one particular case, in which only
## lymphocytes were tested. This happened on day 4 after the injury and it was the second test of the day.

## NOTE: in the following analysis, this test will be discarded (because it's kind of weird).
## Remove datasets that aren't useful anymore
rm(df_wbcs_tested_one, df_wbcs_tested_some, df_wbcs_tested_molten)


########################################################################################################## #
#### WBCs analysis: when tests are performed ####       
########################################################################################################## #
# Plot one histogram
# To have day and time ordered in the graph
df_plot <- df_wbcs_tested_all %>%
  mutate(DayTime_ord = factor(DayTime, levels = mixedsort(unique(DayTime)))) %>% 
  mutate(Time = str_sub(DayTime, - 1, - 1))

ggplot(data = df_plot, aes(x = DayTime_ord)) +
  geom_bar(width = 0.7, aes(fill = Time, color = Time)) +
  scale_fill_manual(values = c("dodgerblue", "pink"),
                    name = "Time in day", 
                    labels = c("First", "Second")) +
  scale_color_manual(values = c("blue", "magenta"),
                    name = "Time in day", 
                    labels = c("First", "Second")) +
  geom_text(stat='count', aes(label=..count..), vjust = -0.3, size = 2.5) +
  labs(title = "N of patients tested for WBCs per day", x = "Day", y = "N of patients tested for WBCs") +
  theme(legend.key.size = unit(0.2, "cm"), legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.position="bottom")
ggsave("tests_per_day.png", path = "~/Plots/Plots_10_LM")

# We can see that on day 0 there are not many patients that are tested. This is likely due to the fact that
# surgeons are more focused on keeping patients alive, rather than testing if they have/could develop
# infections.
# The highest number of tests is at day 4: infections had the time to develop and have effect on the body.
# Also, many patients have undergone surgeries and are staying in a hospital, thus they could develop
# infections for this reason.


########################################################################################################## #
#### WBCs analysis: when tests are performed in each patient ####       
########################################################################################################## #
df_plot <- df_wbcs_tested_final %>% dplyr::select(Patient_number, DayTime, a_lymphocyten) %>% 
  mutate(a_lymphocyten = ifelse(is.na(a_lymphocyten), 0, 1))
setnames(df_plot, "Patient_number", "patientennummer")
df_plot <- merge(x = df_plot, y = dplyr::select(df_demog, c("patientennummer", "Grade")), 
                 by = "patientennummer", all.x = T)

df_plot$DayTime <- factor(df_plot$DayTime, levels = unique(mixedsort(df_plot$DayTime)))
df_plot %>% filter(Grade %in% c("A", "B", "C", "D")) %>%
  ggplot(aes(x = DayTime, y = a_lymphocyten, group = patientennummer)) +
  geom_line(size = 0.5, alpha = 0.3, position = position_dodge(width=0.1)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  facet_wrap(~ Grade, ncol = 2) +
  labs(title = "N° of tests per day and per patient", x = "Day", y = "N° of tests")
ggsave(filename = "counts_patients_day2.png", path = "~/Plots/Plots_10_LM")


########################################################################################################## #
#### WBCs analysis: which patients ####       
########################################################################################################## #
# Start with a demographic analysis: characteristics of patients tested n times
df_patients_wbcs <- plyr::count(df_wbcs_tested_all$Patient_number)
# Change first column name to patientennumber to enable merging in the following step
names(df_patients_wbcs) <- c("patientennummer", "N_tests_wbcs")

# Perform a left outer join between df_patients_wbcs and df_demog: all rows of df_patients_wbcs are kept
# and any rows from the df_demog that match the key. In this case, the left join is the same as the inner
# join.
df_patients_wbcs_demog <- merge(x = df_patients_wbcs, y = df_demog, by = "patientennummer", all.x = TRUE)

df_patients_wbcs_demog %>% filter(Grade %in% c("A", "B", "C", "D")) %>% 
  ggplot(aes(x = cut(AgeAtDOI, breaks = seq.int(0, 100, 5), right = FALSE),
                                          y = N_tests_wbcs)) + 
  geom_bin_2d(stat = "bin2d", position = "identity", na.rm = FALSE,) +
  stat_bin2d(geom = "text", aes(label = ..count..)) +
  scale_fill_gradient(low = "turquoise", high = "blue") +
  facet_wrap(~ Grade, ncol=2) +
  labs(x = "Age")
ggsave("age_grade_heatmap.png", path = "~/Plots/Plots_10_LM", width = 8, height = 6,
       units = "in")

df_patients_wbcs_demog %>% filter(Grade %in% c("A", "B", "C", "D")) %>% 
  ggplot(aes(x = cut(AgeAtDOI, breaks = seq.int(0, 100, 5), right = FALSE),
             y = N_tests_wbcs)) + 
  geom_bin_2d(stat = "bin2d", position = "identity", na.rm = FALSE,) +
  stat_bin2d(geom = "text", aes(label = ..count..)) +
  scale_fill_gradient(low = "turquoise", high = "blue") +
  facet_wrap(~ Sex, ncol=2) +
  labs(x = "Age")
ggsave("age_sex_heatmap.png", path = "~/Plots/Plots_10_LM", width = 8, height = 6,
       units = "in")

df_patients_wbcs_demog %>% filter(Grade %in% c("A", "B", "C", "D")) %>% filter(!is.na(Plegia)) %>% 
  ggplot(aes(x = cut(AgeAtDOI, breaks = seq.int(0, 100, 5), right = FALSE),
             y = N_tests_wbcs)) + 
  geom_bin_2d(stat = "bin2d", position = "identity", na.rm = FALSE,) +
  stat_bin2d(geom = "text", aes(label = ..count..)) +
  scale_fill_gradient(low = "turquoise", high = "blue") +
  facet_wrap(~ Plegia, ncol=2) +
  labs(x = "Age")
ggsave("age_plegia_heatmap.png", path = "~/Plots/Plots_10_LM", width = 8, height = 6,
       units = "in")

df_patients_wbcs_demog %>% filter(Grade %in% c("A", "B", "C", "D")) %>% filter(!is.na(Plegia)) %>% 
  ggplot(aes(x = Grade, y = N_tests_wbcs)) + 
  geom_bin_2d(stat = "bin2d", position = "identity") +
  stat_bin2d(geom = "text", aes(label = ..count..)) +
  scale_fill_gradient(low = "turquoise", high = "blue") +
  facet_wrap(~ Plegia, ncol=2)
ggsave("plegia_heatmap.png", path = "~/Plots/Plots_10_LM")


########################################################################################################## #
#### Check that there is no significant difference between tested and not tested patients for WBCs ####       
########################################################################################################## #
# Define the two datasets to work with:
## Tested: patients that were tested at least once for WBCs
## Not_tested: patients that were never tested for WBCs
tested_df <- df_patients_wbcs_demog
temp_df <- df_wbcs[rowSums(is.na(df_wbcs)) == (ncol(df_wbcs)-1), ] %>% select(Patient_number)
names(temp_df) <- "patientennummer"
not_tested_df <- merge(x = temp_df, y = df_demog, by = "patientennummer", all.x = TRUE)

# --
# 1. In the age distribution: t-test
## First, we need to check that the two distributions are roughy normally distributed
ggplot(data = tested_df, aes(x = cut(AgeAtDOI, breaks = seq.int(0, 100, 5), right = FALSE))) +
  geom_bar(width = 0.7, fill = "turquoise", na.rm = T) +
  geom_text(stat = 'count', aes(label=..count..), vjust = -0.3, size = 2.5) +
  labs(title = "Age distribution of patients that were tested for WBCs at least once",
       x = "Age",
       y = "N patients")

ggplot(data = not_tested_df, aes(x = cut(AgeAtDOI, breaks = seq.int(0, 100, 5), right = FALSE))) +
  geom_bar(width = 0.7, fill = "turquoise", na.rm = T) +
  geom_text(stat = 'count', aes(label=..count..), vjust = -0.3, size = 2.5) +
  labs(title = "Age distribution of patients that were tested for WBCs at least once",
       x = "Age",
       y = "N patients")

df_plot <- tested_df %>% select(AgeAtDOI) %>% mutate(Tested = rep("yes", nrow(.))) %>% 
  rbind(., (not_tested_df %>% select(AgeAtDOI) %>% mutate(Tested = rep("no", nrow(.)))))

ggplot(data = df_plot, aes(x = cut(AgeAtDOI, breaks = seq.int(0, 100, 5), right = FALSE), fill = Tested)) +
  geom_bar(width = 0.7, position = "dodge", na.rm = T) +
  geom_text(stat = 'count', aes(label=..count..), position = position_dodge(width = 0.7),
            vjust=-0.25, size = 2.5) +
  labs(title = "Age distribution distinguishing patients that were never tested for WBCs",
       x = "Age",
       y = "N patients")
ggsave("age_distributions.png", path = "~/Plots/Plots_10_LM",
       width = 8, height = 4, units = "in")

## T-test assuming normality, which is not really the case. However, in the paper they used it for age
## distributions comparison.
wilcox.test(AgeAtDOI ~ Tested, data = df_plot, alternative = "two.sided")
t.test(not_tested_df$AgeAtDOI, tested_df$AgeAtDOI, alternative = c("two.sided"))
## The p-value is > 0.05, thus we can conclude that there is no significant difference between the age
## distribution of patients that were tested for WBCs and patients that were not tested for WBCs.

# --
# 2. In the sex distribution: chi-squared
sex_tested <- tested_df %>% plyr::count("Sex") %>% column_to_rownames(var = "Sex")
names(sex_tested) <- "tested"
sex_not_tested <- not_tested_df %>% plyr::count("Sex") %>% column_to_rownames(var = "Sex")
names(sex_not_tested) <- "not_tested"

## Contingency table for the chi-squared test
df_chi <- cbind(sex_tested, sex_not_tested)

## Perform the actual chi-squared test.
## Null hypothesis: the two variables, i.e. sex and test for WBCs are independent. This means that the
## proportion of patients that were tested for WBCs at least once is independent from their sex.
chisq_sex <- chisq.test(df_chi)
chisq_sex$p.value
## The p-value = 0.0007901214, thus there is a dependence between test for WBCs and sex.

chisq_sex$residuals
corrplot(chisq_sex$residuals, is.cor = F)
# Females are more prone to be tested than males??

# --
# 3. In the plegia distribution: same test as for the sex distribution
plegia_tested <- tested_df %>% filter(!is.na(Plegia)) %>%
  plyr::count("Plegia") %>%
  column_to_rownames(var = "Plegia")
names(plegia_tested) <- "tested"

plegia_not_tested <- not_tested_df %>% filter(!is.na(Plegia)) %>%
  plyr::count("Plegia") %>%
  column_to_rownames(var = "Plegia")
names(plegia_not_tested) <- "not_tested"

## Contingency table for the chi-squared test
df_chi_plegia <- cbind(plegia_tested, plegia_not_tested)

## Perform the actual chi-squared test.
## Null hypothesis: the two variables, i.e. plegia and test for WBCs are independent. This means that the
## proportion of patients that were tested for WBCs at least once is independent from whether they are
##  para/tetraplegic.
chisq_plegia <- chisq.test(df_chi_plegia, correct = F)
chisq_plegia$p.value
## When Yates correction is not used, the p-value = 0.0461277, thus we can conclude that plegia and test for
## WBCs are dependent variables.

chisq_plegia$residuals
corrplot(chisq_plegia$residuals, is.cor = F)
# Tetraplegic patients are more prone to be tested than paraplegic patients??

# --
# 4. In the grade distribution: ASIA grade is a categorical variable: chi-squared test can be used
grade_tested <- tested_df %>% filter(Grade %in% c("A", "B", "C", "D")) %>% 
  plyr::count("Grade") %>% 
  column_to_rownames(var = "Grade")
names(grade_tested) <- "tested"

grade_not_tested <- not_tested_df %>% filter(Grade %in% c("A", "B", "C", "D")) %>% 
  plyr::count("Grade") %>% 
  column_to_rownames(var = "Grade")
names(grade_not_tested) <- "not_tested"

## Contingency table for the chi-squared test
df_chi_grade <- cbind(grade_tested, grade_not_tested)

## Perform the actual chi-squared test.
## Null hypothesis: the two variables, i.e. grade and test for WBCs are independent. This means that the
## proportion of patients that were tested for WBCs at least once is independent from their grade.
chisq_grade <- chisq.test(df_chi_grade, correct = F)
chisq_grade$p.value

chisq_grade$residuals
corrplot(chisq_grade$residuals, is.cor = F)
