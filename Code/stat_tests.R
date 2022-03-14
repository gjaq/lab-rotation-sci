# Load packages ----
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)

# Load the data ----
# Murnau dataset
df_murnau <- read.csv("~/Data/HematologicalBiomark_DATA_2021-01-07_0729.csv",
                      na.strings = c(""))

# Murnau dictionary
path <- "~/Data/HematologicalBiomarkersForSCI_DataDictionary_2019-10-01.csv"
dictionary_murnau <- read.csv(path)
# Change columns names of the dictionary (remove dots)
names(dictionary_murnau) <- gsub("\\.", "", names(dictionary_murnau))

# Load the 2 weeks dataset
df_2_weeks <- read.csv("~/Data/Lab_tests_two_weeks.csv")


########################################################################################################## #
#### Useful functions ####                                                                                   
########################################################################################################## #

rm_time_in_day = function(test_names) {
  # Remove the last character of each name, indicating whether this is the first, second, ...
  # test done on one specific day
  test_names <- gsub(".{1}$", "", test_names)
  return(test_names)
}

extract_unique_times = function(col_names) {
  times <- sub(".*[0-9]", "", col_names)
  return(unique(times))
}

return_grade = function(df_murnau_ = df_murnau) {
  df <- data.frame(Grade_va = df_murnau_$va_ais, Grade_ai = df_murnau_$ai_ais)
  
  df <- df %>% mutate(Grade = ifelse(Grade_va == "ND" | is.na(Grade_va), Grade_ai, Grade_va))
  return(df$Grade)
}

return_level = function(df_murnau_ = df_murnau) {
  df <- data.frame(Level_va = df_murnau_$va_nli, Level_ai = df_murnau_$ai_nli)
  
  df <- df %>% mutate(Level = ifelse(Level_va == "ND" | is.na(Level_va), Level_ai, Level_va))
  return(df$Level)
}


########################################################################################################## #
#### Dataset: count number of blood draws per day and per patient ####       
########################################################################################################## #
# Create a binary dataset and remove column "ausg_ek_1h" (since no "ausg_ek_1g")
df_2_weeks_bin <- df_2_weeks %>% select(-ausg_ek_1h) %>% 
  mutate_all(function(x) ifelse(is.na(x), 0, 1))

# Make one dataset per day
make_df_day = function(day, df = df_2_weeks_bin) {
  return(df[, as.numeric(gsub(".*_(\\d+).*", "\\1", names(df))) == day])
}

# Make one dataset per time
make_df_time = function(df_day, time) {
  cols_to_select <- names(df_day)[sub(".*[0-9]", "", names(df_day)) == time]
  return(df_day %>% select(all_of(cols_to_select)))
}

# Function that checks if all variables in the second, third test, ... are also in the previous ones
check_if_in_previous = function(df) {
  times <- extract_unique_times(names(df))
  n <- length(times)
  was_tested <- rep(F, (n-1))
  
  for(i in 1:(n-1)) {
    was_tested[i] <- (length(setdiff(rm_time_in_day(names(make_df_time(df, times[i+1]))),
                             rm_time_in_day(names(make_df_time(df, times[i]))))) == 0)
  }
  return(was_tested)
}

# Count blood draws per day
count_blood_draws_per_day <- function(df) {
  times <- extract_unique_times(names(df))
  n <- length(times)
  tested <- matrix(0, 363, n)
  
  # For each times of the day, create a separate dataset to then fill the corresponding column
  for (i in 1:n) {
    time <- times[i]
    df_time <- make_df_time(df, time) %>% 
      mutate(Total = rowSums(.)) %>% 
      mutate(Test = ifelse(Total == 0, 0, 1))
    tested[, i] <- df_time$Test
  }
  return(rowSums(tested))
}

# Make a dataframe for each day
days <- 1:14
list_df_per_day <- mapply(make_df_day, days)

# Check if a marker tested at time n was also tested at times 1:n-1
checks <- sapply(list_df_per_day, check_if_in_previous)
# Yes, they are

# Make a dataframe that, for each day, records the number of blood draws per patient
n_test_day <- as.data.frame(sapply(list_df_per_day, count_blood_draws_per_day))
names(n_test_day) <- as.character(1:14)

# Add columns patient, grade and NLI and remove patients that don't have a grade (or have grade E) or NLI
n_test_day <- n_test_day %>% add_column(Patient = 1:363, .before = "1") %>% 
  add_column(Grade = return_grade(), .before = "1") %>% 
  add_column(NLI = return_level(), .before = "1") %>% 
  filter(!Grade %in% c("ND", "E")) %>% 
  filter(!is.na(Grade)) %>% 
  filter(!is.na(NLI)) %>% 
  filter(!NLI %in% c("NA", "INT"))

# Add a column plegia
tetra <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "T1")
n_test_day <- n_test_day %>% 
  add_column(Plegia = (ifelse(n_test_day$NLI %in% tetra, "tetra", "para")), .before = "1")


########################################################################################################## #
#### Plot histograms of count each day ####       
########################################################################################################## #
n_test_day_per_grade <- n_test_day %>% select(-c("Patient", "NLI")) %>% 
  group_by(Grade, Plegia) %>% 
  mutate(across("1":"14", mean)) %>% 
  distinct() %>%
  pivot_longer(cols = "1":"14", names_to = "Day", values_to = "N_tests")

# Plot the histogram
ggplot(data = n_test_day_per_grade, aes(x = as.numeric(Day), y = N_tests, fill = Plegia)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ Grade, ncol = 2) +
  labs(title = "Mean counts per grade and per day")
ggsave(filename = "mean_counts.png", path = "~/Plots/Plots_stat_tests")


########################################################################################################## #
#### Visual check of normality ####       
########################################################################################################## #
n_test_tot <- n_test_day %>% mutate(Total = rowSums(select(., "1":"14")))

# Plot the boxplots per each grade
ggplot(data = n_test_tot, aes(x = Grade, y = Total, fill = Grade)) +
  geom_violin() +
  geom_point(aes(fill = Grade), size = 1, shape = 21, position = position_jitterdodge()) +
  labs(title = "Distribution of the n° of blood draws per patient \n in the first 2 weeks according to their AIS grade", 
       x = "", y = "N° blood draws in the first 2 weeks")
ggsave(filename = "violin_tot.png", path = "~/Plots/Plots_stat_tests")

# Plot the boxplots dividing patient with tetra and paraplegia
ggplot(data = n_test_tot, aes(x = Grade, y = Total, fill = Grade)) +
  geom_violin() +
  geom_point(aes(fill = Grade), size = 1, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~ Plegia, nrow = 2) +
  labs(title = "Distribution of the n° of blood draws per patient \n in the first 2 weeks according to their AIS grade", 
       x = "", y = "N° blood draws in the first 2 weeks")
ggsave(filename = "violin_plegia_tot.png", path = "~/Plots/Plots_stat_tests")

# Plot the histogram for each group, just to have a sense of the distribution
# Number of bins: 24, like the maximum number of blood draws that a patient had in two weeks
ggplot(data = n_test_tot, aes(x = Total, fill = Grade)) +
  geom_histogram(bins = 24, color = "white") +
  facet_wrap(~ Grade, ncol = 2) +
  labs(title = "Distribution of the n° of blood draws per patient \n in the first 2 weeks according to their AIS grade", 
       x = "N° blood draws in the first 2 weeks", y = "N° patients")
ggsave(filename = "histogram_tot.png", path = "~/Plots/Plots_stat_tests")

# Plot the histogram for each group, just to have a sense of the distribution
# Number of bins: 24, like the maximum number of blood draws that a patient had in two weeks
ggplot(data = n_test_tot, aes(x = Total, fill = Grade)) +
  geom_histogram(bins = 24, color = "white") +
  facet_grid(Grade ~ Plegia) +
  labs(title = "Distribution of the n° of blood draws per patient \n in the first 2 weeks according to their AIS grade", 
       x = "N° blood draws in the first 2 weeks", y = "N° patients")
ggsave(filename = "histogram_plegia_tot.png", path = "~/Plots/Plots_stat_tests")


########################################################################################################## #
#### Shapiro check of normality ####       
########################################################################################################## #
df_counts_a <- n_test_tot %>% filter(Grade == "A") %>% select(c("Patient", "Grade", "Plegia", "Total"))
df_counts_b <- n_test_tot %>% filter(Grade == "B") %>% select(c("Patient", "Grade", "Plegia", "Total"))
df_counts_c <- n_test_tot %>% filter(Grade == "C") %>% select(c("Patient", "Grade", "Plegia", "Total"))
df_counts_d <- n_test_tot %>% filter(Grade == "D") %>% select(c("Patient", "Grade", "Plegia", "Total"))

df_counts_a_p <- df_counts_a %>% filter(Plegia == "para")
df_counts_a_t <- df_counts_a %>% filter(Plegia == "tetra")
df_counts_b_p <- df_counts_b %>% filter(Plegia == "para")
df_counts_b_t <- df_counts_b %>% filter(Plegia == "tetra")
df_counts_c_p <- df_counts_c %>% filter(Plegia == "para")
df_counts_c_t <- df_counts_c %>% filter(Plegia == "tetra")
df_counts_d_p <- df_counts_d %>% filter(Plegia == "para")
df_counts_d_t <- df_counts_d %>% filter(Plegia == "tetra")

# Statistical test: Shapiro-Wilk's method
# Tests if the four distributions are normally distributed
shapiro.test(df_counts_a$Total)
shapiro.test(df_counts_b$Total)
shapiro.test(df_counts_c$Total)
shapiro.test(df_counts_d$Total)

shapiro.test(df_counts_a_p$Total)
shapiro.test(df_counts_a_t$Total)
shapiro.test(df_counts_b_p$Total)
shapiro.test(df_counts_b_t$Total)
shapiro.test(df_counts_c_p$Total)
shapiro.test(df_counts_c_t$Total)
shapiro.test(df_counts_d_p$Total)
shapiro.test(df_counts_d_t$Total)


########################################################################################################## #
#### Kruskal-Wallis test ####       
########################################################################################################## #
kruskal.test(Total ~ Grade, data = n_test_tot)
kruskal.test(Total ~ Plegia, data = n_test_tot)

# Kruskal test only in the paraplegia group
n_test_tot_p <- n_test_tot %>% filter(Plegia == "para")
kruskal.test(Total ~ Grade, data = n_test_tot_p)

# Kruskal test only in the tetraplegia group
n_test_tot_t <- n_test_tot %>% filter(Plegia == "tetra")
kruskal.test(Total ~ Grade, data = n_test_tot_t)


# Check which groups are causing this difference
pairwise.wilcox.test(n_test_tot$Total, n_test_tot$Grade, p.adjust.method = "BH")
pairwise.wilcox.test(n_test_tot_p$Total, n_test_tot_p$Grade, p.adjust.method = "BH")
pairwise.wilcox.test(n_test_tot_t$Total, n_test_tot_t$Grade, p.adjust.method = "BH")


########################################################################################################## #
#### Poisson regression ####       
########################################################################################################## #
n_test_tot <- n_test_tot %>% mutate(Grade = as.factor(Grade)) %>% mutate(Plegia = as.factor(Plegia))
counts_pois <- glm(Total ~ Plegia * Grade, family = "poisson", data = n_test_tot)
counts_pois_t <- glm(Total ~ Grade, family = "quasipoisson", data = n_test_tot_t)
counts_pois_p <- glm(Total ~ Grade, family = "quasipoisson", data = n_test_tot_p)
plot(counts_pois)

library(car)
Anova(counts_pois)
Anova(counts_pois_t)
Anova(counts_pois_p)

library(emmeans)
marginal = emmeans(counts_pois, ~ Grade)
marginal_p = emmeans(counts_pois_p, ~ Grade)
marginal_t = emmeans(counts_pois_t, ~ Grade)
pairs(marginal)
pairs(marginal_p)
pairs(marginal_t)