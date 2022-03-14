# Load packages ----
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(glmmTMB)
library(DHARMa)
library(MASS)
library(emmeans)
library(car)
library(bbmle)

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

# Load the dataset with all useful demographic information
demog_info <- read.csv("~/Data/demographic_info.csv")


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
df_2_weeks_bin <- df_2_weeks %>% dplyr::select(-ausg_ek_1h) %>%
  mutate_all(function(x) ifelse(is.na(x), 0, 1))

# Make one dataset per day
make_df_day = function(day, df = df_2_weeks_bin) {
  return(df[, as.numeric(gsub(".*_(\\d+).*", "\\1", names(df))) == day])
}

# Make one dataset per time
make_df_time = function(df_day, time) {
  cols_to_select <- names(df_day)[sub(".*[0-9]", "", names(df_day)) == time]
  return(df_day %>% dplyr::select(all_of(cols_to_select)))
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
  add_column(patientennummer = df_murnau$patientennummer, .before = "1") %>% 
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

# Final dataset: add a column recording the total number of tests per each patient in the first two weeks and
# ensure that both the Grade and the Plegia columns are seen as factor variables.
n_test_tot <- n_test_day %>% mutate(Total = rowSums(dplyr::select(., "1":"14"))) %>% 
  mutate(Grade = as.factor(Grade)) %>%
  mutate(Plegia = as.factor(Plegia))


########################################################################################################## #
#### 1. Test the effect of the grade (alone) ####       
########################################################################################################## #
# 1.1 Poisson distribution: less likely to fit to real biological data well because thes are often
#     overdispersed.
gp <- glmmTMB(Total ~ Grade, data = n_test_tot, family = poisson)

# 1.2 Negative binomial fit
gnb1 <- glmmTMB(Total ~ Grade, data = n_test_tot, family = nbinom1)
gnb2 <- glmmTMB(Total ~ Grade, data = n_test_tot, family = nbinom2)

# 1.3 Zero-inflated negative binomial
gpzi <- glmmTMB(Total~Grade, data = n_test_tot, ziformula = ~1, family = poisson)
gnb1zi <- glmmTMB(Total~Grade, data = n_test_tot, ziformula = ~1, family = nbinom1)
gnb2zi <- glmmTMB(Total~Grade, data = n_test_tot, ziformula = ~1, family = nbinom2)

# Model comparison using AIC
bbmle::AICtab(gp, gnb1, gnb2, gpzi, gnb1zi, gnb2zi)

# Residual analysis using DHARMa
n_sim <- 250
simulationOutput <- simulateResiduals(fittedModel = gp, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gnb1, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gnb2, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpzi, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gnb1zi, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gnb2zi, n = n_sim, seed = 123, plot = T)

## COMMENT: both AIC and DHARMa analysis of residuals seem to point out that the zero-inflated model with
## negative binomial distribution is the best among the 6 models tested.


########################################################################################################## #
#### 2. Test the effect of the plegia variable (alone) ####       
########################################################################################################## #
# 2.1 Poisson distribution: less likely to fit to real biological data well because thes are often
#     overdispersed.
pp <- glmmTMB(Total ~ Plegia, data = n_test_tot, family = poisson)

# 2.2 Negative binomial fit
pnb1 <- glmmTMB(Total ~ Plegia, data = n_test_tot, family = nbinom1)
pnb2 <- glmmTMB(Total ~ Plegia, data = n_test_tot, family = nbinom2)

# 2.3 Zero-inflated negative binomial
ppzi <- glmmTMB(Total~Plegia, data = n_test_tot, ziformula = ~1, family = poisson)
pnb1zi <- glmmTMB(Total~Plegia, data = n_test_tot, ziformula = ~1, family = nbinom1)
pnb2zi <- glmmTMB(Total~Plegia, data = n_test_tot, ziformula = ~1, family = nbinom2)

# Model comparison using AIC
bbmle::AICtab(pp, pnb1, pnb2, ppzi, pnb1zi, pnb2zi)

# Residual analysis using DHARMa
n_sim <- 250
simulationOutput <- simulateResiduals(fittedModel = pp, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = pnb1, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = pnb2, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = ppzi, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = pnb1zi, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = pnb2zi, n = n_sim, seed = 123, plot = T)

## COMMENT: from the AIC, the zero-inflated model with nbinom2 seems slightly better than the model with nbinom1
## and the residuals analysis with DHARMa shows the same thing. Nbinom2 has a higher p-value for the KS test, 
## while nbinom1 has higher p-values for outliers and dispersion tests, although all p-values are not 
## significant.


########################################################################################################## #
#### 3. Test the effect of the grade and plegia ####       
########################################################################################################## #
# 3.1 Poisson distribution: less likely to fit to real biological data well because thes are often
#     overdispersed.
gpp <- glmmTMB(Total ~ Grade * Plegia, data = n_test_tot, family = poisson)

# 3.2 Negative binomial fit
gpnb1 <- glmmTMB(Total ~ Grade * Plegia, data = n_test_tot, family = nbinom1)
gpnb2 <- glmmTMB(Total ~ Grade * Plegia, data = n_test_tot, family = nbinom2)

# 3.3 Zero-inflated negative binomial
gppzi <- glmmTMB(Total ~ Grade * Plegia, data = n_test_tot, ziformula = ~1, family = poisson)
gpnb1zi <- glmmTMB(Total ~ Grade * Plegia, data = n_test_tot, ziformula = ~1, family = nbinom1)
gpnb2zi <- glmmTMB(Total ~ Grade * Plegia, data = n_test_tot, ziformula = ~1, family = nbinom2)

# Model comparison using AIC
bbmle::AICtab(gpp, gpnb1, gpnb2, gppzi, gpnb1zi, gpnb2zi)

# Residual analysis using DHARMa
n_sim <- 250
simulationOutput <- simulateResiduals(fittedModel = gpp, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb1, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb2, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gppzi, n = n_sim, seed = 123, plot = T)
simulationOutput_c <- simulateResiduals(fittedModel = gpnb1zi, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb2zi, n = n_sim, seed = 123, plot = T)

plotQQunif(simulationOutput_c)

## COMMENT: the AIC indicates that the best model is the zero-inflated model with nbinom1, while DHARMa has 
## higher p-values for zero-inflated model with nbinom2. However, none of the models passes the Leven test for
## homogeneity of variances.
## Interestingly, in this case, the two negative binomial models without the zero-inflation already pass the 
## KS test.


########################################################################################################## #
#### 4. Model adjusted for age and sex as confounding variables ####       
########################################################################################################## #
# Join all information needed in one dataset
demog_info_red <- demog_info %>% dplyr::select(c("patientennummer", "AgeAtDOI", "Sex", "days_diff"))
n_test_tot_red <- n_test_tot %>% dplyr::select(c("patientennummer", "Grade", "Plegia", "Total"))
model_df <- merge(x = n_test_tot_red, y = demog_info_red, by = "patientennummer")

# 4.1 Poisson distribution: less likely to fit to real biological data well because thes are often
#     overdispersed.
gpp_adj <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex, data = model_df, family = poisson)

# 4.2 Negative binomial fit
gpnb1_adj <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex, data = model_df, family = nbinom1)
gpnb2_adj <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex, data = model_df, family = nbinom2)

# 4.3 Zero-inflated negative binomial
gppzi_adj <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex, model_df, ziformula = ~1, family = poisson)
gpnb1zi_adj <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex, model_df, ziformula = ~1, family = nbinom1)
gpnb2zi_adj <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex, model_df, ziformula = ~1, family = nbinom2)

# Model comparison using AIC
bbmle::AICtab(gpp_adj, gpnb1_adj, gpnb2_adj, gppzi_adj, gpnb1zi_adj, gpnb2zi_adj)

# Residual analysis using DHARMa
n_sim <- 250
simulationOutput <- simulateResiduals(fittedModel = gpp_adj, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb1_adj, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb2_adj, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gppzi_adj, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb1zi_adj, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb2zi_adj, n = n_sim, seed = 123, plot = T)

## COMMENT: again, the zero-inflated model with nbinom1 family is the one that seems to be the best, according
## to both the AIC and the residual analysis with the DHARMa package.


########################################################################################################## #
#### 5. Model adjusted for also days_diff as confounding variable ####       
########################################################################################################## #
# 5.1 Poisson distribution: less likely to fit to real biological data well because thes are often
#     overdispersed.
gpp_adj_d <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex + days_diff, data = model_df, family = poisson)

# 5.2 Negative binomial fit
gpnb1_adj_d <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex + days_diff, data = model_df, family = nbinom1)
gpnb2_adj_d <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex + days_diff, data = model_df, family = nbinom2)

# 5.3 Zero-inflated negative binomial
gppzi_adj_d <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex + days_diff,
                       model_df, ziformula = ~1, family = poisson)
gpnb1zi_adj_d <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex + days_diff,
                         model_df, ziformula = ~1, family = nbinom1)
gpnb2zi_adj_d <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex + days_diff,
                         model_df, ziformula = ~1, family = nbinom2)

# 5.4 Zero-inflated model depending on the days_diff variable
gppzid_adj_d <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex + days_diff,
                        model_df, ziformula = ~days_diff, family = poisson)
gpnb1zid_adj_d <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex + days_diff,
                          model_df, ziformula = ~days_diff, family = nbinom1)
gpnb2zid_adj_d <- glmmTMB(Total ~ Grade * Plegia + AgeAtDOI + Sex + days_diff,
                          model_df, ziformula = ~days_diff, family = nbinom2)

# Model comparison using AIC
bbmle::AICtab(gpp_adj, gpnb1_adj, gpnb2_adj, gppzi_adj, gpnb1zi_adj, gpnb2zi_adj, gppzid_adj_d, gpnb1zid_adj_d,
              gpnb2zid_adj_d)

# Residual analysis using DHARMa
n_sim <- 250
simulationOutput <- simulateResiduals(fittedModel = gpp_adj_d, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb1_adj_d, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb2_adj_d, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gppzi_adj_d, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb1zi_adj_d, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb2zi_adj_d, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gppzid_adj_d, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb1zid_adj_d, n = n_sim, seed = 123, plot = T)
simulationOutput <- simulateResiduals(fittedModel = gpnb2zid_adj_d, n = n_sim, seed = 123, plot = T)


########################################################################################################## #
#### 5. Emmeans analysis with the best model ####       
########################################################################################################## #
emmeans(gpnb1zi_adj, specs = pairwise ~ Grade | Plegia, type = "response")