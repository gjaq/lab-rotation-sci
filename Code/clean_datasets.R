# Load packages ----
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(visdat)
library(naniar)
library(reshape2)
library(RColorBrewer)
library(UpSetR)

#---
# Load the data ----
df_murnau <- read.csv("~/Data/HematologicalBiomark_DATA_2021-01-07_0729.csv",
                      na.strings = c(""))
path <- "~/Data/HematologicalBiomarkersForSCI_DataDictionary_2019-10-01.csv"
dictionary_murnau <- read.csv(path)

# Change columns names of the dictionary (remove dots)
names(dictionary_murnau) <- gsub("\\.", "", names(dictionary_murnau))

#---
# Create the 28 datasets with the 28 routinely assessed blood markers
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

# Amylase ----
# Good, only need to trasform integers to numeric
df_amylase <- df_amylase %>% mutate_all(as.numeric)


# Alkaline phosphatase ----
# Patient 325 in columns 0a is "zu haemolytic"
df_ap$ap_0a[325] <- NA
df_ap$ap_97a[299] <- NA
# Replace "< ..." by "0"
replace_by_zero <- function(x) {
  replace(x, grep("[<]", x), "0")
}
df_ap <- df_ap %>% mutate_all(replace_by_zero)

# Put everything as numeric
df_ap <- df_ap %>% mutate_all(as.numeric)


# Calcium ----
# Patient "zu haemolytic"
df_calcium$calcium_0a[325] <- NA
df_calcium$calcium_97a[299] <- NA

# Some number have comma instead of dots
df_calcium <- df_calcium %>% mutate_all(funs(str_replace(., ",", "\\.")))

# Put everything as numeric
df_calcium <- df_calcium %>% mutate_all(as.numeric)


# Cholinesterase ----
# Patient "zu haemolytic"
df_che$che_97a[299] <- NA

# One value is "<2000". Setting it to 1500 (outside the normal range 5320-12920)
replace_by_ <- function(x) {
  replace(x, grep("[<]", x), "1500")
}
df_che <- df_che %>% mutate_all(replace_by_)

# Put everything as numeric
df_che <- df_che %>% mutate_all(as.numeric)


# Creatinin ----
# Patient "zu haemolytic"
df_creatinin$creatinin_0a[325] <- NA
df_creatinin$creatinin_97a[299] <- NA

df_creatinin <- df_creatinin %>% mutate_all(as.numeric)


# CRP ----
# Patient "zu haemolytic"
df_crp$crp_0a[325] <- NA
df_crp$crp_97a[299] <- NA

# Many values have "<". Replace by zero, but not sure it's correct
df_crp <- df_crp %>% mutate_all(replace_by_zero)

# Put everything as numeric
df_crp <- df_crp %>% mutate_all(as.numeric)


# Gamma GT ----
# Patient "zu haemolytic"
df_gamma_gt$gamma_gt_0a[325] <- NA
df_gamma_gt$gamma_gt_97a[299] <- NA

# Put everything as numeric
df_gamma_gt <- df_gamma_gt %>% mutate_all(as.numeric)


# Total bilirubin ----
# Patient "zu haemolytic"
df_bilirubin$gesamt_bilirubin_97a[299] <- NA

df_bilirubin <- df_bilirubin %>% mutate_all(replace_by_zero)

# Put everything as numeric
df_bilirubin <- df_bilirubin %>% mutate_all(as.numeric)


# Total proteins ----
# Everything's fine
df_gesamteiweiss <- df_gesamteiweiss %>% mutate_all(as.numeric)


# Glucose ----
# Don't know how to deal with "negativ", "++", ...


# Got ----
# Patient "zu haemolytic"
df_got$got_97a[299] <- NA

# One "<3" in one column: convert as zero
df_got <- df_got %>% mutate_all(replace_by_zero)

df_got <- df_got %>% mutate_all(as.numeric)






