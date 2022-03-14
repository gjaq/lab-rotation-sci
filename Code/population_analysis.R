# Load packages ----
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)

#---
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
df_demog <- read_excel('~/Data/qry_CATHERINE_Demogr_60000.xlsx')
df_age <- read_excel('~/Data/qry_CATHERINE_Age_60000.xlsx')
df_age_up <- read.csv("~/Data/df_age_updated_murnau.csv")


########################################################################################################## #
#### Build the dataset for the demographic analysis ####  
########################################################################################################## #
# Use the dictionary to select interesting columns in the Murnau dataset
# The ones containing "very_acute" or "personel_information"
to_select <- dictionary_murnau %>% 
  filter(grepl('very_acute|information', FormName))

murnau_demog <- df_murnau %>% dplyr::select(all_of(to_select$VariableFieldName)) %>% 
  dplyr::select(!contains("date"))

# Verify that all patients in the dataset with population information are the same as the patients in the
# Murnau dataset of hematological markers
setequal(murnau_demog$patientennummer, df_demog$Patientennummer)
# Yes, they are
# We can thus add the Age and Sex columns from the demographic datset to the murnau dataset
murnau_demog <- murnau_demog %>% mutate(Age = df_demog$AgeAtDOI, Sex = df_demog$Sex)


########################################################################################################## #
#### Check general characteristics ####  
########################################################################################################## #
# Sex distribution
sex_counts <- table(murnau_demog$Sex)

# Age distribution
age_mean <- mean(murnau_demog$Age, na.rm = T)
age_sd <- sd(murnau_demog$Age, na.rm = T)
age_min <- min(murnau_demog$Age, na.rm = T)
age_max <- max(murnau_demog$Age, na.rm = T)

# Make a plot showing how many people per age
ggplot(data = df_demog, aes(x = cut(AgeAtDOI, breaks = seq.int(0, 100, 5), right = FALSE))) +
  geom_bar(width = 0.7, fill = "turquoise") +
  geom_text(stat = 'count', aes(label=..count..), vjust = -0.3, size = 2.5) +
  labs(title = "Age distribution among patients", x = "Age", y = "N patients")
ggsave("age_distribution.png", path = "~/Plots/Plots_population_analysis",
       width = 8, height = 4, units = "in")

# Injury level
level_na <- murnau_demog %>% filter(is.na(va_nli)) %>% nrow()
level_counts <- table(murnau_demog$va_nli)

# Motor score
upper_na <- murnau_demog %>% filter(is.na(va_uems)) %>% nrow()
upper_nd <- murnau_demog %>% filter(va_uems == "ND") %>% nrow()
lower_na <- murnau_demog %>% filter(is.na(va_lems)) %>% nrow()
lower_nd <- murnau_demog %>% filter(va_lems == "ND") %>% nrow()

upper <- murnau_demog %>% filter(va_uems != "ND")
upper_mean <- mean(as.numeric(upper$va_uems), na.rm = T)
upper_sd <- sd(as.numeric(upper$va_uems), na.rm = T)
lower <- murnau_demog %>% filter(va_lems != "ND")
lower_mean <- mean(as.numeric(lower$va_lems), na.rm = T)
lower_sd <- sd(as.numeric(lower$va_lems), na.rm = T)

return_level = function(df_murnau_ = df_murnau) {
  df <- data.frame(Level_va = df_murnau_$va_nli, Level_ai = df_murnau_$ai_nli)
  
  df <- df %>% mutate(Level = ifelse(Level_va == "ND" | is.na(Level_va), Level_ai, Level_va))
  return(df$Level)
}

level_counts <- table(return_level())


########################################################################################################## #
#### Construct a general demographic dataset ####  
########################################################################################################## #
# Columns:
#   - Patient number
#   - Age
#   - Sex
#   - Date of injury
#   - Admission and demission (two separate columns)
#   - Grade and stage at which the grade was taken: very acute when possible, acute I when grade at the very
#     acute stage was missing or not determined
#   - Level and stage at which the level was taken: very acute when possible, acute I when level at the very 
#     acute stage was missing or not determined
tetra <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "T1")

# Prepare the two datasets for merging
df1 <- df_murnau %>% dplyr::select(c("patientennummer", "va_ais", "ai_ais", "va_nli", "ai_nli",
                                     "date_of_injury", "admission", "demission", "geburtstag")) %>%
  mutate(year_of_injury = as.integer(gsub("-.*", "", date_of_injury))) %>% 
  mutate(AgeAtDOI = year_of_injury - geburtstag) %>% 
  mutate(days_diff = as.Date(admission, format = "%Y-%m-%d") - as.Date(date_of_injury, format = "%Y-%m-%d"))

df2 <- df_demog %>% dplyr::select(c("Patientennummer", "Sex"))
names(df2) <- c("patientennummer", "Sex")

# Merge the two datasets
murnau_demog2 <- merge(x = df1, y = df2, by = "patientennummer")
murnau_demog_final <- murnau_demog2 %>% 
  mutate(Grade = ifelse(va_ais == "ND" | is.na(va_ais), ai_ais, va_ais)) %>% 
  mutate(Grade_given_at = ifelse(va_ais == "ND" | is.na(va_ais), "acute_i", "very_acute")) %>%
  mutate(NLI = ifelse(va_nli == "ND" | is.na(va_nli), ai_nli, va_nli)) %>% 
  mutate(NLI_given_at = ifelse(va_nli == "ND" | is.na(va_nli), "acute_i", "very_acute")) %>% 
  mutate(Plegia = (ifelse(is.na(NLI), NA, ifelse(NLI %in% tetra, "tetra", "para"))))
murnau_demog_final$Sex[murnau_demog_final$Sex == "w"] <- "f"

# Save the dataframe as a csv table
write.csv(murnau_demog_final, "~/Data/demographic_info.csv", row.names = FALSE)






