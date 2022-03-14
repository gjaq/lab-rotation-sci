# Load packages ----
library(plyr)
library(dplyr)
library(gtools)
library(tibble)


# Load the data ----
# Table containing patients and the number of test they did in the first two weeks
n_test_tot <- read.csv("~/Data/n_test.csv")

# Load the table containing all the useful demographic information
demo_info <- read.csv("~/Data/demographic_info.csv")

# Load the Murnau dataset
df_murnau <- read.csv("~/Data/HematologicalBiomark_DATA_2021-01-07_0729.csv",
                      na.strings = c(""))

# Murnau dictionary
path <- "~/Data/HematologicalBiomarkersForSCI_DataDictionary_2019-10-01.csv"
dictionary_murnau <- read.csv(path)
names(dictionary_murnau) <- gsub("\\.", "", names(dictionary_murnau))

########################################################################################################## #
#### Make the dataset of never tested patients, with all necassary information ####       
########################################################################################################## #
never_tested <- n_test_tot %>% 
  dplyr::select(c("patientennummer", "Total")) %>% 
  merge(x = ., y = demo_info, by = "patientennummer") %>% 
  filter(Total == 0)


########################################################################################################## #
#### Check if those patients are tested after the first two weeks ####       
########################################################################################################## #
test_names <- dictionary_murnau %>% filter(grepl("lab", FormName)) %>% 
  filter(!grepl("neg", FormName))

never_tested_all <- df_murnau %>% 
  filter(patientennummer %in% never_tested$patientennummer) %>%
  dplyr::select(all_of(append(test_names$VariableFieldName, "patientennummer"))) %>% 
  mutate(N_test = rowSums(!is.na(.)) - 1)

never_tested_all <- Filter(function(x) !all(is.na(x)), never_tested_all)
never_tested_all <- never_tested_all[, mixedorder(names(never_tested_all))] %>% 
  dplyr::select(patientennummer, N_test, everything()) %>% 
  add_column(Grade = never_tested$Grade, .after = "patientennummer") %>% 
  add_column(Plegia = never_tested$Plegia, .after = "patientennummer") %>% 
  add_column(Age_at_DOI = never_tested$Age_at_DOI, .after = "patientennummer") %>% 
  add_column(Sex = never_tested$Sex, .after = "patientennummer") %>% 
  add_column(days_diff = never_tested$days_diff, .after = "N_test")

# Among the 14 patients (that had a value for the NLI and the grade, that was not E) that were never tested 
# in the first two weeks, there are five patients that were never tested at all (B:1, C: 1, D: 2).
# Among the other, there are three patients that were tested less than 100 times (A: 1, B: 1, D: 1), four
# patients that were tested between 310 and 468 times (A: 2, D: 2), one patient with grade D that was tested 
# 936 times and one patient that was tested 1526 times.

library(xlsx)
never_tested <- never_tested %>% mutate(N_test = never_tested_all$N_test)
write.xlsx(never_tested, "~/Data/patients_never_tested.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)



########################################################################################################## #
#### See which markers are tested in the most tested patient ####       
########################################################################################################## #
most_tested <- never_tested_all %>% filter(N_test > 1000) %>% select_if(!is.na(.))
markers <- sub("_[^_]+$", "", names(most_tested)[6:1531])
markers_table <- plyr::count(markers)

# The markers tested in this patient are the following:
#   - WBCs: about 30 times
#   - 28 markers: ranging from 28 (amylase) to blood urea (65)
#   - Other markers: tested a couple of times



