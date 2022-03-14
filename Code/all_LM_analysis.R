# Load packages ----
library(dplyr)
# library(plyr)
library(sjmisc)
library(ggplot2)
library(gtools)
library(reshape2)
library(UpSetR)
library(tidyr)
library(stringr)


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
add_grade_to_df = function(df, df_murnau_ = df_murnau) {
  df$Grade_va <- df_murnau_$va_ais
  df$Grade_ai <- df_murnau_$ai_ais
  
  df <- df %>% mutate(Grade = ifelse(Grade_va == "ND" | is.na(Grade_va), Grade_ai, Grade_va)) %>%
    dplyr::select(-c("Grade_va", "Grade_ai"))
  return(df)
}


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
#### Extract all the unique lab tests made from the Murnau dataset ####       
########################################################################################################## #

# Use the dictionary to do it
lab_test_df <- filter(dictionary_murnau, grepl("lab", FormName))
lab_test <- lab_test_df$VariableFieldName
lab_test_df <- dplyr::select(df_murnau, any_of(lab_test))

# Remove the columns that only contain NAs
lab_test_df_no_na_col <- Filter(function(x) !all(is.na(x)), lab_test_df)

# Replace all values that are empty strings with NAs
lab_test_df_no_na_col <- mutate_all(lab_test_df_no_na_col, list(~na_if(.,"")))

# Only select the columns that contains tests done after the injury
lab_tests_tot <- lab_test_df_no_na_col[, str_contains("neg", 
                                                      names(lab_test_df_no_na_col), switch = TRUE) == F]
#write.csv(lab_test_df_no_na_col,"~/Data/Lab_tests_tot.csv", row.names = FALSE)

# Find all unique tests 
test_names <- sub("[0-9]+.*", "", names(lab_test_df_no_na_col))
test_names_unique <- unique(test_names)
# Remove the last character: _
test_names_unique <- substr(test_names_unique, 1, nchar(test_names_unique)-1)
cat("Number of unique tests:", length(test_names_unique))


########################################################################################################## #
#### Total number of tests per marker (counting columns) ####
########################################################################################################## #

test_names_df <- as.data.frame(plyr::count(test_names))
test_names_df$is_neg <- str_contains("neg", test_names_df$x, switch = TRUE)

# Plot the numbers
ggplot(test_names_df, aes(x=x, y=freq)) +
  geom_bar(stat="identity", aes(fill=is_neg), width = 0.7) +
  geom_text(aes(label=freq), vjust=-0.3, size=1.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 2)) +
  labs(title = "Number of tests per laboratory marker", x = "Markers", y = "Number of tests") +
  scale_fill_manual(values = c("turquoise", "magenta"),
                    name = "Time of test", 
                    labels = c("After injury", "Before injury"))
ggsave("n_tests.png", path = "~/Plots/Plots_all_LM")

# Plot only the tests before injury
ggplot(subset(test_names_df, is_neg == TRUE), aes(x=x, y=freq)) +
  geom_bar(stat="identity", fill = "magenta", width = 0.7) +
  geom_text(aes(label=freq), vjust=-0.3, size=1.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3.5)) +
  labs(title = "Number of tests per laboratory marker taken before injury", x = "Markers",
       y = "Number of tests")
ggsave("n_tests_before.png", path = "~/Plots/Plots_all_LM")

cat("Number of unique tests done before injury:", sum(test_names_df$is_neg == T))

# Plot only the tests before injury
ggplot(subset(test_names_df, is_neg == FALSE), aes(x=x, y=freq)) +
  geom_bar(stat="identity", fill = "turquoise", width = 0.7) +
  geom_text(aes(label=freq), vjust=-0.3, size=1.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
  labs(title = "Number of tests per laboratory marker taken after injury", x = "Markers",
       y = "Number of tests")
ggsave("n_tests_after.png", path = "~/Plots/Plots_all_LM")

cat("Number of unique tests done after injury:", sum(test_names_df$is_neg == F))


########################################################################################################## #
#### First two weeks analysis: subset the dataset lab_test_df ####
########################################################################################################## #
# Extract the days from the dataset without the time of the day
days <- extract_days(names(lab_tests_tot), T)
# Select the columns whose day is less than 14
inds <- as.numeric(days) <= 14
df_2_weeks <- lab_tests_tot[, inds]
# Order the columns by mixed order (numbers treat as numbers)
df_2_weeks <- df_2_weeks[, mixedorder(names(df_2_weeks))]
#write.csv(df_2_weeks,"~/Data/Lab_tests_two_weeks.csv", row.names = FALSE)

cat("Number of unique tests done in the first two weeks:",
    length(unique(sub("[0-9]+.*", "", names(df_2_weeks)))))

########################################################################################################## #
#### Check the markers that are not tested in the first two weeks ####
########################################################################################################## #
markers_not_2_weeks <- unique(sub("[0-9]+.*", "", names(df_2_weeks)))
markers_not_2_weeks <- substr(markers_not_2_weeks, 1, nchar(markers_not_2_weeks)-1)
markers_not_2_weeks <- setdiff(test_names_unique, markers_not_2_weeks)

# Extract the days in which those markers were tested
N <- length(markers_not_2_weeks)
list_markers_not_2_weeks <- vector("list", N)
i <- 1
while(i <= N) {
  list_markers_not_2_weeks[[i]] <- names(select(lab_test_df_no_na_col, starts_with(markers_not_2_weeks[i])))
  i <- i+1
}
days_markers_not_2_weeks <- unlist(list_markers_not_2_weeks)
days_markers_not_2_weeks <- extract_days(days_markers_not_2_weeks, F)


########################################################################################################## #
#### Count the total number of tests that was done for each marker ####                                 
########################################################################################################## #
# Add a patient columns to make the melting easier
df_2_weeks$Patient <- 1:363
# Add the grade columns
df_2_weeks <- add_grade_to_df(df_2_weeks)
df_2_weeks <- mutate_all(df_2_weeks, as.character)

# Melt the dataset
df_melted <- pivot_longer(df_2_weeks, !c("Grade", "Patient"), names_to = "Day", values_to = "Value")
# Add columns for day and marker using the column variable, then remove the latter
df_melted <- df_melted %>% 
  mutate(DayTime = gsub(".*_", "", Day)) %>% 
  mutate(Marker = sub("_[^_]+$", "", Day)) %>% 
  select(-Day)
# Reorder the columns
df_melted <- df_melted[, c("Patient", "Grade", "Marker", "DayTime", "Value")]
# Reshape the dataframe into a long format, so that the name of the markers are in columns
df_cast_marker <- df_melted %>% select(-c("Grade"))
df_cast_marker <- pivot_wider(df_cast_marker, id_cols = c("Patient", "DayTime"), 
                              names_from = Marker, values_from = Value)
# Reorder the rows in DayTime
df_cast_marker <- df_cast_marker[mixedorder(df_cast_marker$DayTime), ]

# Store the DayTime column
# It will be useful when using upset each day
day_col <- gsub("[^[:digit:] ]", "", df_cast_marker$DayTime)

# Transform the dataframe into a binary matrix
df_upset <- as.data.frame(df_cast_marker)
df_upset <- select(df_upset, -c("Patient", "DayTime"))
df_upset <- df_upset %>% mutate_all(function(x) ifelse(is.na(x), 0, 1))


# Before using upset: plot the number of non-NA values per marker
# Fill 28 routinely assessed markers with another color
n_no_na <- colSums(df_upset)
list_df_names_murnau = c("amylase", "ap", "calcium", "che", "creatinin", "crp", "ery", "gamma_gt", 
                         "gesamt_bilirubin", "gesamteiweiss", "glucose", "harnstoff", "hb", "hbe",
                         "hk", "inr", "kalium", "ldh", "leuco_nl", "lipase","mchc", "mcv", "natrium", 
                         "ptt", "quick", "thrombo", "got", "gpt")
is_28 <- names(df_upset) %in% list_df_names_murnau
n_no_na_df <- data.frame(Marker = names(df_upset), N_tests = n_no_na, Assessed = is_28)

ggplot(data = n_no_na_df, aes(x = reorder(Marker, -N_tests), y = N_tests, fill = Assessed)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
  labs(title = "Number of non-Na values per laboratory marker", x = "Markers", y = "Number of tests") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("turquoise", "magenta"),
                    name = "Routinely assessed", 
                    labels = c("No", "Yes"))
ggsave("n_no_na_ordered.png", path = "~/Plots/Plots_all_LM")

# Isolate markers that were tested less than 10 times, to see whether they contain useful information
less_10 <- rownames(filter(n_no_na_df, N_tests <= 10))
less_10_df <- df_2_weeks
less_10_df <- less_10_df[, sub("_[0-9]+.*", "", names(less_10_df)) %in% less_10] 
less_10_df <- less_10_df %>% mutate(Patient = df_2_weeks$Patient) %>% 
  mutate(Grade = df_2_weeks$Grade)
less_10_df <- less_10_df[rowSums(is.na(less_10_df)) != ncol(less_10_df)-2, ]
less_10_count_grades <- table(less_10_df$Grade, useNA = "always")


########################################################################################################## #
#### Identify which markers are often tested together ####                                 
########################################################################################################## #
# Use the upset function
upset(df_upset, nsets = 107, nintersects = 12, order.by = "freq")
upset(df_upset, nsets = 38, order.by = "freq")
upset(df_upset, nsets = 39, order.by = "freq")

# Remove the columns referring to the 28 routinely tested markers
df_upset_no_28 <- select(df_upset, -all_of(list_df_names_murnau))
upset(df_upset_no_28, nsets = 10, nintersects = 10, order.by = "freq")
upset(df_upset_no_28, nsets = 13, order.by = "freq")
upset(df_upset_no_28, nsets = 48, order.by = "freq")

# Make upset plots for each day
for (i in 0:14) {
  df <- df_upset %>% mutate(Day = day_col) %>% 
    filter(Day == as.character(i)) %>% 
    select(-Day) %>% 
    select(-all_of(list_df_names_murnau))
  
  p <- upset(df, nsets = 13, order.by = "freq")
  print(p)
}


########################################################################################################## #
#### Analysis of middle markers ####                                 
########################################################################################################## #
list_names_middle_markers <- c("a_luc", "ckmb_prozent", "ck_mb", "a_basophile", "a_eosinophile", "a_monocyten",
                               "a_lymphocyten", "atrophile", "nrbc", "gesamt_ck")
list_middle_markers <- c("LUC", "Creatin kinase MB", "Creatin kinase MB %", "Basophile", "Eosinophile",
                         "Monocytes", "Lymphocytes", "Atrophile", "NRBC", "Total CK")
df_middle_markers_melt <- filter(df_melted, Marker %in% list_names_middle_markers)
df_middle_markers_melt <- mutate(df_middle_markers_melt, Day = as.numeric(gsub(".{1}$", "", DayTime)))

# For each patient and for each marker, count how many tests were done per day
df_middle_markers_count <- filter(df_middle_markers_melt, !is.na(Value))
df_middle_markers_count <- as.data.frame(dplyr::count(df_middle_markers_count, Patient, Grade, Marker, Day))

df_middle_markers_count1 <- df_middle_markers_melt %>% 
  filter(!str_detect(DayTime, "b"))  %>% 
  filter(is.na(Value))
df_middle_markers_count1 <- as.data.frame(dplyr::count(df_middle_markers_count1, Patient, Grade, Marker, Day))
df_middle_markers_count1 <- mutate(df_middle_markers_count1, n = rep(0, nrow(df_middle_markers_count1)))

df_middle_markers_count <- rbind(df_middle_markers_count, df_middle_markers_count1)

# Plot counts per day and per patient using geom_bar
# Define the palette that will be used
cols <- c("#c10315", "#ff82c3", "#6d3d8e", "#4c9eff", "#00a34a", "#edb700")

for (marker in list_names_middle_markers) {
  df <- filter(df_middle_markers_count, Marker == marker) %>% 
    mutate(Cases = df_summary_grades[Grade,]$Grade_given) %>% 
    mutate(Weights = 1/Cases) %>%
    select(-Cases)
  
  p1 <- ggplot(data = df, aes(x = n, weights = Weights, fill = Grade)) + 
    geom_bar() + 
    stat_count(aes(label = round(..count.., digits = 2)), 
               geom = "text", 
               position = position_stack(vjust=0.5), 
               size=1.5) +
    facet_wrap(~ Day, ncol = 8) +
    scale_fill_manual(values = cols, na.value = "grey",
                      name = "AIS grade",
                      labels = c("A", "B", "C", "D", "E", "ND", "NA")) +
    labs(title = "Patients divided according to the n° of tests/day", x = "N° of tests", y = "Patients count",
         subtitle = marker)
  
  f1 <- paste0(marker, "_bar.png")
  dir1 <- "~/Plots/Plots_all_LM/Patients_test_day"
  ggsave(filename = f1, plot = p1, path = dir1)
}

for (marker in list_names_middle_markers) {
  df <- filter(df_middle_markers_count, Marker == marker) %>% 
    mutate(Cases = df_summary_grades[Grade,]$Grade_given) %>% 
    mutate(Weights = 1/Cases) %>%
    select(-Cases) %>%
    filter(Grade != "ND") %>% 
    filter(!is.na(Grade))
  
  p1 <- ggplot(data = df, aes(x = n, weights = Weights, fill = Grade)) + 
    geom_bar() + 
    stat_count(aes(label = round(..count.., digits = 2)), 
               geom = "text", 
               position = position_stack(vjust=0.5), 
               size=1.5) +
    facet_wrap(~ Day, ncol = 8) +
    scale_fill_manual(values = cols,
                      name = "AIS grade",
                      labels = c("A", "B", "C", "D", "E")) +
    labs(title = "Patients divided according to the n° of tests/day", x = "N° of tests", y = "Patients count",
         subtitle = marker)
  
  f1 <- paste0(marker, "_bar_nonand.png")
  dir1 <- "~/Plots/Plots_all_LM/Patients_test_day"
  ggsave(filename = f1, plot = p1, path = dir1)
}

