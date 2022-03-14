# Load packages ----
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(visdat)
library(naniar)
library(tidyr)
library(RColorBrewer)
library(UpSetR)
library(gtools)
library(data.table)

#---
# Load the data ----
df_murnau <- read.csv("~/Data/HematologicalBiomark_DATA_2021-01-07_0729.csv",
                      na.strings = c(""))
path <- "~/HematologicalBiomarkersForSCI_DataDictionary_2019-10-01.csv"
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
#write.csv(df_summary_grades,"~/Data/Summary_grades.csv", row.names = T)


########################################################################################################## #
#### Extract all unique times a test was performed ####
########################################################################################################## #

times_in_day <- rep(F, 28)
test_times <- mapply(extract_days, list_murnau_sorted, list_df_names_murnau, times_in_day)
unique_times <- unique(unlist(test_times))

unique_days <- gsub(".{1}$", "", unique_times)
unique_days_df <- as.data.frame(plyr::count(unique_days))


########################################################################################################## #
#### Plot the size of each dataset, in order to verify that is all good ####
########################################################################################################## #

sizes <- sapply(list_murnau_sorted, ncol)
df_plot <- data.frame(lab_test = list_names_murnau, n_cols = sizes)

ggplot(data = df_plot, aes(x=lab_test, y=n_cols)) +
  geom_bar(stat="identity", color="blue", fill="dodgerblue", width = 0.7) +
  geom_text(aes(label=n_cols), vjust=-0.3, size=1.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title = "Dataset size per laboratory marker", x = "Markers", y = "Number of columns (test times)")
ggsave("dataset_sizes.png", path = "~/Plots/Plots_28_LM")


########################################################################################################## #
#### Plot how many times per day a specific marker is tested ####
########################################################################################################## #

# Function that will be applied to the list of sorted datasets
plot_tests_per_day = function(df, marker, name_df, s = T) {
  test_days <- extract_days(df, name_df)
  test_days_df <- as.data.frame(plyr::count(test_days))
  names(test_days_df) <- c("Day", "Freq")
  test_days_df$Day <- as.numeric(test_days_df$Day)
  test_days_df <- test_days_df[order(test_days_df$Day), ]
  
  # Create a dataset with all days from 0 to the last day of the dataset
  # Initialise all test counts per day to zero
  last_day <- tail(test_days_df$Day, n=1)
  days <- seq(0, last_day, by=1)
  n_test <- rep(0, last_day+1)
  days_df <- data.frame(Day = days, N_test = n_test)
  
  # For the days in which there was a test, replace with the true value from test_days_df
  ind_test <- days_df$Day %in% test_days_df$Day
  days_df$N_test[ind_test] <- test_days_df$Freq
  
  # Make a dataframe useful for plotting the different phases
  rects <- data.frame(xstart = c(0, 14, 28, 84, 182),
                      xend = c(14, 28, 84, 182, last_day), 
                      phase = c("Very acute", "Acute I", "Acute II", "Acute III", "Chronic"))
  
  # Make the title
  title <- paste(marker, ": Number of tests per day")
  
  # Plot the data
  p1 <- ggplot() +
    geom_line(data = days_df, 
              aes(x = Day, y = N_test, group = 1),
              size = 0.2) +
    geom_point(data = test_days_df,
               aes(x = Day, y = Freq),
               size = 0.2) +
    geom_rect(data = rects,
              aes(xmin = xstart, xmax = xend,
                  ymin = -Inf, ymax = Inf,
                  fill = phase),
              alpha = 0.5) +
    scale_fill_brewer(palette = 'Pastel1', name = 'Phase') +
    theme(legend.position = c(0.85, 0.8), legend.direction = "vertical") +
    labs(title = title, y = "N° tests") +
    ylim(0, 4)
    
  
  if(s) {
    # Save plot
    f1 <- paste0(name_df, ".png")
    dir1 <- "~/Plots/Plots_28_LM/Test_per_day"
    save_plot(filename = f1, plot = p1, path = dir1)
  } else {
    return(p1)
  }
}

mapply(plot_tests_per_day, list_murnau_sorted, list_names_murnau, list_df_names_murnau)

# Now plot the 28 plots in only two images, to have a quick comparison
## First build the list of plots
myplots <- list()
for (i in 1:28) {
  myplots[[i]] <- local({
    i <- i
    p1 <- plot_tests_per_day(list_murnau_sorted[[i]], list_names_murnau[[i]], list_df_names_murnau[[i]], F)
  })
}

## Second, extract the legend from one of the plots
legend <- get_legend(myplots[[1]] + theme(legend.position="bottom",
                                          legend.title = element_text(size = 30),
                                          legend.text = element_text(size = 20)))

## Third, define a function to remove the legend and apply it to all plots
remove_legend = function(plot) {
  return(plot + theme(legend.position="none"))
}

myplots <- lapply(myplots, remove_legend)

# Finally, make the plot grid
p1 <- plot_grid(plotlist = myplots[1:14], nrow = 5, ncol = 3) + draw_grob(legend, 2/3, 0, 1/3, 1/5)
save_plot(filename = "first_half.png", plot = p1, path = "~/Plots/Plots_28_LM/Test_per_day",
          nrow = 5, ncol = 3)
p2 <- plot_grid(plotlist = myplots[15:28], ncol=3) + draw_grob(legend, 2/3, 0, 1/3, 1/5)
save_plot(filename = "second_half.png", plot = p2, path = "~/Plots/Plots_28_LM/Test_per_day",
          nrow = 5, ncol = 3)


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
#### Plot the size of each dataset, in order to verify that is all good ####
########################################################################################################## #

sizes <- sapply(list_murnau_2_weeks, ncol)
df_plot <- data.frame(lab_test = list_names_murnau, n_cols = sizes)

ggplot(data = df_plot, aes(x=lab_test, y=n_cols)) +
  geom_bar(stat="identity", color="blue", fill="dodgerblue", width = 0.7) +
  geom_text(aes(label=n_cols), vjust=-0.3, size=1.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title = "Dataset sizes considering the first 2 weeks",
       x = "Markers", y = "Number of columns (test times)")
ggsave("dataset_sizes_2_weeks.png", path = "~/Plots/Plots_28_LM/2_weeks_analysis")


########################################################################################################## #
#### Glance at NAs data ####
########################################################################################################## #

plot_glance_at_na = function(df, marker, name_df, s = T) {
  p_dat <- vis_dat(df,palette = "qual", sort_type = F)
  p_miss <- vis_miss(df)
  
  if(s) {
    # Save plot
    f1 <- paste0(name_df, ".png")
    f2 <- paste0(name_df, "_miss.png")
    dir <- "~/Plots/Plots_28_LM/NA_glance"
    save_plot(filename = f1, plot = p_dat, path = dir)
    save_plot(filename = f2, plot = p_miss, path = dir)
  } else {
    return(list(p_dat, p_miss))
  }
}

mapply(plot_glance_at_na, list_murnau_2_weeks, list_names_murnau, list_df_names_murnau)

# Group all plots into two images
## Define the two lists of plots
myplots <- list()
myplots_miss <- list()

for (i in 1:28) {
  plots <- plot_glance_at_na(list_murnau_2_weeks[[i]], list_names_murnau[[i]], list_df_names_murnau[[i]], F)
  myplots[[i]] <- plots[[1]]
  myplots_miss[[i]] <- plots[[2]]
}

## Second, extract the legend from one of the plots
legend <- get_legend(myplots[[1]])

## Third, remove the legend
myplots <- lapply(myplots, remove_legend)

## Finally, make the plot grid for vis_dat plots
path1 <- "~/Plots/Plots_28_LM/NA_glance"
p1 <- plot_grid(plotlist = myplots[1:14], nrow = 3, ncol = 5) + draw_grob(legend, 2/3, 0, 1/3, 1/5)
save_plot(filename = "first_half.png", plot = p1, path = path1,
          nrow = 3, ncol = 5)
p2 <- plot_grid(plotlist = myplots[15:28], nrow = 3, ncol = 5) + draw_grob(legend, 2/3, 0, 1/3, 1/5)
save_plot(filename = "second_half.png", plot = p2, path = path1,
          nrow = 3, ncol = 5)

## Make the plot grid for vis_miss plots
p1 <- plot_grid(plotlist = myplots_miss[1:14], nrow = 3, ncol = 5)
save_plot(filename = "first_half_miss.png", plot = p1, path = path1,
          nrow = 3, ncol = 5)
p2 <- plot_grid(plotlist = myplots_miss[15:28], nrow = 3, ncol = 5)
save_plot(filename = "second_half_miss.png", plot = p2, path = path1,
          nrow = 3, ncol = 5)


########################################################################################################## #
#### Plot no NA values per patient for the first 2 weeks ####
########################################################################################################## #

# This means, plotting how many patients have n no-NA in each of the 14 days
# Build the function to create the dataframes containing the count of test per day per each patient
make_df_no_na_per_day = function(df, name_df) {
  # Extract the days
  unique_days <- as.integer(extract_days(df, name_df, T))
  
  # Count how many times in a day a test was performed
  test_number <- plyr::count(unique_days)$freq
  
  # Build the matrix that counts the number of tests per day per each patient
  l <- length(test_number)
  no_na_mat <- matrix(0, nrow = nrow(df), ncol = l)
  i = 1
  start = 1
  
  while(i <= l) {
    end = start + test_number[i] - 1
    if (start == end) {
      no_na_mat[, i] <- as.integer(!is.na(df[, start]))
    } else {
      no_na_mat[, i] <- rowSums(!is.na(df[, start:end]))
    }
    start = end + 1
    i = i + 1
  }
  
  # Transform matrix into a dataframe
  n_no_na_df <- as.data.frame(no_na_mat)
  # Change the names of the dataframe
  names(n_no_na_df) <- as.character(0:(l-1))
  
  return(n_no_na_df)
}


# Apply the function to the list of dataframes
# tpd: test per day
list_patients_tpd <- mapply(make_df_no_na_per_day, list_murnau_2_weeks, list_df_names_murnau,
                            SIMPLIFY = F)
list_patients_tpd <- lapply(list_patients_tpd, add_grade_to_df)

# Plot patients according to the number of test per day using geom_jitter
plot_patients_test_day_jitter = function(df, marker, name_df, dir1) {
  df_melt <- pivot_longer(df, !Grade, names_to = "Day", values_to = "N_test")
  df_melt <- df_melt[order(df_melt$Day), ]
  df_melt <- df_melt %>% mutate(Day = as.numeric(Day))
  
  # Plot counts per day and per patient using geom_jitter
  cols <- c("#c10315", "#ff82c3", "#6d3d8e",
            "#4c9eff", "#00a34a", "#edb700")
  p1 <- ggplot(data = df_melt, aes(x = Day, y = N_test)) +
    geom_jitter(aes(color = Grade), size=0.5, alpha=0.9, height = 0.2) +
    scale_color_manual(values = cols, na.value = "grey",
                       name = "AIS grade",
                       labels = c("A", "B", "C", "D", "E", "ND", "NA")) +
    labs(title = "Patients divided according to the n° of tests/day",
         subtitle = marker) +
    theme(legend.position = "bottom", legend.box = "horizontal")
  
  f1 <- paste0(name_df, ".png")
  save_plot(filename = f1, plot = p1, path = dir1)
}
dirs <- rep("~/Plots/Plots_28_LM/2_weeks_analysis/Patient_tests_per_day", 28)
mapply(plot_patients_test_day_jitter, list_patients_tpd, list_names_murnau, list_df_names_murnau, dirs)

# Plot patients according to the number of test per day using geom_bar
plot_patients_test_day_bar = function(df, marker, name_df, dir1) {
  df_melt <- pivot_longer(df, !Grade, names_to = "Day", values_to = "N_test")
  df_melt <- df_melt[order(df_melt$Day), ]
  df_melt <- df_melt %>% mutate(Day = as.numeric(Day))
  
  # Plot counts per day and per patient using geom_bar
  cols <- c("#c10315", "#ff82c3", "#6d3d8e",
            "#4c9eff", "#00a34a", "#edb700")
  p1 <- ggplot(data = df_melt, aes(x = N_test, fill = Grade)) + 
    geom_bar() + 
    stat_count(aes(label = ..count..), geom = "text", position = position_stack(vjust=0.5), size=1.5) +
    facet_wrap(~ Day, ncol = 8) +
    scale_fill_manual(values = cols, na.value = "grey",
                      name = "AIS grade",
                      labels = c("A", "B", "C", "D", "E", "ND", "NA")) +
    labs(title = "Patients divided according to the n° of tests/day", x = "N° of tests", y = "Patients count",
         subtitle = marker)
  
  f1 <- paste0(name_df, "_bar.png")
  save_plot(filename = f1, plot = p1, path = dir1)
}
mapply(plot_patients_test_day_bar, list_patients_tpd, list_names_murnau, list_df_names_murnau, dirs)


# Plot patients according to the number of test per day using geom_bar and position = fill to obtain proportions
plot_patients_test_day_prop = function(df, marker, name_df, dir1) {
  df_melt <- pivot_longer(df, !Grade, names_to = "Day", values_to = "N_test")
  df_melt <- df_melt[order(df_melt$Day), ]
  df_melt <- df_melt %>% mutate(Day = as.numeric(Day))
  
  # Plot counts per day and per patient using geom_bar and position = fill
  cols <- c("#c10315", "#ff82c3", "#6d3d8e",
            "#4c9eff", "#00a34a", "#edb700")
  p1 <- ggplot(data = df_melt, aes(x = N_test, fill = Grade)) + 
    geom_bar(position = "fill") +
    facet_wrap(~ Day, ncol = 8) +
    scale_fill_manual(values = cols, na.value = "grey",
                      name = "AIS grade",
                      labels = c("A", "B", "C", "D", "E", "ND", "NA")) +
    labs(title = "Patients proportions according to the n° of tests/day", 
         x = "N° of tests", 
         y = "Patients proportion",
         subtitle = marker)
  
  f1 <- paste0(name_df, "_prop.png")
  save_plot(filename = f1, plot = p1, path = dir1)
}
mapply(plot_patients_test_day_prop, list_patients_tpd, list_names_murnau, list_df_names_murnau, dirs)


########################################################################################################## #
#### Re-do the same plots but divide by the total number of patients having the same grade ####
########################################################################################################## #
# Function that divides the number of tests per day per patient by the total number of patients having 
# the same grade as the patient considered. Then plots the results
plot_patients_test_day_bar_freq = function(df, marker, name_df, dir1) {
  df <- df %>% mutate(Cases = df_summary_grades[Grade,]$Grade_given) %>% 
    mutate(Weights = 1/Cases) %>%
    select(-Cases)
  df_melt <- pivot_longer(df, !c("Grade", "Weights"), names_to = "Day", values_to = "N_test")
  df_melt <- df_melt[order(df_melt$Day), ]
  df_melt <- df_melt %>% mutate(Day = as.numeric(Day))
  
  # Plot counts per day and per patient using geom_bar
  cols <- c("#c10315", "#ff82c3", "#6d3d8e",
            "#4c9eff", "#00a34a", "#edb700")
  p1 <- ggplot(data = df_melt, aes(x = N_test, weights = Weights, fill = Grade)) + 
    geom_bar() + 
    stat_count(aes(label = round(..count.., digits = 2)), geom = "text", position = position_stack(vjust=0.5),
               size=1.5) +
    facet_wrap(~ Day, ncol = 8) +
    scale_fill_manual(values = cols, na.value = "grey",
                      name = "AIS grade",
                      labels = c("A", "B", "C", "D", "E", "ND", "NA")) +
    labs(title = "Patients divided according to the n° of tests/day",
         x = "N° of tests", y = "Patients frequency",
         subtitle = marker)
  
  f1 <- paste0(name_df, "_bar.png")
  save_plot(filename = f1, plot = p1, path = dir1)
}
dirs <- rep("~/Plots/Plots_28_LM/2_weeks_analysis/Patient_tests_per_day_freq", 28)
mapply(plot_patients_test_day_bar_freq, list_patients_tpd, list_names_murnau, list_df_names_murnau, dirs)

# Same plot as above, but with bars filling the whole space
plot_patients_test_day_prop_freq = function(df, marker, name_df, dir1) {
  df <- df %>% mutate(Cases = df_summary_grades[Grade,]$Grade_given) %>% 
    mutate(Weights = 1/Cases) %>%
    select(-Cases)
  df_melt <- pivot_longer(df, !c("Grade", "Weights"), names_to = "Day", values_to = "N_test")
  df_melt <- df_melt[order(df_melt$Day), ]
  df_melt <- df_melt %>% mutate(Day = as.numeric(Day))
  
  # Plot counts per day and per patient using geom_bar and position = fill
  cols <- c("#c10315", "#ff82c3", "#6d3d8e",
            "#4c9eff", "#00a34a", "#edb700")
  p1 <- ggplot(data = df_melt, aes(x = N_test, weights = Weights, fill = Grade)) + 
    geom_bar(position = "fill") +
    facet_wrap(~ Day, ncol = 8) +
    scale_fill_manual(values = cols, na.value = "grey",
                      name = "AIS grade",
                      labels = c("A", "B", "C", "D", "E", "ND", "NA")) +
    labs(title = "Patients divided according to the n° of tests/day",
         x = "N° of tests", y = "Patients frequency",
         subtitle = marker)
  
  f1 <- paste0(name_df, "_prop.png")
  save_plot(filename = f1, plot = p1, path = dir1)
}
mapply(plot_patients_test_day_prop_freq, list_patients_tpd, list_names_murnau, list_df_names_murnau, dirs)

########################################################################################################## #
#### Summarise plot using a heatmap ####
########################################################################################################## #
# For each dataset, contruct a 7 rows dataset that records, for each day, the frequency of each grade, weighted
# by the total number of patients having that grade

# Helper functions needed
sum_tested = function(col) {
  return(sum(col != 0))
}

divide_by_tot = function(col) {
  return(col/df_summary_grades$Grade_given)
}

# Main function
make_grades_frequencies_per_day = function(df, marker) {
  tested <- df %>% group_by(Grade) %>% 
    summarise(across(c(0:15), ~sum_tested(.))) %>% 
    mutate(across(c(2:16), ~divide_by_tot(.))) %>% 
    mutate(Marker = rep(marker, 7))
  return(tested)
}

# Apply function: fdp is frequency per day
list_grades_fpd <- mapply(make_grades_frequencies_per_day, list_patients_tpd, list_names_murnau, SIMPLIFY = F)

# Fuse the 28 datasets
df_grades_fpd <- bind_rows(list_grades_fpd)

# Get the datasets for each grade and plot a heatmap
grades <- df_summary_grades %>% filter(AIS_grade %in% c("A", "B", "C", "D")) %>% rownames

plot_heatmap = function(df_grades_fpd, grade) {
  grade <- grades[i]
  if (grade != "NA") {
    df_grade <- df_grades_fpd %>% 
      filter(Grade == grade) %>% 
      pivot_longer(cols = 2:16, names_to = "Day", values_to = "Freq") %>% 
      mutate(Freq = cut(Freq, breaks = c(0, 0.47, 0.53, 1), right = FALSE)) %>% 
      mutate(Day = as.numeric(Day))
  } else {
    df_grade <- df_grades_fpd %>% 
      filter(is.na(Grade)) %>% 
      pivot_longer(cols = 2:16, names_to = "Day", values_to = "Freq") %>% 
      mutate(Freq = cut(Freq, breaks = c(0, 0.47, 0.53, 1), right = FALSE)) %>% 
      mutate(Day = as.numeric(Day))
  }
  
  p1 <- ggplot(df_grade, aes(x = Day, y = Marker)) + 
    geom_tile(aes(fill = Freq), colour = "white") +
    scale_fill_manual(values = c("yellow", "orange", "red")) +
    labs(title = paste("Test frequencies for grade", grade))
  
  return(p1)
}

myplots <- list()
for (i in 1:4) {
  myplots[[i]] <- local({
    i <- i
    p1 <- plot_heatmap(df_grades_fpd, grades[i])
  })
}

## Make the plot grid for the 4 heatmaps (for grades A, B, C and D)
path1 <- "~/Plots/Plots_28_LM/2_weeks_analysis/Patient_tests_per_day_freq"
p1 <- plot_grid(plotlist = myplots, nrow = 2, ncol = 2)
save_plot(filename = "heatmaps.png", plot = p1, path = path1,
          nrow = 2, ncol = 2)


########################################################################################################## #
#### Plot the same information: each line is a patient, is facet a different grade ####
########################################################################################################## #
plot_patients_test_day_lines = function(df, marker) {
  df <- df %>% mutate(Patient = 1:363) %>%
    filter(Grade %in% c("A", "B", "C", "D")) %>% 
    pivot_longer(cols = 1:15, names_to = "Day", values_to = "N_tests")
  
  p <- ggplot(df, aes(x = as.numeric(Day), y = as.integer(N_tests), group = Patient)) +
    geom_line(size = 0.5, alpha = 0.1, position = position_dodge(width=0.1)) +
    theme_bw() +
    theme(panel.grid=element_blank()) +
    facet_wrap(~ Grade, ncol = 2) +
    labs(title = paste0(marker, ": N° of tests per day and per patient"),
         x = "Day", y = "N° of tests",
         subtitle = marker)
  
  ggsave(filename = paste0(marker, "_lines.png"),
            path = "~/Plots/Plots_28_LM/2_weeks_analysis/Patient_tests_per_day")
}

mapply(plot_patients_test_day_lines, list_patients_tpd, list_names_murnau)


########################################################################################################## #
#### Melt all 28 dataset ####
########################################################################################################## #

# Construct a function that melt all values of a dataframe
# So that we obtain a two-columns dataframe, with the first column recording the day and the patient
# (i.e. 0a_1) and the first column the test value
melt_df_day_patient_value = function(df, name_df) {
  # Change the column names so that we only have the day
  days <- extract_days(df, name_df, F)
  names(df) <- days
  
  # Add columns "Grade" and "Patient" to the dataframe
  df$Patient <- 1:363
  df <- add_grade_to_df(df)
  
  # Add columns Cases to the dataframe
  df <- df %>% mutate(Cases = df_summary_grades[Grade,]$Grade_given)
  
  # Melt the dataframe: keep columns grade, patient and cases
  df_melt <- data.table::melt(setDT(df), id.vars = c("Patient", "Grade", "Cases"), variable.factor = FALSE)
  
  # Change the name of the columns, so that all dataframes will have four columns, the first three named
  # "Patient", "Grade", "Cases" and "DayTime" respectively, while the fourth will have the name of the marker
  names(df_melt) <- c("Patient", "Grade", "Cases", "DayTime", name_df)
  
  return(df_melt)
}

list_df_melted <- mapply(melt_df_day_patient_value, list_murnau_2_weeks, list_df_names_murnau,
                         SIMPLIFY = F)

# Merge all the melted datasets together
joined_melted_dfs <- Reduce(function(x, y) merge(x, y, by = c("Patient", "DayTime", "Grade", "Cases"),
                                                 all.x = T, all.y = T), list_df_melted)
# Sort the rows in ascending days
joined_melted_dfs <- joined_melted_dfs[mixedorder(joined_melted_dfs$DayTime), ]
joined_melted_dfs <- as.data.frame(joined_melted_dfs)


########################################################################################################## #
#### Markers often tested together ####
########################################################################################################## #

# Prepare the dataset so that it can be used with upset function of package UpSetR
# Remove the Day_Patient column from the joined dataframe
joined_melted_dfs_upset <- joined_melted_dfs %>% select(-c("Patient", "DayTime", "Grade", "Cases"))

# Transform the melted dataframe in a binary dataframe
# In this way, we can use the upset function: if, for example, P1 has a 1 under the amylase column
# at day 0a, the upset function sees it as belonging to the "Tested amylase" subset
joined_melted_dfs_upset <- joined_melted_dfs_upset %>% mutate_all(function(x) ifelse(is.na(x), 0, 1))
head(joined_melted_dfs_upset)

# Plot of the total number of tests per each laborytory marker
n_test_tot <- colSums(joined_melted_dfs_upset)
df_plot <- data.frame(Marker = names(joined_melted_dfs_upset), N_test = n_test_tot)
ggplot(data = df_plot, aes(x = Marker, y = N_test)) +
  geom_bar(stat="identity", color="blue", fill="dodgerblue", width = 0.7) +
  geom_text(aes(label = N_test), vjust=-0.3, size=1.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title = "Total n° of non NA values per laboratory marker", x = "Markers", y = "Non NA values")
ggsave("non_na_number.png", path = "~/Plots/Plots_28_LM/2_weeks_analysis")

# Use upset
upset(joined_melted_dfs_upset, nsets = 28, nintersects = 10, order.by = "freq")
# I saved it inPlots_28_LM by making a screenshot


########################################################################################################## #
#### Analysis of data with all 28 markers tested ####
########################################################################################################## #

joined_melted_rs <- joined_melted_dfs
joined_melted_rs$Rowsum <- rowSums(joined_melted_dfs_upset)

# Subset with all the patients that have been tested for all 28 markers at least once
joined_melted_allm <- joined_melted_rs %>% filter(Rowsum == 28) %>% 
  select(-names(joined_melted_rs[5:33]))

# Add a column Weights
joined_melted_allm$Weights <- 1/joined_melted_allm$Cases

# For each patient, count how many times he was tested for all 28 markers
joined_melted_allm_patient <- ddply(joined_melted_allm, .(Patient, Grade, Weights), nrow)

## Make an histogram showing how many times and how many patients 28 markers
### Use identity
ggplot(joined_melted_allm_patient, aes(x = V1, fill = Grade)) + 
  geom_histogram(binwidth = 1, color = "white") +
  stat_count(aes(label = round(..count.., digits = 2)), geom = "text", position = position_stack(vjust=0.5),
             size=1.5) +
  labs(title = "Number of patients that were tested for all \n 28 markers n times", x = "Number of times", 
       y = "Number of patients")
ggsave("hist_allm_patients.png", path = "~/Plots_28_LM/2_weeks_analysis")

### Use identity and frequency values
#### Plot
ggplot(joined_melted_allm_patient, aes(x = V1, weights = Weights, fill = Grade)) + 
  geom_histogram(binwidth = 1, color = "white") +
  stat_count(aes(label = round(..count.., digits = 2)), geom = "text", position = position_stack(vjust=0.5),
             size=1.5) +
  labs(title = "Number of patients that were tested for all 28 markers n times",
       x = "Number of times", y = "Number of patients",
       subtitle = "Count of each grade divided by the total number of patients \n having that grade")
ggsave("hist_allm_patients_freq.png", path = "~/Plots/Plots_28_LM/2_weeks_analysis")

### Use proportions
ggplot(joined_melted_allm_patient, aes(x = V1, weights = Weights, fill = Grade)) + 
  geom_bar(position = "fill") +
  labs(title = "Proportion of patients that were tested \n for all 28 markers n times", x = "Number of times", 
       y = "Patient proportions")
ggsave("hist_allm_patients_prop.png", path = "~/Plots_28_LM/2_weeks_analysis")

### Use proportions and frequency values
ggplot(joined_melted_allm_patient, aes(x = V1, fill = Grade)) + 
  geom_bar(position = "fill") +
  labs(title = "Proportion of patients that were tested \n for all 28 markers n times",
       x = "Number of times", y = "Number of patients",
       subtitle = "Count of each grade divided by the total number of patients \n having that grade")
ggsave("hist_allm_patients_prop_freq.png", path = "~/Plots/Plots_28_LM/2_weeks_analysis")

## For each day, count how many times there was a test for all 28 markers
joined_melted_allm_day <- ddply(joined_melted_allm, .(DayTime, Grade, Weights), nrow)
joined_melted_allm_day <- joined_melted_allm_day[mixedorder(joined_melted_allm_day$DayTime), ]

## Make an histogram showing how many patients were tested with all 28 markers for each day in the first 2 weeks
joined_melted_allm_day$DayTime <- factor(joined_melted_allm_day$DayTime,
                                         levels = unique(joined_melted_allm_day$DayTime))

### Add a column that divides each V1 by the total number of patients having the same grade
joined_melted_allm_day <- joined_melted_allm_day %>% mutate(Freq = V1/df_summary_grades[Grade, ]$Grade_given)

### Use identity
ggplot(data = joined_melted_allm_day, aes(x = DayTime, y = V1, fill = Grade)) +
  geom_bar(stat="identity", width = 0.7) +
  geom_text(aes(label = V1), position = position_stack(vjust=0.5), size = 2.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title = "Number of patients that were tested with all 28 markers \n in the first two weeks",
       x = "Days", y = "Number of patients")
ggsave("hist_allm_day.png", path = "~/Plots_28_LM/2_weeks_analysis")

### Use identity and frequencies
ggplot(data = joined_melted_allm_day, aes(x = DayTime, y = Freq, weights = Weights, fill = Grade)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(Freq, digits = 2)), position = position_stack(vjust=0.5), size = 2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title = "Number of patients that were tested with all 28 markers",
       x = "Days", y = "Number of patients",
       subtitle = "Count of each grade divided by the total number of patients \n having that grade")
ggsave("hist_allm_day_freq.png", path = "~/Plots/Plots_28_LM/2_weeks_analysis")

## Use position = fill
ggplot(data = joined_melted_allm_day, aes(x = DayTime, y = V1, fill = Grade)) +
  geom_col(position = "fill", width = 0.7) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title = "Proportion of patients that were tested with all 28 markers",
       x = "Days", y = "Patient proportions")
ggsave("hist_allm_day_prop.png", path = "~/Plots/Plots_28_LM/2_weeks_analysis")

## Use position = fill and frequencies
ggplot(data = joined_melted_allm_day, aes(x = DayTime, y = Freq, fill = Grade)) +
  geom_col(position = "fill", width = 0.7) +
  geom_text(aes(label = round(Freq, digits = 2)), size = 2, position = "fill", hjust = 0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title = "Proportion of patients that were tested with all 28 markers",
       x = "Days", y = "Patient proportions",
       subtitle = "Count of each grade divided by the total number of patients \n having that grade")
ggsave("hist_allm_day_prop_freq.png", path = "~/Plots/Plots_28_LM/2_weeks_analysis")

