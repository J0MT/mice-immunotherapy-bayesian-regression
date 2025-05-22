# Load necessary libraries
library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)

# Load and prepare dataset1
dataset1 <- read_excel("C:/Users/jomkr/OneDrive/Documents/Imperial_Y3/UROP/Data/CBD-IL-12 myeloid.xlsx", sheet = "Myeloid", range = "A15:Y81")
names(dataset1)[1] <- "mice" 
puredata <- dataset1 %>%
  separate(col = mice, into = c('analysis', 'number'), sep = 2) %>%
  mutate(number = as.integer(number)) %>%
  mutate(days = case_when(
    number >= 1 & number <= 4 ~ 0,
    number >= 5 & number <= 9 ~ 1,
    number >= 10 & number <= 14 ~ 2,
    number >= 15 & number <= 19 ~ 4,
    number >= 20 & number <= 24 ~ 7,
    number >= 25 & number <= 28 ~ 10,
    number >= 29 & number <= 34 ~ 13,
    number >= 35 ~ 14
  ))

puredata_long <- puredata %>%
  pivot_longer(cols = c(3:26), names_to = 'metric', values_to = 'value')

# Filter for TM analysis
TM_data <- puredata_long %>%
  filter(analysis == "TM") %>%
  mutate(metric = factor(metric, levels = c('CD45+', 'CD11b+CD45+', 'M-MDSC', 'G-MDSC', 'Macrophages', 'M2 macrophages', 'M1 macrophages',
                                            'B cells', 'Endothelial cells', 'Dendritic cells', '% M-MDSC/CD45', '% G-MDSC/CD45', '%macrophage/CD45',
                                            '%M2 macrophage/CD45', '%M1 macrophage/CD45', '% B cells/CD45', '% Endothelial cells/CD45', '%Dendritic cells/CD45',
                                            'tumor mg', '# macrophages/tumor', '#M2 macrophages/tumor mg', '# M1 macrophage/tumor mg', '#M-MDSC/tumor mg',
                                            '# B cells/tumor')))

# Load and prepare dataset2
dataset2 <- read_excel("C:/Users/jomkr/OneDrive/Documents/Imperial_Y3/UROP/Data/CBD-IL-12 single injection B16F10.xlsx", 
                       range = "A1:AN10", col_types = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", "numeric"))
names(dataset2)[1] <- "days"

# Ensure 'number' column is included for dataset2
data_long <- dataset2 %>%
  pivot_longer(cols = c(2:40), names_to = "measurement", values_to = "value") %>%
  mutate(treatment = case_when(
    measurement %in% c("PBS", "...3", "...4", "5", "...6", "...7") ~ "PBS",
    measurement %in% c("1 μg CBD-IL-12", "...11", "...12", "...13", "...14", "...15") ~ "1ug CBD-IL-12",
    measurement %in% c("10 μg CBD-IL-12", "...19", "...20", "...21", "...22", "...23", "...24") ~ "10ug CBD-IL-12",
    measurement %in% c("25 μg CBD-IL-12", "...27", "...28", "...29", "...30", "...31", "...32") ~ "25ug CBD-IL-12",
    measurement %in% c("8 μg CBD-IL-2", "...35", "...36", "...37", "...38", "...39", "...40") ~ "8ug CBD-IL-2",
    TRUE ~ NA_character_)) %>%
  na.omit()

# Add 'number' column to data_long and combine with TM data
data_long <- data_long %>%
  group_by(measurement) %>%
  mutate(number = row_number()) %>%
  ungroup()

# Combine TM data and dataset2 for plotting
combined_data <- TM_data %>%
  mutate(source = "TM Analysis") %>%
  rename(treatment = analysis) %>%
  bind_rows(data_long %>%
              mutate(metric = treatment, source = "Tumor Volume")) %>%
  mutate(treatment = factor(treatment),
         metric = factor(metric, levels = c(levels(TM_data$metric), "PBS", "1ug CBD-IL-12", "10ug CBD-IL-12", "25ug CBD-IL-12", "8ug CBD-IL-2")))

# Plot combined data with sample numbers represented as different colors
ggplot(combined_data, aes(x = days, y = value, color = factor(number))) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ metric, nrow = 5, scales = "free_y") +
  labs(title = "Combined Metrics and Tumor Volume over Time",
       x = "Days",
       y = "Value",
       color = "Sample Number") +
  theme_bw() +
  theme(legend.position="right")

