library(dplyr)

setwd("C:/Users/mds376/Desktop/Manu_Peter_metaanalise2024")

new_test <- read.csv("ant_review_spreadsheet_sanger_ngs_revisions.csv") 

# Package needed to manipulate data in long format
library(tidyr)

# Clean the data: Remove extra spaces and standardize to lowercase
new_test$Bacterial_order1 <- trimws(tolower(new_test$Bacterial_order1))
new_test$Bacterial_order2 <- trimws(tolower(new_test$Bacterial_order2))
new_test$Bacterial_order3 <- trimws(tolower(new_test$Bacterial_order3))

# Combine the three bacterial order columns into a single column, keeping the diet
long_data <- new_test %>%
  gather(key = "Bacterial_order_type", value = "Bacterial_order", Bacterial_order1, Bacterial_order2, Bacterial_order3) %>%
  filter(!is.na(Bacterial_order))  # Remove rows with missing bacterial orders only

# Check if all bacteria were preserved correctly
unique(long_data$Bacterial_order)  # This will show all bacteria present after cleaning

# Convert bacterial names to title case
long_data$Bacterial_order <- tools::toTitleCase(long_data$Bacterial_order)

# Create a contingency table between 'Diet' and each 'Bacterial_order'
tab_Diet_Bacterial_order <- table(long_data$Diet, long_data$Bacterial_order)

# Subfamily
# Create a contingency table between 'Subfamily' and each 'Bacterial_order'
tab_Subfamily_Bacterial_order <- table(long_data$Subfamily, long_data$Bacterial_order)

tab_Subfamily_Bacterial_order <- table(long_data$Bacterial_order, long_data$Subfamily)

# Continent
tab_Continent_Bacterial_order <- table(long_data$Continent, long_data$Bacterial_order)

# Continent (again)
tab_Continent_Bacterial_order <- table(long_data$Bacterial_order, long_data$Continent)

# View the contingency tables to ensure bacteria are present
View(tab_Diet_Bacterial_order)

# Continent
View(tab_Continent_Bacterial_order)

# Subfamily
View(tab_Subfamily_Bacterial_order)

View(long_data)

# Perform Chi-square test to check for association
chisq.test(tab_Diet_Bacterial_order, simulate.p.value = TRUE, B = 5000)

# Perform Chi-square test to check for association
chisq.test(tab_Continent_Bacterial_order, simulate.p.value = TRUE, B = 5000)

# Perform Chi-square test to check for association
chisq.test(tab_Subfamily_Bacterial_order, simulate.p.value = TRUE, B = 5000)

# Install the package if necessary
install.packages("metafor")

# Load the package
library(metafor)

# Example preparation (will depend on your actual file)
# Assuming we have: Bacterial_order, Diet, and Study_ID

# Here we count how many times each bacterium appeared by Diet and Study
library(dplyr)
meta_data <- long_data %>%
  group_by(Paper.Title, Diet, Bacterial_order) %>%
  summarise(Count = n()) %>%
  ungroup()

# Transform into proportions or compute log-odds (if it's presence/absence)
# Meta-analysis model (simplified example using effect size and variance)
res <- rma.mv(yi = log(Count), V = 1/Count, mods = ~ Diet,
              random = ~ 1 | Paper.Title/Bacterial_order,
              data = meta_data)
summary(res)

# Here we count how many times each bacterium appeared by Continent and Study
library(dplyr)
meta_data_1 <- long_data %>%
  group_by(Paper.Title, Continent, Bacterial_order) %>%
  summarise(Count = n()) %>%
  ungroup()

# Transform into proportions or compute log-odds (if it's presence/absence)
# Meta-analysis model (simplified example using effect size and variance)
res <- rma.mv(yi = log(Count), V = 1/Count, mods = ~ Continent,
              random = ~ 1 | Paper.Title/Bacterial_order,
              data = meta_data_1)
summary(res)

# Here we count how many times each bacterium appeared by Subfamily and Study
library(dplyr)
meta_data_2 <- long_data %>%
  group_by(Paper.Title, Subfamily, Bacterial_order) %>%
  summarise(Count = n()) %>%
  ungroup()

# Transform into proportions or compute log-odds (if it's presence/absence)
# Meta-analysis model (simplified example using effect size and variance)
res <- rma.mv(yi = log(Count), V = 1/Count, mods = ~ Subfamily,
              random = ~ 1 | Paper.Title/Bacterial_order,
              data = meta_data_2)
summary(res)
