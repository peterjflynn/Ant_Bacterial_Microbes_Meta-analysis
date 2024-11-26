# Set working directory
setwd("/Users/peterflynn/Documents/Ant_Bacterial_Review_Myrmeco/All_Files/")
setwd("~/Documents/Ant_Bacterial_Review_Myrmeco/")
#### richness estimates from Kass 2022####
library(terra)
library(sf)

# Step 1: Load the raster using the correct file path
rf_richness <- rast("main_analysis_data/overlays/for_figs3-6/ants.rf.plus.irs.hs.tif")

# Step 2: Reproject the raster to WGS84 (latitude/longitude) if necessary
rf_richness_wgs84 <- project(rf_richness, "EPSG:4326")

# Step 3: Define continent boundaries (in WGS84, EPSG:4326)
north_america <- ext(-168, -52, 15, 72)   # Long, Lat range for North America
central_america <- ext(-92, -77, 7, 23)   # Long, Lat range for Central America
south_america <- ext(-81, -34, -56, 12)   # South America
europe <- ext(-10, 40, 35, 72)  # Adjusted eastern boundary to 40°E to restrict Eastern Europe
africa <- ext(-17, 54, -35, 37)  # Expand the eastern boundary to 54°E to include Madagascar
asia <- ext(60, 180, -10, 81)  # Expanding to the south to include islands like the Philippines and Indonesia
australia <- ext(113, 153, -44, -10)      # Australia & Oceania

# Step 4: Crop the raster for each continent using the defined boundaries
north_america_richness <- crop(rf_richness_wgs84, north_america)
central_america_richness <- crop(rf_richness_wgs84, central_america)
south_america_richness <- crop(rf_richness_wgs84, south_america)
europe_richness <- crop(rf_richness_wgs84, europe)
africa_richness <- crop(rf_richness_wgs84, africa)
asia_richness <- crop(rf_richness_wgs84, asia)
australia_richness <- crop(rf_richness_wgs84, australia)

# Step 5: Summarize richness for each continent
continents <- list(north_america_richness, central_america_richness, south_america_richness,
                   europe_richness, africa_richness, asia_richness, australia_richness)

continent_names <- c("North America", "Central America", "South America", "Europe", "Africa", "Asia", "Australia & Oceania")

# Inspect raster values for each cropped continent to check if they contain any data
for (i in 1:length(continents)) {
  print(paste("Checking data for", continent_names[i]))
  print(global(continents[[i]], fun = 'sum', na.rm = TRUE))  # This will help see if there is any valid data
}

# Step 6: Print the summaries for each continent
for (i in 1:length(continents)) {
  # Check if the raster has valid values before summarizing
  count_valid <- global(continents[[i]], fun = 'sum', na.rm = TRUE)[[1]]
  if (!is.null(count_valid) && count_valid > 0) {
    summary_value <- global(continents[[i]], fun = 'mean', na.rm = TRUE)[[1]]  # Extract the numeric value
    cat(continent_names[i], ": ", summary_value, "\n")
  } else {
    cat(continent_names[i], ": No valid data\n")
  }
}

#total bacterial count per continent
#Africa: 12 Asia: 27 Central America: 98 Europe: 18 North America: 105 Oceania: 9 South America: 193
## mean = Africa: 12.48 #Asia: 28.24 #Central America: 11.89 #Europe: 104.94 #North America: 355.09 #Oceania: 6.00 #South America: 99.00

##### bar chart map just NGS#####
# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(rnaturalearth)
library(sf)
library(cowplot)

# Load the original data
ant_review_data <- read_csv("ant_review_spreadsheet_nosanger_Chansonadded.csv")

# Clean the data by converting "NA" strings to actual NA values
ant_review_clean <- ant_review_data %>%
  mutate(across(where(is.character), ~na_if(., "NA")))

# Transform the data by unpivoting bacterial orders
ant_review_long <- ant_review_clean %>%
  select(Genus, Continent, Bacterial_order1, Bacterial_order2, Bacterial_order3) %>%
  pivot_longer(cols = c(Bacterial_order1, Bacterial_order2, Bacterial_order3),
               names_to = "OrderType", values_to = "Bacterial_Order") %>%
  filter(!is.na(Bacterial_Order))

# Aggregate data for genus-level bacterial counts
genus_bacterial_counts <- ant_review_long %>%
  group_by(Genus, Bacterial_Order) %>%
  summarise(Count = n(), .groups = "drop")

# Aggregate data for continent-level bacterial counts
continent_bacterial_counts <- ant_review_long %>%
  group_by(Continent, Bacterial_Order) %>%
  summarise(Count = n(), .groups = "drop")

# Save the cleaned results to CSV files
write_csv(genus_bacterial_counts, "Cleaned_Full_Bacterial_Order_Counts_nosanger.csv")
write_csv(continent_bacterial_counts, "Cleaned_Microbial_Summary_by_Continent_nosanger.csv")

# Define the bacterial orders color palette
color_palette <- c(
  "Rickettsiales" = "#FF0000", "Bacillales" = "#222222", "Enterobacterales" = "#f3c300",
  "Entomoplasmatales" = "#875692", "Hyphomicrobiales" = "#f38400", "Kitasatosporales" = "#a1caf1",
  "Pseudomonadales" = "#be0032", "Pseudonocardiales" = "#c2b280", "Unclassified Firmicutes" = "#848482",
  "Actinomycetales" = "#008856", "Burkholderiales" = "#e68fac", "Lactobacillales" = "#0067a5",
  "Moraxellales" = "#f99379", "Sphingobacteriales" = "#604e97", "Xanthomonadales" = "#f6a600",
  "Micrococcales" = "#b3446c", "Mycobacteriales" = "#dcd300", "Rhodospirillales" = "#882d17",
  "Terriglobales" = "#8db600", "Campylobacterales" = "#654522", "Flavobacteriales" = "#e25822",
  "Opitutales" = "#2b3d26", "Rhizobiales" = "#8DD3C7", "Pasteurellales" = "#DCF1B9",
  "Acetobacterales" = "#E4E3C2", "Bacteroidales" = "#C4B4CF", "Nitrospinales" = "#EE8B86",
  "Caulobacteriales" = "#BD98A2", "Sphingomonadales" = "#99B1BC", "Erysipelotrichales" = "#F0B36D",
  "Solirubrobacterales" = "#D0CD66", "Lysobacterales" = "#C8D88E", "Propionibacteriales" = "#FCCDE5"
)

# Add coordinates for continents (adjusted positions)
continent_coords <- data.frame(
  Continent = c("Africa", "Asia", "Europe", "North America", "South America", "Oceania", "Central America"),
  Lon = c(10, 80, 0, -130, -65, 110, -100),
  Lat = c(0, 30, 50, 40, -10, -25, 15)
)

# Merge the coordinates
continent_bacterial_counts <- continent_bacterial_counts %>%
  left_join(continent_coords, by = "Continent")

# Create a function for bar charts without count labels and with custom y-axis colors
create_bar_chart_no_legend <- function(data_row, y_axis_color = "black") {
  data_row <- data_row %>%
    mutate(Bacterial_Order = fct_reorder(Bacterial_Order, Count, .desc = TRUE))  # Order by abundance
  
  ggplot(data_row, aes(x = Bacterial_Order, y = Count, fill = Bacterial_Order)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = color_palette) +
    theme_minimal() +
    labs(x = NULL, y = "Count") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),  # Remove x-axis text
          axis.title.y = element_text(size = 10, color = y_axis_color, face = "bold"),
          axis.text.y = element_text(size = 8, color = y_axis_color, face = "bold"),
          panel.grid = element_blank())                                           # Remove grid lines
}

# Create bar charts for each continent with custom y-axis colors
bar_charts_updated <- list()
for (i in unique(continent_bacterial_counts$Continent)) {
  continent_data <- continent_bacterial_counts %>% filter(Continent == i)
  
  if (i == "Africa" || i == "Asia") {
    bar_chart <- create_bar_chart_no_legend(continent_data, y_axis_color = "white")
  } else if (i == "Oceania" || i == "Central America") {
    bar_chart <- create_bar_chart_no_legend(continent_data, y_axis_color = "darkgrey")
  } else {
    bar_chart <- create_bar_chart_no_legend(continent_data, y_axis_color = "black")
  }
  
  bar_charts_updated[[i]] <- bar_chart
}

# Create the bacterial diversity map
bacterial_diversity <- data.frame(
  Continent = c("Africa", "Asia", "Central America", "Europe", "North America", "Oceania", "South America"),
  Category = c("Undersampled", "Undersampled", "Undersampled", "Sampled", "Oversampled", "Undersampled", "Sampled")
)

# Define the colors for each category
category_colors <- c(
  "Oversampled" = "#ffcccc",    # Light red (North America)
  "Sampled" = "#ff9999",        # Medium red (Europe, South America)
  "Undersampled" = "#660000"    # Darkest red (Oceania, Africa, Asia, Central America)
)

# Load the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Assign continent categories to the world map
world$Continent <- ifelse(world$region_un == "Africa", "Africa",
                          ifelse(world$region_un == "Asia", "Asia",
                                 ifelse(world$region_un == "Europe", "Europe",
                                        ifelse(world$region_un == "Oceania", "Oceania",
                                               ifelse(world$region_un == "Americas" & world$subregion == "South America", "South America",
                                                      ifelse(world$region_un == "Americas" & world$subregion == "Central America", "Central America",
                                                             ifelse(world$region_un == "Americas", "North America", NA)))))))

# Merge the bacterial diversity categories with the world map data
world <- merge(world, bacterial_diversity, by.x = "Continent", by.y = "Continent", all.x = TRUE)

# Create the map with custom colors for each sampling category
world_map <- ggplot(data = world) +
  geom_sf(aes(fill = Category), color = "black") +
  scale_fill_manual(values = category_colors) +
  theme_void() +
  labs(fill = "Sampling Level")

# Combine the bacterial diversity map with bar charts
combined_bar_chart_map <- ggdraw() +
  draw_plot(world_map, 0, 0, 1, 1) +
  draw_plot(bar_charts_updated[['Africa']], 0.45, 0.50, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Asia']], 0.58, 0.56, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Europe']], 0.44, 0.65, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['North America']], 0.18, 0.64, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['South America']], 0.27, 0.4, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Oceania']], 0.71, 0.38, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Central America']], 0.18, 0.53, 0.10, 0.10)

# Save the combined bar chart map
ggsave("final_combined_bar_chart_map_justNGS.pdf", plot = combined_bar_chart_map, width = 15, height = 10)

##### Pie Chart map just NGS #####
# Load libraries
library(ggplot2)
library(rnaturalearth)
library(sf)
library(cowplot)

# Step 2: Load the data
ant_review_data <- read_csv("ant_review_spreadsheet_nosanger_Chansonadded.csv")

# Step 3: Clean the data
ant_review_data_clean <- ant_review_data %>%
  mutate(across(where(is.character), ~na_if(., "NA"))) 

# Step 4: Transform the data by unpivoting bacterial orders
ant_review_long <- ant_review_data_clean %>%
  select(Genus, Bacterial_order1, Bacterial_order2, Bacterial_order3, Continent) %>%
  pivot_longer(cols = c(Bacterial_order1, Bacterial_order2, Bacterial_order3),
               names_to = "OrderType", values_to = "Bacterial_Order") %>%
  filter(!is.na(Bacterial_Order) & !is.na(Continent))

# Step 5: Aggregate data by continent and bacterial order
continent_bacterial_counts <- ant_review_long %>%
  group_by(Continent, Bacterial_Order) %>%
  summarise(Count = n(), .groups = "drop")

# Step 6: Define the bacterial orders color palette
color_palette <- c("Rickettsiales" = "#FF0000", "Bacillales" = "#222222", "Enterobacterales" = "#f3c300", 
                   "Entomoplasmatales" = "#875692", "Hyphomicrobiales" = "#f38400", "Kitasatosporales" = "#a1caf1",
                   "Pseudomonadales" = "#be0032", "Pseudonocardiales" = "#c2b280", "Unclassified Firmicutes" = "#848482",
                   "Actinomycetales" = "#008856", "Burkholderiales" = "#e68fac", "Lactobacillales" = "#0067a5",
                   "Moraxellales" = "#f99379", "Sphingobacteriales" = "#604e97", "Xanthomonadales" = "#f6a600",
                   "Micrococcales" = "#b3446c", "Mycobacteriales" = "#dcd300", "Rhodospirillales" = "#882d17",
                   "Terriglobales" = "#8db600", "Campylobacterales" = "#654522", "Flavobacteriales" = "#e25822",
                   "Opitutales" = "#2b3d26", "Rhizobiales" = "#8DD3C7", "Pasteurellales" = "#DCF1B9",
                   "Acetobacterales" = "#E4E3C2", "Bacteroidales" = "#C4B4CF", "Nitrospinales" = "#EE8B86",
                   "Caulobacteriales" = "#BD98A2", "Sphingomonadales" = "#99B1BC", "Erysipelotrichales" = "#F0B36D",
                   "Solirubrobacterales" = "#D0CD66", "Lysobacterales" = "#C8D88E", "Propionibacteriales" = "#FCCDE5")

# Step 7: Add coordinates for continents (adjusted positions)
continent_coords <- data.frame(
  Continent = c("Africa", "Asia", "Europe", "North America", "South America", "Oceania", "Central America"),
  Lon = c(10, 80, 0, -130, -65, 110, -100),  # North America, Australia, Central America moved left
  Lat = c(0, 30, 50, 40, -10, -25, 15)
)

# Merge the coordinates
continent_bacterial_counts <- continent_bacterial_counts %>%
  left_join(continent_coords, by = "Continent")

# Step 8: Create a function for pie charts
create_pie_chart_no_legend <- function(data_row) {
  ggplot(data_row, aes(x = "", y = Count, fill = Bacterial_Order)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = color_palette) +
    theme_void() +
    theme(legend.position = "none")
}

# Create pie charts for each continent
pie_charts_updated <- list()
for (i in unique(continent_bacterial_counts$Continent)) {
  continent_data <- continent_bacterial_counts %>% filter(Continent == i)
  pie_chart <- create_pie_chart_no_legend(continent_data)
  pie_charts_updated[[i]] <- pie_chart
}

# Step 9: Create the bacterial diversity map
bacterial_diversity <- data.frame(
  Continent = c("Africa", "Asia", "Central America", "Europe", "North America", "Oceania", "South America"),
  Category = c("Undersampled", "Undersampled", "Undersampled", "Sampled", "Oversampled", "Undersampled", "Sampled")
)
# Define the colors for each category
category_colors <- c(
  "Oversampled" = "#ffcccc",    # Light red (North America)
  "Sampled" = "#ff9999",        # Medium red (Europe, South America)
  "Undersampled" = "#660000"    # Darkest red (Oceania, Africa, Asia, Central America)
)

world <- ne_countries(scale = "medium", returnclass = "sf")
world$Continent <- ifelse(world$region_un == "Africa", "Africa",
                          ifelse(world$region_un == "Asia", "Asia",
                                 ifelse(world$region_un == "Europe", "Europe",
                                        ifelse(world$region_un == "Oceania", "Oceania",
                                               ifelse(world$region_un == "Americas" & world$subregion == "South America", "South America",
                                                      ifelse(world$region_un == "Americas" & world$subregion == "Central America", "Central America",
                                                             ifelse(world$region_un == "Americas", "North America", NA)))))))

world <- merge(world, bacterial_diversity, by.x = "Continent", by.y = "Continent", all.x = TRUE)

# Combine the bacterial diversity map with pie charts
combined_pie_chart_map <- ggdraw() +
  draw_plot(world_map, 0, 0, 1, 1) +  
  draw_plot(pie_charts_updated[['Africa']], 0.45, 0.50, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Asia']], 0.58, 0.56, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Europe']], 0.44, 0.65, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['North America']], 0.18, 0.64, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['South America']], 0.27, 0.4, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Oceania']], 0.71, 0.38, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Central America']], 0.18, 0.53, 0.10, 0.10)

# Save the combined pie chart map
ggsave("final_combined_pie_chart_map_justNGS.pdf", plot = combined_pie_chart_map, width = 15, height = 10)


##### legend #####
# Step 6: Define the bacterial orders color palette
color_palette <- c("Rickettsiales" = "#FF0000", "Bacillales" = "#222222", "Enterobacterales" = "#f3c300", 
                   "Entomoplasmatales" = "#875692", "Hyphomicrobiales" = "#f38400", "Kitasatosporales" = "#a1caf1",
                   "Pseudomonadales" = "#be0032", "Pseudonocardiales" = "#c2b280", "Unclassified Firmicutes" = "#848482",
                   "Actinomycetales" = "#008856", "Burkholderiales" = "#e68fac", "Lactobacillales" = "#0067a5",
                   "Moraxellales" = "#f99379", "Sphingobacteriales" = "#604e97", "Xanthomonadales" = "#f6a600",
                   "Micrococcales" = "#b3446c", "Mycobacteriales" = "#dcd300", "Rhodospirillales" = "#882d17",
                   "Terriglobales" = "#8db600", "Campylobacterales" = "#654522", "Flavobacteriales" = "#e25822",
                   "Opitutales" = "#2b3d26", "Pasteurellales" = "#DCF1B9",
                   "Acetobacterales" = "#E4E3C2", "Bacteroidales" = "#C4B4CF", "Nitrospinales" = "#EE8B86",
                   "Caulobacteriales" = "#BD98A2", "Sphingomonadales" = "#99B1BC", "Erysipelotrichales" = "#F0B36D",
                   "Solirubrobacterales" = "#D0CD66", "Lysobacterales" = "#C8D88E", "Propionibacteriales" = "#FCCDE5")

# Step 8: Create a dummy plot just to generate the legend
dummy_plot <- ggplot(data.frame(Bacterial_Order = names(color_palette)), aes(x = Bacterial_Order, fill = Bacterial_Order)) +
  geom_bar() +
  scale_fill_manual(values = color_palette) +
  theme_void() +
  guides(fill = guide_legend(title = "Bacterial Orders"))

# Extract the legend
legend <- get_legend(dummy_plot)

# Display the legend on its own
plot_grid(legend, ncol = 1)

# Optionally, save the legend to a PDF or PNG file
ggsave("bacterial_order_legend.pdf", plot = legend, width = 8, height = 8)


# bar chart #
# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(rnaturalearth)
library(sf)
library(cowplot)

# Load the original data
ant_review_data <- read_csv("ant_review_spreadsheet_sanger_ngs_Chansonadded.csv")

# Clean the data by converting "NA" strings to actual NA values
ant_review_clean <- ant_review_data %>%
  mutate(across(where(is.character), ~na_if(., "NA")))

# Transform the data by unpivoting bacterial orders
ant_review_long <- ant_review_clean %>%
  select(Genus, Continent, Bacterial_order1, Bacterial_order2, Bacterial_order3) %>%
  pivot_longer(cols = c(Bacterial_order1, Bacterial_order2, Bacterial_order3),
               names_to = "OrderType", values_to = "Bacterial_Order") %>%
  filter(!is.na(Bacterial_Order))

# Aggregate data for genus-level bacterial counts
genus_bacterial_counts <- ant_review_long %>%
  group_by(Genus, Bacterial_Order) %>%
  summarise(Count = n(), .groups = "drop")

# Aggregate data for continent-level bacterial counts
continent_bacterial_counts <- ant_review_long %>%
  group_by(Continent, Bacterial_Order) %>%
  summarise(Count = n(), .groups = "drop")

# Save the cleaned results to CSV files
write_csv(genus_bacterial_counts, "Cleaned_Full_Bacterial_Order_Counts_sanger_ngs.csv")
write_csv(continent_bacterial_counts, "Cleaned_Microbial_Summary_by_Continent_sanger_ngs.csv")

# Define the bacterial orders color palette
color_palette <- c(
  "Rickettsiales" = "#FF0000", "Bacillales" = "#222222", "Enterobacterales" = "#f3c300",
  "Entomoplasmatales" = "#875692", "Hyphomicrobiales" = "#f38400", "Kitasatosporales" = "#a1caf1",
  "Pseudomonadales" = "#be0032", "Pseudonocardiales" = "#c2b280", "Unclassified Firmicutes" = "#848482",
  "Actinomycetales" = "#008856", "Burkholderiales" = "#e68fac", "Lactobacillales" = "#0067a5",
  "Moraxellales" = "#f99379", "Sphingobacteriales" = "#604e97", "Xanthomonadales" = "#f6a600",
  "Micrococcales" = "#b3446c", "Mycobacteriales" = "#dcd300", "Rhodospirillales" = "#882d17",
  "Terriglobales" = "#8db600", "Campylobacterales" = "#654522", "Flavobacteriales" = "#e25822",
  "Opitutales" = "#2b3d26", "Rhizobiales" = "#8DD3C7", "Pasteurellales" = "#DCF1B9",
  "Acetobacterales" = "#E4E3C2", "Bacteroidales" = "#C4B4CF", "Nitrospinales" = "#EE8B86",
  "Caulobacteriales" = "#BD98A2", "Sphingomonadales" = "#99B1BC", "Erysipelotrichales" = "#F0B36D",
  "Solirubrobacterales" = "#D0CD66", "Lysobacterales" = "#C8D88E", "Propionibacteriales" = "#FCCDE5"
)

# Add coordinates for continents (adjusted positions)
continent_coords <- data.frame(
  Continent = c("Africa", "Asia", "Europe", "North America", "South America", "Oceania", "Central America"),
  Lon = c(10, 80, 0, -130, -65, 110, -100),
  Lat = c(0, 30, 50, 40, -10, -25, 15)
)

# Merge the coordinates
continent_bacterial_counts <- continent_bacterial_counts %>%
  left_join(continent_coords, by = "Continent")

# Create a function for bar charts without count labels and with custom y-axis colors
create_bar_chart_no_legend <- function(data_row, y_axis_color = "black") {
  data_row <- data_row %>%
    mutate(Bacterial_Order = fct_reorder(Bacterial_Order, Count, .desc = TRUE))  # Order by abundance
  
  ggplot(data_row, aes(x = Bacterial_Order, y = Count, fill = Bacterial_Order)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = color_palette) +
    theme_minimal() +
    labs(x = NULL, y = "Count") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),  # Remove x-axis text
          axis.title.y = element_text(size = 10, color = y_axis_color, face = "bold"),
          axis.text.y = element_text(size = 8, color = y_axis_color, face = "bold"),
          panel.grid = element_blank())                                           # Remove grid lines
}

# Create bar charts for each continent with custom y-axis colors
bar_charts_updated <- list()
for (i in unique(continent_bacterial_counts$Continent)) {
  continent_data <- continent_bacterial_counts %>% filter(Continent == i)
  
  if (i == "Africa" || i == "Asia") {
    bar_chart <- create_bar_chart_no_legend(continent_data, y_axis_color = "white")
  } else if (i == "Oceania" || i == "Central America") {
    bar_chart <- create_bar_chart_no_legend(continent_data, y_axis_color = "darkgrey")
  } else {
    bar_chart <- create_bar_chart_no_legend(continent_data, y_axis_color = "black")
  }
  
  bar_charts_updated[[i]] <- bar_chart
}

# Create the bacterial diversity map
bacterial_diversity <- data.frame(
  Continent = c("Africa", "Asia", "Central America", "Europe", "North America", "Oceania", "South America"),
  Category = c("Undersampled", "Undersampled", "Undersampled", "Sampled", "Oversampled", "Undersampled", "Sampled")
)

# Define the colors for each category
category_colors <- c(
  "Oversampled" = "#ffcccc",    # Light red (North America)
  "Sampled" = "#ff9999",        # Medium red (Europe, South America)
  "Undersampled" = "#660000"    # Darkest red (Oceania, Africa, Asia, Central America)
)

# Load the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Assign continent categories to the world map
world$Continent <- ifelse(world$region_un == "Africa", "Africa",
                          ifelse(world$region_un == "Asia", "Asia",
                                 ifelse(world$region_un == "Europe", "Europe",
                                        ifelse(world$region_un == "Oceania", "Oceania",
                                               ifelse(world$region_un == "Americas" & world$subregion == "South America", "South America",
                                                      ifelse(world$region_un == "Americas" & world$subregion == "Central America", "Central America",
                                                             ifelse(world$region_un == "Americas", "North America", NA)))))))

# Merge the bacterial diversity categories with the world map data
world <- merge(world, bacterial_diversity, by.x = "Continent", by.y = "Continent", all.x = TRUE)

# Create the map with custom colors for each sampling category
world_map <- ggplot(data = world) +
  geom_sf(aes(fill = Category), color = "black") +
  scale_fill_manual(values = category_colors) +
  theme_void() +
  labs(fill = "Sampling Level")

# Combine the bacterial diversity map with bar charts
combined_bar_chart_map <- ggdraw() +
  draw_plot(world_map, 0, 0, 1, 1) +
  draw_plot(bar_charts_updated[['Africa']], 0.45, 0.50, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Asia']], 0.58, 0.56, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Europe']], 0.44, 0.65, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['North America']], 0.18, 0.64, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['South America']], 0.27, 0.4, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Oceania']], 0.71, 0.38, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Central America']], 0.18, 0.53, 0.10, 0.10)

# Save the combined bar chart map
ggsave("final_combined_bar_chart_map_sanger_ngs.pdf", plot = combined_bar_chart_map, width = 15, height = 10)

##### Bar Chart map sanger and NGS #####
# Load libraries
library(ggplot2)
library(rnaturalearth)
library(sf)
library(cowplot)

# Step 2: Load the data
ant_review_data <- read_csv("ant_review_spreadsheet_sanger_ngs_Chansonadded.csv")

# Step 3: Clean the data
ant_review_data_clean <- ant_review_data %>%
  mutate(across(where(is.character), ~na_if(., "NA"))) 

# Step 4: Transform the data by unpivoting bacterial orders
ant_review_long <- ant_review_data_clean %>%
  select(Genus, Bacterial_order1, Bacterial_order2, Bacterial_order3, Continent) %>%
  pivot_longer(cols = c(Bacterial_order1, Bacterial_order2, Bacterial_order3),
               names_to = "OrderType", values_to = "Bacterial_Order") %>%
  filter(!is.na(Bacterial_Order) & !is.na(Continent))

# Step 5: Aggregate data by continent and bacterial order
continent_bacterial_counts <- ant_review_long %>%
  group_by(Continent, Bacterial_Order) %>%
  summarise(Count = n(), .groups = "drop")

# Step 6: Define the bacterial orders color palette
color_palette <- c("Rickettsiales" = "#FF0000", "Bacillales" = "#222222", "Enterobacterales" = "#f3c300", 
                   "Entomoplasmatales" = "#875692", "Hyphomicrobiales" = "#f38400", "Kitasatosporales" = "#a1caf1",
                   "Pseudomonadales" = "#be0032", "Pseudonocardiales" = "#c2b280", "Unclassified Firmicutes" = "#848482",
                   "Actinomycetales" = "#008856", "Burkholderiales" = "#e68fac", "Lactobacillales" = "#0067a5",
                   "Moraxellales" = "#f99379", "Sphingobacteriales" = "#604e97", "Xanthomonadales" = "#f6a600",
                   "Micrococcales" = "#b3446c", "Mycobacteriales" = "#dcd300", "Rhodospirillales" = "#882d17",
                   "Terriglobales" = "#8db600", "Campylobacterales" = "#654522", "Flavobacteriales" = "#e25822",
                   "Opitutales" = "#2b3d26", "Rhizobiales" = "#8DD3C7", "Pasteurellales" = "#DCF1B9",
                   "Acetobacterales" = "#E4E3C2", "Bacteroidales" = "#C4B4CF", "Nitrospinales" = "#EE8B86",
                   "Caulobacteriales" = "#BD98A2", "Sphingomonadales" = "#99B1BC", "Erysipelotrichales" = "#F0B36D",
                   "Solirubrobacterales" = "#D0CD66", "Lysobacterales" = "#C8D88E", "Propionibacteriales" = "#FCCDE5")

# Step 7: Add coordinates for continents (adjusted positions)
continent_coords <- data.frame(
  Continent = c("Africa", "Asia", "Europe", "North America", "South America", "Oceania", "Central America"),
  Lon = c(10, 80, 0, -130, -65, 110, -100),  # North America, Australia, Central America moved left
  Lat = c(0, 30, 50, 40, -10, -25, 15)
)

# Merge the coordinates
continent_bacterial_counts <- continent_bacterial_counts %>%
  left_join(continent_coords, by = "Continent")

# Step 8: Create a function for pie charts
create_pie_chart_no_legend <- function(data_row) {
  ggplot(data_row, aes(x = "", y = Count, fill = Bacterial_Order)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = color_palette) +
    theme_void() +
    theme(legend.position = "none")
}

# Create pie charts for each continent
pie_charts_updated <- list()
for (i in unique(continent_bacterial_counts$Continent)) {
  continent_data <- continent_bacterial_counts %>% filter(Continent == i)
  pie_chart <- create_pie_chart_no_legend(continent_data)
  pie_charts_updated[[i]] <- pie_chart
}

# Step 9: Create the bacterial diversity map
bacterial_diversity <- data.frame(
  Continent = c("Africa", "Asia", "Central America", "Europe", "North America", "Oceania", "South America"),
  Category = c("Undersampled", "Undersampled", "Undersampled", "Sampled", "Oversampled", "Undersampled", "Sampled")
)
# Define the colors for each category
category_colors <- c(
  "Oversampled" = "#ffcccc",    # Light red (North America)
  "Sampled" = "#ff9999",        # Medium red (Europe, South America)
  "Undersampled" = "#660000"    # Darkest red (Oceania, Africa, Asia, Central America)
)

world <- ne_countries(scale = "medium", returnclass = "sf")
world$Continent <- ifelse(world$region_un == "Africa", "Africa",
                          ifelse(world$region_un == "Asia", "Asia",
                                 ifelse(world$region_un == "Europe", "Europe",
                                        ifelse(world$region_un == "Oceania", "Oceania",
                                               ifelse(world$region_un == "Americas" & world$subregion == "South America", "South America",
                                                      ifelse(world$region_un == "Americas" & world$subregion == "Central America", "Central America",
                                                             ifelse(world$region_un == "Americas", "North America", NA)))))))

world <- merge(world, bacterial_diversity, by.x = "Continent", by.y = "Continent", all.x = TRUE)

# Combine the bacterial diversity map with pie charts
combined_pie_chart_map <- ggdraw() +
  draw_plot(world_map, 0, 0, 1, 1) +  
  draw_plot(pie_charts_updated[['Africa']], 0.45, 0.50, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Asia']], 0.58, 0.56, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Europe']], 0.44, 0.65, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['North America']], 0.18, 0.64, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['South America']], 0.27, 0.4, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Oceania']], 0.71, 0.38, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Central America']], 0.18, 0.53, 0.10, 0.10)

# Save the combined pie chart map
ggsave("final_combined_pie_chart_map_sanger_ngs.pdf", plot = combined_pie_chart_map, width = 15, height = 10)

# bar chart for SANGER AND NGS #
# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(rnaturalearth)
library(sf)
library(cowplot)

# Load the original data
ant_review_data <- read_csv("ant_review_spreadsheet_sanger_ngs_Chansonadded.csv")

# Clean the data by converting "NA" strings to actual NA values
ant_review_clean <- ant_review_data %>%
  mutate(across(where(is.character), ~na_if(., "NA")))

# Transform the data by unpivoting bacterial orders
ant_review_long <- ant_review_clean %>%
  select(Genus, Continent, Bacterial_order1, Bacterial_order2, Bacterial_order3) %>%
  pivot_longer(cols = c(Bacterial_order1, Bacterial_order2, Bacterial_order3),
               names_to = "OrderType", values_to = "Bacterial_Order") %>%
  filter(!is.na(Bacterial_Order))

# Aggregate data for genus-level bacterial counts
genus_bacterial_counts <- ant_review_long %>%
  group_by(Genus, Bacterial_Order) %>%
  summarise(Count = n(), .groups = "drop")

# Aggregate data for continent-level bacterial counts
continent_bacterial_counts <- ant_review_long %>%
  group_by(Continent, Bacterial_Order) %>%
  summarise(Count = n(), .groups = "drop")

# Save the cleaned results to CSV files
write_csv(genus_bacterial_counts, "Cleaned_Full_Bacterial_Order_Counts_sanger_ngs.csv")
write_csv(continent_bacterial_counts, "Cleaned_Microbial_Summary_by_Continent_sanger_ngs.csv")

# Define the bacterial orders color palette
color_palette <- c(
  "Rickettsiales" = "#FF0000", "Bacillales" = "#222222", "Enterobacterales" = "#f3c300",
  "Entomoplasmatales" = "#875692", "Hyphomicrobiales" = "#f38400", "Kitasatosporales" = "#a1caf1",
  "Pseudomonadales" = "#be0032", "Pseudonocardiales" = "#c2b280", "Unclassified Firmicutes" = "#848482",
  "Actinomycetales" = "#008856", "Burkholderiales" = "#e68fac", "Lactobacillales" = "#0067a5",
  "Moraxellales" = "#f99379", "Sphingobacteriales" = "#604e97", "Xanthomonadales" = "#f6a600",
  "Micrococcales" = "#b3446c", "Mycobacteriales" = "#dcd300", "Rhodospirillales" = "#882d17",
  "Terriglobales" = "#8db600", "Campylobacterales" = "#654522", "Flavobacteriales" = "#e25822",
  "Opitutales" = "#2b3d26", "Rhizobiales" = "#8DD3C7", "Pasteurellales" = "#DCF1B9",
  "Acetobacterales" = "#E4E3C2", "Bacteroidales" = "#C4B4CF", "Nitrospinales" = "#EE8B86",
  "Caulobacteriales" = "#BD98A2", "Sphingomonadales" = "#99B1BC", "Erysipelotrichales" = "#F0B36D",
  "Solirubrobacterales" = "#D0CD66", "Lysobacterales" = "#C8D88E", "Propionibacteriales" = "#FCCDE5"
)

# Add coordinates for continents (adjusted positions)
continent_coords <- data.frame(
  Continent = c("Africa", "Asia", "Europe", "North America", "South America", "Oceania", "Central America"),
  Lon = c(10, 80, 0, -130, -65, 110, -100),
  Lat = c(0, 30, 50, 40, -10, -25, 15)
)

# Merge the coordinates
continent_bacterial_counts <- continent_bacterial_counts %>%
  left_join(continent_coords, by = "Continent")

# Create a function for bar charts without count labels and with custom y-axis colors
create_bar_chart_no_legend <- function(data_row, y_axis_color = "black") {
  data_row <- data_row %>%
    mutate(Bacterial_Order = fct_reorder(Bacterial_Order, Count, .desc = TRUE))  # Order by abundance
  
  ggplot(data_row, aes(x = Bacterial_Order, y = Count, fill = Bacterial_Order)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = color_palette) +
    theme_minimal() +
    labs(x = NULL, y = "Count") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),  # Remove x-axis text
          axis.title.y = element_text(size = 10, color = y_axis_color, face = "bold"),
          axis.text.y = element_text(size = 8, color = y_axis_color, face = "bold"),
          panel.grid = element_blank())                                           # Remove grid lines
}

# Create bar charts for each continent with custom y-axis colors
bar_charts_updated <- list()
for (i in unique(continent_bacterial_counts$Continent)) {
  continent_data <- continent_bacterial_counts %>% filter(Continent == i)
  
  if (i == "Africa" || i == "Asia") {
    bar_chart <- create_bar_chart_no_legend(continent_data, y_axis_color = "white")
  } else if (i == "Oceania" || i == "Central America") {
    bar_chart <- create_bar_chart_no_legend(continent_data, y_axis_color = "darkgrey")
  } else {
    bar_chart <- create_bar_chart_no_legend(continent_data, y_axis_color = "black")
  }
  
  bar_charts_updated[[i]] <- bar_chart
}

# Create the bacterial diversity map
bacterial_diversity <- data.frame(
  Continent = c("Africa", "Asia", "Central America", "Europe", "North America", "Oceania", "South America"),
  Category = c("Undersampled", "Undersampled", "Undersampled", "Sampled", "Oversampled", "Undersampled", "Sampled")
)

# Define the colors for each category
category_colors <- c(
  "Oversampled" = "#ffcccc",    # Light red (North America)
  "Sampled" = "#ff9999",        # Medium red (Europe, South America)
  "Undersampled" = "#660000"    # Darkest red (Oceania, Africa, Asia, Central America)
)

# Load the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Assign continent categories to the world map
world$Continent <- ifelse(world$region_un == "Africa", "Africa",
                          ifelse(world$region_un == "Asia", "Asia",
                                 ifelse(world$region_un == "Europe", "Europe",
                                        ifelse(world$region_un == "Oceania", "Oceania",
                                               ifelse(world$region_un == "Americas" & world$subregion == "South America", "South America",
                                                      ifelse(world$region_un == "Americas" & world$subregion == "Central America", "Central America",
                                                             ifelse(world$region_un == "Americas", "North America", NA)))))))

# Merge the bacterial diversity categories with the world map data
world <- merge(world, bacterial_diversity, by.x = "Continent", by.y = "Continent", all.x = TRUE)

# Create the map with custom colors for each sampling category
world_map <- ggplot(data = world) +
  geom_sf(aes(fill = Category), color = "black") +
  scale_fill_manual(values = category_colors) +
  theme_void() +
  labs(fill = "Sampling Level")

# Combine the bacterial diversity map with bar charts
combined_bar_chart_map <- ggdraw() +
  draw_plot(world_map, 0, 0, 1, 1) +
  draw_plot(bar_charts_updated[['Africa']], 0.45, 0.50, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Asia']], 0.58, 0.56, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Europe']], 0.44, 0.65, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['North America']], 0.18, 0.64, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['South America']], 0.27, 0.4, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Oceania']], 0.71, 0.38, 0.10, 0.10) +
  draw_plot(bar_charts_updated[['Central America']], 0.18, 0.53, 0.10, 0.10)

# Save the combined bar chart map
ggsave("final_combined_bar_chart_map_sanger_ngs.pdf", plot = combined_bar_chart_map, width = 15, height = 10)

##### Pie Chart map for SANGER AND NGS#####
# Load libraries
library(ggplot2)
library(rnaturalearth)
library(sf)
library(cowplot)

# Step 2: Load the data
ant_review_data <- read_csv("ant_review_spreadsheet_sanger_ngs_Chansonadded.csv")

# Step 3: Clean the data
ant_review_data_clean <- ant_review_data %>%
  mutate(across(where(is.character), ~na_if(., "NA"))) 

# Step 4: Transform the data by unpivoting bacterial orders
ant_review_long <- ant_review_data_clean %>%
  select(Genus, Bacterial_order1, Bacterial_order2, Bacterial_order3, Continent) %>%
  pivot_longer(cols = c(Bacterial_order1, Bacterial_order2, Bacterial_order3),
               names_to = "OrderType", values_to = "Bacterial_Order") %>%
  filter(!is.na(Bacterial_Order) & !is.na(Continent))

# Step 5: Aggregate data by continent and bacterial order
continent_bacterial_counts <- ant_review_long %>%
  group_by(Continent, Bacterial_Order) %>%
  summarise(Count = n(), .groups = "drop")

# Step 6: Define the bacterial orders color palette
color_palette <- c("Rickettsiales" = "#FF0000", "Bacillales" = "#222222", "Enterobacterales" = "#f3c300", 
                   "Entomoplasmatales" = "#875692", "Hyphomicrobiales" = "#f38400", "Kitasatosporales" = "#a1caf1",
                   "Pseudomonadales" = "#be0032", "Pseudonocardiales" = "#c2b280", "Unclassified Firmicutes" = "#848482",
                   "Actinomycetales" = "#008856", "Burkholderiales" = "#e68fac", "Lactobacillales" = "#0067a5",
                   "Moraxellales" = "#f99379", "Sphingobacteriales" = "#604e97", "Xanthomonadales" = "#f6a600",
                   "Micrococcales" = "#b3446c", "Mycobacteriales" = "#dcd300", "Rhodospirillales" = "#882d17",
                   "Terriglobales" = "#8db600", "Campylobacterales" = "#654522", "Flavobacteriales" = "#e25822",
                   "Opitutales" = "#2b3d26", "Rhizobiales" = "#8DD3C7", "Pasteurellales" = "#DCF1B9",
                   "Acetobacterales" = "#E4E3C2", "Bacteroidales" = "#C4B4CF", "Nitrospinales" = "#EE8B86",
                   "Caulobacteriales" = "#BD98A2", "Sphingomonadales" = "#99B1BC", "Erysipelotrichales" = "#F0B36D",
                   "Solirubrobacterales" = "#D0CD66", "Lysobacterales" = "#C8D88E", "Propionibacteriales" = "#FCCDE5")

# Step 7: Add coordinates for continents (adjusted positions)
continent_coords <- data.frame(
  Continent = c("Africa", "Asia", "Europe", "North America", "South America", "Oceania", "Central America"),
  Lon = c(10, 80, 0, -130, -65, 110, -100),  # North America, Australia, Central America moved left
  Lat = c(0, 30, 50, 40, -10, -25, 15)
)

# Merge the coordinates
continent_bacterial_counts <- continent_bacterial_counts %>%
  left_join(continent_coords, by = "Continent")

# Step 8: Create a function for pie charts
create_pie_chart_no_legend <- function(data_row) {
  ggplot(data_row, aes(x = "", y = Count, fill = Bacterial_Order)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = color_palette) +
    theme_void() +
    theme(legend.position = "none")
}

# Create pie charts for each continent
pie_charts_updated <- list()
for (i in unique(continent_bacterial_counts$Continent)) {
  continent_data <- continent_bacterial_counts %>% filter(Continent == i)
  pie_chart <- create_pie_chart_no_legend(continent_data)
  pie_charts_updated[[i]] <- pie_chart
}

# Step 9: Create the bacterial diversity map
bacterial_diversity <- data.frame(
  Continent = c("Africa", "Asia", "Central America", "Europe", "North America", "Oceania", "South America"),
  Category = c("Undersampled", "Undersampled", "Undersampled", "Sampled", "Oversampled", "Undersampled", "Sampled")
)
# Define the colors for each category
category_colors <- c(
  "Oversampled" = "#ffcccc",    # Light red (North America)
  "Sampled" = "#ff9999",        # Medium red (Europe, South America)
  "Undersampled" = "#660000"    # Darkest red (Oceania, Africa, Asia, Central America)
)

world <- ne_countries(scale = "medium", returnclass = "sf")
world$Continent <- ifelse(world$region_un == "Africa", "Africa",
                          ifelse(world$region_un == "Asia", "Asia",
                                 ifelse(world$region_un == "Europe", "Europe",
                                        ifelse(world$region_un == "Oceania", "Oceania",
                                               ifelse(world$region_un == "Americas" & world$subregion == "South America", "South America",
                                                      ifelse(world$region_un == "Americas" & world$subregion == "Central America", "Central America",
                                                             ifelse(world$region_un == "Americas", "North America", NA)))))))

world <- merge(world, bacterial_diversity, by.x = "Continent", by.y = "Continent", all.x = TRUE)

# Combine the bacterial diversity map with pie charts
combined_pie_chart_map <- ggdraw() +
  draw_plot(world_map, 0, 0, 1, 1) +  
  draw_plot(pie_charts_updated[['Africa']], 0.45, 0.50, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Asia']], 0.58, 0.56, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Europe']], 0.44, 0.65, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['North America']], 0.18, 0.64, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['South America']], 0.27, 0.4, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Oceania']], 0.71, 0.38, 0.10, 0.10) +  
  draw_plot(pie_charts_updated[['Central America']], 0.18, 0.53, 0.10, 0.10)

# Save the combined pie chart map
ggsave("final_combined_pie_chart_map_sanger_ngs.pdf", plot = combined_pie_chart_map, width = 15, height = 10)


##### Association Figures Sanger and NGS ####

# Load necessary libraries
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(tools)

# Define custom color palette for bacterial orders
color_palette <- c(
  "Rickettsiales" = "#FF0000", "Bacillales" = "#222222", "Enterobacterales" = "#f3c300",
  "Entomoplasmatales" = "#875692", "Hyphomicrobiales" = "#f38400", "Kitasatosporales" = "#a1caf1",
  "Pseudomonadales" = "#be0032", "Pseudonocardiales" = "#c2b280", "Unclassified Firmicutes" = "#848482",
  "Actinomycetales" = "#008856", "Burkholderiales" = "#e68fac", "Lactobacillales" = "#0067a5",
  "Moraxellales" = "#f99379", "Sphingobacteriales" = "#604e97", "Xanthomonadales" = "#f6a600",
  "Micrococcales" = "#b3446c", "Mycobacteriales" = "#dcd300", "Rhodospirillales" = "#882d17",
  "Terriglobales" = "#8db600", "Campylobacterales" = "#654522", "Flavobacteriales" = "#e25822",
  "Opitutales" = "#2b3d26", "Pasteurellales" = "#DCF1B9",
  "Acetobacterales" = "#E4E3C2", "Bacteroidales" = "#C4B4CF", "Nitrospinales" = "#EE8B86",
  "Caulobacteriales" = "#BD98A2", "Sphingomonadales" = "#99B1BC", "Erysipelotrichales" = "#F0B36D",
  "Solirubrobacterales" = "#D0CD66", "Lysobacterales" = "#C8D88E", "Propionibacteriales" = "#FCCDE5"
)

# Load and clean the dataset
test <- read.csv("ant_review_spreadsheet_sanger_ngs_Chansonadded.csv")
test$Bacterial_order1 <- trimws(tolower(test$Bacterial_order1))
test$Bacterial_order2 <- trimws(tolower(test$Bacterial_order2))
test$Bacterial_order3 <- trimws(tolower(test$Bacterial_order3))
test$Diet <- trimws(test$Diet)
test$Continent <- trimws(test$Continent)
test$Subfamily <- trimws(test$Subfamily)

# Create a cleaned subset of the data, allowing at least one bacterial order to be non-empty
clean_test <- test %>%
  filter(
    (Bacterial_order1 != "" | Bacterial_order2 != "" | Bacterial_order3 != "") &
      Diet != "" & !is.na(Diet) &
      Continent != "" & !is.na(Continent) &
      Subfamily != "" & !is.na(Subfamily)
  )

# Combine bacterial order columns into long format for plotting
long_data_clean <- clean_test %>%
  gather(key = "Bacterial_order_type", value = "Bacterial_order", Bacterial_order1, Bacterial_order2, Bacterial_order3) %>%
  filter(!is.na(Bacterial_order) & Bacterial_order != "")

# Standardize bacterial order names to Title Case to match the color_palette
long_data_clean$Bacterial_order <- tools::toTitleCase(long_data_clean$Bacterial_order)

# Helper function to reorder both rows and columns based on abundance
reorder_table_colors <- function(tab, palette) {
  # Calculate the total abundance of each bacterial order (rows)
  total_abundance_row <- rowSums(tab)
  
  # Calculate the total abundance for each category (columns)
  total_abundance_col <- colSums(tab)
  
  # Reorder the table by decreasing abundance of rows (bacterial orders)
  tab <- tab[order(total_abundance_row, decreasing = TRUE), ]
  
  # Reorder the columns by decreasing abundance of categories (Diet, Continent, Subfamily)
  tab <- tab[, order(total_abundance_col, decreasing = TRUE)]
  
  # Reorder the colors to match the new bacterial order ordering
  ordered_bacterial_orders <- rownames(tab)
  colors <- palette[ordered_bacterial_orders]
  
  return(list(tab = tab, colors = colors))
}

# Helper function to calculate the maximum value for x-axis limits using colSums
get_max_value <- function(tab) {
  max(colSums(tab))
}

### Diet and Bacterial Orders (using clean subset)
tab_Diet_Bacterial_order <- table(long_data_clean$Bacterial_order, long_data_clean$Diet)

# Only keep valid bacterial orders that exist in the custom color palette
valid_bacterial_orders <- rownames(tab_Diet_Bacterial_order)[rownames(tab_Diet_Bacterial_order) %in% names(color_palette)]

# Check if valid bacterial orders exist
if (length(valid_bacterial_orders) > 0) {
  tab_Diet_Bacterial_order <- tab_Diet_Bacterial_order[valid_bacterial_orders, ]
  
  # Reorder the table and colors by abundance
  reordered_diet <- reorder_table_colors(tab_Diet_Bacterial_order, color_palette)
  tab_Diet_Bacterial_order <- reordered_diet$tab
  colors <- reordered_diet$colors
  
  # Set x-axis limit based on the actual data with a larger buffer
  x_max <- get_max_value(tab_Diet_Bacterial_order)
  xlim_diet_adjusted <- c(0, x_max * 1.3)  # Increase the buffer to 1.3 times the max value
  
  # Save the Diet vs. Bacterial Orders plot as PNG
  png("Diet_vs_Bacterial_Orders_sanger_NGS.png", width = 12, height = 10, units = "in", res = 300)
  par(mar = c(5, 20, 5, 2))  # Adjust margins for better fit
  barplot(tab_Diet_Bacterial_order,
          main = "Association between Diet and Bacterial Orders",
          col = colors,
          horiz = TRUE,
          cex.axis = 0.8,
          cex.lab = 0.8,
          cex.names = 1.2,  # Larger text on the Y-axis
          las = 1,
          xlim = xlim_diet_adjusted)
  dev.off()  # Save the file
} else {
  print("No valid bacterial orders found for Diet after filtering.")
}

### Continent and Bacterial Orders (using clean subset)
tab_Continent_Bacterial_order <- table(long_data_clean$Bacterial_order, long_data_clean$Continent)

# Only keep valid bacterial orders that exist in the custom color palette
valid_bacterial_orders <- rownames(tab_Continent_Bacterial_order)[rownames(tab_Continent_Bacterial_order) %in% names(color_palette)]

# Check if valid bacterial orders exist
if (length(valid_bacterial_orders) > 0) {
  tab_Continent_Bacterial_order <- tab_Continent_Bacterial_order[valid_bacterial_orders, ]
  
  # Reorder the table and colors by abundance
  reordered_continent <- reorder_table_colors(tab_Continent_Bacterial_order, color_palette)
  tab_Continent_Bacterial_order <- reordered_continent$tab
  colors <- reordered_continent$colors
  
  # Set x-axis limit based on the actual data with a larger buffer
  x_max <- get_max_value(tab_Continent_Bacterial_order)
  xlim_continent_adjusted <- c(0, x_max * 1.5)  # Increase the buffer to 1.5 times the max value
  
  # Save the Continent vs. Bacterial Orders plot as PNG
  png("Continent_vs_Bacterial_Orders_sanger_NGS.png", width = 12, height = 10, units = "in", res = 300)
  par(mar = c(5, 20, 5, 2))  # Adjust margins for better fit
  barplot(tab_Continent_Bacterial_order,
          main = "Association between Continent and Bacterial Orders",
          col = colors,
          horiz = TRUE,
          cex.axis = 0.8,
          cex.lab = 0.8,
          cex.names = 1.2,  # Larger text on the Y-axis
          las = 1,
          xlim = xlim_continent_adjusted)
  dev.off()  # Save the file
} else {
  print("No valid bacterial orders found for Continent after filtering.")
}

### Subfamily and Bacterial Orders (using clean subset)
tab_Subfamily_Bacterial_order <- table(long_data_clean$Bacterial_order, long_data_clean$Subfamily)

# Only keep valid bacterial orders that exist in the custom color palette
valid_bacterial_orders <- rownames(tab_Subfamily_Bacterial_order)[rownames(tab_Subfamily_Bacterial_order) %in% names(color_palette)]

# Check if valid bacterial orders exist
if (length(valid_bacterial_orders) > 0) {
  tab_Subfamily_Bacterial_order <- tab_Subfamily_Bacterial_order[valid_bacterial_orders, ]
  
  # Reorder the table and colors by abundance
  reordered_subfamily <- reorder_table_colors(tab_Subfamily_Bacterial_order, color_palette)
  tab_Subfamily_Bacterial_order <- reordered_subfamily$tab
  colors <- reordered_subfamily$colors
  
  # Set x-axis limit based on the actual data with a larger buffer
  x_max <- get_max_value(tab_Subfamily_Bacterial_order)
  xlim_subfamily_adjusted <- c(0, x_max * 1.3)
  
  # Save the Subfamily vs. Bacterial Orders plot as PNG
  png("Subfamily_vs_Bacterial_Orders_sanger_NGS.png", width = 12, height = 10, units = "in", res = 300)
  par(mar = c(5, 20, 5, 2))  # Adjust margins for better fit
  barplot(tab_Subfamily_Bacterial_order,
          main = "Association between Subfamily and Bacterial Orders",
          col = colors,
          horiz = TRUE,
          cex.axis = 0.8,
          cex.lab = 0.8,
          cex.names = 1.2,  # Larger text on the Y-axis
          las = 1,
          xlim = xlim_subfamily_adjusted)
  dev.off()  # Save the file
} else {
  print("No valid bacterial orders found for Subfamily after filtering.")
}

# Step 4: Chi-square tests using the clean subset
chisq.test(tab_Diet_Bacterial_order, simulate.p.value = TRUE, B = 5000)
chisq.test(tab_Continent_Bacterial_order, simulate.p.value = TRUE, B = 5000)
chisq.test(tab_Subfamily_Bacterial_order, simulate.p.value = TRUE, B = 5000)

# Create and save the legend for bacterial orders, sorted by abundance
# Calculate the total abundance for each bacterial order
long_data_clean$Bacterial_order <- tools::toTitleCase(long_data_clean$Bacterial_order)
bacterial_abundance <- table(long_data_clean$Bacterial_order)
sorted_bacterial_orders_by_abundance <- names(sort(bacterial_abundance, decreasing = TRUE))

# Create a blank plot for the legend
png("Bacterial_Orders_Legend_by_assocation.png", width = 8, height = 8, units = "in", res = 300)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Bacterial Orders Legend (by Abundance)")

# Add the legend to the blank plot, with two columns, sorted by abundance
legend("center", 
       legend = sorted_bacterial_orders_by_abundance, 
       fill = color_palette[sorted_bacterial_orders_by_abundance], 
       cex = 0.8,  # Adjust the size of the text
       bty = "n",  # Remove the box around the legend
       ncol = 2,   # Set number of columns
       title = "Bacterial Orders")

# Save the legend as a PNG
dev.off()


##### Association Figures just NGS ####

# Load necessary libraries
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(tools)

# Define custom color palette for bacterial orders
color_palette <- c(
  "Rickettsiales" = "#FF0000", "Bacillales" = "#222222", "Enterobacterales" = "#f3c300",
  "Entomoplasmatales" = "#875692", "Hyphomicrobiales" = "#f38400", "Kitasatosporales" = "#a1caf1",
  "Pseudomonadales" = "#be0032", "Pseudonocardiales" = "#c2b280", "Unclassified Firmicutes" = "#848482",
  "Actinomycetales" = "#008856", "Burkholderiales" = "#e68fac", "Lactobacillales" = "#0067a5",
  "Moraxellales" = "#f99379", "Sphingobacteriales" = "#604e97", "Xanthomonadales" = "#f6a600",
  "Micrococcales" = "#b3446c", "Mycobacteriales" = "#dcd300", "Rhodospirillales" = "#882d17",
  "Terriglobales" = "#8db600", "Campylobacterales" = "#654522", "Flavobacteriales" = "#e25822",
  "Opitutales" = "#2b3d26", "Pasteurellales" = "#DCF1B9",
  "Acetobacterales" = "#E4E3C2", "Bacteroidales" = "#C4B4CF", "Nitrospinales" = "#EE8B86",
  "Caulobacteriales" = "#BD98A2", "Sphingomonadales" = "#99B1BC", "Erysipelotrichales" = "#F0B36D",
  "Solirubrobacterales" = "#D0CD66", "Lysobacterales" = "#C8D88E", "Propionibacteriales" = "#FCCDE5"
)

# Load and clean the dataset
test <- read.csv("ant_review_spreadsheet_nosanger_Chansonadded.csv")
test$Bacterial_order1 <- trimws(tolower(test$Bacterial_order1))
test$Bacterial_order2 <- trimws(tolower(test$Bacterial_order2))
test$Bacterial_order3 <- trimws(tolower(test$Bacterial_order3))
test$Diet <- trimws(test$Diet)
test$Continent <- trimws(test$Continent)
test$Subfamily <- trimws(test$Subfamily)

# Create a cleaned subset of the data, allowing at least one bacterial order to be non-empty
clean_test <- test %>%
  filter(
    (Bacterial_order1 != "" | Bacterial_order2 != "" | Bacterial_order3 != "") &
      Diet != "" & !is.na(Diet) &
      Continent != "" & !is.na(Continent) &
      Subfamily != "" & !is.na(Subfamily)
  )

# Combine bacterial order columns into long format for plotting
long_data_clean <- clean_test %>%
  gather(key = "Bacterial_order_type", value = "Bacterial_order", Bacterial_order1, Bacterial_order2, Bacterial_order3) %>%
  filter(!is.na(Bacterial_order) & Bacterial_order != "")

# Standardize bacterial order names to Title Case to match the color_palette
long_data_clean$Bacterial_order <- tools::toTitleCase(long_data_clean$Bacterial_order)

# Helper function to reorder both rows and columns based on abundance
reorder_table_colors <- function(tab, palette) {
  # Calculate the total abundance of each bacterial order (rows)
  total_abundance_row <- rowSums(tab)
  
  # Calculate the total abundance for each category (columns)
  total_abundance_col <- colSums(tab)
  
  # Reorder the table by decreasing abundance of rows (bacterial orders)
  tab <- tab[order(total_abundance_row, decreasing = TRUE), ]
  
  # Reorder the columns by decreasing abundance of categories (Diet, Continent, Subfamily)
  tab <- tab[, order(total_abundance_col, decreasing = TRUE)]
  
  # Reorder the colors to match the new bacterial order ordering
  ordered_bacterial_orders <- rownames(tab)
  colors <- palette[ordered_bacterial_orders]
  
  return(list(tab = tab, colors = colors))
}

# Helper function to calculate the maximum value for x-axis limits using colSums
get_max_value <- function(tab) {
  max(colSums(tab))
}

### Diet and Bacterial Orders (using clean subset)
tab_Diet_Bacterial_order <- table(long_data_clean$Bacterial_order, long_data_clean$Diet)

# Only keep valid bacterial orders that exist in the custom color palette
valid_bacterial_orders <- rownames(tab_Diet_Bacterial_order)[rownames(tab_Diet_Bacterial_order) %in% names(color_palette)]

# Check if valid bacterial orders exist
if (length(valid_bacterial_orders) > 0) {
  tab_Diet_Bacterial_order <- tab_Diet_Bacterial_order[valid_bacterial_orders, ]
  
  # Reorder the table and colors by abundance
  reordered_diet <- reorder_table_colors(tab_Diet_Bacterial_order, color_palette)
  tab_Diet_Bacterial_order <- reordered_diet$tab
  colors <- reordered_diet$colors
  
  # Set x-axis limit based on the actual data with a larger buffer
  x_max <- get_max_value(tab_Diet_Bacterial_order)
  xlim_diet_adjusted <- c(0, x_max * 1.3)  # Increase the buffer to 1.3 times the max value
  
  # Save the Diet vs. Bacterial Orders plot as PNG
  png("Diet_vs_Bacterial_Orders_justNGS.png", width = 12, height = 10, units = "in", res = 300)
  par(mar = c(5, 20, 5, 2))  # Adjust margins for better fit
  barplot(tab_Diet_Bacterial_order,
          main = "Association between Diet and Bacterial Orders",
          col = colors,
          horiz = TRUE,
          cex.axis = 0.8,
          cex.lab = 0.8,
          cex.names = 1.2,  # Larger text on the Y-axis
          las = 1,
          xlim = xlim_diet_adjusted)
  dev.off()  # Save the file
} else {
  print("No valid bacterial orders found for Diet after filtering.")
}

### Continent and Bacterial Orders (using clean subset)
tab_Continent_Bacterial_order <- table(long_data_clean$Bacterial_order, long_data_clean$Continent)

# Only keep valid bacterial orders that exist in the custom color palette
valid_bacterial_orders <- rownames(tab_Continent_Bacterial_order)[rownames(tab_Continent_Bacterial_order) %in% names(color_palette)]

# Check if valid bacterial orders exist
if (length(valid_bacterial_orders) > 0) {
  tab_Continent_Bacterial_order <- tab_Continent_Bacterial_order[valid_bacterial_orders, ]
  
  # Reorder the table and colors by abundance
  reordered_continent <- reorder_table_colors(tab_Continent_Bacterial_order, color_palette)
  tab_Continent_Bacterial_order <- reordered_continent$tab
  colors <- reordered_continent$colors
  
  # Set x-axis limit based on the actual data with a larger buffer
  x_max <- get_max_value(tab_Continent_Bacterial_order)
  xlim_continent_adjusted <- c(0, x_max * 1.5)  # Increase the buffer to 1.5 times the max value
  
  # Save the Continent vs. Bacterial Orders plot as PNG
  png("Continent_vs_Bacterial_Orders_justNGS.png", width = 12, height = 10, units = "in", res = 300)
  par(mar = c(5, 20, 5, 2))  # Adjust margins for better fit
  barplot(tab_Continent_Bacterial_order,
          main = "Association between Continent and Bacterial Orders",
          col = colors,
          horiz = TRUE,
          cex.axis = 0.8,
          cex.lab = 0.8,
          cex.names = 1.2,  # Larger text on the Y-axis
          las = 1,
          xlim = xlim_continent_adjusted)
  dev.off()  # Save the file
} else {
  print("No valid bacterial orders found for Continent after filtering.")
}

### Subfamily and Bacterial Orders (using clean subset)
tab_Subfamily_Bacterial_order <- table(long_data_clean$Bacterial_order, long_data_clean$Subfamily)

# Only keep valid bacterial orders that exist in the custom color palette
valid_bacterial_orders <- rownames(tab_Subfamily_Bacterial_order)[rownames(tab_Subfamily_Bacterial_order) %in% names(color_palette)]

# Check if valid bacterial orders exist
if (length(valid_bacterial_orders) > 0) {
  tab_Subfamily_Bacterial_order <- tab_Subfamily_Bacterial_order[valid_bacterial_orders, ]
  
  # Reorder the table and colors by abundance
  reordered_subfamily <- reorder_table_colors(tab_Subfamily_Bacterial_order, color_palette)
  tab_Subfamily_Bacterial_order <- reordered_subfamily$tab
  colors <- reordered_subfamily$colors
  
  # Set x-axis limit based on the actual data with a larger buffer
  x_max <- get_max_value(tab_Subfamily_Bacterial_order)
  xlim_subfamily_adjusted <- c(0, x_max * 1.3)
  
  # Save the Subfamily vs. Bacterial Orders plot as PNG
  png("Subfamily_vs_Bacterial_Orders_justNGS.png", width = 12, height = 10, units = "in", res = 300)
  par(mar = c(5, 20, 5, 2))  # Adjust margins for better fit
  barplot(tab_Subfamily_Bacterial_order,
          main = "Association between Subfamily and Bacterial Orders",
          col = colors,
          horiz = TRUE,
          cex.axis = 0.8,
          cex.lab = 0.8,
          cex.names = 1.2,  # Larger text on the Y-axis
          las = 1,
          xlim = xlim_subfamily_adjusted)
  dev.off()  # Save the file
} else {
  print("No valid bacterial orders found for Subfamily after filtering.")
}

# Step 4: Chi-square tests using the clean subset
chisq.test(tab_Diet_Bacterial_order, simulate.p.value = TRUE, B = 5000)
chisq.test(tab_Continent_Bacterial_order, simulate.p.value = TRUE, B = 5000)
chisq.test(tab_Subfamily_Bacterial_order, simulate.p.value = TRUE, B = 5000)

dev.off()
####Association Figures by Bacterial Order Sanger and NGS ####

# Load necessary libraries
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(tools)

# Define the custom color palette for each subfamily
dark_palette <- c(
  "Myrmicinae" = "#4C72B0", 
  "Formicinae" = "#55A868", 
  "Dorylinae" = "#C44E52", 
  "Amblyoponinae" = "#8172B3", 
  "Dolichoderinae" = "#CCB974", 
  "Ectatomminae" = "#64B5CD", 
  "Ponerinae" = "#8C8C8C", 
  "Paraponerinae" = "#E58606", 
  "Pseudomyrmecinae" = "#B40F0F"
)

# Load and clean the dataset
test <- read.csv("ant_review_spreadsheet_sanger_ngs_Chansonadded.csv")
test$Bacterial_order1 <- trimws(tolower(test$Bacterial_order1))
test$Bacterial_order2 <- trimws(tolower(test$Bacterial_order2))
test$Bacterial_order3 <- trimws(tolower(test$Bacterial_order3))
test$Diet <- trimws(test$Diet)
test$Continent <- trimws(test$Continent)
test$Subfamily <- trimws(test$Subfamily)

# Combine bacterial order columns into long format
long_data_clean <- test %>%
  gather(key = "Bacterial_order_type", value = "Bacterial_order", Bacterial_order1, Bacterial_order2, Bacterial_order3) %>%
  filter(!is.na(Bacterial_order) & Bacterial_order != "")

# Standardize bacterial order names to Title Case
long_data_clean$Bacterial_order <- tools::toTitleCase(long_data_clean$Bacterial_order)

# Function to plot bacterial orders on Y-axis with custom colors for different groupings
plot_bacterial_orders_stacked <- function(data, group_var, title, output_file, palette = NULL, order_by_group = FALSE) {
  
  # Create a contingency table with Bacterial Orders as rows and the selected grouping variable (e.g., Subfamily, Country, Diet) as columns
  tab <- table(data$Bacterial_order, data[[group_var]])
  
  # Only keep valid bacterial orders that exist in the dataset
  valid_bacterial_orders <- rownames(tab)[rownames(tab) %in% unique(data$Bacterial_order)]
  
  if (length(valid_bacterial_orders) > 0) {
    tab <- tab[valid_bacterial_orders, ]
    
    # Reorder the grouping variable (e.g., Subfamily, Country, Diet) by abundance if requested
    if (order_by_group) {
      # Calculate the total abundance for each group
      total_abundance <- colSums(tab)
      ordered_groups <- names(sort(total_abundance, decreasing = TRUE))
      tab <- tab[, ordered_groups]  # Reorder by group (subfamilies, countries, or diets)
    }
    
    # Reorder bacterial orders by total abundance (row sums)
    bacterial_order_totals <- rowSums(tab)
    tab <- tab[order(bacterial_order_totals, decreasing = TRUE), ]
    
    # Get the group levels in the correct order of abundance
    group_levels <- colnames(tab)
    
    # Ensure the colors are assigned correctly to group levels
    group_colors <- palette[group_levels]
    
    # Dynamic x-axis limit based on the maximum value in the table + buffer (20% buffer)
    x_max <- max(rowSums(tab)) * 1.2  # Increase buffer to 20%
    
    # Plot using custom colors for the grouping variable
    png(output_file, width = 12, height = 10, units = "in", res = 300)
    par(mar = c(5, 20, 5, 2))  # Adjust margins for better readability
    
    # Plot using consistent colors for both the barplot and the legend
    barplot(t(tab),  # Use t() to stack across Bacterial Orders
            main = title,
            col = group_colors,  # Use the custom palette for the group
            horiz = TRUE,  # Horizontal bars
            cex.names = 1.2,  # Adjust Y-axis label size
            las = 1,  # Rotate Y-axis labels
            legend.text = group_levels,  # Legend for the group
            args.legend = list(x = "topright", bty = "n", fill = group_colors),  # Legend with matching colors
            xlab = "Count",
            xlim = c(0, x_max))  # Dynamic x-axis limit with 20% buffer
    
    dev.off()
  } else {
    print(paste("No valid bacterial orders found for", group_var))
  }
}

# Plot for Subfamily using manual color assignments and ordered by subfamily and bacterial order abundance
plot_bacterial_orders_stacked(long_data_clean, "Subfamily", 
                              "Bacterial Orders vs Subfamily (Stacked with Custom Colors, Ordered by Subfamily and Bacterial Orders)", 
                              "Bacterial_Orders_vs_Subfamily_Stacked_sanger_ngs.png", 
                              dark_palette, order_by_group = TRUE)

# Custom pastel palette for Country (more saturated pastel shades)
country_palette <- c("#FFB347", "#D8BFD8", "#C3B091", "#99CCFF", "#B0E57C","#FF9999", "#FF6961", "#FDFD96")
country_palette <- setNames(country_palette[1:length(unique(long_data_clean$Continent))], unique(long_data_clean$Continent))

# For Diet, use the Dark2 palette for more distinct colors
diet_palette <- setNames(brewer.pal(n = min(8, length(unique(long_data_clean$Diet))), name = "Dark2"), unique(long_data_clean$Diet))

# Plot for Country with a custom more saturated pastel palette and bacterial orders ordered by abundance
plot_bacterial_orders_stacked(long_data_clean, "Continent", 
                              "Bacterial Orders vs Country (Stacked, Ordered by Country and Bacterial Orders)", 
                              "Bacterial_Orders_vs_Country_Stacked_sanger_ngs.png", 
                              country_palette, 
                              order_by_group = TRUE)

# Plot for Diet with a distinct color palette (Dark2) and bacterial orders ordered by abundance
plot_bacterial_orders_stacked(long_data_clean, "Diet", 
                              "Bacterial Orders vs Diet (Stacked, Ordered by Diet and Bacterial Orders)", 
                              "Bacterial_Orders_vs_Diet_Stacked_sanger_ngs.png", 
                              diet_palette, 
                              order_by_group = TRUE)

####Association Figures by Bacterial Order just NGS ####

# Load necessary libraries
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(tools)

# Define the custom color palette for each subfamily
dark_palette <- c(
  "Myrmicinae" = "#4C72B0", 
  "Formicinae" = "#55A868", 
  "Dorylinae" = "#C44E52", 
  "Amblyoponinae" = "#8172B3", 
  "Dolichoderinae" = "#CCB974", 
  "Ectatomminae" = "#64B5CD", 
  "Ponerinae" = "#8C8C8C", 
  "Paraponerinae" = "#E58606", 
  "Pseudomyrmecinae" = "#B40F0F"
)

# Load and clean the dataset
test <- read.csv("ant_review_spreadsheet_nosanger_Chansonadded.csv")
test$Bacterial_order1 <- trimws(tolower(test$Bacterial_order1))
test$Bacterial_order2 <- trimws(tolower(test$Bacterial_order2))
test$Bacterial_order3 <- trimws(tolower(test$Bacterial_order3))
test$Diet <- trimws(test$Diet)
test$Continent <- trimws(test$Continent)
test$Subfamily <- trimws(test$Subfamily)

# Combine bacterial order columns into long format
long_data_clean <- test %>%
  gather(key = "Bacterial_order_type", value = "Bacterial_order", Bacterial_order1, Bacterial_order2, Bacterial_order3) %>%
  filter(!is.na(Bacterial_order) & Bacterial_order != "")

# Standardize bacterial order names to Title Case
long_data_clean$Bacterial_order <- tools::toTitleCase(long_data_clean$Bacterial_order)

# Function to plot bacterial orders on Y-axis with custom colors for different groupings
plot_bacterial_orders_stacked <- function(data, group_var, title, output_file, palette = NULL, order_by_group = FALSE) {
  
  # Create a contingency table with Bacterial Orders as rows and the selected grouping variable (e.g., Subfamily, Country, Diet) as columns
  tab <- table(data$Bacterial_order, data[[group_var]])
  
  # Only keep valid bacterial orders that exist in the dataset
  valid_bacterial_orders <- rownames(tab)[rownames(tab) %in% unique(data$Bacterial_order)]
  
  if (length(valid_bacterial_orders) > 0) {
    tab <- tab[valid_bacterial_orders, ]
    
    # Reorder the grouping variable (e.g., Subfamily, Country, Diet) by abundance if requested
    if (order_by_group) {
      # Calculate the total abundance for each group
      total_abundance <- colSums(tab)
      ordered_groups <- names(sort(total_abundance, decreasing = TRUE))
      tab <- tab[, ordered_groups]  # Reorder by group (subfamilies, countries, or diets)
    }
    
    # Reorder bacterial orders by total abundance (row sums)
    bacterial_order_totals <- rowSums(tab)
    tab <- tab[order(bacterial_order_totals, decreasing = TRUE), ]
    
    # Get the group levels in the correct order of abundance
    group_levels <- colnames(tab)
    
    # Ensure the colors are assigned correctly to group levels
    group_colors <- palette[group_levels]
    
    # Dynamic x-axis limit based on the maximum value in the table + buffer (20% buffer)
    x_max <- max(rowSums(tab)) * 1.2  # Increase buffer to 20%
    
    # Plot using custom colors for the grouping variable
    png(output_file, width = 12, height = 10, units = "in", res = 300)
    par(mar = c(5, 20, 5, 2))  # Adjust margins for better readability
    
    # Plot using consistent colors for both the barplot and the legend
    barplot(t(tab),  # Use t() to stack across Bacterial Orders
            main = title,
            col = group_colors,  # Use the custom palette for the group
            horiz = TRUE,  # Horizontal bars
            cex.names = 1.2,  # Adjust Y-axis label size
            las = 1,  # Rotate Y-axis labels
            legend.text = group_levels,  # Legend for the group
            args.legend = list(x = "topright", bty = "n", fill = group_colors),  # Legend with matching colors
            xlab = "Count",
            xlim = c(0, x_max))  # Dynamic x-axis limit with 20% buffer
    
    dev.off()
  } else {
    print(paste("No valid bacterial orders found for", group_var))
  }
}

# Plot for Subfamily using manual color assignments and ordered by subfamily and bacterial order abundance
plot_bacterial_orders_stacked(long_data_clean, "Subfamily", 
                              "Bacterial Orders vs Subfamily (Stacked with Custom Colors, Ordered by Subfamily and Bacterial Orders)", 
                              "Bacterial_Orders_vs_Subfamily_justNGS.png", 
                              dark_palette, order_by_group = TRUE)

# Custom pastel palette for Country (more saturated pastel shades)
country_palette <- c("#FFB347", "#D8BFD8", "#C3B091", "#99CCFF", "#B0E57C","#FF9999", "#FF6961", "#FDFD96")
country_palette <- setNames(country_palette[1:length(unique(long_data_clean$Continent))], unique(long_data_clean$Continent))

# For Diet, use the Dark2 palette for more distinct colors
diet_palette <- setNames(brewer.pal(n = min(8, length(unique(long_data_clean$Diet))), name = "Dark2"), unique(long_data_clean$Diet))

# Plot for Country with a custom more saturated pastel palette and bacterial orders ordered by abundance
plot_bacterial_orders_stacked(long_data_clean, "Continent", 
                              "Bacterial Orders vs Country (Stacked, Ordered by Country and Bacterial Orders)", 
                              "Bacterial_Orders_vs_Country_Stacked_justNGS.png", 
                              country_palette, 
                              order_by_group = TRUE)

# Plot for Diet with a distinct color palette (Dark2) and bacterial orders ordered by abundance
plot_bacterial_orders_stacked(long_data_clean, "Diet", 
                              "Bacterial Orders vs Diet (Stacked, Ordered by Diet and Bacterial Orders)", 
                              "Bacterial_Orders_vs_Diet_Stacked_justNGS.png", 
                              diet_palette, 
                              order_by_group = TRUE)


##### Phylogeny Figures just NGS ####
# Load necessary libraries
library(ggtree)
library(tidytree)
library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)  # For placing bar graphs on the map
library(forcats)  # For reordering factors

# Read the Newick tree (rerooted from iTOL)
tree <- read.tree("rerooted_tree_justNGS.txt")

# Remove branch lengths if necessary
tree$edge.length <- NULL

# Read the subfamily assignment file
subfamily_data <- read.table("genera_subfamily_ngs.txt", header = TRUE, sep = "\t")

# Merge the subfamily data with the tree data
tree_data <- as_tibble(tree)
tree_data <- left_join(tree_data, subfamily_data, by = c("label" = "Genus"))

# Assign darker colors to subfamilies
dark_palette <- c("#4C72B0", "#55A868", "#C44E52", "#CCB974", "#8C8C8C", "#E58606", "#B40F0F", "#80b1d3")
names(dark_palette) <- unique(subfamily_data$Subfamily)
subfamily_colors <- scale_color_manual(values = dark_palette)

# Plot the tree with bold labels and more spread out branches
p <- ggtree(tree) %<+% tree_data +
  geom_tiplab(aes(color = Subfamily), size = 4, offset = 0.1, fontface = "bold") +  # Bold labels with color
  subfamily_colors +
  theme_tree2() +    # Automatically manage space for branches and labels
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +  # Remove scale
  scale_x_continuous(expand = c(0.25, 0.25))  # Spread out branches more

# Create a separate data frame for the legend
legend_data <- data.frame(Subfamily = unique(subfamily_data$Subfamily))

# Add custom legend using geom_point()
p <- p + geom_point(data = legend_data, aes(x = Inf, y = Inf, color = Subfamily), size = 5, shape = 15, show.legend = TRUE) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5)))

# Save the plot to a file with larger dimensions
ggsave("final_tree_custom_legend_justNGS.pdf", plot = p, width = 18, height = 12)

### Ordered Bacterial Clades #
# Load bacterial data
file_path <- "Cleaned_Full_Bacterial_Order_Counts_nosanger.csv"
bacterial_counts_full <- read.csv(file_path)

# Assuming 'tree' is already loaded, with genera ordered according to the phylogenetic tree
ordered_genera <- rev(tree$tip.label)  # Reverse the order of the genera according to the phylogenetic tree

# Filter the bacterial data to only include genera present in the tree
bacterial_counts_filtered <- bacterial_counts_full %>%
  filter(Genus %in% ordered_genera)

# Ensure that the Genus column is ordered according to the reversed phylogenetic tree
bacterial_counts_filtered$Genus <- factor(bacterial_counts_filtered$Genus, levels = ordered_genera)

# Reorder Bacterial Orders by their total abundance (most abundant first)
bacterial_order_totals <- bacterial_counts_filtered %>%
  group_by(Bacterial_Order) %>%
  summarise(Total_Count = sum(Count)) %>%
  arrange(desc(Total_Count))

# Reverse the order so that the most abundant bacterial orders appear at the bottom
bacterial_counts_filtered$Bacterial_Order <- factor(bacterial_counts_filtered$Bacterial_Order, 
                                                    levels = rev(bacterial_order_totals$Bacterial_Order))

# Use the provided manual color palette for bacterial orders
color_palette <- c(
  "Rickettsiales" = "#FF0000",
  "Bacillales" = "#222222",
  "Enterobacterales" = "#f3c300",
  "Entomoplasmatales" = "#875692",
  "Hyphomicrobiales" = "#f38400",
  "Kitasatosporales" = "#a1caf1",
  "Pseudomonadales" = "#be0032",
  "Pseudonocardiales" = "#c2b280",
  "Unclassified Firmicutes" = "#848482",
  "Actinomycetales" = "#008856",
  "Burkholderiales" = "#e68fac",
  "Lactobacillales" = "#0067a5",
  "Moraxellales" = "#f99379",
  "Sphingobacteriales" = "#604e97",
  "Xanthomonadales" = "#f6a600",
  "Micrococcales" = "#b3446c",
  "Mycobacteriales" = "#dcd300",
  "Rhodospirillales" = "#882d17",
  "Terriglobales" = "#8db600",
  "Campylobacterales" = "#654522",
  "Flavobacteriales" = "#e25822",
  "Opitutales" = "#2b3d26",
  "Rhizobiales" = "#8DD3C7",
  "Pasteurellales" = "#DCF1B9",
  "Acetobacterales" = "#E4E3C2",
  "Bacteroidales" = "#C4B4CF",
  "Nitrospinales" = "#EE8B86",
  "Caulobacteriales" = "#BD98A2",
  "Sphingomonadales" = "#99B1BC",
  "Erysipelotrichales" = "#F0B36D",
  "Solirubrobacterales" = "#D0CD66",
  "Lysobacterales" = "#C8D88E",
  "Propionibacteriales" = "#FCCDE5"
)

# Ensure that the color palette matches the bacterial orders in the filtered dataset
bacterial_orders_unique <- unique(bacterial_counts_filtered$Bacterial_Order)
color_palette <- color_palette[names(color_palette) %in% bacterial_orders_unique]

# Apply the reordered color palette to the plot and reorder the legend
ggplot(bacterial_counts_filtered, aes(x = Genus, y = Count, fill = Bacterial_Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Bacterial Orders per Genus (Reversed Phylogeny Order)", x = "Genus", y = "Count", fill = "Bacterial Order") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +  # Adjust label size and angle
  scale_x_discrete(limits = ordered_genera) +  # Use the reversed order
  scale_y_continuous(expand = c(0, 0)) +  # Ensure y-axis starts at 0
  scale_fill_manual(values = color_palette, guide = guide_legend(reverse = TRUE)) +  # Apply reversed legend order
  coord_cartesian(clip = "off")  # Prevent text from being cut off

# Save the bar plot to a file
ggsave("bacterial_orders_by_reversed_phylogeny_justNGS.pdf", width = 14, height = 10)

##### Phylogeny Figures- Sanger and NGS ####
# Load necessary libraries
library(ggtree)
library(tidytree)
library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)  # For placing bar graphs on the map
library(forcats)  # For reordering factors

# Read the Newick tree (rerooted from iTOL)
tree <- read.tree("pruned_tree_sanger_ngs.txt")  # Use the rerooted tree from iTOL

# Remove branch lengths if necessary
tree$edge.length <- NULL

# Read the subfamily assignment file
subfamily_data <- read.table("genera_subfamily.txt", header = TRUE, sep = "\t")

# Merge the subfamily data with the tree data
tree_data <- as_tibble(tree)
tree_data <- left_join(tree_data, subfamily_data, by = c("label" = "Genus"))

# Assign darker colors to subfamilies
dark_palette <- c("#4C72B0", "#55A868", "#C44E52", "#8172B3", "#CCB974", "#64B5CD", "#8C8C8C", "#E58606", "#B40F0F")
names(dark_palette) <- unique(subfamily_data$Subfamily)
subfamily_colors <- scale_color_manual(values = dark_palette)

# Plot the tree with bold labels and more spread out branches
p <- ggtree(tree) %<+% tree_data +
  geom_tiplab(aes(color = Subfamily), size = 4, offset = 0.1, fontface = "bold") +  # Bold labels with color
  subfamily_colors +
  theme_tree2() +    # Automatically manage space for branches and labels
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +  # Remove scale
  scale_x_continuous(expand = c(0.25, 0.25))  # Spread out branches more

# Create a separate data frame for the legend
legend_data <- data.frame(Subfamily = unique(subfamily_data$Subfamily))

# Add custom legend using geom_point()
p <- p + geom_point(data = legend_data, aes(x = Inf, y = Inf, color = Subfamily), size = 5, shape = 15, show.legend = TRUE) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5)))

# Save the plot to a file with larger dimensions
ggsave("final_tree_custom_legend_sanger_ngs.pdf", plot = p, width = 18, height = 12)

### Ordered Bacterial Clades #
# Load bacterial data
file_path <- "Cleaned_Full_Bacterial_Order_Counts_sanger_ngs.csv"
bacterial_counts_full <- read.csv(file_path)

# Assuming 'tree' is already loaded, with genera ordered according to the phylogenetic tree
ordered_genera <- rev(tree$tip.label)  # Reverse the order of the genera according to the phylogenetic tree

# Filter the bacterial data to only include genera present in the tree
bacterial_counts_filtered <- bacterial_counts_full %>%
  filter(Genus %in% ordered_genera)

# Ensure that the Genus column is ordered according to the reversed phylogenetic tree
bacterial_counts_filtered$Genus <- factor(bacterial_counts_filtered$Genus, levels = ordered_genera)

# Reorder Bacterial Orders by their total abundance (most abundant first)
bacterial_order_totals <- bacterial_counts_filtered %>%
  group_by(Bacterial_Order) %>%
  summarise(Total_Count = sum(Count)) %>%
  arrange(desc(Total_Count))

# Reverse the order so that the most abundant bacterial orders appear at the bottom
bacterial_counts_filtered$Bacterial_Order <- factor(bacterial_counts_filtered$Bacterial_Order, 
                                                    levels = rev(bacterial_order_totals$Bacterial_Order))

# Use the provided manual color palette for bacterial orders
color_palette <- c(
  "Rickettsiales" = "#FF0000",
  "Bacillales" = "#222222",
  "Enterobacterales" = "#f3c300",
  "Entomoplasmatales" = "#875692",
  "Hyphomicrobiales" = "#f38400",
  "Kitasatosporales" = "#a1caf1",
  "Pseudomonadales" = "#be0032",
  "Pseudonocardiales" = "#c2b280",
  "Unclassified Firmicutes" = "#848482",
  "Actinomycetales" = "#008856",
  "Burkholderiales" = "#e68fac",
  "Lactobacillales" = "#0067a5",
  "Moraxellales" = "#f99379",
  "Sphingobacteriales" = "#604e97",
  "Xanthomonadales" = "#f6a600",
  "Micrococcales" = "#b3446c",
  "Mycobacteriales" = "#dcd300",
  "Rhodospirillales" = "#882d17",
  "Terriglobales" = "#8db600",
  "Campylobacterales" = "#654522",
  "Flavobacteriales" = "#e25822",
  "Opitutales" = "#2b3d26",
  "Rhizobiales" = "#8DD3C7",
  "Pasteurellales" = "#DCF1B9",
  "Acetobacterales" = "#E4E3C2",
  "Bacteroidales" = "#C4B4CF",
  "Nitrospinales" = "#EE8B86",
  "Caulobacteriales" = "#BD98A2",
  "Sphingomonadales" = "#99B1BC",
  "Erysipelotrichales" = "#F0B36D",
  "Solirubrobacterales" = "#D0CD66",
  "Lysobacterales" = "#C8D88E",
  "Propionibacteriales" = "#FCCDE5"
)

# Ensure that the color palette matches the bacterial orders in the filtered dataset
bacterial_orders_unique <- unique(bacterial_counts_filtered$Bacterial_Order)
color_palette <- color_palette[names(color_palette) %in% bacterial_orders_unique]

# Apply the reordered color palette to the plot and reorder the legend
ggplot(bacterial_counts_filtered, aes(x = Genus, y = Count, fill = Bacterial_Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Bacterial Orders per Genus (Reversed Phylogeny Order)", x = "Genus", y = "Count", fill = "Bacterial Order") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +  # Adjust label size and angle
  scale_x_discrete(limits = ordered_genera) +  # Use the reversed order
  scale_y_continuous(expand = c(0, 0)) +  # Ensure y-axis starts at 0
  scale_fill_manual(values = color_palette, guide = guide_legend(reverse = TRUE)) +  # Apply reversed legend order
  coord_cartesian(clip = "off")  # Prevent text from being cut off

# Save the bar plot to a file
ggsave("bacterial_orders_by_reversed_phylogeny_sanger_ngs.pdf", width = 14, height = 10)
