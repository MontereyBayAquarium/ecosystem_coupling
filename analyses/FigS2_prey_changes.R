

rm(list=ls())

################################################################################
#Prep workspace

librarian::shelf(tidyverse, gtools, ggplot2, googledrive)

#Set directories and load data
datin <- here::here("output","processed")
figdir <- here::here("figures")

# Define the Google Drive file ID for sofa output
# This file exceeds 100MB GitHub limit so is stored externally
file_id <- "1yC8yrS69uOqbvlzQJakL0z-DC-yMhNel"

temp_file <- tempfile(fileext = ".rdata")
drive_download(as_id(file_id), path = temp_file, overwrite = TRUE)
SOFA_rslts <- mget(load(temp_file))

################################################################################
#Prep workspace

# Check structure of Fdat to identify prey-related fields
str(SOFA_rslts$Fdat)

# Extract relevant prey data from Fdat and exclude 'Sho'
df_prey <- SOFA_rslts$Fdat %>%
  select(Date, Prey, N_items) %>%  # Select date, prey type, and number of items (or another effort metric)
  filter(!is.na(Prey) & Prey != "" & Prey != "sho") %>%  # Filter out rows with no prey data and exclude 'Sho'
  mutate(Year = as.numeric(format(Date, "%Y")),
         Period = ifelse(Year < 2013, "Pre-2013", "Post-2013"))  # Classify periods

# Calculate pre- and post-2013 averages for each prey type
pre_post_averages <- df_prey %>%
  group_by(Prey, Period) %>%
  summarise(
    mean_items = mean(N_items, na.rm = TRUE),  # Calculate mean number of items captured (or other metric)
    .groups = 'drop'
  )

# Calculate percent change for each prey type from pre-2013 to post-2013
percent_change_prey <- pre_post_averages %>%
  pivot_wider(names_from = Period, values_from = mean_items) %>%
  mutate(
    percent_change = ((`Post-2013` - `Pre-2013`) / `Pre-2013`) * 100
  ) %>%
  select(Prey, percent_change) %>%
  arrange(desc(percent_change)) %>%
  filter(!is.na(percent_change))

# Print percent change summary for all prey types
View(percent_change_prey)

################################################################################
#Plot
base_theme <-  theme(axis.text=element_text(size=8, color = "black"),
                     axis.title=element_text(size=9,color = "black"),
                     plot.tag=element_text(size=8,color = "black"),
                     plot.title=element_text(size=9,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.4, "cm"), 
                     #legend.key = element_rect(fill = "white"), # Set it to transparent
                     legend.spacing.y = unit(0.4, "cm"),  
                     legend.text=element_text(size=8,color = "black"),
                     legend.title=element_blank(),
                     #legend.key.height = unit(0.1, "cm"),
                     #legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=8, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())


# Create the tornado plot for all prey items excluding 'Sho'
p1 <- ggplot(percent_change_prey, aes(x = percent_change, y = reorder(Prey, percent_change))) +
  geom_bar(stat = "identity", aes(fill = percent_change > 0), width = 0.6) +
  scale_fill_manual(values = c("TRUE" = "indianred", "FALSE" = "navyblue"), guide = FALSE) +
  labs(x = "Percent Change in Foraging Effort (%)", y = "Prey Type") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )+
  theme_bw() + base_theme
p1



