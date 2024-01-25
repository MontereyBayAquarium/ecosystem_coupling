#Joshua G. Smith; jossmith@mbayaq.org

#This processing script reads the raw census count data obtained from USGS
#and converts total counts per ATOS to density estimates. Census data were obtained
#by request from USGS. 

rm(list=ls())

######
#required packages
librarian::shelf(tidyverse,sf, janitor)

#set directories 
basedir <- here::here("output")

#read data
census_orig <- read_csv(file.path(basedir,"raw/census_data/ATOS_362-449_modifed_census_results.csv"))
atos_orig <- st_read(file.path(basedir,"raw/gis_data/atos/ATOS_polygon_teale83.shp"))


################################################################################
#Step 1 - process ATOS geometries

atos_build1 <- atos_orig %>%
               #merge offshore / onshore geometrues
               group_by(ATOS_ID) %>%
               summarize(geometry = st_union(geometry))%>%
               mutate(area_m2 = st_area(geometry),
                      area_km2 = area_m2 / 1e6,
                      #format
                      area_km2 = as.numeric(sprintf("%.7f", area_km2))
                      ) %>%
                clean_names()


################################################################################
#Step 2 - inspect and process CENSUS data

census_build1 <- census_orig %>%
                  #drop index column
                  select(-`...1`) %>%
                  clean_names() %>%
                  select(-yr) %>%
                  rename(atos_id = atos)%>%
                  #join atos geometry
                  left_join(atos_build1, by = "atos_id", relationship = "many-to-many") %>%
                  select(-geometry, -trend5yr)%>%
                  #convert density estimates to counts
                  mutate(n_indep = lin_dens*area_km2,
                         n_pup = n_indep*pupratio,
                         #set fields
                         atos_id = as.character(atos_id),
                         area_m2 = as.character(area_m2)
                         ) 

#check trends to see if the numbers make sense

summarized_data <- census_build1 %>%
  group_by(year) %>%
  summarize(total_n_indep = sum(n_indep))

base_theme <-  theme(axis.text=element_text(size=10, color = "black"),
                     axis.title=element_text(size=10,color = "black"),
                     plot.tag=element_text(size=9,color = "black"),
                     plot.title=element_text(size=12,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     #legend.key = element_rect(fill = "white"), # Set it to transparent
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=8,color = "black"),
                     legend.title=element_blank(),
                     #legend.key.height = unit(0.1, "cm"),
                     #legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=10, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())

# Plotting the summarized data
g <- ggplot(summarized_data, aes(x = year, y = total_n_indep)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "Year",
       y = "No. independents")+
  theme_bw() + base_theme

#ggsave(g, file = file.path("/Users/jossmith/Downloads/indep_plot.png"), width = 6,
 #      height = 5, units = "in")

################################################################################
#Step 3 - export

write.csv(census_build1, file = file.path(basedir, "processed/census/census_data_processed.csv"))












