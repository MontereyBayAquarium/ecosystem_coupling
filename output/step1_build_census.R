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
                  #convert density estimates to counts
                  mutate(num_indep = lin_dens*area_orig,
                         num_pup = num_indep*pupratio
                         )

#check trends to see if the numbers make sense

summarized_data <- census_build1 %>%
  group_by(year) %>%
  summarize(total_num_indep = sum(num_indep))

# Plotting the summarized data
ggplot(summarized_data, aes(x = year, y = total_num_indep)) +
  geom_line() +
  geom_point() +
  labs(title = "Total num_indep by Year",
       x = "Year",
       y = "Total num_indep")















