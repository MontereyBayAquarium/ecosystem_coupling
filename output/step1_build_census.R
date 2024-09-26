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
census_orig <- read_csv(file.path(basedir,"raw/census_data/ATOS_362-449_modifed_census_results (1985-2023).csv")) %>% clean_names()
atos_orig <- st_read(file.path(basedir,"raw/gis_data/atos/ATOS_polygon_teale83.shp"))%>% clean_names()
atos_df <- load(file.path(basedir, "raw/census_data/ATOS_DF.rdata")) 


################################################################################
#Step 1 - process ATOS geometries

atos_areas <- ATOS_DF %>% clean_names() %>% rename(
  atos_id = atos,
  area_bay = bay,
  area_near = nearshore,
  area_off = offshore,
  area_wayoff = wayoffsh
)

atos_build1 <- atos_orig %>%
               #merge offshore / onshore geometrues
               group_by(atos_id) %>%
               summarize(geometry = st_union(geometry))%>%
               mutate(area_m2 = st_area(geometry),
                      area_km2 = area_m2 / 1e6,
                      #format
                      area_km2 = as.numeric(sprintf("%.7f", area_km2))
                      ) %>%
                clean_names() %>%
              left_join(atos_areas, by = "atos_id")


################################################################################
#Step 2 - inspect and process CENSUS data

census_build1 <- census_orig %>%
                  #drop index column
                  #select(-`...1`) %>%
                  clean_names() %>%
                  select(-yr) %>%
                  rename(atos_id = atos)%>%
                  #join atos geometry
                  left_join(atos_build1, by = "atos_id", relationship = "many-to-many") %>%
                  select(-geometry, -trend5yr)%>%
                  #convert density estimates to counts
                    mutate(total_counts = ((dns_bay * area_bay) + (dns_near * area_near) +
                                             (dns_off * area_off) + (dns_wayoff * area_wayoff)), #calculate total counts
                           pup_prop = pupratio/(1+pupratio), #convert pup:indep ratio into a proportion of pups out of the total
                           pup_counts = total_counts * pup_prop, #determine pup counts
                           n_indep = total_counts - pup_counts, #calculate number of independents
                           atos_id = as.character(atos_id),
                           area_m2 = as.character(area_m2)
                    )


summarized_data <- census_build1 %>%
  group_by(year) %>%
  summarize(total_n_indep = sum(n_indep),
            total_pup = sum(pup_counts))

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

g <- ggplot(summarized_data, aes(x = year)) +
  geom_line(aes(y = total_n_indep, color = "Independents")) +
  geom_point(aes(y = total_n_indep, color = "Independents")) +
  geom_line(aes(y = total_pup, color = "Pups")) +
  geom_point(aes(y = total_pup, color = "Pups")) +
  labs(title = "",
       x = "Year",
       y = "Counts along the Monterey Peninsula")+
  theme_bw() + base_theme +
  scale_color_manual(values = c(Independents = "forestgreen", Pups = "orange"))

g


ggsave(g, file = file.path("/Users/jossmith/Downloads/indep_plot.png"), width = 6,
       height = 5, units = "in")

################################################################################
#Step 3 - export

write.csv(census_build1, file = file.path(basedir, "processed/census/census_data_processed.csv")) #last write 26 September 2024












