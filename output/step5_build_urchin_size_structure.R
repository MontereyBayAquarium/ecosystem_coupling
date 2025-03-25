#Sea urchin abundance data processing
#Joshua G. Smith
#August 1, 2024

rm(list=ls())

require(librarian)
librarian::shelf(tidyverse, here)


###subtidal monitoring data were accessed at
###https://opc.dataone.org/view/doi%3A10.25494%2FP6%2FMLPA_kelpforest.7

##require:

#MLPA_kelpforest_swath.6.csv
#MLPA_kelpforest_site_table.6.csv
#MLPA_kelpforest_taxon_table.6.csv

################################################################################
#set directories and load data
basedir <- here::here("output","raw","urchin_data")
output <- here::here("output","processed","urchin_data")

kelp_swath_raw <- read.csv(file.path(basedir, "MLPA_kelpforest_swath.6.csv")) %>%
  janitor::clean_names()

kelp_taxon <- read.csv(file.path(basedir, "MLPA_kelpforest_taxon_table.6.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(classcode, species_definition) %>%
  distinct()

site_table <- read.csv(file.path(basedir, "MLPA_kelpforest_site_table.6.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(site, latitude, longitude, ca_mpa_name_short, mpa_class=site_designation, mpa_designation=site_status, 
                baseline_region)%>%
  distinct() #remove duplicates

################################################################################
#process kelp swath

#select vars and clean
kelp_swath_build1 <- kelp_swath_raw %>%
  #drop juvenile sea urchins, since these were not consistently recorded 
  dplyr::filter(classcode == "STRPURAD" |
                    classcode == "MESFRAAD") %>%
  #join species names by class code 
  left_join(., kelp_taxon, by="classcode") %>%
  #add affiliated_mpa and lat/long
  left_join(., site_table, by="site") %>%
  #filter to focal study area
  filter(latitude >= 36.47986 & latitude <= 36.64640) %>%
  ungroup() %>%
  dplyr::select(year, baseline_region, site,latitude, longitude, affiliated_mpa=ca_mpa_name_short,mpa_class, mpa_designation, zone, transect, species = species_definition,
                count, size_cm = size) %>%
  mutate(mpa_designation = ifelse(mpa_designation=="reference","ref",mpa_class)) %>%
  #note: size fq did not begin until 2011
  filter(year > 2010) %>% 
  drop_na(size_cm) %>%
  #calculate size frequency at site-year level
  group_by(year, latitude, longitude, site, species, size_cm) %>%
  summarize(total_counts = sum(count)) %>%
  #set column types
  mutate(site = as.factor(site),
         species = as.factor(species))



################################################################################
#export summarized data

#write_csv(kelp_swath_build1, file = file.path(output, "urchin_sizefq.csv")) #last write 01 August 2024









