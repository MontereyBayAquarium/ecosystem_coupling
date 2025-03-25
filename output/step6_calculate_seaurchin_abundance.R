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
output <- here::here("output","processed","annual_urc_mus_data")

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

mus_dat <- load(file = file.path("output/processed/rocky_intertidal/pisaster_mytilus_processed.rdata"))


################################################################################
#process urchin data

#select vars and clean
kelp_swath_build1 <- kelp_swath_raw %>%
  #drop juvenile sea urchins, since these were not consistently recorded 
  dplyr::filter(!(classcode == "STRPURREC" |
                    classcode == "MESFRAREC"))%>%
  dplyr::select(year, site, zone, transect, classcode, count, size)%>%
  group_by(year, site, zone, transect, classcode)%>%
  dplyr::summarize(total_count = sum(count), #counts in raw data are grouped by size class. Take summary across all sizes
                   total_size = sum(size)) %>% #this is for kelp stipe counts only 
  mutate(total_count = ifelse(classcode == "MACPYRAD",total_size, total_count)) %>% #replace num plants with total stipes
  dplyr::select(!(total_size)) %>%
  #select urchins only
  filter(classcode == "STRPURAD" | classcode == "MESFRAAD") 


#join species names by class code 
kelp_swath_build2 <- left_join(kelp_swath_build1, kelp_taxon, by="classcode")

kelp_swath_build3 <- kelp_swath_build2 %>%
  dplyr::select(year, site, zone, transect, species=species_definition, total_count)


#add affiliated_mpa and lat/long
kelp_swath_build4 <- left_join(kelp_swath_build3, site_table, by="site")

kelp_swath_build5 <- kelp_swath_build4 %>%
  ungroup() %>%
  dplyr::select(year, baseline_region, site,latitude, longitude, affiliated_mpa=ca_mpa_name_short,mpa_class, mpa_designation, zone, transect, species,
                total_count) %>%
  mutate(mpa_designation = ifelse(mpa_designation=="reference","ref",mpa_class)) %>%
  #make sure species are not duplicates by summarizing at the transect level (total counts)
  group_by(year, baseline_region, latitude, longitude, site, affiliated_mpa, mpa_class, mpa_designation, zone,  transect, species) %>%
  dplyr::summarise(counts = sum(total_count))


#filter to sites within focal stusy area
kelp_swath_build6 <- kelp_swath_build5 %>%
                      filter(latitude >= 36.47986 & latitude <= 36.64640)


#set column types
kelp_swath_build7 <- kelp_swath_build6 %>%
                      data.frame() %>%
                      mutate(year = as.numeric(year),
                             baseline_region = factor(baseline_region),
                             latitude = as.numeric(latitude),
                             longitude = as.numeric(longitude),
                             site = as.factor(site),
                             affiliated_mpa = as.factor(affiliated_mpa),
                             mpa_class = as.factor(mpa_class),
                             mpa_designation = as.factor(mpa_designation),
                             zone = as.factor(zone),
                             transect = as.factor(transect),
                             species = as.factor(species),
                             counts = as.numeric(counts))


################################################################################
#aggregate data to site-level means
kelp_swath_build8 <- kelp_swath_build7 %>%
                     #calcualte density
                     mutate(density_m2 = counts / 60) %>%
                     group_by(year, baseline_region, site, species) %>%
                     dplyr::summarize(
                               density_m2_mean = mean(density_m2, na.rm = TRUE),
                               #density_m2_sd = sd(density_m2, na.rm = TRUE)
                               ) %>%
                     #drop 1999
                     filter(year > 1999)

################################################################################
#aggregate data to annual means with sd calculate among sites

kelp_swath_build9 <- kelp_swath_build8 %>%
  #calcualte density
  group_by(year, baseline_region, species) %>%
  dplyr::summarize(mean_density_m2 = mean(density_m2_mean, na.rm = TRUE),
                   sd_density_m2 = sd(density_m2_mean, na.rm=TRUE)) %>%
  #drop 1999
  filter(year > 1999) %>%
  data.frame()%>%
  dplyr::select(-baseline_region)

################################################################################
#check trends

ggplot(kelp_swath_build9, aes(x = year, y = mean_density_m2, color = species, group = species)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_density_m2 - sd_density_m2, ymax = mean_density_m2 + sd_density_m2), width = 0.2) +
  geom_line() +
  labs(
    title = "Overall Mean Density Across All Sites for Each Year by Species",
    x = "Year",
    y = "Mean Density (m^2)",
    color = "Species"
  ) +
  theme_minimal()

################################################################################
#process mussel data

mus_build <- mus_build1 %>%
              data.frame()%>%
              group_by(year, species = species_lump) %>%
              summarize(mean_cov = mean(percent_cover, na.rm=TRUE),
                        sd_cov = sd(percent_cover, na.rm=TRUE))

################################################################################
#check trends

ggplot(mus_build, aes(x = year, y = mean_cov)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_cov - sd_cov, ymax = mean_cov + sd_cov), width = 0.2) +
  geom_line() +
  labs(
    title = "Overall Mean Density Across All Sites for Each Year by Species",
    x = "Year",
    y = "Mean Density (m^2)",
    color = "Species"
  ) +
  theme_minimal()

################################################################################
#export summarized data

#write_csv(kelp_swath_build8, file = file.path(output, "urchin_site_level_abundances.csv")) #last write 01 August 2024
#write_csv(kelp_swath_build9, file = file.path(output, "urchin_annual_density.csv")) #last write 09 August 2024
#write_csv(mus_build, file = file.path(output, "mus_annual_cov.csv")) #last write 09 August 2024







