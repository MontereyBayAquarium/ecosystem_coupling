#Joshua G. Smith; jossmith@mbayaq.org

#This processing script reads the raw rocky intertidal data obtained by request
#from the Multi-Agency Rocky Intertidal Network (MARINe)


rm(list=ls())

######
#required packages
require(librarian)
librarian::shelf(tidyverse, sf, janitor, here)

#set directories 
localdir <- here::here("output")

# Get rocky intertidal position data
meta_dat <- readxl::read_xlsx(file.path(localdir,"raw/rocky_intertidal/mytilus_pisaster_position_data.xlsx"),sheet = 1)
mus_pos_raw <- readxl::read_xlsx(file.path(localdir,"raw/rocky_intertidal/mytilus_pisaster_position_data.xlsx"),sheet = 2)
star_pos_raw <- readxl::read_xlsx(file.path(localdir,"raw/rocky_intertidal/mytilus_pisaster_position_data.xlsx"),sheet = 3)


mus_el <- readxl::read_xlsx(file.path(localdir,"raw/rocky_intertidal/mussel_elevation_data_20240909.xlsx"),sheet = 2)%>%
  janitor::clean_names() 

#read mussel size fq. 
#mus_size_orig <- readxl::read_xlsx(file.path(localdir,"raw/rocky_intertidal/mussel_size_fq.xlsx"),sheet = 1) #original file
mus_size_orig <- read_csv(file.path(localdir,"raw/rocky_intertidal/mussel_sizes_20240603.csv"))

#read mussel perc cov
mus_cov_orig <- readxl::read_xlsx(file.path(localdir,"raw/rocky_intertidal/mytilus_cov_pis_den_20240924.xlsx"),sheet = 2)


################################################################################
#Step 1 - filter to focal study area

mus_pos_build1 <- mus_pos_raw %>% filter(latitude >= 36.47986 & latitude <= 36.64640) %>%
  mutate(ssw_period = ifelse(year <=2013, "Before","After")) %>%
  #drop asilomar since there is only one year
  filter(!(intertidal_sitename %in% c("Asilomar","China Rocks")))

#check
unique(mus_pos_build1$intertidal_sitename)

mus_size_build1 <- mus_size_orig %>% 
  filter(marine_site_name %in% c("Point Lobos","Hopkins","Stillwater","Point Pinos"))%>%
  rename(year = marine_common_year)%>%
  mutate(ssw_period = ifelse(year <=2013, "Before","After")) 


mus_cov_period <- mus_cov_orig %>% filter(latitude >= 36.47986 & latitude <= 36.64640) %>%
  mutate(ssw_period = ifelse(year <=2013, "Before","After")) %>%
  #drop asilomar since there is only one year
  filter(!(marine_site_name %in% c("Asilomar","China Rocks"))) 


mus_elev <- mus_el %>% filter(latitude >= 36.47986 & latitude <= 36.64640) %>%
  mutate(ssw_period = ifelse(year <=2013, "Before","After")) %>%
  #drop asilomar since there is only one year
  filter(!(intertidal_sitename %in% c("Asilomar","China Rocks")))

################################################################################
#export
#save(mus_pos_build1, mus_size_build1, mus_cov_period, mus_elev, file = 
 #      file.path(localdir, "processed/rocky_intertidal/position_data.rdata")) #last write 26 Sept 2024



