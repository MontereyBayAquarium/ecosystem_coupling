#Joshua G. Smith; jossmith@mbayaq.org

#This processing script reads the raw rocky intertidal data obtained by request
#from the Multi-Agency Rocky Intertidal Network (MARINe)


rm(list=ls())

######
#required packages
require(librarian)
librarian::shelf(tidyverse,sf, janitor, here)

#set directories 
basedir <- here::here("output")

#read data
meta_dat <- readxl::read_xlsx(file.path(basedir,"/raw/rocky_intertidal/mytilus_cov_pis_den_20240924.xlsx"),sheet = 1) 
mus_orig <- readxl::read_xlsx(file.path(basedir,"/raw/rocky_intertidal/mytilus_cov_pis_den_20240924.xlsx"),sheet = 2) 
pis_orig <- readxl::read_xlsx(file.path(basedir,"/raw/rocky_intertidal/mytilus_cov_pis_den_20240924.xlsx"),sheet = 3)


################################################################################
##Step 1 - filter rocky intertidal to focal study area

mus_build1 <- mus_orig %>% filter(latitude >= 36.47986 & latitude <= 36.64640) %>%
  #add sea star wasting period
  mutate(ssw_period = ifelse(year <=2013, "before","after")) %>%
  #drop asilomar since there is only one year
  filter(!(marine_site_name %in% c("Asilomar","China Rocks"))) 

#check
unique(mus_build1$marine_site_name)

pis_build1 <- pis_orig %>% filter(latitude >= 36.47986 & latitude <= 36.64640) %>%
  #add sea star wasting period
  mutate(ssw_period = ifelse(year <=2013, "before","after")) %>%
  #drop asilomar since there is only one year
  filter(!(marine_site_name %in% c("Asilomar","China Rocks"))) 


################################################################################
##Step 2 - save


#save(mus_build1, pis_build1, file = file.path(basedir, "processed/rocky_intertidal/pisaster_mytilus_processed.rdata")) #last write 26 Sept 2024






