#Joshua G. Smith; jossmith@mbayaq.org

#This processing script reads the raw rocky intertidal data obtained by request
#from the Multi-Agency Rocky Intertidal Network (MARINe)


rm(list=ls())

######
#required packages
librarian::shelf(tidyverse,sf, janitor)

#set directories 
basedir <- here::here("output")

#read data
meta_dat <- readxl::read_xlsx(file.path(basedir,"/raw/rocky_intertidal/mytilus_cov_pis_den.xlsx"),sheet = 1)
mus_orig <- readxl::read_xlsx(file.path(basedir,"/raw/rocky_intertidal/mytilus_cov_pis_den.xlsx"),sheet = 2)
pis_orig <- readxl::read_xlsx(file.path(basedir,"/raw/rocky_intertidal/mytilus_cov_pis_den.xlsx"),sheet = 3)