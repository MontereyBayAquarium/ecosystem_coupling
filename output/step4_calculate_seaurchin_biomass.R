

################################################################################
# About
# This is the processing script used to derive purple sea urchin length-weight
# conversion parameters for soft tissue biomass (exclusing the hard shell / test). 
#The conversion parameters were obtained using 111 purple sea urchins (strongylocentrotus purpuratus) collected 
#from a kelp forest in May, 2024. Sea urchins were brought to the lab, measured, weighed, and dissected. 


################################################################################

rm(list=ls())

librarian::shelf(tidyverse,here, janitor)

#set paths
datin <- here::here("output","raw")
figdir <- here::here("figures")

#read urchin data
urch_dat_orig <- read_csv(file.path(datin,"urchin_data/purple_urchin_biomass.csv")) %>% clean_names() %>%
                    #drop missing data
                    drop_na(test_diameter_mm) %>%
                    #filter
                    filter(test_diameter_mm < 90)


################################################################################
# derive parameters

#take a look
plot(urch_dat_orig$test_diameter_mm, urch_dat_orig$soft_tissue_mass_g)

# Initial estimates based on data
a_init <- -20
b_init <- 10
c_init <- 0.03

# Fit the biomass_fun model to the data with initial estimates
set.seed(1985)
fit <- nls(soft_tissue_mass_g ~ a + b * exp(c * test_diameter_mm), 
           data = urch_dat_orig,
           start = list(a = a_init, b = b_init, c = c_init))


# Extract the estimated parameters
a_est <- coef(fit)["a"]
b_est <- coef(fit)["b"]
c_est <- coef(fit)["c"]

# Print the estimated parameters
cat("Estimated Parameters:\n")
cat("a:", a_est, "\n")
cat("b:", b_est, "\n")
cat("c:", c_est, "\n")


################################################################################
#determine fit
# Calculate the predicted values from the model
predicted_values <- fitted(fit)

# Calculate the residuals
residuals <- urch_dat_orig$soft_tissue_mass_g - predicted_values

# Calculate the RSS (Residual Sum of Squares)
rss <- sum(residuals^2)

# Calculate the TSS (Total Sum of Squares)
mean_soft_mass <- mean(urch_dat_orig$soft_tissue_mass_g)
tss <- sum((urch_dat_orig$soft_tissue_mass_g - mean_soft_mass)^2)

# Calculate R-squared
r_squared <- 1 - (rss / tss)

# Print the R-squared value
cat("R-squared:", r_squared, "\n")


################################################################################
#plot

base_theme <-  theme(axis.text=element_text(size=12, color = "black"),
                     axis.title=element_text(size=12,color = "black"),
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


#generate equation expression
equation_text <- substitute(italic(y) == a + b %*% e^(c * italic(x)) * "," ~ italic(R)^2 ~ "=" ~ r2,
                            list(a = round(a_est, 2), 
                                 b = round(b_est, 2), 
                                 c = round(c_est, 2),
                                 r2 = round(r_squared, 2)))

# Create the plot
g <- ggplot(urch_dat_orig, aes(x = test_diameter_mm, y = soft_tissue_mass_g)) +
  geom_point() +
  geom_line(aes(y = a_est + b_est * exp(c_est * test_diameter_mm)), color = "purple", size = 1) +
  labs(x = "Test Diameter (mm)", y = "Consumable wet biomass (g)") +
  theme_bw() + 
  base_theme +
  annotate("text", x = min(urch_dat_orig$test_diameter_mm), y = max(urch_dat_orig$soft_tissue_mass_g), 
           label = as.expression(equation_text), hjust = 0, vjust = 1, size = 5, color = "black")

ggsave(g, file = file.path(figdir, "purple_urchin_biomass.png"), width = 6,
       height = 5, units = "in")


