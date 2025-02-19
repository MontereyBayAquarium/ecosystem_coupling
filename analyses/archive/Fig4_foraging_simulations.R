

rm(list=ls())

################################################################################
#Prep workspace

#install.packages("dirichlet", repos="http://R-Forge.R-project.org")
librarian::shelf(dplyr, gtools, ggplot2, dirichlet, here, googledrive, ggsignif, rstan)

rstan::rstan_options(javascript=FALSE)

################################################################################
#Set directories and load data
datin <- here::here("output","processed")
figdir <- here::here("figures")

# Define the Google Drive file ID for sofa output
# This file exceeds 100MB GitHub limit so is stored externally
file_id <- "1yC8yrS69uOqbvlzQJakL0z-DC-yMhNel"

temp_file <- tempfile(fileext = ".rdata")
drive_download(as_id(file_id), path = temp_file, overwrite = TRUE)
SOFA_rslts <- mget(load(temp_file))

df_mus = read.csv(file.path(datin, "annual_urc_mus_data/mus_annual_cov.csv"))
df_urc = read.csv(file.path(datin, "annual_urc_mus_data/urchin_annual_density.csv"))

################################################################################
#Process data
df_urc$sp = NA #character(length = nrow(df_urc))
df_urc$sp[df_urc$species=="Strongylocentrotus purpuratus"] = "pur"
df_urc$sp[is.na( df_urc$sp)] = "red"
df_red = df_urc[df_urc$sp=="red",]
df_pur = df_urc[df_urc$sp=="pur",]
#
# NOTE: need to adjust purple urchin density to account for nutritional changes
#  (reduced gonad index) at higher densities: get "edible" urchin density 
#  (*** use Josh's data sets to better anchor this adjustment to data?)
# x = log(df_pur$mean_density_m2)
#df_pur$dens_adj =  exp(x - (x-min(x)+.01)*.37)
df_pur$lg_dns = log(df_pur$mean_density_m2)
df_red$lg_dns = log(df_red$mean_density_m2)
# plot(df_pur$year,df_pur$lg_dns); points(df_red$year,df_red$lg_dns+.6,col="red")
df_mus$lam = inv.logit(log(df_mus$mean_cov)-3)
df_mus$lg_dns = log(df_mus$mean_cov)
#
sumstats = SOFA_rslts$sumstats; mcmc = SOFA_rslts$mcmc; vn = SOFA_rslts$vn; 
Nsims = dim(mcmc)[1]
Years = as.integer(SOFA_rslts$Groupname)
N = length(Years)
Y0 = min(Years)-1
K = SOFA_rslts$Nptps
dfPtp = SOFA_rslts$dfPtp
ii = which(startsWith(vn,"ER["))
E_pr = as.matrix(mcmc[,ii])
################################################################################
#Run simulations

# Assume constant E intake per prey type
# Use moving average to smooth eta... 3-year MAvg span?  
# Alternative scenario... eta for urc and mus after yr 6 are avg of years 1:6,
#   and the difference is made up by increasing eta[K] ("all other" prey type)
# Ebar is calculated as eta * E_pr
eta_sims = array(0,dim = c(Nsims, K, N))
eta_sims_alt = array(0,dim = c(Nsims, K, N))
ones = array(1,dim=Nsims)
tau = numeric(length = N)
eta_obs = matrix(0,nrow = N, ncol = K)
eta_lo = matrix(0,nrow = N, ncol = K)
eta_hi = matrix(0,nrow = N, ncol = K)
eta_alt = matrix(0,nrow = N, ncol = K)
eta_alt_lo = matrix(0,nrow = N, ncol = K)
eta_alt_hi = matrix(0,nrow = N, ncol = K)
Ebar = numeric(length = N)
Ebar_lo = numeric(length = N)
Ebar_hi = numeric(length = N)
Ebar_alt = numeric(length = N)
Ebar_alt_lo = numeric(length = N)
Ebar_alt_hi = numeric(length = N)
urcP_obs = numeric(length = N); urcR_obs = numeric(length = N); mus_obs= numeric(length = N)
DP_flag = numeric(length = N); DR_flag = numeric(length = N); DM_flag= numeric(length = N)
for(t in 1:N){
  # Invert data
  yr = Years[t]
  ii = which(df_pur$year==yr)
  if(length(ii)>0){DP_flag[t] = 1; urcP_obs[t] = df_pur$lg_dns[ii] }
  ii = which(df_red$year==yr)
  if(length(ii)>0){DR_flag[t] = 1; urcR_obs[t] = df_red$lg_dns[ii] }
  ii = which(df_mus$year==yr)
  if(length(ii)>0){DM_flag[t] = 1; mus_obs[t] = df_mus$lg_dns[ii] }
  # Foraging data
  ii = which(startsWith(vn,paste0("etaG[",t,",")))
  eta_sims[,,t] = mcmc[,ii]
  ft = dirichlet::fit.dirichlet(mcmc[,ii])
  tau[t] = ft$most.likely.k
  if (t == 1 ){
    ii = which(startsWith(vn,paste0("etaG[",t+1,",")))
    eta_sims[,,t] = eta_sims[,,t] + mcmc[,ii]
    eta_sims[,,t] = eta_sims[,,t] / 2
    eta_sims[,,t] = sweep(eta_sims[,,t],1,apply(eta_sims[,,t],1,sum),FUN = '/')
  }else if(t == N){
    ii = which(startsWith(vn,paste0("etaG[",t-1,",")))
    eta_sims[,,t] = eta_sims[,,t] + mcmc[,ii]
    eta_sims[,,t] = eta_sims[,,t] / 2
    eta_sims[,,t] = sweep(eta_sims[,,t],1,apply(eta_sims[,,t],1,sum),FUN = '/')
  }else{
    ii = which(startsWith(vn,paste0("etaG[",t+1,",")))
    eta_sims[,,t] = eta_sims[,,t] + mcmc[,ii]
    ii = which(startsWith(vn,paste0("etaG[",t-1,",")))
    eta_sims[,,t] = eta_sims[,,t] + mcmc[,ii]
    eta_sims[,,t] = eta_sims[,,t] / 3
    eta_sims[,,t] = sweep(eta_sims[,,t],1,apply(eta_sims[,,t],1,sum),FUN = '/')
  }
  eta_sims_alt[,,t] = eta_sims[,,t]
  if(t > 6){
    eta_sims_alt[,1,t] = apply(eta_sims[,1,1:6],1,mean)
    eta_sims_alt[,2,t] = apply(eta_sims[,2,1:6],1,mean)
    dif = ones - apply(eta_sims_alt[,,t],1,sum)
    eta_sims_alt[,K,t] = eta_sims_alt[,K,t] + dif
  }
  eta_obs[t,] = apply(eta_sims[,,t],2,mean)
  eta_lo[t,] = apply(eta_sims[,,t],2,quantile,probs=0.1)
  eta_hi[t,] = apply(eta_sims[,,t],2,quantile,probs=0.9)
  eta_alt[t,] = apply(eta_sims_alt[,,t],2,mean)
  eta_alt_lo[t,] = apply(eta_sims_alt[,,t],2,quantile,probs=0.1)
  eta_alt_hi[t,] = apply(eta_sims_alt[,,t],2,quantile,probs=0.9)
  Ebar[t] = mean(rowSums(eta_sims[,,t] * E_pr))  
  Ebar_lo[t] = quantile(rowSums(eta_sims[,,t] * E_pr),0.1)  
  Ebar_hi[t] = quantile(rowSums(eta_sims[,,t] * E_pr),0.9)  
  Ebar_alt[t] = mean(rowSums(eta_sims_alt[,,t] * E_pr))  
  Ebar_alt_lo[t] = quantile(rowSums(eta_sims_alt[,,t] * E_pr),0.1)  
  Ebar_alt_hi[t] = quantile(rowSums(eta_sims_alt[,,t] * E_pr),0.9)  
}

df_Erate_est = data.frame(Year=Years,Scenario=rep("Actual",N),
                          Erate=Ebar,
                          Erate_lo=Ebar_lo,Erate_hi=Ebar_hi)
df_Erate_est = rbind(df_Erate_est, 
                     data.frame(Year=Years,Scenario=rep("Alternate",N),
                                Erate=Ebar_alt,
                                Erate_lo=Ebar_alt_lo,Erate_hi=Ebar_alt_hi) )
#
plt_Erate = ggplot(df_Erate_est,aes(x=Year,y=Erate,group=Scenario)) +
  geom_ribbon(aes(ymin=Erate_lo,ymax=Erate_hi,fill=Scenario),alpha=.25) +
  geom_line(aes(color=Scenario)) + 
  labs(y="Energy Intake Rate") +
  ggtitle("Effect of urchin and mussel increase on Energy Intake") +
  theme_classic()
print(plt_Erate)

df_diet = data.frame(Year=Years,Prey_type = rep("urchins",N),
                     Scenario=rep("Actual",N),
                     eta = eta_obs[,1],eta_lo=eta_lo[,1],eta_hi=eta_hi[,1])
df_diet = rbind(df_diet,
                data.frame(Year=Years,Prey_type = rep("mussels",N),
                           Scenario=rep("Actual",N),
                           eta = eta_obs[,2],eta_lo=eta_lo[,2],eta_hi=eta_hi[,2]) )

df_diet = rbind(df_diet,
                data.frame(Year=Years,Prey_type = rep("urchins",N),
                           Scenario=rep("Alternate",N),
                           eta = eta_alt[,1],eta_lo=eta_alt_lo[,1],eta_hi=eta_alt_hi[,1]) )

df_diet = rbind(df_diet,
                data.frame(Year=Years,Prey_type = rep("mussels",N),
                           Scenario=rep("Alternate",N),
                           eta = eta_alt[,2],eta_lo=eta_alt_lo[,2],eta_hi=eta_alt_hi[,2]) )

ii = which(df_diet$Scenario=="Actual")
plt_effort = ggplot(data = df_diet[ii,],aes(x=Year,group=Prey_type)) +
  geom_ribbon(aes(ymin=eta_lo,ymax=eta_hi,fill=Prey_type),alpha=0.25) +
  geom_line(aes(y=eta, color=Prey_type)) +
  labs(x = "Year", y = "Proportion of foraging effort") +
  theme_classic() + #theme(legend.position = "inside") +
  scale_x_continuous(breaks = as.integer(Years),
                     guide = guide_axis(angle = 45))
print(plt_effort)

################################################################################
#plot Figure 5 with matched theme


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

# Determine the common range of years for both plots
year_range <- range(df_Erate_est$Year, df_diet$Year)


p1 <- ggplot(df_Erate_est, aes(x = Year, y = Erate, group = Scenario)) +
  geom_ribbon(aes(ymin = Erate_lo, ymax = Erate_hi, fill = Scenario), alpha = 0.25) +
  geom_line(aes(color = Scenario)) + 
  labs(y = "Energy intake rate", x = "Year", tag = "B") +
  ggtitle("Effect of urchin and mussel increase on energy intake") +
  theme_bw() + base_theme +
  scale_fill_manual(values = c("Actual" = "indianred", "Alternate" = "navyblue")) +
  scale_color_manual(values = c("Actual" = "indianred", "Alternate" = "navyblue")) +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  annotate(geom = "text", label = "SSW", x = 2010.5, y = 9.1, size = 2.5) +
  annotate("segment", x = 2010.8, y = 9.0, xend = 2012.7, yend = 8.6,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  scale_x_continuous(limits = year_range, breaks = seq(min(year_range), max(year_range), by = 2)) 
 
p1

# Second plot (p2)
p2 <- ggplot(data = df_diet[ii,], aes(x = Year, group = Prey_type)) +
  geom_ribbon(aes(ymin = eta_lo, ymax = eta_hi, fill = Prey_type), alpha = 0.25) +
  geom_line(aes(y = eta, color = Prey_type)) +
  labs(x = "", y = "Proportion of foraging effort", tag = "A") +
  ggtitle("Foraging effort over time") +
  theme_classic() +
  scale_x_continuous(limits = year_range, breaks = seq(min(year_range), max(year_range), by = 2),  
                     guide = guide_axis(angle = 45)) +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  annotate(geom = "text", label = "SSW", x = 2010.5, y = 0.39, size = 2.5) +
  annotate("segment", x = 2010.8, y = 0.375, xend = 2012.7, yend = 0.32,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  theme_bw() + base_theme +
  scale_fill_manual(values = c("mussels" = "orange", "urchins" = "purple"),
                    labels = c("mussels" = "Mussels", "urchins" = "Sea urchins")) +
  scale_color_manual(values = c("mussels" = "orange", "urchins" = "purple"),
                     labels = c("mussels" = "Mussels", "urchins" = "Sea urchins"))+
  theme(axis.text.x = element_blank())  
p2

p <- ggpubr::ggarrange(p2, p1, nrow = 2, align = "v")
p

#ggsave(p, filename = file.path(figdir, "Fig4_lv_model.png"), 
 #      width =5, height = 5, units = "in", dpi = 600, bg = "white") #last write 26 Sept 2024


################################################################################
#statistical comparisons

df_diet <- df_diet %>%
  data.frame()%>%
  mutate(Prey_type = recode(Prey_type,
                            "mussels" = "Mussels",
                            "urchins" = "Sea urchins"),
         Period = ifelse(Year > 2013, "Post-2013","Pre-2013"))

# Create the boxplot 
b1 <- ggplot(df_diet, aes(x = Prey_type, y = eta, fill = Period)) +
  geom_boxplot(position = position_dodge(width = 0.75)) + 
  labs(y = "Proportion of foraging effort", x = "Prey type", fill = "Period") +
  scale_fill_manual(values = c("Pre-2013" = "navyblue", "Post-2013" = "indianred")) +  
  theme(legend.position = "right",  
        plot.title = element_blank()) +
  theme_bw() + base_theme
b1

#ggsave(b1, filename = file.path(figdir, "FigS2_boxplots.png"), 
 #      width =5, height = 5, units = "in", dpi = 600, bg = "white") #last write 26 Sept 2024

# Calculate pre- and post-2013 averages for each species
pre_post_averages <- df_diet %>%
  group_by(Prey_type, Period) %>%
  summarise(
    mean_eta = mean(eta, na.rm = TRUE),  # Calculate the mean proportion of foraging effort
    se_eta = sd(eta, na.rm = TRUE) / sqrt(n()),  # Calculate standard error
    n = n(),  # Count of observations
    .groups = 'drop'
  )

# Display the pre- and post-2013 averages summary
print(pre_post_averages)

# Calculate fold increase for each prey type
fold_increase <- pre_post_averages %>%
  group_by(Prey_type) %>%
  summarise(
    fold_increase = mean_eta[Period == "Post-2013"] / mean_eta[Period == "Pre-2013"],
    .groups = 'drop'
  )

# Display the fold increase results
print(fold_increase)


# Perform t-tests to statistically compare pre- vs. post-2013 proportions for each prey type using raw data
t_test_results_raw <- df_diet %>%
  group_by(Prey_type) %>%
  summarise(
    t_test = list(t.test(eta ~ Period)),
    .groups = 'drop'
  )

# Extract p-values and other test details
t_test_summary_raw <- t_test_results_raw %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    statistic = sapply(t_test, function(x) x$statistic),
    conf_int = sapply(t_test, function(x) paste(round(x$conf.int, 3), collapse = " to ")),
    mean_diff = sapply(t_test, function(x) diff(x$estimate))
  ) %>%
  select(Prey_type, mean_diff, statistic, p_value, conf_int)

# Display the t-test summary results
print(t_test_summary_raw)


