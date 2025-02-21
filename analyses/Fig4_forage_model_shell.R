
rm(list=ls())

################################################################################
#Prep workspace

#install.packages("dirichlet", repos="http://R-Forge.R-project.org")
librarian::shelf(dplyr, parallel, cmdstanr, posterior, gtools, ggplot2, bayesplot,
                 readxl, RColorBrewer, dirichlet, googledrive)

rstan::rstan_options(javascript=FALSE)

################################################################################
#Set directories and load data
datin <- here::here("output","processed")
figdir <- here::here("figures")
# Load spreadsheets with summary statistics from SOFA model fit
df_Prey = read_excel(file.path(datin,"/sofa/Preytypes.xlsx"))
df_SOFA_rslts = read_excel(file.path(datin,"/sofa/SOFA_posterior_sums.xlsx"))

################################################################################
# Process data 
N_base = c(3,6) 
Years = seq(min(df_SOFA_rslts$Year,na.rm = T),max(df_SOFA_rslts$Year,na.rm = T))
N = length(Years)
Y0 = min(Years)-1
K = nrow(df_Prey)
PreyTypes = df_Prey$PreyType
log_G_pri = array(0, dim=c(2,K))
log_G_pri[1,1:K] = df_Prey$log_G_mn
log_G_pri[2,1:K] = df_Prey$log_G_sd
mu_obs = matrix(df_SOFA_rslts$Value[df_SOFA_rslts$Param=="mu_est"],nrow = N, ncol=K, byrow = F)
sig_mu = df_SOFA_rslts$Value[df_SOFA_rslts$Param=="sig_mu"]
V_mu = sig_mu^2
eta_obs = matrix(df_SOFA_rslts$Value[df_SOFA_rslts$Param=="eta"],nrow = N, ncol=K, byrow = F)
tau =  df_SOFA_rslts$Value[df_SOFA_rslts$Param=="tau"]
E_obs = exp( colMeans(mu_obs) + V_mu/2)
sd_E = sqrt( (exp(V_mu) - 1)*exp(2*colMeans(mu_obs) + V_mu) )
#
# Set up Stan fitting --------------------------------------
nsamples = 3000
nburnin = 2000
cores = detectCores()
ncore = max(3,min(5,cores-2))
Niter = round(nsamples/ncore)
fitmodel = here::here("analyses","Foraging_fit.stan")
#
stan.data = list(N=N,K=K,eta_obs=eta_obs,tau=tau,mu_obs=mu_obs,sig_mu=sig_mu,
                 N_base=N_base,log_G_pri=log_G_pri)
#
parms = c("sig_L","sig_E","Ebar","Ebar_alt","pi","delta","mu",
          "U_prd","M_prd","U_prd_alt","M_prd_alt")
#
mod <- cmdstan_model(fitmodel)

init_fun <- function() {list(sig_E=runif(K, .04, .06), 
                             sig_L=runif(3, .2, .25), 
                             mu = runif(K, .9*colMeans(mu_obs), 1.1*colMeans(mu_obs)) ) }
#
suppressMessages(
  suppressWarnings ( 
    fit <- mod$sample(
      data = stan.data,
      init <- init_fun ,
      seed = 123,
      chains = ncore,
      parallel_chains = ncore,
      refresh = 100,
      # thin=thinval,
      iter_warmup = nburnin,
      iter_sampling = Niter )
  )
)
source(here::here("analyses","cmdstan_sumstats.r"))
#

################################################################################
#Plot Figure 4 - proportion foraging effort

ii_urchins = which(startsWith(rownames(sumstats), "pi[") & endsWith(rownames(sumstats), ",1]"))
ii_mussels = which(startsWith(rownames(sumstats), "pi[") & endsWith(rownames(sumstats), ",2]"))

df_diet = data.frame(
  Year = Years, 
  Prey_type = rep("urchins", N),
  Scenario = rep("Actual", N),
  eta = sumstats$mean[ii_urchins],  # Mean value from sumstats
  eta_lo = sumstats$`5%`[ii_urchins],  
  eta_hi = sumstats$`95%`[ii_urchins]  
)

df_diet = rbind(df_diet, 
                data.frame(
                  Year = Years, 
                  Prey_type = rep("mussels", N),
                  Scenario = rep("Actual", N),
                  eta = sumstats$mean[ii_mussels],  
                  eta_lo = sumstats$`5%`[ii_mussels], 
                  eta_hi = sumstats$`95%`[ii_mussels]
                ))


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

# Define the common range of years
year_range <- range(df_diet$Year)


# Generate the p2 plot
p <- ggplot(data = df_diet, aes(x = Year, group = Prey_type)) +
  geom_ribbon(aes(ymin = eta_lo, ymax = eta_hi, fill = Prey_type), alpha = 0.25) +
  geom_line(aes(y = eta, color = Prey_type)) +
  labs(x = "", y = "Proportion of foraging effort") +
  ggtitle("Foraging effort over time") +
  theme_classic() +
  scale_x_continuous(limits = year_range, breaks = seq(min(year_range), max(year_range), by = 2),  
                     guide = guide_axis(angle = 45)) +
  geom_vline(xintercept = 2013, linetype = "dotted", linewidth = 0.6) +
  annotate(geom = "text", label = "SSW", x = 2010.5, y = 0.39, size = 2.5) +
  annotate("segment", x = 2010.8, y = 0.375, xend = 2012.7, yend = 0.32,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  theme_bw() + base_theme +
  scale_fill_manual(values = c("mussels" = "orange", "urchins" = "purple"),
                    labels = c("mussels" = "Mussels", "urchins" = "Sea urchins")) +
  scale_color_manual(values = c("mussels" = "orange", "urchins" = "purple"),
                     labels = c("mussels" = "Mussels", "urchins" = "Sea urchins")) 

p

#ggsave(p, filename = file.path(figdir, "Fig4_proportion_effort.png"), 
 #      width =5, height = 3, units = "in", dpi = 600, bg = "white")

################################################################################
#Plot posteriors and diagnostics

color_scheme_set("mix-viridis-orange-purple")
mcmc_trace(mcmc_array,regex_pars=("sig_L"))
mcmc_trace(mcmc_array,regex_pars=("sig_E"))
# mcmc_trace(mcmc_array,pars="alpha")
# mcmc_trace(mcmc_array,regex_pars="delta")
#
color_scheme_set("blue")
mcmc_areas(mcmc, pars= paste0("Ebar[",seq(1,N),"]"),
           area_method="equal height",
           prob = 0.8) + scale_y_discrete(labels=as.character(Years)) + labs(x="Energy gain")

Preynames = df_Prey$PreyType[1:K]
plt_trends = list()
df_prey_est = data.frame(Prey=character(),
                         Year=numeric(),
                         Dns=numeric(),
                         Dns_lo=numeric(),
                         Dns_hi=numeric() )


#
for(i in 1:(K-1)){
  ii = which(startsWith(vn,"delta[") & endsWith(vn,paste0(",",i,"]")))
  Dns = sumstats$mean[ii]
  Dns_lo = sumstats$`5%`[ii] 
  Dns_hi = sumstats$`95%`[ii]
  df_pr_est = data.frame(Prey = rep(Preynames[i],N), Year=Years,
                         Dns=Dns,Dns_lo=Dns_lo,Dns_hi=Dns_hi)
  plt_trends[[i]] = ggplot(filter(df_pr_est,Year>2007),aes(x=Year,y=Dns)) +
    geom_ribbon(aes(ymin=Dns_lo,ymax=Dns_hi),alpha=.25) +
    # ylab(expression(paste("logit (", lambda,")"))) +
    ylab("Density") +
    ggtitle(paste0("Trends in ",Preynames[i])) +
    geom_line() + theme_classic()
  
  df_prey_est = rbind(df_prey_est, df_pr_est)
}
gridExtra::grid.arrange(grobs = plt_trends)

df_prey_est$Prey = factor(df_prey_est$Prey, levels = PreyTypes)

# Add pretty text version of prey types, convert to factor to maintain same order in plots
Prey_txt = stringr::str_replace_all(df_prey_est$Prey, "_", " ") 
Prey_txt = stringr::str_to_sentence(Prey_txt)
df_prey_est$Prey_txt = factor(Prey_txt, levels = unique(Prey_txt))

# Compute diagnostics
summary_stats <- summarize_draws(mcmc)

# 2️⃣ R-hat Distribution
rhat_plot <- mcmc_rhat(summary_stats$rhat) +
  ggtitle("R-hat Convergence")

# 3️⃣ Effective Sample Size (ESS)
ess_plot <- mcmc_neff(summary_stats$ess_bulk) +
  ggtitle("Effective Sample Size")




################################################################################
#FigS2 Estimated trends in prey encounter rates

plt_trends2 = ggplot(filter(df_prey_est,Year>2007 & (Prey == "urchin" | Prey == "mussel")),aes(x=Year,y=Dns)) +
  geom_ribbon(aes(ymin=Dns_lo,ymax=Dns_hi,fill=Prey),alpha=.3) +
  geom_line(aes(color=Prey),linewidth=1.1) +
  geom_vline(xintercept = 2013, linetype="dashed") +
  # ylab(expression(paste("logit (", lambda,")"))) +
  ylab("Density") +
  ggtitle(paste0("Trends in prey encounter rates, urchins and mussels")) +
  theme_classic() + 
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  facet_wrap(vars(Prey),nrow = 2)
print(plt_trends2)

new_palette <- brewer.pal(name="Paired",n=10)[c(4,10,1,6,8,3,5,7,2,9)]
base_theme <- theme(axis.text=element_text(size=12,color = "black"),
                     axis.title=element_text(size=13,color = "black"),
                     plot.tag=element_text(size=11,color = "black"),
                     plot.title=element_text(size=11,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     #legend.key = element_rect(fill = "white"), # Set it to transparent
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=9,color = "black"),
                     legend.title=element_blank(),
                     #legend.key.height = unit(0.1, "cm"),
                     #legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=10, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())

plt_trends3 = ggplot(filter(df_prey_est, Year>2007 & Prey != "other"),aes(x=Year,y=Dns)) + # & Prey != "urchin" & Prey != "mussel" 
  geom_ribbon(aes(ymin=Dns_lo,ymax=Dns_hi),alpha=.25) +
  geom_line(linewidth=1.1) +
  # ylab(expression(paste("logit (", lambda,")"))) +
  geom_vline(xintercept = 2013, linetype="dashed") +
  ylab("Density") +
  ggtitle(paste0("")) +
  scale_color_manual(values = new_palette) +
  scale_fill_manual(values = new_palette) +
  facet_wrap(vars(Prey_txt),nrow = 5, scales = "free")+
  theme_classic() + base_theme

print(plt_trends3)

ggsave(plt_trends3, filename = file.path(figdir, "FigS2_prey_density.png"), 
       width =7, height = 8, units = "in", dpi = 600, bg = "white") #last write Feb 19, 2025


################################################################################
#Plot Figure 5 - energetic intake

ii = which(startsWith(vn,"Ebar["))
iiu = which(startsWith(vn,"U_prd[")) 
iim = which(startsWith(vn,"M_prd[")) 

df_Erate_est = data.frame(Year=Years,Scenario=rep("Actual",N),
                          Erate=sumstats$mean[ii],
                          Erate_lo=sumstats$`5%`[ii], 
                          Erate_hi=sumstats$`95%`[ii],
                          U_prd=sumstats$mean[iiu],
                          U_prd_lo=sumstats$`5%`[iiu], 
                          U_prd_hi=sumstats$`95%`[iiu],
                          M_prd=sumstats$mean[iim],
                          M_prd_lo=sumstats$`5%`[iim], 
                          M_prd_hi=sumstats$`95%`[iim] )

for(j in 1:3){
  ii = which(startsWith(vn,paste0("Ebar_alt[",j)))
  iiu = which(startsWith(vn,paste0("U_prd_alt[",j)))
  iim = which(startsWith(vn,paste0("M_prd_alt[",j)))
  df_Erate_est = rbind(df_Erate_est, 
                       data.frame(Year=Years,
                                  Scenario=rep(paste0("Scenario ",j),N),
                                  Erate=sumstats$mean[ii],
                                  Erate_lo=sumstats$`5%`[ii], 
                                  Erate_hi=sumstats$`95%`[ii],
                                  U_prd=sumstats$mean[iiu],
                                  U_prd_lo=sumstats$`5%`[iiu], 
                                  U_prd_hi=sumstats$`95%`[iiu],
                                  M_prd=sumstats$mean[iim],
                                  M_prd_lo=sumstats$`5%`[iim], 
                                  M_prd_hi=sumstats$`95%`[iim] ) )
  
}

base_theme2 <-  theme(axis.text=element_text(size=10, color = "black"),
                     axis.title=element_text(size=11,color = "black"),
                     plot.tag=element_text(size=10,color = "black"),
                     plot.title=element_text(size=11,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     #panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.4, "cm"), 
                     #legend.key = element_rect(fill = "white"), # Set it to transparent
                     legend.spacing.y = unit(0.4, "cm"),  
                     legend.text=element_text(size=10,color = "black"),
                     legend.title=element_blank(),
                     #legend.key.height = unit(0.1, "cm"),
                     #legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=9, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())

ii = which(df_Erate_est$Scenario == "Actual")
tmp = df_Erate_est[ii,c(1,3,4,5)]


df_Erate_est$Scenario2 <- factor(df_Erate_est$Scenario, levels = c("Actual", "Scenario 1", "Scenario 2", "Scenario 3"),
                                labels = c("Observed energetic intake",
                                           "Scenario 1: no sea urchin increase",
                                           "Scenario 2: no mussel increase",
                                           "Scenario 3: no sea urchin or mussel increase"))



dark2_colors <- RColorBrewer::brewer.pal(8, "Dark2")  
dark2_colors  

plt_Erate = ggplot(df_Erate_est, aes(x = Year, y = Erate)) +
  geom_ribbon(aes(ymin = Erate_lo, ymax = Erate_hi, fill = Scenario), alpha = 0.3) +
  geom_line(aes(color = Scenario), lwd=0.6) + 
  geom_ribbon(data = tmp, aes(ymin = Erate_lo, ymax = Erate_hi), alpha = 0.1) +
  geom_line(data = tmp, aes(x = Year, y = Erate), linetype = "dotdash") + 
  geom_line(data = tmp, aes(x = Year, y = Erate_lo), linetype = "dotted") + 
  geom_line(data = tmp, aes(x = Year, y = Erate_hi), linetype = "dotted") + 
  geom_vline(xintercept = 2013, linetype = "dashed") +
  labs(y = "Energy intake \n(Kcal / min)", 
       title = "") +
  scale_fill_manual(values = c("Actual" = "gray90", 
                               "Scenario 1" = dark2_colors[2],  # Orange
                               "Scenario 2" = dark2_colors[3],  # Purple
                               "Scenario 3" = dark2_colors[4])) +  # Pink
  scale_color_manual(values = c("Actual" = "black", 
                                "Scenario 1" = dark2_colors[2],  
                                "Scenario 2" = dark2_colors[3],  
                                "Scenario 3" = dark2_colors[4])) +
  facet_wrap(~Scenario2, ncol = 1, strip.position = "top") +  
  scale_x_continuous(breaks = seq(2006, max(df_Erate_est$Year), by = 2),
                     guide = guide_axis(angle = 45)) + 
  theme_classic() + base_theme2 +  
  theme(legend.position = "none")

print(plt_Erate)


ggsave(plt_Erate, filename = file.path(figdir, "Fig5_energetic_intake.png"), 
      width =4, height = 7, units = "in", dpi = 600, bg = "white") #last write Feb 19, 2025


################################################################################
#Plot Figure S3 urchin predation


ii = which(df_Erate_est$Scenario == "Actual")
tmp = df_Erate_est[ii,c(1,6,7,8)]

# Ensure colors match plt_Erate
plt_U_pred = ggplot(df_Erate_est, aes(x = Year, y = U_prd)) +
  geom_ribbon(aes(ymin = U_prd_lo, ymax = U_prd_hi, fill = Scenario), alpha = 0.3) +  # Match transparency
  geom_line(aes(color = Scenario), size = 0.6) +  # Match line width
  geom_ribbon(data = tmp, aes(ymin = U_prd_lo, ymax = U_prd_hi), alpha = 0.1) +
  geom_line(data = tmp, aes(x = Year, y = U_prd), linetype = "dotdash", size = 0.6) +  
  geom_line(data = tmp, aes(x = Year, y = U_prd_lo), linetype = "dotted", size = 0.6) +  
  geom_line(data = tmp, aes(x = Year, y = U_prd_hi), linetype = "dotted", size = 0.6) +  
  geom_vline(xintercept = 2013, linetype = "dashed") +
  labs(y = "Urchins consumed per minute", 
       title = "") +  # Remove title to match plt_Erate
  scale_fill_manual(values = c("Actual" = "gray90", 
                               "Scenario 1" = dark2_colors[2],  # Orange
                               "Scenario 2" = dark2_colors[3],  # Purple
                               "Scenario 3" = dark2_colors[4])) +  # Pink
  scale_color_manual(values = c("Actual" = "black", 
                                "Scenario 1" = dark2_colors[2],  
                                "Scenario 2" = dark2_colors[3],  
                                "Scenario 3" = dark2_colors[4])) +
  facet_wrap(~Scenario2, ncol = 1, strip.position = "top") +   # Use facet_wrap() to match plt_Erate
  scale_x_continuous(breaks = seq(2006, max(df_Erate_est$Year), by = 2),
                     guide = guide_axis(angle = 45)) + 
  theme_classic() + base_theme2 +  
  theme(legend.position = "none")

print(plt_U_pred)


ggsave(plt_U_pred, filename = file.path(figdir, "FigS3_urchin_consumption.png"), 
       width =4, height = 7, units = "in", dpi = 600, bg = "white") #last write Feb 19, 2025


################################################################################
#Update tables

# Define indices for pre- and post-2013
pre_index <- which(Years < 2013)
post_index <- which(Years >= 2013)

#verfiy dim
cat("Pre-2013 Years: ", Years[pre_index], "\n")
cat("Post-2013 Years: ", Years[post_index], "\n")

#Get proportion effort
pre_effort_mean <- colMeans(eta_obs[pre_index, ]) * 100
pre_effort_sd <- apply(eta_obs[pre_index, ], 2, sd) * 100
post_effort_mean <- colMeans(eta_obs[post_index, ]) * 100
post_effort_sd <- apply(eta_obs[post_index, ], 2, sd) * 100

cat("Pre-2013 Effort Means: ", pre_effort_mean, "\n")
cat("Pre-2013 Effort SDs: ", pre_effort_sd, "\n")
cat("Post-2013 Effort Means: ", post_effort_mean, "\n")
cat("Post-2013 Effort SDs: ", post_effort_sd, "\n")

#Get proportion effort
E_mean <- E_obs
E_sd <- sd_E

df_prey_table <- data.frame(
  PreyType = PreyTypes,
  E_kcal_m = E_mean,
  E_SD = E_sd,
  Pre2013_Effort = pre_effort_mean,
  Pre2013_Effort_SD = pre_effort_sd,
  Post2013_Effort = post_effort_mean,
  Post2013_Effort_SD = post_effort_sd
)

write.csv(df_prey_table, file.path(here::here("tables","prey_table.csv"))) #last write 19 Feb 2025






