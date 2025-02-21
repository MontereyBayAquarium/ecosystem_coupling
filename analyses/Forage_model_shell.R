require(dplyr)
require(parallel)
library(cmdstanr)
library(posterior)
library(gtools)
library(ggplot2)
library(bayesplot)
library(readxl)
library(RColorBrewer)
#devtools::install_github("dkahle/dirichlet")
library(dirichlet)
rstan::rstan_options(javascript=FALSE)
# Load data ----------------------------------------------
# SOFA_rslts <- mget(load("Rslt_Grp-Period_2024_Sep_04_14hr.rdata", envir=(NE. <- new.env())), envir=NE.)
# rm(NE.)
df_Prey = read_excel("Preytypes.xlsx")
df_SOFA_rslts = read_excel("SOFA_posterior_sums.xlsx")

# Process data -------------------------------------------
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
V_mu = df_SOFA_rslts$Value[df_SOFA_rslts$Param=="V_mu"]
pi_obs = matrix(df_SOFA_rslts$Value[df_SOFA_rslts$Param=="pi"],nrow = N, ncol=K, byrow = F)
tau =  df_SOFA_rslts$Value[df_SOFA_rslts$Param=="tau"]
E_obs = exp( colMeans(mu_obs) +V_mu/2)
sd_E = sqrt( (exp(V_mu) - 1)*exp(2*colMeans(mu_obs) + V_mu) )
#
# Set up Stan fitting --------------------------------------
nsamples = 3000
nburnin = 2000
cores = detectCores()
ncore = max(3,min(5,cores-2))
Niter = round(nsamples/ncore)
fitmodel = "Foraging_fit.stan"
#
stan.data = list(N=N,K=K,pi_obs=pi_obs,tau=tau,mu_obs=mu_obs,V_mu=V_mu,
                 N_base=N_base,log_G_pri=log_G_pri)
#
parms = c("sig_L","sig_E","Ebar","Ebar_alt","pi","delta",
          "U_prd","M_prd","U_prd_alt","M_prd_alt")
#
mod <- cmdstan_model(fitmodel)

init_fun <- function() {list(sig_E=runif(K, .04, .06), 
                             sig_L=runif(1, .2, .25), 
                             mu = runif(K, .9*colMeans(mu_obs), 1.1*colMeans(mu_obs)) ) }

#
suppressMessages(
  suppressWarnings ( 
    fit <- mod$sample(
      data = stan.data,
      init <- init_fun ,
      seed = 111,
      chains = ncore,
      parallel_chains = ncore,
      refresh = 100,
      # thin=thinval,
      iter_warmup = nburnin,
      iter_sampling = Niter )
  )
)
source("cmdstan_sumstats.r")
#
# Examine results -----------------------------------------------
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

df_prey_est$Prey = factor(df_prey_est$Prey, levels = Preynames)
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

plt_trends2 = ggplot(filter(df_prey_est,Year>2007 & (Prey == "urchin" | Prey == "mussel")),aes(x=Year,y=Dns)) +
  geom_ribbon(aes(ymin=Dns_lo,ymax=Dns_hi,fill=Prey),alpha=.3) +
  geom_line(aes(color=Prey),linewidth=1.1) +
  geom_vline(xintercept = 2013, linetype="dashed") +
  # ylab(expression(paste("logit (", lambda,")"))) +
  ylab("Density") +
  ggtitle(paste0("Trends in prey density, urchins and mussels")) +
  theme_classic() + 
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  facet_wrap(vars(Prey),nrow = 2)
print(plt_trends2)

new_palette <- brewer.pal(name="Paired",n=10)[c(4,10,1,6,8,3,5,7,2,9)]
plt_trends3 = ggplot(filter(df_prey_est,Year>2007 & Prey != "other"),aes(x=Year,y=Dns)) + # & Prey != "urchin" & Prey != "mussel" 
  geom_ribbon(aes(ymin=Dns_lo,ymax=Dns_hi,fill=Prey),alpha=.25) +
  geom_line(aes(color=Prey),linewidth=1.1) +
  # ylab(expression(paste("logit (", lambda,")"))) +
  geom_vline(xintercept = 2013, linetype="dashed") +
  ylab("Density") +
  ggtitle(paste0("Trends in Prey density, alternative prey")) +
  theme_classic() + 
  scale_color_manual(values = new_palette) +
  scale_fill_manual(values = new_palette) +
  facet_wrap(vars(Prey),nrow = 5, scales = "free")
print(plt_trends3)

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

ii = which(df_Erate_est$Scenario == "Actual")
tmp = df_Erate_est[ii,c(1,3,4,5)]
plt_Erate = ggplot(df_Erate_est,aes(x=Year,y=Erate)) +
  geom_ribbon(aes(ymin=Erate_lo,ymax=Erate_hi,fill=Scenario),alpha=.3) +
  geom_line(aes(color=Scenario)) + 
  geom_ribbon(data=tmp,aes(ymin=Erate_lo,ymax=Erate_hi),alpha=.1) +
  geom_line(data=tmp,aes(x=Year,y=Erate),linetype="dotdash") + 
  geom_line(data=tmp,aes(x=Year,y=Erate_lo),linetype="dotted") + 
  geom_line(data=tmp,aes(x=Year,y=Erate_hi),linetype="dotted") + 
  geom_vline(xintercept = 2013, linetype="dashed") +
  labs(y="Energy Intake") +
  ggtitle("Effect of urchin and mussel increase on Energy Intake",
          subtitle = paste0("Scenario 1 = no urchin increase, ",
                            "Scenario 2 = no mussel increase, ",
                            "Scenario 3 = no urchin or mussel increase")) +
  theme_classic() + 
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(rows = vars(Scenario))
print(plt_Erate)

ii = which(df_Erate_est$Scenario == "Actual")
tmp = df_Erate_est[ii,c(1,6,7,8)]
plt_U_pred = ggplot(df_Erate_est,aes(x=Year,y=U_prd)) +
  geom_ribbon(aes(ymin=U_prd_lo,ymax=U_prd_hi,fill=Scenario),alpha=.25) +
  geom_line(aes(color=Scenario)) + 
  geom_ribbon(data=tmp,aes(ymin=U_prd_lo,ymax=U_prd_hi),alpha=.1) +
  geom_line(data=tmp,aes(x=Year,y=U_prd),linetype="dashed") + 
  geom_line(data=tmp,aes(x=Year,y=U_prd_lo),linetype="dotted") + 
  geom_line(data=tmp,aes(x=Year,y=U_prd_hi),linetype="dotted") + 
  geom_vline(xintercept = 2013, linetype="dashed") +
  labs(y="Urchins consumed per minute") +
  ggtitle("Effect of urchin and mussel increase on urchin predation rate",
          subtitle = paste0("Scenario 1 = no urchin increase, ",
                            "Scenario 2 = no mussel increase, ",
                            "Scenario 3 = no urchin or mussel increase")) +
  theme_classic() + 
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  facet_grid(rows = vars(Scenario))
print(plt_U_pred)
