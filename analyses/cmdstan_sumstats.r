# Script to generate sumstats and mcmc matrix similar to rstan,
#  after fitting a model using cmdstan 
#  NOTE: following lines show sample of cmdstan fitting code
# ------------------------------------------------------
# require(parallel)
# require(cmdstanr)
# require(posterior)
# nburnin = 500                    # number warm-up (burn-in) samples 
# nsamples = 10000                 # desired total number samples
# fitmodel = c("filename.stan")    # name of file with stan code
# parms = c("Par1","Par2")         # vector of parameters to save
# stan.data <- list(N=N,X=X,Y=Y)   # list of input data variables
# cores = detectCores()
# ncore = min(40,cores-1)
# Niter = round(nsamples/ncore)
# mod <- cmdstan_model(fitmodel)   # compiles model (if necessary) 
# suppressMessages(                # Suppress messages/warnings (if desired)
#   suppressWarnings ( 
#     fit <- mod$sample(
#       data = stan.data,
#       seed = 1234,
#       chains = ncore,
#       parallel_chains = ncore,
#       refresh = 100,
#       iter_warmup = nburnin,
#       iter_sampling = Niter,
#       max_treedepth = 12,
#       adapt_delta = 0.8
#     )
#   )
# )
# ------------------------------------------------------

if (!exists("select_chains")) {
  fit_draws = fit$draws(parms) 
}else{
  fit_draws = fit$draws(parms) %>% subset_draws(chain=select_chains)
}

# sumstats = fit$draws(parms) %>% subset_draws(chain=select_chains) %>%  
sumstats = fit_draws %>% 
  summarise_draws(mean, mcse = mcse_mean, sd, 
                  ~quantile(.x, probs = c(0.025, 0.05, .5, .95, .975), na.rm=T),
                  N_eff = ess_bulk, rhat)
sumstats = as.data.frame(sumstats)
row.names(sumstats) = sumstats$variable; sumstats = sumstats[,-1] 
#
mcmc = as_draws_matrix(fit_draws,variables = parms)
Nsims = nrow(mcmc)
mcmc_array = as_draws_array(fit_draws,variables = parms)
vn = colnames(mcmc); vns = row.names(sumstats)
#
# sumstats = as.data.frame(fit$summary(variables = parms))
# row.names(sumstats) = sumstats$variable; sumstats = sumstats[,-1] 
# tmp = as.data.frame(fit$summary(variables = parms, mcse = mcse_mean, 
#                                 ~quantile(.x, probs = c(0.025, 0.975), na.rm=T)))
# sumstats$mcse = tmp$mcse 
# sumstats$q2.5 = tmp$`2.5%` 
# sumstats$q97.5 = tmp$`97.5%`
# sumstats$q50 = sumstats$median 
# sumstats$N_eff = sumstats$ess_bulk
# col_order = c("mean","mcse","sd","q2.5","q10","q50","q90","q97.5","N_eff","rhat")
# sumstats = sumstats[, col_order]
# mcmc = as_draws_matrix(fit$draws(variables = parms))
# vn = colnames(mcmc); vns = row.names(sumstats)
# Nsims = nrow(mcmc)
# rm(tmp,col_order)   # remove temporary variables
# rm(mod,fit)       # can also remove mod & fit from memory if no longer needed 
