//
// This Stan program defines a foraging model for sea otters,
// to predict changes in diet and energy gain reflecting changes in 
// relative abundance of alternate prey types. Assumes optimal 
// allocation of effort among patches of different prey types 
// (i.e. to maximize overall rate of energy return)
//
// The input data
data {
  int<lower=0> N ;                 // Number years/periods of observation
  int<lower=0> K ;                 // Number of prey types
  array[N] simplex[K] eta_obs ;    // observation-based estimate of foraging effort allocation, by year
  array[N] real<lower=0> tau ;     // relative precision (dirichlet) of observed foraging effort allocations, by year
  array[N] vector[K] mu_obs ;      // prey-specific log rates of energy return, from obs foraging data (SOFA)  
  vector[K] sig_mu ;               // variance in log rate of energy return estimates from obs foraging data (SOFA)  
  array[2] int<lower=0> N_base ;   // year range for "base" dynamics
  array[2] vector[K] log_G_pri ;   // prior params for log caloric intake per item (mean and sd)
}
//
transformed data {
  int Lng_base = N_base[2] - N_base[1] + 1 ;
  array[Lng_base] int iib = linspaced_int_array(Lng_base, N_base[1], N_base[2]) ;
}
//
// The base parameters accepted by the model. 
parameters {
  vector[K] log_G ;                // mean log Caloric intake per item of prey consumed for prey type i  
  vector<lower=0>[K] mu ;          // mean log rate of energy return during dives allocated to foraging on prey type i  
  vector<lower=0>[K] sig_E;        // sd in log rate of energy return 
  vector<lower=0>[3] sig_D ;       // year-to-year variation in relative abundance of prey (urchins, musseles, all other taxa) 
  matrix[N,K] D ;                  // log relative abundance of prey types by year (log-instantaneous encounter rates)   
}
//
// Derived/transformed parameters
transformed parameters {    
  array[N] vector[K] delta ;       // relative density for each prey type by year (instantaneous encounter rates)
  array[N] vector[K] lambda ;      // finite patch encounter probability for each prey type, each year
  array[N] vector[K] pi ;          // expected alocation of foraging effort among patches of different prey types
  // For mean energy gain, log-normal variance adjustment accounts for estimation uncertainty in observed log mu:
  vector[K] E = exp(mu + square(sig_mu) ./ 2) ; // rate of energy return during dives allocated to foraging on prey type i  
  // ()
  array[K-1] vector[K-1] W ;       // switch variable: identifies prey types more profitable than focal prey   
  for(i in 1:(K-1)){
    W[i] = inv_logit(10 * (E[1:(K-1)] - E[i]) ) ;
    W[i][i] = 0 ;
  }
  for(t in 1:N){    
    vector[K-1] log_pi ;    
    delta[t] = exp(D[t]') ;
    lambda[t] = 1 - exp(-delta[t]) ;
    for(i in 1:(K-1)){
      // Allocation effort: proportional to joint probability of encountering patch of
      // focal prey type & not encountering patches of more profitable prey, on any dive
      // (most profitable prey types always selected when encountered)
      log_pi[i] = log(lambda[t][i]) + sum(log(1 - lambda[t][1:(K-1)] .* W[i] ) ) ; 
    }
    pi[t][1:(K-1)] = exp(log_pi) ;
    pi[t][K] = 1 - sum(exp(log_pi)) ;
  }
}
//
// The model to be estimated. 
model {
  // Observed data: 
  for(t in 1:N){  
    mu_obs[t] ~ normal(mu,sig_E) ;           // Observed mean log rate of energy return 
    eta_obs[t] ~ dirichlet(pi[t] * tau[t]) ; // Observed foraging effort allocation by year
  }
  // Priors:
  log_G ~ normal(log_G_pri[1],log_G_pri[2]) ;
  mu ~ cauchy(0,1) ;
  sig_E ~ cauchy(0,.1) ;  
  sig_D ~ cauchy(0,.1) ;
  // Auto-regressive model (AR[1]) describing log instantaneous prey patch encounter rates over time
  for(i in 1:(K-1)){
    int j = i < 3 ? i : 3 ;
    D[1,i] ~ normal(-2.5, 2.5) ;
    D[2:N,i] ~ normal( D[1:(N-1),i], sig_D[j]) ;
  }  
  D[1:N,K] ~ normal(0,.05) ; // "other" prey category assumed to remain constant  
}
//
generated quantities {
  vector[N] Ebar ;                   // estimated energy intake rate, by year, base scenario     
  vector[N] U_prd ;                  // Urchin predation rate, by year, base scenario     
  vector[N] M_prd ;                  // Mussel predation rate, by year, base scenario     
  array[3] vector[N] Ebar_alt ;      // estimated energy intake, alternative scenarios
  array[3] vector[N] U_prd_alt ;     // Urchin predation rate, alternative scenarios
  array[3] vector[N] M_prd_alt ;     // Mussel predation rate, alternative scenarios  
  vector[K] G = exp(log_G) ;         // mean Caloric intake per item of prey consumed for prey type i  
  for(t in 1:N){
    array[3] vector[K] lam = rep_array(lambda[t], 3) ;
    vector[K-1] log_pi ;    
    array[3] vector[K] pi_alt ;
    // Calculate mean energy gain and urchin/mussel consumotion rate for observed (default) scenario:  
    Ebar[t] = sum(pi[t] .* E ) ;
    U_prd[t] = pi[t][1] * (E[1] / G[1]) ;
    M_prd[t] = pi[t][2] * (E[2] / G[2]) ;
    if(t > N_base[2]){
      // Alternative scenario 1: urchins remain fixed at avg of base years
      lam[1][1] = 1 - exp(-exp( normal_rng(mean(D[iib,1]), sd(D[iib,1])) )) ;
      // Alternative scenario 2: mussels remain fixed at avg of base years
      lam[2][2] =  1 - exp(-exp( normal_rng(mean(D[iib,2]), sd(D[iib,2])) )) ;
      // Alternative scenario 3: urchins & mussels remain fixed at avg of base years
      lam[3][1] = lam[1][1] ;
      lam[3][2] = lam[2][2] ;
    }
    // Calculate mean energy gain and urchin/mussel consumotion rate for alternative scenarios:  
    for(j in 1:3){
      for(i in 1:(K-1)){
        log_pi[i] = log(lam[j][i]) + sum(log(1 - lam[j][1:(K-1)] .* W[i] ) ) ; 
      }
      pi_alt[j][1:(K-1)] = exp(log_pi) ;
      pi_alt[j][K] = 1 - sum(exp(log_pi)) ;
      Ebar_alt[j][t] = sum(pi_alt[j] .* E ) ;
      U_prd_alt[j][t] = pi_alt[j][1] * (E[1] / G[1]) ;
      M_prd_alt[j][t] = pi_alt[j][2] * (E[2] / G[2]) ;
    }
  }
}
