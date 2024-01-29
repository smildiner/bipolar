#---------------------------------------------------
# Fit model for M hidden states and save results
#---------------------------------------------------

# devtools::install_github("smildiner/mHMMbayes", ref = "develop")
# devtools::install_github("smildiner/simHMM", ref = "dev")

# Preferred parallel library
library(tidyverse)
library(mHMMbayes)
library(simHMM)
library(pbmcapply)

#------------------------------------------------------------------------------#

# Def
m <- 4
n_dep <- 12


# Get stating values
out <- readRDS(paste0("results/out_cont_gamma_prior_emiss_prior_m4_12dv_it4000_c2_night.rds"))


#=============================
# Def starting values gamma
gamma_map <- int_to_prob(matrix(apply(out$gamma_int_bar[2001:4000,], 2, median),
                                nrow = m, ncol = m-1, byrow = TRUE))

gamma_map[1,1] <- gamma_map[1,1] + 1-apply(gamma_map, 1, sum)[1]
gamma_map[2,2] <- gamma_map[2,2] + 1-apply(gamma_map, 1, sum)[2]
gamma_map[3,3] <- gamma_map[3,3] + 1-apply(gamma_map, 1, sum)[3]
gamma_map[4,4] <- gamma_map[4,4] + 1-apply(gamma_map, 1, sum)[4]

# Def starting values emiss
emiss_map <- lapply(1:12, function(q) matrix(c(apply(out$emiss_mu_bar[[q]][2001:4000,],2,median),
                                               apply(out$emiss_var_bar[[q]][2001:4000,],2,median)),
                                             nrow = m, ncol = 2, byrow = FALSE))

## Random effect: between subject variance
# Transitions
gamma_var_map <- apply(out$gamma_V_int_bar[2001:4000,] %>%
                         as.data.frame() %>%
                         dplyr::select(paste0("var_int_S",rep(1:4,each=3),"toS",2:4,"_with_int_S",rep(1:4,each=3),"toS",2:4)),
                       2, median) %>%
  as.numeric()
gamma_var_map <- matrix(gamma_var_map, nrow = m, ncol = m-1, byrow = TRUE)

# Emissions
emiss_varmu_map <- lapply(out$emiss_varmu_bar, function(emiss) {
  emiss %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter > 2000) %>%
    dplyr::select(-iter) %>%
    summarise(across(.cols = everything(), median)) %>%
    gather(mu, value) %>%
    dplyr::select(value) %>%
    as.matrix()
})

#=============================
# Add hyperpriors
# using the function prob_to_int to obtain intercept values for the above specified
# transition probability matrix gamma
gamma_mu0 <-  prob_to_int(matrix(c(0.7, 0.1, 0.1, 0.1,
                                   0.1, 0.7, 0.1, 0.1,
                                   0.1, 0.1, 0.7, 0.1,
                                   0.1, 0.1, 0.1, 0.7), nrow = m, ncol = m, byrow = TRUE))
gamma_mu0 <- list(matrix(gamma_mu0[1,], nrow = 1),
                  matrix(gamma_mu0[2,], nrow = 1),
                  matrix(gamma_mu0[3,], nrow = 1),
                  matrix(gamma_mu0[4,], nrow = 1))
gamma_K0 <- 0.1
gamma_nu <- 5
gamma_V <- diag(8, m - 1)

gamma_hyp_prior <- list("gamma_mu0" = gamma_mu0,
                        "gamma_K0" = gamma_K0,
                        "gamma_nu" = gamma_nu,
                        "gamma_V" = gamma_V)

# Specify hyper-prior for the emission distribution
emiss_hyp_prior <- list(
  emiss_mu0 = list(matrix(c(5, 25, 50, 50), nrow = 1), # down
                   matrix(c(5, 25, 50, 50), nrow = 1), # dread rest of day
                   matrix(c(5, 25, 50, 50), nrow = 1), # worry
                   matrix(c(5, 25, 50, 50), nrow = 1), # inadequate
                   matrix(c(25, 25, 50, 50), nrow = 1), # tired
                   matrix(c(50, 25, 5, 5), nrow = 1), # content
                   
                   matrix(c(5, 25, 25, 10), nrow = 1), # agitated
                   matrix(c(5, 25, 25, 10), nrow = 1), # irritated
                   matrix(c(25, 25, 25, 10), nrow = 1), # switch
                   matrix(c(25, 25, 25, 10), nrow = 1), # extremely well
                   matrix(c(25, 25, 25, 10), nrow = 1), # ideas
                   matrix(c(25, 25, 25, 10), nrow = 1)  # thoughts racing
  ),
  emiss_K0  = rep(list(0.1),n_dep),
  emiss_nu  = rep(list(1),n_dep),
  emiss_V   = rep(list(rep(400, m)),n_dep),
  emiss_a0  = rep(list(rep(0.001, m)),n_dep),
  emiss_b0  = rep(list(rep(0.001, m)),n_dep))

rm(out)

#==================================================================================================#
# Outstart clean:
out_list <- pbmcapply::pbmclapply(1:100, function(s) {
  
  # Gen
  case_out <- NULL
  case_out_miss <- NULL
  case_out_clean <- NULL
  sim_data <- NULL
  
  # Specify simulation parameters  
  m       <- 4
  n_dep   <- 12
  n_subj  <- 20
  n_t     <- 30 * 4 * 8 * 2  # 30 (days per month) * 4 (months) * 8 (occasions per day)
  n_iter  <- 4000
  burn_in <- 2000
  
  # 1) Simulate data
  while(is.null(sim_data)){
    try({
      sim_data <- mHMMbayes::sim_mHMM_plnorm(n_t = n_t, n = n_subj,
                                             data_distr = "continuous",
                                             m = m, n_dep = n_dep,
                                             gamma = gamma_map,
                                             emiss_distr = emiss_map,
                                             var_gamma = gamma_var_map,
                                             var_emiss = emiss_varmu_map,
                                             return_ind_par = TRUE)
    })
  }
  
  # Add floor and ceilling effects
  sim_data$obs[sim_data$obs < 0] <- 0
  sim_data$obs[sim_data$obs >= 100] <- 100
  sim_data$obs[sample(1:nrow(sim_data$obs), size = 0.24*nrow(sim_data$obs),replace = FALSE),-1] <- NA
  
  # Delete
  night_occasions <- as.numeric(t(sapply(6:8, function(o) seq(0, 120*8*20-8, 8)+o)))
  
  train_idx <- rep(c(rep(TRUE, 960), rep(FALSE, 960)), 20)
  test_idx <- !train_idx
  
  train_data_raw <- sim_data$obs[train_idx,]
  test_data_raw <- sim_data$obs[test_idx,]
  
  test_data_night <- test_data_raw
  test_data_night[night_occasions, -1] <- NA
  test_data_clean <- test_data_raw[-night_occasions,]
  
  train_data_night <- train_data_raw
  train_data_night[night_occasions, -1] <- NA
  train_data_clean <- train_data_raw[-night_occasions,]
  
  # 4) Fit model on complete data
  try({
    ti <- Sys.time()
    case_out <- mHMMbayes::mHMM_cont(s_data = train_data_raw,
                                     gen = list(m = m, n_dep = n_dep),
                                     start_val = c(list(gamma_map), emiss_map),
                                     gamma_hyp_prior = gamma_hyp_prior,
                                     emiss_hyp_prior = emiss_hyp_prior,
                                     mcmc = list(J = n_iter, burn_in = burn_in),
                                     show_progress = TRUE)
    case_out[["time"]] <- Sys.time() - ti
  })
  
  # 8) Save output
  if(!is.null(case_out)){
    
    # Get MAP estimates
    map_out <- simHMM::MAP(case_out)
    
    # Get credibility intervals
    cci_out <- simHMM::get_cci(case_out)
    
    # Define new output object
    case_out <- list(
      truth = sim_data,
      train = train_data_raw,
      test = test_data_raw,
      map = map_out,
      cci = cci_out)
    
    # Save output
    saveRDS(case_out, paste0("outputs/res/out",
                             "_m",m,
                             "_it",n_iter,
                             "_burnin",burn_in,
                             "_rep",s,"_raw.rds"))
    
  }
  
  # Fit model on data imputing overnight measurements
  try({
    ti <- Sys.time()
    case_out_miss <- mHMMbayes::mHMM_cont(s_data = train_data_night,
                                          gen = list(m = m, n_dep = n_dep),
                                          start_val = c(list(gamma_map), emiss_map),
                                          gamma_hyp_prior = gamma_hyp_prior,
                                          emiss_hyp_prior = emiss_hyp_prior,
                                          mcmc = list(J = n_iter, burn_in = burn_in),
                                          show_progress = TRUE)
    case_out_miss[["time"]] <- Sys.time() - ti
  })
  
  # 8) Save output
  if(!is.null(case_out_miss)){
    
    # Get MAP estimates
    map_out <- simHMM::MAP(case_out_miss)
    
    # Get credibility intervals
    cci_out <- simHMM::get_cci(case_out_miss)
    
    # Define new output object
    case_out_miss <- list(
      truth = sim_data,
      train = train_data_night,
      test = test_data_night,
      map = map_out,
      cci = cci_out)
    
    # Save output
    saveRDS(case_out_miss, paste0("outputs/res/out",
                                  "_m",m,
                                  "_it",n_iter,
                                  "_burnin",burn_in,
                                  "_rep",s,"_imputed.rds"))
    
  }
  
  # Fit model on data omitting overnight measurements
  try({
    ti <- Sys.time()
    case_out_clean <- mHMMbayes::mHMM_cont(s_data = train_data_clean,
                                           gen = list(m = m, n_dep = n_dep),
                                           start_val = c(list(gamma_map), emiss_map),
                                           gamma_hyp_prior = gamma_hyp_prior,
                                           emiss_hyp_prior = emiss_hyp_prior,
                                           mcmc = list(J = n_iter, burn_in = burn_in),
                                           show_progress = TRUE)
    case_out_clean[["time"]] <- Sys.time() - ti
  })
  
  # 8) Save output
  if(!is.null(case_out_clean)){
    
    # Get MAP estimates
    map_out <- simHMM::MAP(case_out_clean)
    
    # Get credibility intervals
    cci_out <- simHMM::get_cci(case_out_clean)
    
    # Define new output object
    case_out_clean <- list(
      truth = sim_data,
      train = train_data_clean,
      test = test_data_clean,
      map = map_out,
      cci = cci_out)
    
    # Save output
    saveRDS(case_out_clean, paste0("outputs/res/out",
                                   "_m",m,
                                   "_it",n_iter,
                                   "_burnin",burn_in,
                                   "_rep",s,"_omitted.rds"))
    
  }
  
}, mc.cores = 100, mc.set.seed = 1L)
