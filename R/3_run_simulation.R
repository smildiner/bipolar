#---------------------------------------------------
# Fit model for M hidden states and save results
#---------------------------------------------------

# devtools::install_github("smildiner/simHMM", ref = "dev")
# devtools::install_github("smildiner/mHMMbayes", ref = "develop")

# Preferred parallel library
library(tidyverse)
library(mHMMbayes)
library(simHMM)
library(pbmcapply)

#------------------------------------------------------------------------------#

# Def
m <- 4
n_dep <- 12

# Get starting values
out <- readRDS(paste0("results/out_cont_gamma_prior_emiss_prior_m4_12dv_it4000_c2.rds"))


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
gamma_K0 <- 1
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
  emiss_K0  = rep(list(1),n_dep),
  emiss_nu  = rep(list(1),n_dep),
  emiss_V   = rep(list(rep(400, m)),n_dep),
  emiss_a0  = rep(list(rep(0.001, m)),n_dep),
  emiss_b0  = rep(list(rep(0.001, m)),n_dep))


#==================================================================================================#
# Simulate:

pbmcapply::pbmclapply(1:200, function(s) {
  
  # Gen
  case_out <- NULL
  
  # Specify simulation parameters  
  m       <- 4
  n_dep   <- 12
  n_subj  <- 20
  n_t     <- 642
  n_iter  <- 4000
  burn_in <- 2000
  
  # 1) Simulate data
  sim_data <- mHMMbayes::sim_mHMM_plnorm(n_t = n_t, n = n_subj,
                                         data_distr = "continuous",
                                         m = m, n_dep = n_dep,
                                         gamma = gamma_map,
                                         emiss_distr = emiss_map,
                                         var_gamma = gamma_var_map,
                                         var_emiss = emiss_varmu_map,
                                         return_ind_par = TRUE)
  
  sim_data$obs[sim_data$obs < 0] <- 0
  sim_data$obs[sim_data$obs >= 100] <- 100
  sim_data$obs[sample(1:nrow(sim_data$obs), replace = FALSE, size = nrow(sim_data$obs)*0.24), -1] <- NA
  
  colnames(sim_data$obs)[-1] <- out$input$dep_labels
  
  # 2) Define hyper priors
  
  # 4) Fit MHMM
  ti <- Sys.time()
  case_out <- mHMMbayes::mHMM_cont(s_data = sim_data$obs,
                                   gen = list(m = m, n_dep = n_dep),
                                   start_val = c(list(gamma_map), emiss_map),
                                   gamma_hyp_prior = gamma_hyp_prior,
                                   emiss_hyp_prior = emiss_hyp_prior,
                                   mcmc = list(J = n_iter, burn_in = burn_in),
                                   show_progress = TRUE)
  case_out[["time"]] <- Sys.time() - ti
  
  # 8) Save output
  if(!is.null(case_out)){
    
    # Get MAP estimates
    map_out <- simHMM::MAP(case_out)
    
    # Get credibility intervals
    cci_out <- simHMM::get_cci(case_out)
    
    # Define new output object
    case_out <- list(
      truth = sim_data,
      map = map_out,
      cci = cci_out)
    
    # Save output
    saveRDS(case_out, paste0("/gpfs/home3/moragas/umcg/outputs/res/out_weakly_inf",
                             "_m",m,
                             "_it",n_iter,
                             "_burnin",burn_in,
                             "_rep",s,"_outstart_clean.rds"))
    
    # Return out
    return(case_out)
    
  }
  
}, mc.cores = 128, mc.set.seed = 1L)
