#---------------------------------------------------
# Fit model for 4 hidden states and save results
#---------------------------------------------------

# library(devtools)
# devtools::install_github("smildiner/mHMMbayes", ref = "develop")

# Preferred parallel library
library(tidyverse)
library(future.apply) # Works in both Windows and Unix
library(progressr)
library(mHMMbayes)
library(depmixS4)

# # Load utility functions
source("R/utils.R")

# load data
data <- foreign::read.spss("data/ESM bipolar cleaned data_UU project.sav", to.data.frame = TRUE)

#------------------------------------------------------------------------------#

# Put data in right format and elongate data with "missing" night occasions
extended_data <- data %>%
  group_by(patient_id) %>%
  summarise(maxdays = max(dayno)) %>%
  group_by(patient_id) %>%
  summarise(dayno = rep(1:maxdays,each=8), beepno = rep(1:8, maxdays))

train_df <- train_df_mar <- left_join(extended_data, data) %>%
  dplyr::select(patient_id, time,
                bs_diary_5, bs_diary_13, bs_diary_22,
                bs_diary_15, bs_diary_9, bs_diary_10,
                bs_diary_7, bs_diary_11, bs_diary_17,
                bs_diary_8, bs_diary_14, bs_diary_16) %>%
  group_by(patient_id) %>%
  mutate(patient_id = cur_group_id(), time = row_number()) %>%
  ungroup() %>%
  arrange(patient_id, time) %>%
  dplyr::select(-time) %>%
  as.matrix()

# Def
m <- 4
n_dep <- 12

# Get starting values from single level HMM
mod <- depmix(list(bs_diary_5 ~ 1, bs_diary_13 ~ 1, bs_diary_22 ~ 1,
                   bs_diary_15 ~ 1, bs_diary_9 ~ 1, bs_diary_10 ~ 1,
                   bs_diary_7 ~ 1, bs_diary_11 ~ 1, bs_diary_17 ~ 1,
                   bs_diary_8 ~ 1, bs_diary_14 ~ 1, bs_diary_16 ~ 1),
              data = as.data.frame(train_df_mar),
              nstates = 4,
              family = list(gaussian(), gaussian(), gaussian(),
                            gaussian(), gaussian(), gaussian(),
                            gaussian(), gaussian(), gaussian(),
                            gaussian(), gaussian(), gaussian()))
em_out <- NULL
while(is.null(em_out)){
  try({em_out <- fit(mod, verbose = 1)})
}

# Starting values for the group level transition distribution
em_start_gamma <- matrix(as.numeric(getpars(em_out)[(m+1):(m+m^2)]), nrow = m, byrow = TRUE)


# Starting values for the group level emission distribution
em_start_val_mean <- as.numeric(getpars(em_out)[(m+m^2+1):(m+m^2+m*n_dep*2)][seq(1,(m*n_dep*2),2)])
em_start_val_sd <- as.numeric(getpars(em_out)[(m+m^2+1):(m+m^2+m*n_dep*2)][seq(2,(m*n_dep*2),2)])*15

# Specify hyper-priors
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
emiss_hyp_pr <- list(
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



# Train model

# Set up cluster
plan(multisession, workers = 4)

# Fit models in parallel
with_progress({
  
  # Define progressor
  p <- progressor(along = 1:4)
  
  out_objs <- future_lapply(1:4, function(s) {
    
    # Gen
    m       <- 4
    n_dep   <- 12           # Number of dependent variables
    attempt <- 1
    out <- NULL
    
    while(attempt <= 1 & is.null(out)){
      
      cat("\nCurrently fitting model for m =",m,", this is the attempt no. ",attempt,"\n")
      
      attempt <- attempt + 1
      
      # Starting values for the group level transition distribution
      start_gamma <- int_to_prob(matrix(ifelse(as.numeric(prob_to_int(em_start_gamma)) < 0,-1,1)*
                                          runif(m*(m-1),
                                                min = abs(as.numeric(prob_to_int(em_start_gamma)))*0.9,
                                                max = abs(as.numeric(prob_to_int(em_start_gamma)))*1.1),
                                        nrow = m, ncol = m-1, byrow = FALSE))
      
      # Starting values for the group level emission distribution
      start_emiss <- vector("list", n_dep)
      
      for(q in 1:n_dep){
        start_emiss[[q]] <- matrix(c(runif(m,
                                           em_start_val_mean[(1+m*(q-1)):(m+m*(q-1))]*0.9,
                                           em_start_val_mean[(1+m*(q-1)):(m+m*(q-1))]*1.1),
                                     runif(m,
                                           em_start_val_sd[(1+m*(q-1)):(m+m*(q-1))]*0.9,
                                           em_start_val_sd[(1+m*(q-1)):(m+m*(q-1))]*1.1)),
                                   nrow = m, ncol = 2, byrow = FALSE)
      }
      
      # Fit model
      out <- mHMMbayes::mHMM_cont(s_data = train_df_mar,
                                  gen = list(m = m, n_dep = n_dep),
                                  start_val = c(list(start_gamma), start_emiss),
                                  gamma_hyp_prior = gamma_hyp_prior,
                                  emiss_hyp_prior = emiss_hyp_pr,
                                  mcmc = list(J = 4000, burn_in = 2000), # J = number of iterations
                                  show_progress = TRUE,
                                  return_path = TRUE)
      
    }
    
    class(out) <- "mHMM_cont"
    saveRDS(out, paste0("results/out_cont_gamma_prior_emiss_prior_m",m,"_12dv_it4000_c",s,"_night.rds"))
    
    p(sprintf("x=%g", s))
    
    return(out)
    
    
  }, future.seed = 42L
  )
})

# Close cluster
plan(sequential)

class(out) <- "mHMM_cont"
summary(out)


#------------------------------------------------------------------------------#
# Assess convergence:
#   - analytically with PSRF (i.e., Rhat <= 1.05)
#   - visually, looking at the trace plots


library(coda)

out <- lapply(1:4, function(c) readRDS(paste0("results/out_cont_gamma_prior_emiss_prior_m4_12dv_it4000_c",c,".rds")))

out[[1]]
out[[2]] # Best
out[[3]]
out[[4]]

m <- 4
n_dep <- 12

walk(out, summary)

# Transitions: convergence
c1_mcmc <- mcmc(out[[1]]$gamma_int_bar)
c2_mcmc <- mcmc(out[[2]]$gamma_int_bar)
c3_mcmc <- mcmc(out[[3]]$gamma_int_bar)
c4_mcmc <- mcmc(out[[4]]$gamma_int_bar)

gamma_mcmc <- mcmc.list(c1_mcmc, c2_mcmc, c3_mcmc, c4_mcmc)
# c3_mcmc, c4_mcmc)
gelman.diag(gamma_mcmc, multivariate = FALSE)

# Plot
gamma_prob_bar <- do.call(rbind, lapply(1:length(out), function(c){
  out[[c]]$gamma_int_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number(),
           chain = paste0("chain_",c)) %>%
    gather(state, value, -iter, -chain)
}))

gamma_prob_bar %>%
  ggplot(aes(x = iter, y = value, color = chain)) +
  geom_line(alpha = 0.6) +
  facet_wrap(state~., nrow = m)

gamma_prob_bar %>%
  filter(iter > 1000) %>%
  ggplot(aes(x = value, color = chain)) +
  geom_density(alpha = 0.6) +
  facet_wrap(state~., nrow = m)


# Emissions: convergence
lapply(1:n_dep, function(q) {
  c1_mcmc <- mcmc(out[[1]]$emiss_mu_bar[[q]])
  c2_mcmc <- mcmc(out[[2]]$emiss_mu_bar[[q]])
  c3_mcmc <- mcmc(out[[3]]$emiss_mu_bar[[q]])
  c4_mcmc <- mcmc(out[[4]]$emiss_mu_bar[[q]])
  
  gamma_mcmc <- mcmc.list(c1_mcmc, c2_mcmc, c3_mcmc, c4_mcmc)
  # c3_mcmc, c4_mcmc)
  gelman.diag(gamma_mcmc, multivariate = FALSE)
})

# Plot
emiss_mu_bar <- do.call(rbind, lapply(1:length(out), function(c){
  do.call(rbind, lapply(1:length(out[[c]]$emiss_mu_bar), function(q){
    out[[c]]$emiss_mu_bar[[q]] %>%
      as.data.frame() %>%
      mutate(iter = row_number(),
             dep = names(out[[c]]$emiss_mu_bar)[q],
             chain = paste0("chain_",c)) %>%
      gather(mu, value, -iter, -chain, -dep)
  }))
}))

emiss_mu_bar %>%
  ggplot(aes(x = iter, y = value, color = chain)) +
  geom_line(alpha = 0.6) +
  facet_grid(dep~mu)

emiss_mu_bar %>%
  filter(iter > 1000) %>%
  ggplot(aes(x = value, color = chain)) +
  geom_density(alpha = 0.6) +
  facet_grid(dep~mu)


# Extract emission means
emiss_mu_bar <- do.call(rbind, lapply(1:length(out), function(c){
  do.call(rbind, lapply(1:length(out[[c]]$emiss_mu_bar), function(q){
    out[[c]]$emiss_mu_bar[[q]] %>%
      as.data.frame() %>%
      mutate(iter = row_number(),
             dep = names(out[[c]]$emiss_mu_bar)[q],
             chain = paste0("chain_",c)) %>%
      gather(mu, value, -iter, -chain, -dep)
  }))
}))

# Add labels to constructs
dep_vars <- c("bs_diary_5", "bs_diary_9", "bs_diary_10",
              "bs_diary_13", "bs_diary_15", "bs_diary_22")
man_vars <- c("bs_diary_7", "bs_diary_8", "bs_diary_11",
              "bs_diary_14", "bs_diary_16", "bs_diary_17")

# Plot
emiss_mu_bar %>%
  filter(iter > 1000, chain == "chain_1") %>%
  mutate(construct = case_when(dep %in% dep_vars ~ "depression",
                               dep %in% man_vars ~ "manic"),
         dep = factor(dep, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                      "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                      "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                      "bs_diary_8", "bs_diary_14", "bs_diary_16"))) %>%
  ggplot(aes(x = value, fill = mu)) +
  geom_density(alpha = 0.4) +
  facet_wrap(construct~dep, ncol = 6, scales = "free_y")

emiss_mu_bar %>%
  filter(iter > 1000, chain == "chain_2") %>%
  mutate(construct = case_when(dep %in% dep_vars ~ "depression",
                               dep %in% man_vars ~ "manic"),
         dep = factor(dep, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                      "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                      "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                      "bs_diary_8", "bs_diary_14", "bs_diary_16"))) %>%
  ggplot(aes(x = mu, y = value, fill = mu)) +
  geom_boxplot() +
  facet_wrap(construct~dep, ncol = 6, scales = "free_y")

emiss_mu_bar %>%
  filter(iter > 1000, chain == "chain_1") %>%
  mutate(construct = case_when(dep %in% dep_vars ~ "depression",
                               dep %in% man_vars ~ "manic"),
         dep = factor(dep, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                      "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                      "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                      "bs_diary_8", "bs_diary_14", "bs_diary_16"))) %>%
  ggplot(aes(x = forcats::fct_rev(dep), y = mu, fill = value)) +
  geom_tile() +
  coord_flip()