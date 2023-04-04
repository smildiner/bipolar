#==============================================================================#
#   Analysis small decoding simulation
#==============================================================================#

# Preferred parallel library
library(tidyverse)
library(mHMMbayes)
library(simHMM)
library(pbmcapply)
library(future)
library(future.apply)
library(progressr)

source("R/utils_sim.R")

# State decoding
m <- 4
n_dep <- 12

#=====================================================================================#
# Decode
out_files <- dir("outputs/decoding_sim_weak/")[stringr::str_detect(dir("outputs/decoding_sim_weak/"), "outstart_clean.rds")]

decoding_outstart_clean <- do.call(bind_rows, pbmcapply::pbmclapply(1:length(out_files), function(rep){
  out <- readRDS(paste0("outputs/decoding_sim_weak/",out_files[rep]))
  decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = out$truth$obs, data_distr = "continuous")
  
  truth <- dplyr::bind_cols(as.data.frame(out$truth$states), as.data.frame(out$truth$obs)[,-1]) %>%
    gather(emiss, value, -subj, -state) %>%
    group_by(emiss, state) %>%
    summarise(exp_mean = mean(value, na.rm = TRUE)) %>%
    spread(key = state, value = exp_mean) %>%
    ungroup() %>%
    dplyr::select(-emiss)
  
  pred <- dplyr::bind_cols(decoded$state_seq %>%
                             as.data.frame() %>%
                             gather(subj, state), as.data.frame(out$truth$obs)[,-1]) %>%
    gather(emiss, value, -subj, -state) %>%
    group_by(emiss, state) %>%
    summarise(exp_mean = mean(value, na.rm = TRUE)) %>%
    spread(key = state, value = exp_mean) %>%
    ungroup() %>%
    dplyr::select(-emiss)
  
  out$truth$states[,2] <- as.numeric(as.character(factor(out$truth$states[,2], levels = 1:m, labels = perms)))
  
  bind_cols("rep" = rep, do.call(bind_rows, lapply(1:length(decoded$state_probs), function(s) {
    probs <- as.data.frame(decoded$state_probs[[s]])
    names(probs) <- paste0("S",1:m)
    bind_cols("id" = s,
              "true" = out$truth$states[which(out$truth$states[,1] == s),2],
              "predicted" = apply(decoded$state_probs[[s]],1,which.max),
              "fw_prob_true" = sapply(1:nrow(decoded$state_probs[[s]]),function(t) decoded$state_probs[[s]][t,out$truth$states[which(out$truth$states[,1] == s),2][t]]),
              probs)
  })))
}, mc.cores = 4))


# Main metrics
decoding_outstart_clean %>%
  mutate(true = factor(true, levels = 1:4),
         predicted = factor(predicted, levels = 1:4)) %>%
  group_by(rep, id) %>%
  summarise(acc = yardstick::accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
            bal_acc = yardstick::bal_accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
            f1 = yardstick::f_meas_vec(truth = true, estimate = predicted, na_rm = TRUE),
            kap = yardstick::kap_vec(truth = true, estimate = predicted, na_rm = TRUE)
  ) %>%
  mutate(kap = ifelse(is.nan(kap),1,kap),
         bal_acc = ifelse(is.na(bal_acc), acc, bal_acc)) %>%
  dplyr::select(-id) %>%
  group_by(rep) %>%
  summarise(across(.fns = mean)) %>%
  dplyr::select(-rep) %>%
  summarise(across(.fns = list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE), median =  ~median(., na.rm = TRUE), ci_lwr = ~quantile(., 0.025, na.rm = TRUE), ci_upr = ~quantile(., 0.975, na.rm = TRUE)),
                   .names = "{.col}.{.fn}")) %>%
  as.data.frame()

#=====================================================================================#
# Get metrics
out_files <- dir("outputs/decoding_sim_weak/")[stringr::str_detect(dir("outputs/decoding_sim_weak/"), "clean.rds")]

#=============================
# Get stating values
out <- readRDS(paste0("snellius/umcg/inputs/out_cont_gamma_prior_emiss_prior_m4_12dv_it4000_c2.rds"))

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

res <- get_all_metrics(out_files = dir("outputs/decoding_sim_weak/")[stringr::str_detect(dir("outputs/decoding_sim_weak/"), "clean.rds")],
                            path = "outputs/decoding_sim_weak/",
                            model = "cont",
                            gamma_map = gamma_map, gamma_var = gamma_var_map,
                            emiss_map = emiss_map, emiss_var = emiss_vemiss_varmu_maparmu_ppc,
                            gamma_var_idx = which(colnames(out$gamma_V_int_bar[2001:4000,]) %in% paste0("var_int_S",rep(1:4,each=3),"toS",2:4,"_with_int_S",rep(1:4,each=3),"toS",2:4)))


#=====================================================================================#
# Plots

## Transitions group-level means (fixed effects) Prob scale
# Bias
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = (m-1)),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = bias,
                      group = scenario, colour = scenario,
                      ymin=bias-bias_mcmc_se, ymax=bias+bias_mcmc_se),
                  position=position_dodge(width=0.25)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  coord_flip() +
  theme_minimal()

# Relative bias
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = (m-1)),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = rel_bias,
                      group = scenario, colour = scenario,
                      ymin=rel_bias-rel_bias_mcmc_se, ymax=rel_bias+rel_bias_mcmc_se),
                  position=position_dodge(width=0.25)
  ) +
  geom_hline(yintercept = 0, colour = "dark red") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed") +
  coord_flip(ylim = c(-1,2)) +
  theme_minimal()

# Empirical SE (precision)
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = (m-1)),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = empirical_se,
                      group = scenario, colour = scenario,
                      ymin=empirical_se-empirical_se_mcmc_se, ymax=empirical_se+empirical_se_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  theme_minimal()

# MSE
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = (m-1)),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = MSE,
                      group = scenario, colour = scenario,
                      ymin=MSE-MSE_mcmc_se, ymax=MSE+MSE_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  theme_minimal()

# Coverage
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = (m-1)),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = coverage,
                      group = scenario, colour = scenario,
                      ymin=coverage-coverage_mcmc_se, ymax=coverage+coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  coord_flip() +
  theme_minimal()

# Bias-corrected coverage
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = (m-1)),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = bias_corr_coverage,
                      group = scenario, colour = scenario,
                      ymin=bias_corr_coverage-bias_corr_coverage_mcmc_se, ymax=bias_corr_coverage+bias_corr_coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  coord_flip() +
  theme_minimal()


## Transitions group-level means (fixed effects) Logit scale
# Bias
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = bias,
                      group = scenario, colour = scenario,
                      ymin=bias-bias_mcmc_se, ymax=bias+bias_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  coord_flip() +
  theme_minimal()

# Relative bias
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = rel_bias,
                      group = scenario, colour = scenario,
                      ymin=rel_bias-rel_bias_mcmc_se, ymax=rel_bias+rel_bias_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed") +
  coord_flip() +
  theme_minimal()

# Empirical SE (precision)
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = empirical_se,
                      group = scenario, colour = scenario,
                      ymin=empirical_se-empirical_se_mcmc_se, ymax=empirical_se+empirical_se_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  theme_minimal()

# MSE
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = MSE,
                      group = scenario, colour = scenario,
                      ymin=MSE-MSE_mcmc_se, ymax=MSE+MSE_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  theme_minimal()

# Coverage
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = coverage,
                      group = scenario, colour = scenario,
                      ymin=coverage-coverage_mcmc_se, ymax=coverage+coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  coord_flip() +
  theme_minimal()

# Bias-corrected coverage
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = bias_corr_coverage,
                      group = scenario, colour = scenario,
                      ymin=bias_corr_coverage-bias_corr_coverage_mcmc_se, ymax=bias_corr_coverage+bias_corr_coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  coord_flip() +
  theme_minimal()


#-----------------------------------------------------------------------------#
## Transitions random effects (fixed effects) Logit scale
# Bias
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = bias,
                      group = scenario, colour = scenario,
                      ymin=bias-bias_mcmc_se, ymax=bias+bias_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  coord_flip() +
  theme_minimal()

# Relative bias
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = rel_bias,
                      group = scenario, colour = scenario,
                      ymin=rel_bias-rel_bias_mcmc_se, ymax=rel_bias+rel_bias_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed") +
  coord_flip() +
  theme_minimal()

# Empirical SE (precision)
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = empirical_se,
                      group = scenario, colour = scenario,
                      ymin=empirical_se-empirical_se_mcmc_se, ymax=empirical_se+empirical_se_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  theme_minimal()

# MSE
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = MSE,
                      group = scenario, colour = scenario,
                      ymin=MSE-MSE_mcmc_se, ymax=MSE+MSE_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  theme_minimal()

# Coverage
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = coverage,
                      group = scenario, colour = scenario,
                      ymin=coverage-coverage_mcmc_se, ymax=coverage+coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  coord_flip() +
  theme_minimal()

# Bias-corrected coverage
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = parameter, y = bias_corr_coverage,
                      group = scenario, colour = scenario,
                      ymin=bias_corr_coverage-bias_corr_coverage_mcmc_se, ymax=bias_corr_coverage+bias_corr_coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  coord_flip() +
  theme_minimal()


#-----------------------------------------------------------------------------#
## Emissions group-level means (fixed effects)
# Bias
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  # mutate(ndep = factor(ndep, levels = paste0("ndep",1:12))) %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = bias,
                      group = scenario, colour = scenario,
                      ymin=bias-bias_mcmc_se, ymax=bias+bias_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  facet_wrap(ndep~., ncol = 4)  +
  theme_minimal()

# Relative bias
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = rel_bias,
                      group = scenario, colour = scenario,
                      ymin=rel_bias-rel_bias_mcmc_se, ymax=rel_bias+rel_bias_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed") +
  coord_flip() +
  facet_wrap(ndep~., nrow = 2) +
  theme_minimal()

# Empirical SE (precision)
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = empirical_se,
                      group = scenario, colour = scenario,
                      ymin=empirical_se-empirical_se_mcmc_se, ymax=empirical_se+empirical_se_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  facet_wrap(ndep~., nrow = 2) +
  theme_minimal()

# MSE
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = MSE,
                      group = scenario, colour = scenario,
                      ymin=MSE-MSE_mcmc_se, ymax=MSE+MSE_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  facet_wrap(ndep~., nrow = 2) +
  theme_minimal()


# Coverage
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = coverage,
                      group = scenario, colour = scenario,
                      ymin=coverage-coverage_mcmc_se, ymax=coverage+coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  coord_flip() +
  facet_wrap(ndep~., nrow = 2) +
  theme_minimal()

# Bias-corrected coverage
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = bias_corr_coverage,
                      group = scenario, colour = scenario,
                      ymin=bias_corr_coverage-bias_corr_coverage_mcmc_se, ymax=bias_corr_coverage+bias_corr_coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), colour = "grey", linetype = "dashed") +
  coord_flip() +
  facet_wrap(ndep~., nrow = 2) +
  theme_minimal()



