#==============================================================================#
#   Analysis small decoding simulation
#==============================================================================#

# Preferred parallel library
library(tidyverse)
library(mHMMbayes)
library(simHMM)
library(pbmcapply)

source("R/utils_sim.R")

# State decoding
m <- 4
n_dep <- 12



#==============================================================================#
# Decode simulated data

# Raw
out_files <- sample(dir("outputs/res/")[stringr::str_detect(dir("outputs/res/"), "raw")], 100, replace = FALSE)

decoding_raw_full <- do.call(rbind, pbmcapply::pbmclapply(1:length(out_files), function(rep){
  out <- readRDS(paste0("outputs/res/",out_files[rep]))
  m <- nrow(out$truth$subject_gamma[[1]])
  
  sim_data <- out$truth
  
  train_idx <- rep(c(rep(TRUE, 960), rep(FALSE, 960)), 20)
  test_idx <- !train_idx
  
  train_data <- sim_data$obs[train_idx,]
  test_data <- sim_data$obs[test_idx,]
  
  night_occasions <- as.numeric(t(sapply(6:8, function(o) seq(0, 120*8*20-8, 8)+o)))
  
  test_data_night <- test_data
  test_data_night[night_occasions, -1] <- NA
  test_data_clean <- test_data[-night_occasions,]
  
  train_raw_data <- train_data
  train_night_data <- train_raw_data
  train_night_data[night_occasions, -1] <- NA
  train_clean_data <- train_raw_data[-night_occasions,]
  
  decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = train_data, data_distr = "continuous")
  bind_cols("rep" = rep, do.call(bind_rows, lapply(1:length(decoded$state_probs), function(s) {
    probs <- as.data.frame(decoded$state_probs[[s]])
    names(probs) <- paste0("S",1:m)
    bind_cols("id" = s,
              "true" = out$truth$states[which(out$truth$states[,1] == s),2][1:960],
              "predicted" = apply(decoded$state_probs[[s]],1,which.max),
              "fw_prob_true" = sapply(1:nrow(decoded$state_probs[[s]]),function(t) decoded$state_probs[[s]][t,out$truth$states[which(out$truth$states[,1] == s),2][t]]),
              probs)
  })))
  
}, mc.cores = 25)) %>%
  mutate("train" = "full", "test" = "full") %>%
  mutate(true = factor(true, levels = 1:4),
         predicted = factor(predicted, levels = 1:4)) %>%
  group_by(train, test, rep, id) %>%
  summarise(
    bal_acc = yardstick::bal_accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
    f_meas = yardstick::f_meas_vec(truth = true, estimate = predicted, na_rm = TRUE),
    kappa = yardstick::kap_vec(truth = true, estimate = predicted, na_rm = TRUE)
  ) %>%
  group_by(test, train, rep) %>%
  summarise(across(.fns = mean, na.rm=TRUE)) %>%
  group_by(test, train) %>%
  dplyr::select(-rep,-id) %>%
  summarise(across(.fns = list( ~mean(., na.rm=TRUE), ~sd(., na.rm=TRUE), ~quantile(., 0.025, na.rm = TRUE), ~quantile(., 0.975, na.rm = TRUE))))

decoding_raw_night <- do.call(rbind, pbmcapply::pbmclapply(1:length(out_files), function(rep){
  out <- readRDS(paste0("outputs/res/",out_files[rep]))
  m <- nrow(out$truth$subject_gamma[[1]])
  
  sim_data <- out$truth
  
  train_idx <- rep(c(rep(TRUE, 960), rep(FALSE, 960)), 20)
  test_idx <- !train_idx
  
  train_data <- sim_data$obs[train_idx,]
  test_data <- sim_data$obs[test_idx,]
  
  night_occasions <- as.numeric(t(sapply(6:8, function(o) seq(0, 120*8*20-8, 8)+o)))
  
  test_data_night <- test_data
  test_data_night[night_occasions, -1] <- NA
  test_data_clean <- test_data[-night_occasions,]
  
  train_raw_data <- train_data
  train_night_data <- train_raw_data
  train_night_data[night_occasions, -1] <- NA
  train_clean_data <- train_raw_data[-night_occasions,]
  
  decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = train_night_data, data_distr = "continuous")
  bind_cols("rep" = rep, do.call(bind_rows, lapply(1:length(decoded$state_probs), function(s) {
    probs <- as.data.frame(decoded$state_probs[[s]])
    names(probs) <- paste0("S",1:m)
    bind_cols("id" = s,
              "true" = out$truth$states[which(out$truth$states[,1] == s),2][1:960],
              "predicted" = apply(decoded$state_probs[[s]],1,which.max),
              "fw_prob_true" = sapply(1:nrow(decoded$state_probs[[s]]),function(t) decoded$state_probs[[s]][t,out$truth$states[which(out$truth$states[,1] == s),2][t]]),
              probs)
  })))
  
}, mc.cores = 25)) %>%
  mutate("train" = "full", "test" = "night") %>%
  mutate(true = factor(true, levels = 1:4),
         predicted = factor(predicted, levels = 1:4)) %>%
  group_by(train, test, rep, id) %>%
  summarise(
    bal_acc = yardstick::bal_accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
    f_meas = yardstick::f_meas_vec(truth = true, estimate = predicted, na_rm = TRUE),
    kappa = yardstick::kap_vec(truth = true, estimate = predicted, na_rm = TRUE)
  ) %>%
  group_by(test, train, rep) %>%
  summarise(across(.fns = mean, na.rm=TRUE)) %>%
  group_by(test, train) %>%
  dplyr::select(-rep,-id) %>%
  summarise(across(.fns = list( ~mean(., na.rm=TRUE), ~sd(., na.rm=TRUE), ~quantile(., 0.025, na.rm = TRUE), ~quantile(., 0.975, na.rm = TRUE))))

decoding_raw_omit <- do.call(rbind, pbmcapply::pbmclapply(1:length(out_files), function(rep){
  out <- readRDS(paste0("outputs/res/",out_files[rep]))
  m <- nrow(out$truth$subject_gamma[[1]])
  
  sim_data <- out$truth
  
  train_idx <- rep(c(rep(TRUE, 960), rep(FALSE, 960)), 20)
  test_idx <- !train_idx
  
  train_data <- sim_data$obs[train_idx,]
  test_data <- sim_data$obs[test_idx,]
  
  night_occasions <- as.numeric(t(sapply(6:8, function(o) seq(0, 120*8*20-8, 8)+o)))
  
  test_data_night <- test_data
  test_data_night[night_occasions, -1] <- NA
  test_data_clean <- test_data[-night_occasions,]
  
  train_raw_data <- train_data
  train_night_data <- train_raw_data
  train_night_data[night_occasions, -1] <- NA
  train_clean_data <- train_raw_data[-night_occasions,]
  
  decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = train_clean_data, data_distr = "continuous")
  bind_cols("rep" = rep, do.call(bind_rows, lapply(1:length(decoded$state_probs), function(s) {
    probs <- as.data.frame(decoded$state_probs[[s]])
    names(probs) <- paste0("S",1:m)
    bind_cols("id" = s,
              "true" = out$truth$states[which(out$truth$states[,1] == s),2][1:960][-night_occasions],
              "predicted" = apply(decoded$state_probs[[s]],1,which.max),
              "fw_prob_true" = sapply(1:nrow(decoded$state_probs[[s]]),function(t) decoded$state_probs[[s]][t,out$truth$states[which(out$truth$states[,1] == s),2][t]]),
              probs)
  })))
  
}, mc.cores = 25)) %>%
  mutate("train" = "full", "test" = "omit") %>%
  mutate(true = factor(true, levels = 1:4),
         predicted = factor(predicted, levels = 1:4)) %>%
  group_by(train, test, rep, id) %>%
  summarise(
    bal_acc = yardstick::bal_accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
    f_meas = yardstick::f_meas_vec(truth = true, estimate = predicted, na_rm = TRUE),
    kappa = yardstick::kap_vec(truth = true, estimate = predicted, na_rm = TRUE)
  ) %>%
  group_by(test, train, rep) %>%
  summarise(across(.fns = mean, na.rm=TRUE)) %>%
  group_by(test, train) %>%
  dplyr::select(-rep,-id) %>%
  summarise(across(.fns = list( ~mean(., na.rm=TRUE), ~sd(., na.rm=TRUE), ~quantile(., 0.025, na.rm = TRUE), ~quantile(., 0.975, na.rm = TRUE))))

# Imputing night
out_files <- sample(dir("outputs/res/")[stringr::str_detect(dir("outputs/res/"), "imp")], 100, replace = FALSE)

decoding_imp_full <- do.call(rbind, pbmcapply::pbmclapply(1:length(out_files), function(rep){
  out <- readRDS(paste0("outputs/res/",out_files[rep]))
  m <- nrow(out$truth$subject_gamma[[1]])
  
  sim_data <- out$truth
  
  train_idx <- rep(c(rep(TRUE, 960), rep(FALSE, 960)), 20)
  test_idx <- !train_idx
  
  train_data <- sim_data$obs[train_idx,]
  test_data <- sim_data$obs[test_idx,]
  
  night_occasions <- as.numeric(t(sapply(6:8, function(o) seq(0, 120*8*20-8, 8)+o)))
  
  test_data_night <- test_data
  test_data_night[night_occasions, -1] <- NA
  test_data_clean <- test_data[-night_occasions,]
  
  train_raw_data <- train_data
  train_night_data <- train_raw_data
  train_night_data[night_occasions, -1] <- NA
  train_clean_data <- train_raw_data[-night_occasions,]
  
  decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = train_data, data_distr = "continuous")
  bind_cols("rep" = rep, do.call(bind_rows, lapply(1:length(decoded$state_probs), function(s) {
    probs <- as.data.frame(decoded$state_probs[[s]])
    names(probs) <- paste0("S",1:m)
    bind_cols("id" = s,
              "true" = out$truth$states[which(out$truth$states[,1] == s),2][1:960],
              "predicted" = apply(decoded$state_probs[[s]],1,which.max),
              "fw_prob_true" = sapply(1:nrow(decoded$state_probs[[s]]),function(t) decoded$state_probs[[s]][t,out$truth$states[which(out$truth$states[,1] == s),2][t]]),
              probs)
  })))
  
}, mc.cores = 25)) %>%
  mutate("train" = "night", "test" = "full") %>%
  mutate(true = factor(true, levels = 1:4),
         predicted = factor(predicted, levels = 1:4)) %>%
  group_by(train, test, rep, id) %>%
  summarise(
    bal_acc = yardstick::bal_accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
    f_meas = yardstick::f_meas_vec(truth = true, estimate = predicted, na_rm = TRUE),
    kappa = yardstick::kap_vec(truth = true, estimate = predicted, na_rm = TRUE)
  ) %>%
  group_by(test, train, rep) %>%
  summarise(across(.fns = mean, na.rm=TRUE)) %>%
  group_by(test, train) %>%
  dplyr::select(-rep,-id) %>%
  summarise(across(.fns = list( ~mean(., na.rm=TRUE), ~sd(., na.rm=TRUE), ~quantile(., 0.025, na.rm = TRUE), ~quantile(., 0.975, na.rm = TRUE))))

decoding_imp_night <- do.call(rbind, pbmcapply::pbmclapply(1:length(out_files), function(rep){
  out <- readRDS(paste0("outputs/res/",out_files[rep]))
  m <- nrow(out$truth$subject_gamma[[1]])
  
  sim_data <- out$truth
  
  train_idx <- rep(c(rep(TRUE, 960), rep(FALSE, 960)), 20)
  test_idx <- !train_idx
  
  train_data <- sim_data$obs[train_idx,]
  test_data <- sim_data$obs[test_idx,]
  
  night_occasions <- as.numeric(t(sapply(6:8, function(o) seq(0, 120*8*20-8, 8)+o)))
  
  test_data_night <- test_data
  test_data_night[night_occasions, -1] <- NA
  test_data_clean <- test_data[-night_occasions,]
  
  train_raw_data <- train_data
  train_night_data <- train_raw_data
  train_night_data[night_occasions, -1] <- NA
  train_clean_data <- train_raw_data[-night_occasions,]
  
  decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = train_night_data, data_distr = "continuous")
  bind_cols("rep" = rep, do.call(bind_rows, lapply(1:length(decoded$state_probs), function(s) {
    probs <- as.data.frame(decoded$state_probs[[s]])
    names(probs) <- paste0("S",1:m)
    bind_cols("id" = s,
              "true" = out$truth$states[which(out$truth$states[,1] == s),2][1:960],
              "predicted" = apply(decoded$state_probs[[s]],1,which.max),
              "fw_prob_true" = sapply(1:nrow(decoded$state_probs[[s]]),function(t) decoded$state_probs[[s]][t,out$truth$states[which(out$truth$states[,1] == s),2][t]]),
              probs)
  })))
  
}, mc.cores = 25)) %>%
  mutate("train" = "night", "test" = "night") %>%
  mutate(true = factor(true, levels = 1:4),
         predicted = factor(predicted, levels = 1:4)) %>%
  group_by(train, test, rep, id) %>%
  summarise(
    bal_acc = yardstick::bal_accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
    f_meas = yardstick::f_meas_vec(truth = true, estimate = predicted, na_rm = TRUE),
    kappa = yardstick::kap_vec(truth = true, estimate = predicted, na_rm = TRUE)
  ) %>%
  group_by(test, train, rep) %>%
  summarise(across(.fns = mean, na.rm=TRUE)) %>%
  group_by(test, train) %>%
  dplyr::select(-rep,-id) %>%
  summarise(across(.fns = list( ~mean(., na.rm=TRUE), ~sd(., na.rm=TRUE), ~quantile(., 0.025, na.rm = TRUE), ~quantile(., 0.975, na.rm = TRUE))))

decoding_imp_omit <- do.call(rbind, pbmcapply::pbmclapply(1:length(out_files), function(rep){
  out <- readRDS(paste0("outputs/res/",out_files[rep]))
  m <- nrow(out$truth$subject_gamma[[1]])
  
  sim_data <- out$truth
  
  train_idx <- rep(c(rep(TRUE, 960), rep(FALSE, 960)), 20)
  test_idx <- !train_idx
  
  train_data <- sim_data$obs[train_idx,]
  test_data <- sim_data$obs[test_idx,]
  
  night_occasions <- as.numeric(t(sapply(6:8, function(o) seq(0, 120*8*20-8, 8)+o)))
  
  test_data_night <- test_data
  test_data_night[night_occasions, -1] <- NA
  test_data_clean <- test_data[-night_occasions,]
  
  train_raw_data <- train_data
  train_night_data <- train_raw_data
  train_night_data[night_occasions, -1] <- NA
  train_clean_data <- train_raw_data[-night_occasions,]
  
  decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = train_clean_data, data_distr = "continuous")
  bind_cols("rep" = rep, do.call(bind_rows, lapply(1:length(decoded$state_probs), function(s) {
    probs <- as.data.frame(decoded$state_probs[[s]])
    names(probs) <- paste0("S",1:m)
    bind_cols("id" = s,
              "true" = out$truth$states[which(out$truth$states[,1] == s),2][1:960][-night_occasions],
              "predicted" = apply(decoded$state_probs[[s]],1,which.max),
              "fw_prob_true" = sapply(1:nrow(decoded$state_probs[[s]]),function(t) decoded$state_probs[[s]][t,out$truth$states[which(out$truth$states[,1] == s),2][t]]),
              probs)
  })))
  
}, mc.cores = 25)) %>%
  mutate("train" = "night", "test" = "omit") %>%
  mutate(true = factor(true, levels = 1:4),
         predicted = factor(predicted, levels = 1:4)) %>%
  group_by(train, test, rep, id) %>%
  summarise(
    bal_acc = yardstick::bal_accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
    f_meas = yardstick::f_meas_vec(truth = true, estimate = predicted, na_rm = TRUE),
    kappa = yardstick::kap_vec(truth = true, estimate = predicted, na_rm = TRUE)
  ) %>%
  group_by(test, train, rep) %>%
  summarise(across(.fns = mean, na.rm=TRUE)) %>%
  group_by(test, train) %>%
  dplyr::select(-rep,-id) %>%
  summarise(across(.fns = list( ~mean(., na.rm=TRUE), ~sd(., na.rm=TRUE), ~quantile(., 0.025, na.rm = TRUE), ~quantile(., 0.975, na.rm = TRUE))))


# Omitting night
out_files <- sample(dir("outputs/res/")[stringr::str_detect(dir("outputs/res/"), "omit")], 100, replace = FALSE)

decoding_omit_full <- do.call(rbind, pbmcapply::pbmclapply(1:length(out_files), function(rep){
  out <- readRDS(paste0("outputs/res/",out_files[rep]))
  m <- nrow(out$truth$subject_gamma[[1]])
  
  sim_data <- out$truth
  
  train_idx <- rep(c(rep(TRUE, 960), rep(FALSE, 960)), 20)
  test_idx <- !train_idx
  
  train_data <- sim_data$obs[train_idx,]
  test_data <- sim_data$obs[test_idx,]
  
  night_occasions <- as.numeric(t(sapply(6:8, function(o) seq(0, 120*8*20-8, 8)+o)))
  
  test_data_night <- test_data
  test_data_night[night_occasions, -1] <- NA
  test_data_clean <- test_data[-night_occasions,]
  
  train_raw_data <- train_data
  train_night_data <- train_raw_data
  train_night_data[night_occasions, -1] <- NA
  train_clean_data <- train_raw_data[-night_occasions,]
  
  decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = train_data, data_distr = "continuous")
  bind_cols("rep" = rep, do.call(bind_rows, lapply(1:length(decoded$state_probs), function(s) {
    probs <- as.data.frame(decoded$state_probs[[s]])
    names(probs) <- paste0("S",1:m)
    bind_cols("id" = s,
              "true" = out$truth$states[which(out$truth$states[,1] == s),2][1:960],
              "predicted" = apply(decoded$state_probs[[s]],1,which.max),
              "fw_prob_true" = sapply(1:nrow(decoded$state_probs[[s]]),function(t) decoded$state_probs[[s]][t,out$truth$states[which(out$truth$states[,1] == s),2][t]]),
              probs)
  })))
  
}, mc.cores = 25)) %>%
  mutate("train" = "omit", "test" = "full") %>%
  mutate(true = factor(true, levels = 1:4),
         predicted = factor(predicted, levels = 1:4)) %>%
  group_by(train, test, rep, id) %>%
  summarise(
    bal_acc = yardstick::bal_accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
    f_meas = yardstick::f_meas_vec(truth = true, estimate = predicted, na_rm = TRUE),
    kappa = yardstick::kap_vec(truth = true, estimate = predicted, na_rm = TRUE)
  ) %>%
  group_by(test, train, rep) %>%
  summarise(across(.fns = mean, na.rm=TRUE)) %>%
  group_by(test, train) %>%
  dplyr::select(-rep,-id) %>%
  summarise(across(.fns = list( ~mean(., na.rm=TRUE), ~sd(., na.rm=TRUE), ~quantile(., 0.025, na.rm = TRUE), ~quantile(., 0.975, na.rm = TRUE))))

decoding_omit_night <- do.call(rbind, pbmcapply::pbmclapply(1:length(out_files), function(rep){
  out <- readRDS(paste0("outputs/res/",out_files[rep]))
  m <- nrow(out$truth$subject_gamma[[1]])
  
  sim_data <- out$truth
  
  train_idx <- rep(c(rep(TRUE, 960), rep(FALSE, 960)), 20)
  test_idx <- !train_idx
  
  train_data <- sim_data$obs[train_idx,]
  test_data <- sim_data$obs[test_idx,]
  
  night_occasions <- as.numeric(t(sapply(6:8, function(o) seq(0, 120*8*20-8, 8)+o)))
  
  test_data_night <- test_data
  test_data_night[night_occasions, -1] <- NA
  test_data_clean <- test_data[-night_occasions,]
  
  train_raw_data <- train_data
  train_night_data <- train_raw_data
  train_night_data[night_occasions, -1] <- NA
  train_clean_data <- train_raw_data[-night_occasions,]
  
  # decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = out$truth$obs[train_idx,], data_distr = "continuous")
  decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = train_night_data, data_distr = "continuous")
  bind_cols("rep" = rep, do.call(bind_rows, lapply(1:length(decoded$state_probs), function(s) {
    probs <- as.data.frame(decoded$state_probs[[s]])
    names(probs) <- paste0("S",1:m)
    bind_cols("id" = s,
              "true" = out$truth$states[which(out$truth$states[,1] == s),2][1:960],
              "predicted" = apply(decoded$state_probs[[s]],1,which.max),
              "fw_prob_true" = sapply(1:nrow(decoded$state_probs[[s]]),function(t) decoded$state_probs[[s]][t,out$truth$states[which(out$truth$states[,1] == s),2][t]]),
              probs)
  })))
  
}, mc.cores = 25)) %>%
  mutate("train" = "omit", "test" = "night") %>%
  mutate(true = factor(true, levels = 1:4),
         predicted = factor(predicted, levels = 1:4)) %>%
  group_by(train, test, rep, id) %>%
  summarise(
    bal_acc = yardstick::bal_accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
    f_meas = yardstick::f_meas_vec(truth = true, estimate = predicted, na_rm = TRUE),
    kappa = yardstick::kap_vec(truth = true, estimate = predicted, na_rm = TRUE)
  ) %>%
  group_by(test, train, rep) %>%
  summarise(across(.fns = mean, na.rm=TRUE)) %>%
  group_by(test, train) %>%
  dplyr::select(-rep,-id) %>%
  summarise(across(.fns = list( ~mean(., na.rm=TRUE), ~sd(., na.rm=TRUE), ~quantile(., 0.025, na.rm = TRUE), ~quantile(., 0.975, na.rm = TRUE))))

decoding_omit_omit <- do.call(rbind, pbmcapply::pbmclapply(1:length(out_files), function(rep){
  out <- readRDS(paste0("outputs/res/",out_files[rep]))
  m <- nrow(out$truth$subject_gamma[[1]])
  
  sim_data <- out$truth
  
  train_idx <- rep(c(rep(TRUE, 960), rep(FALSE, 960)), 20)
  test_idx <- !train_idx
  
  train_data <- sim_data$obs[train_idx,]
  test_data <- sim_data$obs[test_idx,]
  
  night_occasions <- as.numeric(t(sapply(6:8, function(o) seq(0, 120*8*20-8, 8)+o)))
  
  test_data_night <- test_data
  test_data_night[night_occasions, -1] <- NA
  test_data_clean <- test_data[-night_occasions,]
  
  train_raw_data <- train_data
  train_night_data <- train_raw_data
  train_night_data[night_occasions, -1] <- NA
  train_clean_data <- train_raw_data[-night_occasions,]
  
  decoded <- simHMM::vit_mHMM_map(object = out$map, s_data = train_clean_data, data_distr = "continuous")
  bind_cols("rep" = rep, do.call(bind_rows, lapply(1:length(decoded$state_probs), function(s) {
    probs <- as.data.frame(decoded$state_probs[[s]])
    names(probs) <- paste0("S",1:m)
    bind_cols("id" = s,
              "true" = out$truth$states[which(out$truth$states[,1] == s),2][1:960][-night_occasions],
              "predicted" = apply(decoded$state_probs[[s]],1,which.max),
              "fw_prob_true" = sapply(1:nrow(decoded$state_probs[[s]]),function(t) decoded$state_probs[[s]][t,out$truth$states[which(out$truth$states[,1] == s),2][t]]),
              probs)
  })))
  
}, mc.cores = 25)) %>%
  mutate("train" = "omit", "test" = "omit") %>%
  mutate(true = factor(true, levels = 1:4),
         predicted = factor(predicted, levels = 1:4)) %>%
  group_by(train, test, rep, id) %>%
  summarise(
    bal_acc = yardstick::bal_accuracy_vec(truth = true, estimate = predicted, na_rm = TRUE),
    f_meas = yardstick::f_meas_vec(truth = true, estimate = predicted, na_rm = TRUE),
    kappa = yardstick::kap_vec(truth = true, estimate = predicted, na_rm = TRUE)
  ) %>%
  group_by(test, train, rep) %>%
  summarise(across(.fns = mean, na.rm=TRUE)) %>%
  group_by(test, train) %>%
  dplyr::select(-rep,-id) %>%
  summarise(across(.fns = list( ~mean(., na.rm=TRUE), ~sd(., na.rm=TRUE), ~quantile(., 0.025, na.rm = TRUE), ~quantile(., 0.975, na.rm = TRUE))))

#==============================================================================#
## Put together results from decoding
decoding_results <- rbind(
  decoding_raw_night,
  decoding_imp_night,
  decoding_omit_night,
  decoding_raw_full,
  decoding_imp_full,
  decoding_omit_full,
  decoding_raw_omit,
  decoding_imp_omit,
  decoding_omit_omit)







#=====================================================================================#
# Get parameter estimation metrics

# Define path
path <- "outputs/res/"

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

gamma_var_idx <- which(colnames(out$gamma_V_int_bar) %in% paste0("var_int_S",rep(1:4,each=3),"toS",2:4,"_with_int_S",rep(1:4,each=3),"toS",2:4))

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

# Gamma var index
# out <- readRDS(paste0(path,out_files[[1]]))
# gamma_var_idx <- which(colnames(out$gamma_V_int_bar) %in% paste0("var_int_S",rep(1:m,each=m-1),"toS",2:m,"_with_int_S",rep(1:m,each=m-1),"toS",2:m))
# [1]  1  5  9 10 14 18 19 23 27 28 32 36
# gamma_var_idx <- c(1, 5, 9, 10, 14, 18, 19, 23, 27, 28, 32, 36)

# Get evaluation metrics for MHMM
set.seed(42)
# out_files <- sample(dir("outputs/res/")[stringr::str_detect(dir("outputs/res/"), "raw")], 100, replace = FALSE)

res_raw <- get_all_metrics(out_files = dir(path)[stringr::str_detect(dir(path), "raw")],
                           path = path,
                           # model = "mhmm_plnorm",
                           gamma_map = gamma_map, gamma_var = gamma_var_map,
                           emiss_map = emiss_map, emiss_var = emiss_varmu_map,
                           gamma_var_idx = gamma_var_idx) %>%
  mutate(train = "full")

res_imp <- get_all_metrics(out_files = dir(path)[stringr::str_detect(dir(path), "imp")],
                           path = path,
                           # model = "mhmm_plnorm",
                           gamma_map = gamma_map, gamma_var = gamma_var_map,
                           emiss_map = emiss_map, emiss_var = emiss_varmu_map,
                           gamma_var_idx = gamma_var_idx) %>%
  mutate(train = "impute")

res_omit <- get_all_metrics(out_files = dir(path)[stringr::str_detect(dir(path), "omit")],
                            path = path,
                            # model = "mhmm_plnorm",
                            gamma_map = gamma_map, gamma_var = gamma_var_map,
                            emiss_map = emiss_map, emiss_var = emiss_varmu_map,
                            gamma_var_idx = gamma_var_idx) %>%
  mutate(train = "omit")



# Save files:
res <- bind_rows(res_raw, res_imp, res_omit)
res <- rename(res, "model"='train')


#=====================================================================================#
# Plots


# Transitions probs fixed
metrics_data <- bind_cols(res %>%
                            dplyr::select(parameter, bias, MSE, coverage, model, condition) %>%
                            filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) %>%
                            gather(metric, value, -parameter, -model, -condition),
                          res %>%
                            dplyr::select(parameter, bias_mcmc_se, MSE_mcmc_se, coverage_mcmc_se, model, condition) %>%
                            filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) %>%
                            gather(metric_se, value_se, -parameter, -model, -condition) %>%
                            dplyr::select(metric_se, value_se)) %>%
  mutate(metric = factor(metric, levels = c("bias","MSE","coverage"), labels = c("Bias", "MSE", "Coverage")))

ggplot(data = metrics_data) +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = value,
                      colour = model,
                      shape = model,
                      ymin=value-value_se, ymax=value+value_se),
                  position=position_dodge(width=0.25), size = 0.3) +
  geom_hline(data = filter(metrics_data, metric == "Coverage"), aes(yintercept = 0.95), colour = "dark red") +
  geom_hline(data = filter(metrics_data, metric != "Coverage"), aes(yintercept = 0.0), colour = "dark red") +
  geom_hline(data = filter(metrics_data, metric == "Coverage"), aes(yintercept = 0.925), linetype = "dashed") +
  geom_hline(data = filter(metrics_data, metric == "Coverage"), aes(yintercept = 0.975), linetype = "dashed") +
  facet_grid(condition~metric, scales = "free_x") +
  coord_flip() +
  # ggtitle("Bias in group-level transitions") +
  xlab("") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_rect(fill="white", size = 0))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1))  

# Transitions logit fixed
metrics_data <- bind_cols(res %>%
                            dplyr::select(parameter, bias, MSE, coverage, model, condition) %>%
                            filter(parameter %in% paste0("int_S",rep(1:m, each = m),"toS",1:m)) %>%
                            gather(metric, value, -parameter, -model, -condition),
                          res %>%
                            dplyr::select(parameter, bias_mcmc_se, MSE_mcmc_se, coverage_mcmc_se, model, condition) %>%
                            filter(parameter %in% paste0("int_S",rep(1:m, each = m),"toS",1:m)) %>%
                            gather(metric_se, value_se, -parameter, -model, -condition) %>%
                            dplyr::select(metric_se, value_se)) %>%
  mutate(metric = factor(metric, levels = c("bias","MSE","coverage"), labels = c("Bias", "MSE", "Coverage")))
ggplot(data = metrics_data) +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = value,
                      colour = model,
                      shape = model,
                      ymin=value-value_se, ymax=value+value_se),
                  position=position_dodge(width=0.25), size = 0.3) +
  geom_hline(data = filter(metrics_data, metric == "Coverage"), aes(yintercept = 0.95), colour = "dark red") +
  geom_hline(data = filter(metrics_data, metric != "Coverage"), aes(yintercept = 0.0), colour = "dark red") +
  geom_hline(data = filter(metrics_data, metric == "Coverage"), aes(yintercept = 0.925), linetype = "dashed") +
  geom_hline(data = filter(metrics_data, metric == "Coverage"), aes(yintercept = 0.975), linetype = "dashed") +
  facet_grid(condition~metric, scales = "free_x") +
  coord_flip() +
  xlab("") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_rect(fill="white", size = 0))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1))  


# Emissions fixed
metrics_data <- bind_cols(res %>%
                            dplyr::select(parameter, bias, MSE, coverage, model, condition) %>%
                            filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
                            separate(parameter, c("ndep","mu"),"_") %>%
                            gather(metric, value, -mu, -model, -condition, -n_dep),
                          res %>%
                            dplyr::select(parameter, bias_mcmc_se, MSE_mcmc_se, coverage_mcmc_se, model, condition) %>%
                            filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
                            separate(parameter, c("ndep","mu"),"_") %>%
                            gather(metric_se, value_se, -mu, -model, -condition, -n_dep) %>%
                            dplyr::select(metric_se, value_se)) %>%
  mutate(metric = factor(metric, levels = c("bias","MSE","coverage"), labels = c("Bias", "MSE", "Coverage")))
ggplot(data = metrics_data) +
  geom_pointrange(aes(x = forcats::fct_rev(mu), y = value,
                      colour = model,
                      shape = model,
                      ymin=value-value_se, ymax=value+value_se),
                  position=position_dodge(width=0.25), size = 0.3) +
  geom_hline(data = filter(metrics_data, metric == "Coverage"), aes(yintercept = 0.95), colour = "dark red") +
  geom_hline(data = filter(metrics_data, metric != "Coverage"), aes(yintercept = 0.0), colour = "dark red") +
  geom_hline(data = filter(metrics_data, metric == "Coverage"), aes(yintercept = 0.925), linetype = "dashed") +
  geom_hline(data = filter(metrics_data, metric == "Coverage"), aes(yintercept = 0.975), linetype = "dashed") +
  facet_grid(condition~metric, scales = "free_x") +
  coord_flip() +
  xlab("") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_rect(fill="white", size = 0))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1))  


#=====================================================================================#
# Plots
# Bias
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = bias,
                      # group = condition,
                      colour = model,
                      shape = model,
                      ymin=bias-bias_mcmc_se, ymax=bias+bias_mcmc_se),
                  position=position_dodge(width=0.25), size = 0.5) +
  geom_hline(yintercept = 0, colour = "dark red") +
  coord_flip() +
  xlab("") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.text.y = element_text(angle = 0))  

# Relative bias
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = rel_bias,
                      group = model,
                      colour = model,
                      ymin=rel_bias-rel_bias_mcmc_se, ymax=rel_bias+rel_bias_mcmc_se),
                  position=position_dodge(width=0.25)
  ) +
  geom_hline(yintercept = 0, colour = "dark red") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed") +
  coord_flip()+
  ggtitle("Relative bias in group-level transitions") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Empirical SE (precision)
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = empirical_se,
                      colour = model,
                      ymin=empirical_se-empirical_se_mcmc_se, ymax=empirical_se+empirical_se_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  ggtitle("EmpiricalSE in group-level transitions") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# MSE
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = MSE,
                      colour = model,
                      ymin=MSE-MSE_mcmc_se, ymax=MSE+MSE_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  ggtitle("MSE in group-level transitions") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Coverage
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = coverage,
                      colour = model,
                      ymin=coverage-coverage_mcmc_se, ymax=coverage+coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  coord_flip() +
  ggtitle("Coverage in group-level transitions") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Bias-corrected coverage
res %>%
  filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = bias_corr_coverage,
                      colour = model,
                      ymin=bias_corr_coverage-bias_corr_coverage_mcmc_se, ymax=bias_corr_coverage+bias_corr_coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  coord_flip() +
  ggtitle("Bias-corr. coverage in group-level transitions") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  


## Transitions group-level means (fixed effects) Logit scale
# Bias
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = bias,
                      colour = model,
                      ymin=bias-bias_mcmc_se, ymax=bias+bias_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  coord_flip() +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Relative bias
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = rel_bias,
                      colour = model,
                      ymin=rel_bias-rel_bias_mcmc_se, ymax=rel_bias+rel_bias_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed") +
  coord_flip() +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Empirical SE (precision)
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = empirical_se,
                      colour = model,
                      ymin=empirical_se-empirical_se_mcmc_se, ymax=empirical_se+empirical_se_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# MSE
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = MSE,
                      colour = model,
                      ymin=MSE-MSE_mcmc_se, ymax=MSE+MSE_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Coverage
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = coverage,
                      colour = model,
                      ymin=coverage-coverage_mcmc_se, ymax=coverage+coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  coord_flip() +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Bias-corrected coverage
res %>%
  filter(parameter %in% paste0("int_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = bias_corr_coverage,
                      colour = model,
                      ymin=bias_corr_coverage-bias_corr_coverage_mcmc_se, ymax=bias_corr_coverage+bias_corr_coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  coord_flip() +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  


#-----------------------------------------------------------------------------#
## Transitions random effects (fixed effects) Logit scale
# Bias
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = bias,
                      colour = model,
                      ymin=bias-bias_mcmc_se, ymax=bias+bias_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  coord_flip() +
  ggtitle("Bias in transitions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Relative bias
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = rel_bias,
                      # group = condition, 
                      colour = model,
                      ymin=rel_bias-rel_bias_mcmc_se, ymax=rel_bias+rel_bias_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed") +
  coord_flip() +
  ggtitle("Rel. bias in transitions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Empirical SE (precision)
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = empirical_se,
                      colour = model,
                      ymin=empirical_se-empirical_se_mcmc_se, ymax=empirical_se+empirical_se_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  ggtitle("EmpiricalSE in transitions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# MSE
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = MSE,
                      colour = model,
                      ymin=MSE-MSE_mcmc_se, ymax=MSE+MSE_mcmc_se),
                  position=position_dodge(width=0.75)) +
  coord_flip() +
  ggtitle("MSE in transitions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Coverage
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = coverage,
                      colour = model,
                      ymin=coverage-coverage_mcmc_se, ymax=coverage+coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  coord_flip() +
  ggtitle("Coverage in transitions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Bias-corrected coverage
res %>%
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) %>%
  ggplot() +
  geom_pointrange(aes(x = forcats::fct_rev(parameter), y = bias_corr_coverage,
                      colour = model,
                      ymin=bias_corr_coverage-bias_corr_coverage_mcmc_se, ymax=bias_corr_coverage+bias_corr_coverage_mcmc_se),
                  position=position_dodge(width=0.75)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  coord_flip() +
  ggtitle("Bias-corr. coverage in transitions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  


#-----------------------------------------------------------------------------#
## Emissions group-level means (fixed effects)
# Bias
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = bias,
                      colour = model,
                      ymin=bias-bias_mcmc_se, ymax=bias+bias_mcmc_se),
                  position=position_dodge(width=0.25), size = 0.2) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  facet_wrap(ndep~., ncol = 4)  +
  ggtitle("Bias in emissions fixed eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Relative bias
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = rel_bias,
                      colour = model,
                      ymin=rel_bias-rel_bias_mcmc_se, ymax=rel_bias+rel_bias_mcmc_se),
                  position=position_dodge(width=0.25)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed") +
  coord_flip() +
  facet_wrap(ndep~., nrow = 4) +
  ggtitle("Rel. bias in emissions fixed eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Empirical SE (precision)
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = empirical_se,
                      colour = model,
                      ymin=empirical_se-empirical_se_mcmc_se, ymax=empirical_se+empirical_se_mcmc_se),
                  position=position_dodge(width=0.25)) +
  coord_flip() +
  facet_wrap(ndep~., nrow = 4) +
  ggtitle("EmpiricalSE in emissions fixed eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# MSE
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = MSE,
                      colour = model,
                      ymin=MSE-MSE_mcmc_se, ymax=MSE+MSE_mcmc_se),
                  position=position_dodge(width=0.25)) +
  coord_flip() +
  facet_wrap(ndep~., nrow = 4) +
  ggtitle("MSE in emissions fixed eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  


# Coverage
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = coverage,
                      colour = model,
                      ymin=coverage-coverage_mcmc_se, ymax=coverage+coverage_mcmc_se),
                  position=position_dodge(width=0.25)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  coord_flip() +
  facet_wrap(ndep~., nrow = 4) +
  ggtitle("Coverage in emissions fixed eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Bias-corrected coverage
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = bias_corr_coverage,
                      colour = model,
                      ymin=bias_corr_coverage-bias_corr_coverage_mcmc_se, ymax=bias_corr_coverage+bias_corr_coverage_mcmc_se),
                  position=position_dodge(width=0.25)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  coord_flip() +
  facet_wrap(ndep~., nrow = 4) +
  ggtitle("Bias-corr. coverage in emissions fixed eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  



## Emissions between subject variances (random effects)
# Bias
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = bias,
                      colour = model,
                      ymin=bias-bias_mcmc_se, ymax=bias+bias_mcmc_se),
                  position=position_dodge(width=0.25)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  facet_wrap(ndep~., ncol = 4)  +
  ggtitle("Bias in emissions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Relative bias
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = rel_bias,
                      colour = model,
                      ymin=rel_bias-rel_bias_mcmc_se, ymax=rel_bias+rel_bias_mcmc_se),
                  position=position_dodge(width=0.25)) +
  geom_hline(yintercept = 0, colour = "dark red") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed") +
  coord_flip() +
  facet_wrap(ndep~., nrow = 4) +
  ggtitle("Rel. bias in emissions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Empirical SE (precision)
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = empirical_se,
                      colour = model,
                      ymin=empirical_se-empirical_se_mcmc_se, ymax=empirical_se+empirical_se_mcmc_se),
                  position=position_dodge(width=0.25)) +
  coord_flip() +
  facet_wrap(ndep~., nrow = 4) +
  ggtitle("EmpiricalSE in emissions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# MSE
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = MSE,
                      colour = model,
                      ymin=MSE-MSE_mcmc_se, ymax=MSE+MSE_mcmc_se),
                  position=position_dodge(width=0.25)) +
  coord_flip() +
  facet_wrap(ndep~., nrow = 4) +
  ggtitle("MSE in emissions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  


# Coverage
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = coverage,
                      colour = model,
                      ymin=coverage-coverage_mcmc_se, ymax=coverage+coverage_mcmc_se),
                  position=position_dodge(width=0.25)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  coord_flip() +
  facet_wrap(ndep~., nrow = 4) +
  ggtitle("Coverage in emissions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  

# Bias-corrected coverage
res %>%
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot() +
  geom_pointrange(aes(x = mu, y = bias_corr_coverage,
                      colour = model,
                      ymin=bias_corr_coverage-bias_corr_coverage_mcmc_se, ymax=bias_corr_coverage+bias_corr_coverage_mcmc_se),
                  position=position_dodge(width=0.25)) +
  geom_hline(yintercept = 0.95, colour = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  coord_flip() +
  ggtitle("Bias-corr. coverage in emissions random eff.") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))  






#-----------------------------------------------------------------------------#
# Comparing mean absolute bias

# Transitions fixed effects (prob domain)
res |>
  filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) |>
  ggplot(aes(x = model, y = abs(bias),  fill = model)) +
  geom_boxplot(width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  ggtitle("Transition fixed effects (probs)") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) |>
  ggplot(aes(x = model, y = abs(rel_bias),  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = c(0.1), linetype = "dashed") +
  ggtitle("Transition fixed effects (probs)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) |>
  ggplot(aes(x = model, y = MSE,  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  ggtitle("Transition fixed effects (probs)") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("S",rep(1:m, each = m),"toS",1:m)) |>
  ggplot(aes(x = model, y = coverage,  fill = model)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = 0.95, color = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  ggtitle("Transition fixed effects (probs)") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))



# Transitions fixed effects (logit domain)
res |>
  filter(parameter %in% paste0("int_S",rep(1:m, each = m),"toS",1:m)) |>
  ggplot(aes(x = model, y = abs(bias),  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  ggtitle("Transition fixed effects (probs)") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("int_S",rep(1:m, each = m),"toS",1:m)) |>
  ggplot(aes(x = model, y = abs(rel_bias),  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = c(0.1), linetype = "dashed") +
  ggtitle("Transition fixed effects (probs)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("int_S",rep(1:m, each = m),"toS",1:m)) |>
  ggplot(aes(x = model, y = MSE,  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  ggtitle("Transition fixed effects (probs)") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("int_S",rep(1:m, each = m),"toS",1:m)) |>
  ggplot(aes(x = model, y = coverage,  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = 0.95, color = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  ggtitle("Transition fixed effects (probs)") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))



# Transitions random effects (prob domain)
res |>
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) |>
  ggplot(aes(x = model, y = abs(bias),  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  ggtitle("Transition random effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) |>
  ggplot(aes(x = model, y = abs(rel_bias),  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = c(0.1), linetype = "dashed") +
  ggtitle("Transition random effects") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) |>
  ggplot(aes(x = model, y = MSE,  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  ggtitle("Transition random effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("var_S",rep(1:m, each = (m-1)),"toS",2:m)) |>
  ggplot(aes(x = model, y = coverage,  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = 0.95, color = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  ggtitle("Transition random effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))








# Emission fixed effects
res |>
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) |>
  ggplot(aes(x = model, y = abs(bias),  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  ggtitle("Emission fixed effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) |>
  ggplot(aes(x = model, y = abs(rel_bias),  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = c(0.1), linetype = "dashed") +
  ggtitle("Emission fixed effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))


res |>
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) |>
  ggplot(aes(x = model, y = MSE,  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  ggtitle("Emission fixed effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) |>
  ggplot(aes(x = model, y = coverage,  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = 0.95, color = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  ggtitle("Emission fixed effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) |>
  ggplot(aes(x = model, y = bias_corr_coverage,  fill = model)) +
  geom_boxplot(width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = 0.95, color = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  ggtitle("Emission fixed effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))


# Emission random effects
res |>
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot(aes(x = model, y = abs(bias),  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  ggtitle("Emission random effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot(aes(x = model, y = abs(rel_bias),  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = c(0.1), linetype = "dashed") +
  ggtitle("Emission random effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))


res |>
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot(aes(x = model, y = MSE,  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  ggtitle("Emission random effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot(aes(x = model, y = coverage,  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = 0.95, color = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  ggtitle("Emission random effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))

res |>
  filter(parameter %in% paste0("dep",rep(1:n_dep,each = m),"_var",1:m)) %>%
  separate(parameter, c("ndep","mu"),"_") %>%
  ggplot(aes(x = model, y = bias_corr_coverage,  fill = model)) +
  geom_boxplot(trim = FALSE, width = 0.2, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.3) +
  geom_hline(yintercept = 0.95, color = "dark red") +
  geom_hline(yintercept = c(0.925, 0.975), linetype = "dashed") +
  ggtitle("Emission random effects") +
  theme_minimal() +     theme(plot.title = element_text(hjust = 0.5)) +     theme(strip.text.y = element_text(angle = 0))   +
  theme(plot.title = element_text(hjust = 0.5))



## Statistics

# Transitions
res %>%
  filter(
    model == "impute",
    parameter %in% paste0("S",1:m,"toS",1:m)) %>%
  summarise(mean(abs(rel_bias)),
            median(abs(rel_bias)),
            sd(abs(rel_bias)))
res %>%
  filter(
    model == "impute",
    parameter %in% paste0("S",1:m,"toS",1:m)) %>%
  as.data.frame()

res %>%
  filter(
    model == "impute",
    parameter %in% paste0("int_S",1:m,"toS",1:m)) %>%
  summarise(mean(abs(rel_bias)),
            median(abs(rel_bias)),
            sd(abs(rel_bias)))
res %>%
  filter(
    model == "impute",
    parameter %in% paste0("int_S",1:m,"toS",1:m)) %>%
  as.data.frame()

res %>%
  filter(
    # model == "impute",
    parameter %in% paste0("int_S",1:m,"toS",1:m)) %>%
  summarise(mean(coverage),
            median(coverage),
            sd(coverage),
            min(coverage),
            max(coverage)) %>%
  as.data.frame()

# Emissions
res %>%
  filter(
    model == "impute",
    parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  summarise(mean(abs(rel_bias)),
            median(abs(rel_bias)),
            sd(abs(rel_bias)),
            min(rel_bias),
            max(rel_bias)) %>%
  as.data.frame()

res %>%
  filter(
    parameter %in% paste0("dep",rep(1:n_dep,each = m),"_mu",1:m)) %>%
  summarise(mean(coverage),
            median(coverage),
            sd(coverage),
            min(coverage),
            max(coverage)) %>%
  as.data.frame()
