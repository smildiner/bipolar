#=====================================================================================#
# Get bias, rel_bias, empirical_se, mse, coverage, etc, for group-level and 
#   subject-level parameters.

# Group-level
get_all_metrics <- function(out_files,
                            path = NULL,
                            model = NULL,
                            gamma_map = NULL,
                            gamma_var = NULL,
                            emiss_map = NULL,
                            emiss_var = NULL,
                            gamma_var_idx = NULL){
  
  if(is.null(path)){
    path <- paste0(getwd(),"/")
  }
  
  # True values
  true_gamma_int_bar <- as.numeric(t(prob_to_int(gamma_map)))
  true_gamma_prob_bar <- as.numeric(t(gamma_map))
  true_gamma_V_int_bar <- as.numeric(t(gamma_var))
  true_emiss_mu_bar <- emiss_map
  true_emiss_varmu_bar <- emiss_var
  
  # MAPs
  # Matrices
  gamma_int_bar <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$map$gamma_int_bar[["median"]]))
  gamma_prob_bar <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$map$gamma_prob_bar[["median"]]))
  gamma_V_int_bar <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$map$gamma_V_int_bar[["median"]][gamma_var_idx]))
  
  
  # CCIs
  # Matrices
  gamma_int_bar_cci_lwr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_int_bar[["CCI_95"]][seq(1,2*(m*m),2)]))
  gamma_int_bar_cci_upr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_int_bar[["CCI_95"]][seq(2,2*(m*m),2)]))
  
  gamma_prob_bar_cci_lwr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_prob_bar[["CCI_95"]][seq(1,2*(m*m),2)]))
  gamma_prob_bar_cci_upr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_prob_bar[["CCI_95"]][seq(2,2*(m*m),2)]))
  
  gamma_V_int_bar_cci_lwr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_V_int_bar[["CCI_95"]][seq(1,length(readRDS(paste0(path,out_file))$cci$gamma_V_int_bar[["CCI_95"]]),2)][gamma_var_idx]))
  gamma_V_int_bar_cci_upr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_V_int_bar[["CCI_95"]][seq(2,length(readRDS(paste0(path,out_file))$cci$gamma_V_int_bar[["CCI_95"]]),2)][gamma_var_idx]))
  
  # Lists (Ndep)
  emiss_mu_bar <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$map$emiss_mu_bar[[ndep]][["median"]] )) )
  emiss_varmu_bar <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$map$emiss_varmu_bar[[ndep]][["median"]] )) )
  
  emiss_mu_bar_cci_lwr <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$emiss_mu_bar[[ndep]][["CCI_95"]][seq(1,2*(m*n_dep),2)] )) )
  emiss_mu_bar_cci_upr <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$emiss_mu_bar[[ndep]][["CCI_95"]][seq(2,2*(m*n_dep),2)] )) )
  
  emiss_varmu_bar_cci_lwr <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$emiss_varmu_bar[[ndep]][["CCI_95"]][seq(1,2*(m*n_dep),2)] )) )
  emiss_varmu_bar_cci_upr <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$emiss_varmu_bar[[ndep]][["CCI_95"]][seq(2,2*(m*n_dep),2)] )) )
  
  # Get evaluation metrics
  gamma_int_bar_metrics <- cbind(parameter = paste0("int_S",rep(1:m, each = m-1),"toS",2:m), do.call(rbind, lapply(1:(m*(m-1)), function(i) summarize_simulation_scenario(true_values = true_gamma_int_bar[i], simulated_values = gamma_int_bar[,i], lower_cci = gamma_int_bar_cci_lwr[,i], upper_cci = gamma_int_bar_cci_upr[,i]))))
  gamma_prob_bar_metrics <- cbind(parameter = paste0("S",rep(1:m, each = m),"toS",1:m), do.call(rbind, lapply(1:(m*m), function(i) summarize_simulation_scenario(true_values = true_gamma_prob_bar[i], simulated_values = gamma_prob_bar[,i], lower_cci = gamma_prob_bar_cci_lwr[,i], upper_cci = gamma_prob_bar_cci_upr[,i]))))
  gamma_var_int_bar_metrics <- cbind(parameter = paste0("var_S",rep(1:m, each = m-1),"toS",2:m), do.call(rbind, lapply(1:(m*(m-1)), function(i) summarize_simulation_scenario(true_values = true_gamma_V_int_bar[i], simulated_values = gamma_V_int_bar[,i], lower_cci = gamma_V_int_bar_cci_lwr[,i], upper_cci = gamma_V_int_bar_cci_upr[,i]))))
  emiss_mu_bar_metrics <- do.call(rbind, lapply(1:n_dep, function(ndep) cbind(parameter = paste0("dep",ndep,"_mu",1:m),do.call(rbind, lapply(1:(m), function(s) summarize_simulation_scenario(true_values = true_emiss_mu_bar[[ndep]][s], simulated_values = emiss_mu_bar[[ndep]][,s], lower_cci = emiss_mu_bar_cci_lwr[[ndep]][,s], upper_cci = emiss_mu_bar_cci_upr[[ndep]][,s])
  )))))
  emiss_var_bar_metrics <- do.call(rbind, lapply(1:n_dep, function(ndep) cbind(parameter = paste0("dep",ndep,"_var",1:m),do.call(rbind, lapply(1:(m), function(s) summarize_simulation_scenario(true_values = true_emiss_varmu_bar[[ndep]][s], simulated_values = emiss_varmu_bar[[ndep]][,s], lower_cci = emiss_varmu_bar_cci_lwr[[ndep]][,s], upper_cci = emiss_varmu_bar_cci_upr[[ndep]][,s])
  )))))
  
  #Return
  return(bind_rows(gamma_int_bar_metrics, gamma_prob_bar_metrics, gamma_var_int_bar_metrics, emiss_mu_bar_metrics, emiss_var_bar_metrics))
}

# Subject-level
get_all_metrics_subj <- function(out_files,
                                 path = NULL,
                                 model = NULL,
                                 gamma_map = NULL,
                                 gamma_var = NULL,
                                 emiss_map = NULL,
                                 emiss_var = NULL,
                                 gamma_var_idx = NULL){
  
  if(is.null(path)){
    path <- paste0(getwd(),"/")
  }
  
  # Initialize objects
  n_subj <- length(out$truth$subject_gamma)
  
  # True values
  true_gamma_int_bar <- as.numeric(t(prob_to_int(gamma_map)))
  true_gamma_prob_bar <- as.numeric(t(gamma_map))
  true_gamma_V_int_bar <- as.numeric(t(gamma_var))
  true_emiss_mu_bar <- emiss_map
  true_emiss_varmu_bar <- emiss_var
  
  # MAPs
  # Matrices
  gamma_int_bar <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$map$gamma_int_bar[["median"]]))
  gamma_prob_bar <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$map$gamma_prob_bar[["median"]]))
  gamma_V_int_bar <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$map$gamma_V_int_bar[["median"]][gamma_var_idx]))
  
  
  # CCIs
  # Matrices
  gamma_int_bar_cci_lwr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_int_bar[["CCI_95"]][seq(1,2*(m*m),2)]))
  gamma_int_bar_cci_upr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_int_bar[["CCI_95"]][seq(2,2*(m*m),2)]))
  
  gamma_prob_bar_cci_lwr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_prob_bar[["CCI_95"]][seq(1,2*(m*m),2)]))
  gamma_prob_bar_cci_upr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_prob_bar[["CCI_95"]][seq(2,2*(m*m),2)]))
  
  gamma_V_int_bar_cci_lwr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_V_int_bar[["CCI_95"]][seq(1,length(readRDS(paste0(path,out_file))$cci$gamma_V_int_bar[["CCI_95"]]),2)][gamma_var_idx]))
  gamma_V_int_bar_cci_upr <- do.call(rbind, pbapply::pblapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$gamma_V_int_bar[["CCI_95"]][seq(2,length(readRDS(paste0(path,out_file))$cci$gamma_V_int_bar[["CCI_95"]]),2)][gamma_var_idx]))
  
  # Lists (Ndep)
  emiss_mu_bar <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$map$emiss_mu_bar[[ndep]][["median"]] )) )
  emiss_varmu_bar <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$map$emiss_varmu_bar[[ndep]][["median"]] )) )
  
  emiss_mu_bar_cci_lwr <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$emiss_mu_bar[[ndep]][["CCI_95"]][seq(1,2*(m*n_dep),2)] )) )
  emiss_mu_bar_cci_upr <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$emiss_mu_bar[[ndep]][["CCI_95"]][seq(2,2*(m*n_dep),2)] )) )
  
  emiss_varmu_bar_cci_lwr <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$emiss_varmu_bar[[ndep]][["CCI_95"]][seq(1,2*(m*n_dep),2)] )) )
  emiss_varmu_bar_cci_upr <- pbapply::pblapply(1:n_dep, function(ndep) do.call(rbind, lapply(out_files, function(out_file) readRDS(paste0(path,out_file))$cci$emiss_varmu_bar[[ndep]][["CCI_95"]][seq(2,2*(m*n_dep),2)] )) )
  
  # Get evaluation metrics
  gamma_int_bar_metrics <- cbind(parameter = paste0("int_S",rep(1:m, each = m-1),"toS",2:m), do.call(rbind, lapply(1:(m*(m-1)), function(i) summarize_simulation_scenario(true_values = true_gamma_int_bar[i], simulated_values = gamma_int_bar[,i], lower_cci = gamma_int_bar_cci_lwr[,i], upper_cci = gamma_int_bar_cci_upr[,i]))))
  gamma_prob_bar_metrics <- cbind(parameter = paste0("S",rep(1:m, each = m),"toS",1:m), do.call(rbind, lapply(1:(m*m), function(i) summarize_simulation_scenario(true_values = true_gamma_prob_bar[i], simulated_values = gamma_prob_bar[,i], lower_cci = gamma_prob_bar_cci_lwr[,i], upper_cci = gamma_prob_bar_cci_upr[,i]))))
  gamma_var_int_bar_metrics <- cbind(parameter = paste0("var_S",rep(1:m, each = m-1),"toS",2:m), do.call(rbind, lapply(1:(m*(m-1)), function(i) summarize_simulation_scenario(true_values = true_gamma_V_int_bar[i], simulated_values = gamma_V_int_bar[,i], lower_cci = gamma_V_int_bar_cci_lwr[,i], upper_cci = gamma_V_int_bar_cci_upr[,i]))))
  emiss_mu_bar_metrics <- do.call(rbind, lapply(1:n_dep, function(ndep) cbind(parameter = paste0("dep",ndep,"_mu",1:m),do.call(rbind, lapply(1:(m), function(s) summarize_simulation_scenario(true_values = true_emiss_mu_bar[[ndep]][s], simulated_values = emiss_mu_bar[[ndep]][,s], lower_cci = emiss_mu_bar_cci_lwr[[ndep]][,s], upper_cci = emiss_mu_bar_cci_upr[[ndep]][,s])
  )))))
  emiss_var_bar_metrics <- do.call(rbind, lapply(1:n_dep, function(ndep) cbind(parameter = paste0("dep",ndep,"_var",1:m),do.call(rbind, lapply(1:(m), function(s) summarize_simulation_scenario(true_values = true_emiss_varmu_bar[[ndep]][s], simulated_values = emiss_varmu_bar[[ndep]][,s], lower_cci = emiss_varmu_bar_cci_lwr[[ndep]][,s], upper_cci = emiss_varmu_bar_cci_upr[[ndep]][,s])
  )))))
  
  #Return
  return(bind_rows(gamma_int_bar_metrics, gamma_prob_bar_metrics, gamma_var_int_bar_metrics, emiss_mu_bar_metrics, emiss_var_bar_metrics))
}
