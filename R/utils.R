
#------------------------------------------------#
#               Utility functions                #
#------------------------------------------------#

# computes probabilities from intercepts
int_to_prob <- function(int1) {
  if(is.matrix(int1)){
    prob1 <- matrix(nrow = nrow(int1), ncol = ncol(int1) + 1)
    for(r in 1:nrow(int1)){
      exp_int1 	<- matrix(exp(c(0, int1[r,])), nrow  = 1)
      prob1[r,] <- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
    }
  } else {
    exp_int1 	<- matrix(exp(c(0, int1)), nrow  = 1)
    prob1 		<- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
  }
  return(round(prob1,4))
}

# computes intercepts from probabilities, per row of input matrix
# first catagory is reference catagory
prob_to_int <- function(prob1){
  prob1 <- prob1 + 0.00001
  b0 <- matrix(NA, nrow(prob1), ncol(prob1)-1)
  sum_exp <- numeric(nrow(prob1))
  for(r in 1:nrow(prob1)){
    sum_exp[r] <- (1/prob1[r,1]) - 1
    for(cr in 2:ncol(prob1)){
      #for every b0 except the first collumn (e.g. b012 <- log(y12/y11-y12))
      b0[r,(cr-1)] <- log(prob1[r,cr]*(1+sum_exp[r]))
    }
  }
  return(round(b0,4))
}


# Define functions
get_durations <- function(x, cov_value = NULL, burnin = NULL, m = NULL){
  
  if(is.null(burnin)){
    burnin <- x$input$burn_in
  }
  
  if(is.null(m)){
    m <- x$input$m
  }
  
  if(is.null(cov_value)){
    apply(x$gamma_int_bar,
                    1,
                    function(r) matrix(t(int_to_prob(matrix(r, nrow = m, byrow = TRUE))), nrow = 1)) %>%
                t() %>%
                as.data.frame() %>%
                purrr::set_names(paste0("S",rep(1:m,each = m),"toS",1:m)) %>%
      mutate(iter = row_number()) %>%
      gather(parameter, value, -iter) %>%
      mutate(value = 1/(1-value)) %>%
      filter(iter > burnin,
             parameter %in% paste0("S",1:m,"toS",1:m)) %>%
      group_by(parameter) %>%
      summarise(median_val = median(value), CCI_lwr = quantile(value, 0.025), CCI_upr = quantile(value, 0.975))
  } else {
    do.call(rbind, lapply(1:length(cov_value), function(val) {
      as.data.frame(x$gamma_int_bar + x$gamma_cov_bar * cov_value[val]) %>%
        summarise(apply(.,
                        1,
                        function(r) matrix(t(int_to_prob(matrix(r, nrow = m, byrow = TRUE))), nrow = 1)) %>%
                    t() %>%
                    as.data.frame() %>%
                    purrr::set_names(paste0("S",rep(1:m,each = m),"toS",1:m))) %>%
        mutate(iter = row_number()) %>%
        gather(parameter, value, -iter) %>%
        mutate(value = 1/(1-value)) %>%
        filter(iter > burnin,
               parameter %in% paste0("S",1:m,"toS",1:m)) %>%
        group_by(parameter) %>%
        summarise(median_val = median(value), CCI_lwr = quantile(value, 0.025), CCI_upr = quantile(value, 0.975)) %>%
        mutate(cov_val = cov_value[val])
    }))
  }
  
}

plot_durations <- function(x, cci = TRUE){
  
  p <- ggplot(x, aes(x = cov_val,
                     y = median_val,
                     group = parameter,
                     color = parameter,
                     fill = parameter)) +
    geom_line() +
    # scale_y_log10() +
    theme_minimal()
  
  if(cci == TRUE){
    p <- p + geom_ribbon(aes(ymin = CCI_lwr, ymax = CCI_upr,), alpha = 0.3, linetype = 0)
  }
  
  p
  
}

# # Get expected durations
# exp_durations <- get_durations(x = out, cov_value = seq(1,5,0.1)) # Change values for the desired covariate

# # Get plot
# plot_durations(exp_durations, cci = TRUE)
# plot_durations(exp_durations, cci = FALSE)
