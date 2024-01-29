#------------------------------------------#
#  Missingness utils: auxiliary functions  #
#------------------------------------------#

# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.
cat_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, gamma, delta = NULL){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep)
  out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
  return(out)
}

# Obtain fractional log likelihood for multinomial intercept only model, bassed on P. Rossi 2004
llmnl_int_frac <- function(beta, Obs, n_cat, pooled_likel, w, wgt){
  return((1 - w) * llmnl_int(beta = beta, Obs = Obs, n_cat = n_cat) + w * wgt * pooled_likel)
}

# Obtain mnl -Expected[Hessian]  for intercept only model, bassed on P.Rossi 2004
mnlHess_int <- function(int, Obs, n_cat){
  n_Obs 	<- length(Obs)
  betas   <- matrix(c(0, int), byrow = T, ncol = n_cat)
  prob    <- exp(betas) / sum(exp(betas))
  Hess    <- (diag(x = prob[-1], nrow = n_cat-1) - outer(prob[-1],prob[-1])) * n_Obs
  return(Hess)
}

# one run of the random walk metropolis sampler for an intercept only multinomial distribution
# this means no covariates at the lower/time level
mnl_RW_once <- function(int1, Obs, n_cat, mu_int_bar1, V_int1, scalar, candcov1) {
  # obtain likelihood and transition prob with the parameters sampled in the previous iteration and current sampled state sequence
  oldloglike	 		<- llmnl_int(beta = int1, Obs = Obs, n_cat = n_cat)
  oldpostlike	 		<- oldloglike + dmvnorm(int1, mu_int_bar1, V_int1, log = TRUE)
  probold				  <- int_to_prob(int1)
  
  # obtain new parameters for gamma from proposal distribution plus new likelihood
  int_new		 		  <- int1 + rmvnorm(1, rep(0, (n_cat - 1)), scalar^2 * candcov1, method = "svd")
  newloglike	 		<- llmnl_int(beta = int_new, Obs = Obs, n_cat = n_cat)
  newpostlike	 		<- newloglike + dmvnorm(int_new, mu_int_bar1, V_int1, log = TRUE)
  probnew				  <- int_to_prob(int_new)
  
  # determine to use the updated or current (previous iteration) gamma values of the parameters
  acc 				   <- min(log(1), (newpostlike - oldpostlike))
  if(acc < log(1)) {
    unif         <- log(runif(1))
  } else {
    unif         <- log(1)
  }
  if (unif <= acc) {
    draw_int		<- int_new
    accept			<- 1
    prob			  <- probnew
  } else {
    draw_int		<- int1
    accept			<- 0
    prob			  <- probold
  }
  return(list(draw_int = draw_int, accept = accept, prob = prob))
}


# simple functions used in mHMM
dif_matrix <- function(rows, cols){
  return(matrix(, ncol = cols, nrow = rows))
}

nested_list <- function(n_dep, m){
  return(rep(list(vector("list", n_dep)),m))
}

dif_vector <- function(x){
  return(numeric(x))
}

is.whole <- function(x) {
  return(is.numeric(x) && floor(x) == x)
}

is.mHMM <- function(x) {
  inherits(x, "mHMM")
}

is.mHMM_gamma <- function(x) {
  inherits(x, "mHMM_gamma")
}

hms <- function(t){
  paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
        formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
        formatC(t %% 60, width = 2, format = "d", flag = "0"),
        sep = ":")
}


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

nested_list <- function(n_dep, m){
  return(rep(list(vector("list", n_dep)),m))
}

#------------------------------------------------------------------------------#
library(Rcpp)
cppFunction('List cat_mult_fw_cpp(NumericMatrix allprobs, NumericMatrix gamma, int m, int n, NumericVector delta) {

            int i, t, k;
            NumericVector foo(m);
            NumericVector foo1(m);
            double sumfoo;
            NumericMatrix alpha_prob(m,n);
            double lscale;
            NumericMatrix lalpha(m, n);

            if(std::any_of(allprobs(0,_).cbegin(), allprobs(0,_).cend(), NumericVector::is_na)){
              foo = delta;
            } else {
              foo = delta * allprobs(0,_);
            }
            sumfoo = 0;
            for (i = 0; i < m; i++){
            sumfoo += foo[i];
            }
            for (i = 0; i < m; i++){
            alpha_prob(i,0) = foo[i]/sumfoo;
            }
            lscale = log(sumfoo);
            for(i = 0; i < m; i++){
            lalpha(i, 0) = log(alpha_prob(i, 0)) + lscale;
            }

            for (t = 1; t < n; t++){
            for (i = 0; i < m; i++){
            foo1[i] = 0;
            for (k = 0; k < m; k++){
            foo1[i] += alpha_prob(k, (t - 1)) * gamma(k,i);
            }
            }
            if(std::any_of(allprobs(t,_).cbegin(), allprobs(t,_).cend(), NumericVector::is_na)){
              foo = foo1;
            } else {
              foo = foo1 * allprobs(t,_);
            }
            sumfoo = 0;
            for (i = 0; i < m; i++){
            sumfoo += foo[i];
            }
            for (i = 0; i < m; i++){
            alpha_prob(i,t) = foo[i]/sumfoo;
            }
            lscale = lscale + log(sumfoo);
            for(i = 0; i < m; i++){
            lalpha(i, t) = log(alpha_prob(i, t)) + lscale;
            }
            }

            return List::create(alpha_prob, lalpha);
            }')



cppFunction('double llmnl_int(NumericVector beta, IntegerVector Obs, int n_cat) {

    int n_Obs = Obs.size();
    //2: Calculate log sum only once:
        // double expBetas_log_sum = log(sum(exp(betas)));
        double expBetas_log_sum = 1.0; // std::exp(0)
        for (int i = 1; i < n_cat; i++) {
            expBetas_log_sum += std::exp(beta[i-1]);
        };
        expBetas_log_sum = std::log(expBetas_log_sum);

        double ll_sum = 0;
        //3: Use n_Obs, to avoid calling Xby.size() every time
        for (int i = 0; i < n_Obs; i++) {
            if(Obs[i] == 1L) continue;
            ll_sum += beta[Obs[i]-2L];
        };
        //4: Use that we know denom is the same for all I:
            ll_sum = ll_sum - expBetas_log_sum * n_Obs;
            return ll_sum;
}')

# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.
pois_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, gamma, delta = NULL){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "poisson")
  out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
  return(out)
}

# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.
cont_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, gamma, delta = NULL){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "continuous")
  out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
  return(out)
}

# Calculates the probabilities of observing each state at each point in time given
# the observations of all dependent variables, used for the forward probabilities
# Based on Zuchini 2016.
all1 <- function(x, emiss, n_dep, data_distr){
  inp <- rep(list(NULL), n_dep)
  if(data_distr == "categorical"){
    for(q in 1:n_dep){
      inp[[q]] <- t(emiss[[q]][,x[,q]])
    }
  } else if (data_distr == "continuous"){
    for(q in 1:n_dep){
      inp[[q]] <- outer(x[,q], Y = emiss[[q]][,1], FUN = dnorm, sd = rep(sqrt(emiss[[q]][,2]), each = dim(x)[1]))
    }
  } else if (data_distr == "poisson") {
    for(q in 1:n_dep){
      inp[[q]] <- outer(x[,q], Y = emiss[[q]][,1], FUN = dpois)
    }
  }
  allprobs <- Reduce("*", inp)
  return(allprobs)
}


#------------------------------------------------------------------------------#
# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.
cont_mult_fw <- function(x, m, emiss, n_dep, gamma, delta = NULL){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "continuous")
  out <- mult_fw_mar(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
  return(out)
}

# Could maybe made external
# Calculates the forward probabilities, used for sampling the state sequence
# Based on Zuchini 2016.

mult_fw_mar <- function(allprobs, m, n, emiss, gamma, delta = NULL){
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  lalpha   <- alpha_prob <- matrix(NA_real_, m, n)
  if(any(is.na(allprobs[1,]))){
    foo <- delta
  } else {
    foo             <- delta * allprobs[1, ]
  }
  sumfoo          <- sum(foo)
  alpha_prob[, 1] <- foo/sumfoo
  lscale          <- log(sumfoo)
  lalpha[, 1]     <- log(alpha_prob[, 1]) + lscale
  for (i in 2:n){
    if(any(is.na(allprobs[i,]))){
      foo              <- alpha_prob[, (i - 1)] %*% gamma
    } else {
      foo              <- alpha_prob[, (i - 1)] %*% gamma * allprobs[i, ]
    }
    sumfoo           <- sum(foo)
    alpha_prob[, i]  <- foo / sumfoo
    lscale           <- lscale + log(sumfoo)
    lalpha[, i]       <- log(alpha_prob[, i]) + lscale
  }
  list(forward_p = alpha_prob, la = lalpha)
}



obtain_gamma_cont <- function(object, level = "group", burn_in = NULL){
  # if (!is.mHMM(object) & !is.mHMM_cont(object)){
  #   stop("The input object used should either be from the class mHMM or mHMM_cont, obtained by using the function mHMM or mHMM_cont.")
  # }
  if (level != "group" & level != "subject"){
    stop("The specification at the input variable -level- should be set to either group or subject")
  }
  input   <- object$input
  n_subj  <- input$n_subj
  if (is.null(burn_in)){
    burn_in <- input$burn_in
  }
  J       <- input$J
  if (burn_in >= (J-1)){
    stop(paste("The specified burn in period should be at least 2 points smaller
               compared to the number of iterations J, J =", J))
  }
  m       <- input$m
  if(is.mHMM(object)){
    q_emiss <- input$q_emiss
  }
  n_dep   <- input$n_dep
  est <- matrix(, ncol = m, nrow = m)
  colnames(est) <- paste("To state", 1:m)
  rownames(est) <- paste("From state", 1:m)
  if (level == "group"){
    est[] <- matrix(round(apply(object$gamma_prob_bar[((burn_in + 1): J),], 2, median),3),
                    byrow = TRUE, ncol = m, nrow = m)
    est_gamma <- est
  }
  if (level == "subject"){
    est_gamma <- rep(list(est), n_subj)
    names(est_gamma) <- paste("Subject", 1:n_subj)
    for(i in 1:n_subj){
      if(is.mHMM(object)){
        est_gamma[[i]][] <- matrix(round(apply(object$PD_subj[[i]][burn_in:J, (1 + sum(m * q_emiss)) : (sum(m * q_emiss) + m*m)], 2, median), 3),
                                   byrow = TRUE, ncol = m, nrow = m)
        # } else if (is.mHMM_cont(object)){
      } else {
        est_gamma[[i]][] <- matrix(round(apply(object$PD_subj[[i]][burn_in:J, (1 + m * 2 * n_dep) : (m * 2 * n_dep + m * m)], 2, median), 3),
                                   byrow = TRUE, ncol = m, nrow = m)
      }
    }
  }
  class(est_gamma) <- append(class(est_gamma), "mHMM_gamma")
  return(est_gamma)
}



vit_mHMM_cont_mar <- function(object, s_data, burn_in = NULL){
  
  if (!("mHMM_cont" %in% class(object))){
    stop("The input object used should be from the class mHMM_cont, obtained by using the function mHMM_cont.")
  }
  # if (!is.mHMM_cont(object)){ ## THIS VERSION PRODUCES AN ERROR
  #   stop("The input object used should be from the class mHMM_cont, obtained by using the function mHMM_cont.")
  # }
  id         <- unique(s_data[,1])
  n_subj     <- length(id)
  if(length(object$PD_subj) != n_subj){
    stop("s_data used should be identical to the data used for creating the object in mHMM_cont.
         The number of subjects in the datasets are not the same.")
  }
  n_vary     <- table(s_data[,1])
  max_n      <- max(n_vary)
  state_seq  <- matrix(,ncol = n_subj, nrow = max_n)
  
  input      <- object$input
  n_dep      <- input$n_dep
  m          <- input$m
  q_emiss    <- input$q_emiss
  if(is.null(burn_in)){
    burn_in  <- input$burn_in
  }
  J          <- input$J
  if (burn_in >= (J-1)){
    stop(paste("The specified burn in period should be at least 2 points smaller
               compared to the number of iterations J, J =", J))
  }
  est_emiss  <- rep(list(rep(list(matrix(,nrow = m, ncol = 2)),n_dep)), n_subj)
  for(s in 1:n_subj){
    for(q in 1:n_dep){
      est_emiss[[s]][[q]][] <- matrix(round(apply(object$PD_subj[[s]][burn_in:J, (1+(q-1)*m) : (m+(q-1)*m)], 2, median), 4),
                                      byrow = TRUE, nrow = m, ncol = 1)
      est_emiss[[s]][[q]][] <- matrix(round(c(apply(object$PD_subj[[s]][burn_in:J, ((q-1) * m + 1) : ((q-1) * m + m)], 2, median),
                                              apply(object$PD_subj[[s]][burn_in:J, (n_dep * m + (q-1) * m + 1) : (n_dep * m + (q-1) * m + m)], 2, median)),
                                            3),
                                      ncol = 2,
                                      nrow = m)
    }
  }
  est_gamma <- obtain_gamma_cont(object, level = "subject")
  for(s in 1:n_subj){
    emiss <- est_emiss[[s]]
    gamma    <- est_gamma[[s]]
    probs    <- cont_mult_fw_r_to_cpp(x = as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                      m = m, emiss = emiss, n_dep = n_dep, gamma = gamma)[[1]]
    state_seq[1:n_vary[s], s] <- apply(probs, 2, which.max)
  }
  colnames(state_seq) <- paste("Subj_", id, sep = "")
  return(state_seq)
}




all1 <- function(x, emiss, n_dep, data_distr, left_truncation = NA){
  inp <- rep(list(NULL), n_dep)
  if(data_distr == "categorical"){
    for(q in 1:n_dep){
      inp[[q]] <- t(emiss[[q]][,x[,q]])
    }
  } else if (data_distr == "continuous"){
    for(q in 1:n_dep){
      inp[[q]] <- outer(x[,q], Y = emiss[[q]][,1], FUN = dnorm, sd = rep(sqrt(emiss[[q]][,2]), each = dim(x)[1]))
    }
  } else if (data_distr == "tr_continuous"){
    for(q in 1:n_dep){
      na_idx <- is.na(x[,q])
      inp[[q]] <- outer(x[,q], Y = emiss[[q]][,1], FUN = truncnorm::dtruncnorm, a = left_truncation[q], b = Inf, sd = rep(sqrt(emiss[[q]][,2]), each = dim(x)[1]))
      inp[[q]][na_idx,] <- NA
      # inp[[q]] <- outer(X = x[,q], Y = emiss[[q]][,1], FUN = function(X, Y, a, b, sd) truncnorm::dtruncnorm(x = X, mean = Y, a = left_truncation[q], b = Inf, sd = rep(sqrt(emiss[[q]][,2]), each = dim(x)[1])),a = left_truncation[q], b = Inf, sd = rep(sqrt(emiss[[q]][,2]), each = dim(x)[1]))
    }
  } else if (data_distr == "poisson") {
    for(q in 1:n_dep){
      inp[[q]] <- outer(x[,q], emiss[[q]][,1], FUN = dpois)
    }
  }
  allprobs <- Reduce("*", inp)
  return(allprobs)
}

tr_cont_mult_fw_r_to_cpp <- function(x, m, emiss, n_dep, gamma, delta = NULL, left_truncation = NA){
  if(is.na(left_truncation)){
    left_truncation <- rep(0,n_dep)
  }
  if(is.null(delta)) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  n        <- dim(x)[1]
  allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "tr_continuous", left_truncation = left_truncation)
  out <- cat_mult_fw_cpp(allprobs = allprobs, gamma = gamma, m = m, n = n, delta = delta)
  return(out)
}


vit_mHMM_tr_cont <- function(object, s_data, burn_in = NULL, left_truncation = NA){
  
  if (!("mHMM_cont" %in% class(object))){
    stop("The input object used should be from the class mHMM_cont, obtained by using the function mHMM_cont.")
  }
  # if (!is.mHMM_cont(object)){ ## THIS VERSION PRODUCES AN ERROR
  #   stop("The input object used should be from the class mHMM_cont, obtained by using the function mHMM_cont.")
  # }
  id         <- unique(s_data[,1])
  n_subj     <- length(id)
  if(length(object$PD_subj) != n_subj){
    stop("s_data used should be identical to the data used for creating the object in mHMM_cont.
         The number of subjects in the datasets are not the same.")
  }
  n_vary     <- table(s_data[,1])
  max_n      <- max(n_vary)
  state_seq  <- matrix(NA_integer_,ncol = n_subj, nrow = max_n)
  state_probs <- vector("list", length = n_subj)
  
  input      <- object$input
  n_dep      <- input$n_dep
  m          <- input$m
  q_emiss    <- input$q_emiss
  if(is.null(burn_in)){
    burn_in  <- input$burn_in
  }
  J          <- input$J
  if (burn_in >= (J-1)){
    stop(paste("The specified burn in period should be at least 2 points smaller
               compared to the number of iterations J, J =", J))
  }
  est_emiss  <- rep(list(rep(list(matrix(,nrow = m, ncol = 2)),n_dep)), n_subj)
  for(s in 1:n_subj){
    for(q in 1:n_dep){
      est_emiss[[s]][[q]][] <- matrix(round(apply(object$PD_subj[[s]][burn_in:J, (1+(q-1)*m) : (m+(q-1)*m)], 2, median), 4),
                                      byrow = TRUE, nrow = m, ncol = 1)
      est_emiss[[s]][[q]][] <- matrix(round(c(apply(object$PD_subj[[s]][burn_in:J, ((q-1) * m + 1) : ((q-1) * m + m)], 2, median),
                                              apply(object$PD_subj[[s]][burn_in:J, (n_dep * m + (q-1) * m + 1) : (n_dep * m + (q-1) * m + m)], 2, median)),
                                            3),
                                      ncol = 2,
                                      nrow = m)
    }
  }
  est_gamma <- obtain_gamma_pois(object, level = "subject")
  for(s in 1:n_subj){
    emiss <- est_emiss[[s]]
    gamma    <- est_gamma[[s]]
    probs    <- tr_cont_mult_fw_r_to_cpp(x = as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                         m = m, emiss = emiss, n_dep = n_dep, gamma = gamma, left_truncation = left_truncation)[[1]]
    state_seq[1:n_vary[s], s] <- apply(probs, 2, which.max)
    state_probs[[s]] <- t(probs)
  }
  colnames(state_seq) <- paste0("subj_", id)
  names(state_probs) <- paste0("subj_", id)
  return(list("state_seq" = state_seq, "state_probs" = state_probs))
}




sim_mHMM_tr_cont <- function(n_t, n, data_distr = 'tr_continuous', m, n_dep = 1, left_truncation = 0,
                             start_state = NULL, q_emiss = NULL, gamma, emiss_distr, xx_vec = NULL, beta = NULL,
                             var_gamma = 0.1, var_emiss = NULL, return_ind_par = FALSE){
  
  #############
  # Inbuild checks for correct specification of parameters ---------------------
  #############
  
  if (dim(gamma)[1] != m){
    stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
  }
  if (dim(gamma)[2] != m){
    stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
  }
  if(!isTRUE(all.equal(apply(gamma,1,sum), rep(1,m)))){
    stop("The elements in each row of the transition probability matrix gamma should sum up to 1")
  }
  if(!is.list(emiss_distr)){
    stop("The format of emiss_distr should be a list with", n_dep, "elements.")
  }
  if(length(emiss_distr) != n_dep){
    stop("The number of dependent variables specified in n_dep and the number of elements specified in the list emiss_distr should be equal")
  }
  if(data_distr == "categorical" & length(q_emiss) != n_dep){
    stop("The length of q_emiss specifying the number of output categories for each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
  }
  for(i in 1:n_dep){
    if (dim(emiss_distr[[i]])[1] != m){
      stop(paste("The number of rows of emission distribution matrix in element", i, "should be
                       equal to the number of states, which is", m, "."))
    }
    if(data_distr == 'categorical'){
      if (dim(emiss_distr[[i]])[2] != q_emiss[i]){
        stop(paste("The number of columns of the emission distribution matrix should be
                           equal to the number of observable categories, which is", q_emiss[i], ". See emission distribution in element", i, "."))
      }
      if(!isTRUE(all.equal(apply(emiss_distr[[i]], 1, sum), rep(1, m)))){
        stop("The elements in each row of the emission distribution matrix should sum up to 1, see emission distribution in element", i, ".")
      }
    }
    # if(data_distr == 'continuous'){
    #
    # }
  }
  if(length(left_truncation) != n_dep){
    stop("The input argument left_truncation should be numeric vector with length n_dep")
  }
  if((is.null(xx_vec) & !is.null(beta)) | (!is.null(xx_vec) & is.null(beta))){
    stop("Either only xx_vec or only beta is specified. Please specify both 1) the values for the covariate
             in xx_vec and 2) the values of the regression parameters in beta, to allow correct simulation of the
             data.")
  }
  if(!is.null(xx_vec)){
    if((!is.null(xx_vec[[1]]) & is.null(beta[[1]])) |
       (!is.null(xx_vec[[2]]) & is.null(beta[[2]]))){
      stop("Either only xx_vec or only beta is specified in one of the elements.
                 Please specify both 1) the values for the covariate in xx_vec and 2)
                 the values of the regression parameters in beta if either one is not
                 empty, to allow correct simulation of the data.")
    }
  }
  if(!is.null(beta)){
    # extend to all 1 + n_dep
    if((!is.null(beta[[1]]) & is.null(xx_vec[[1]])) |
       (!is.null(beta[[2]]) & is.null(xx_vec[[2]]))){
      stop("Either only xx_vec or only beta is specified in one of the elements.
                 Please specify both 1) the values for the covariate in xx_vec and 2)
                 the values of the regression parameters in beta if either one is not
                 empty, to allow correct simulation of the data.")
    }
  }
  if(!is.null(xx_vec)){
    # extend to all 1 + n_dep
    if((!is.null(xx_vec[[1]]) & length(xx_vec[[1]]) != n) |
       (!is.null(xx_vec[[2]]) & length(xx_vec[[2]]) != n)){
      stop("The length of the vectors in xx_vec should be equal to the number of subjects to be simulated,
                 set in n, if (the element in) xx_vec is not set to NULL.")
    }
  }
  if (!is.null(beta)){
    if (!is.null(beta[[1]])){
      if ((dim(beta[[1]])[1] != (m)) | (dim(beta[[1]])[2] != (m-1))){
        stop(paste("The first element of beta to predict the transiton probability matrix gamma should be a m (", m, " ) by m - 1 (", m - 1, ") matrix."))
      }
    }
    if (!is.null(beta[[2]]) & data_distr == 'categorical'){
      # extend to all 1 + n_dep and continuous
      if((dim(beta[[2]])[1] != (m)) | (dim(beta[[2]])[2] != (q_emiss[1]-1))){
        stop(paste("The second element of beta to predict the emission distribution should be a m (", m, ") by q_emiss - 1 (", q_emiss[1] - 1, ") matrix."))
      }
    }
  }
  if(is.null(xx_vec)){
    xx_vec <- rep(list(NULL), n_dep + 1)
    for(i in 1:(n_dep + 1)){
      xx_vec[[i]] <- rep(1,n)
    }
  } else {
    for(i in 1:(n_dep + 1)){
      if(is.null(xx_vec[[i]])) {
        xx_vec[[i]] <- rep(1,n)
      }
    }
  }
  if(is.null(beta)){
    beta <- rep(list(NULL), n_dep + 1)
    beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
    for(i in 2:(n_dep + 1)){
      if(data_distr == 'categorical'){
        beta[[i]] <- matrix(0, ncol = q_emiss[i-1] - 1, nrow = m)
      } else if (data_distr %in% c('continuous', 'tr_continuous')){
        beta[[i]] <- matrix(0, ncol = 1, nrow = m)
      } else if (data_distr == 'poisson-lognormal'){
        beta[[i]] <- matrix(0, ncol = 1, nrow = m)
      }
    }
  } else {
    if(is.null(beta[[1]])) {
      beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
    }
    for (i in 2:(n_dep + 1)){
      if (is.null(beta[[i]])) {
        if(data_distr == 'categorical'){
          beta[[i]] <- matrix(0, ncol = q_emiss[i-1] - 1, nrow = m)
        } else if (data_distr %in% c('continuous', 'tr_continuous')){
          beta[[i]] <- matrix(0, ncol = 1, nrow = m)
        } else if (data_distr == 'poisson-lognormal'){
          beta[[i]] <- matrix(0, ncol = 1, nrow = m)
        }
      }
    }
  }
  
  # if(data_distr == 'continuous'){
  #     if(n == 1){
  #         var_gamma <- 0
  #         var_emiss <- rep(0, n_dep)
  #     }
  #
  #     if(is.null(var_emiss)){
  #         var_emiss <- rep(0.1, n_dep)
  #     } else if(length(var_emiss) != n_dep){
  #         stop("The length of var_emiss specifying variance between subjects in each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
  #     }
  # }
  
  if(data_distr %in% c('poisson-lognormal','continuous', 'tr_continuous')){
    
    # If only 1 subject
    if(n == 1){
      var_gamma <- 0
      var_emiss <- rep(0, n_dep)
    }
    
    # If a single value of var_gamma specified, use for all transitions
    if(length(var_gamma) == 1){
      var_gamma <- matrix(rep(var_gamma, m*(m-1)),nrow = m, byrow = TRUE)
    } else if(is.matrix(var_gamma)){
      if (dim(var_gamma)[1] != m){
        stop(paste("The betweem-subject variance matrix for the transition distribution should be a", m, "by", m-1, "matrix."))
      }
      if (dim(var_gamma)[2] != m-1){
        stop(paste("The betweem-subject variance matrix for the transition distribution should be a", m, "by", m-1, "matrix."))
      }
    }
    
    # If a single value of var_emiss specified, use for all states and n_dep
    if(is.null(var_emiss)){
      var_emiss <- rep(list(NULL), n_dep)
      for(i in 1:n_dep){
        var_emiss[[i]] <- matrix(rep(0.1, m), nrow = m, byrow = TRUE)
      }
    } else if(is.numeric(var_emiss) & length(var_emiss) == n_dep){
      arg_var_emiss <- var_emiss
      var_emiss <- rep(list(NULL), n_dep)
      for(i in 1:n_dep){
        var_emiss[[i]] <- matrix(rep(arg_var_emiss[i], m), nrow = m, byrow = TRUE)
      }
    } else if(is.list(var_emiss)){
      for(i in 1:n_dep){
        if(dim(var_emiss[[i]])[1] != m){
          stop(paste0("The number of rows of the between-subject variance for the emission distribution should be
                           equal to ",m,", the number of hidden states chosen."))
        }
        if(dim(var_emiss[[i]])[2] != 1){
          stop(paste0("The number of columns of the between-subject variance for the emission distribution should be
                           equal to one."))
        }
      }
    } else if(length(var_emiss) != n_dep){
      stop("The length of var_emiss specifying variance between subjects in each of the number of dependent variables should equal the number of dependent variables specified in n_dep. Note that var_emiss can either by a list of matrices, or a numeric vector.")
    }
    
  }
  
  # UPDATED CHECKS: else if(length(var_emiss)==){old check}, if(is.list(var_emiss)){new check}
  if(data_distr == 'categorical'){
    
    # If only 1 subject
    if(n == 1){
      var_gamma <- matrix(rep(0, m*(m-1)),nrow = m, byrow = TRUE)
      var_emiss <- rep(list(NULL), n_dep)
      for(i in 1:n_dep){
        var_emiss[[i]] <- matrix(rep(0, m*(q_emiss[i]-1)),nrow = m, byrow = TRUE)
      }
    }
    
    # If a single value of var_gamma specified, use for all categories
    if(length(var_gamma) == 1){
      var_gamma <- matrix(rep(var_gamma, m*(m-1)),nrow = m, byrow = TRUE)
    } else if(is.matrix(var_gamma)){
      if (dim(var_gamma)[1] != m){
        stop(paste("The betweem-subject variance matrix for the transition distribution should be a", m, "by", m-1, "matrix."))
      }
      if (dim(var_gamma)[2] != m-1){
        stop(paste("The betweem-subject variance matrix for the transition distribution should be a", m, "by", m-1, "matrix."))
      }
    }
    
    # If a single value of var_emiss specified, use for all categories and n_dep
    if(is.null(var_emiss)){
      var_emiss <- rep(list(NULL), n_dep)
      for(i in 1:n_dep){
        var_emiss[[i]] <- matrix(rep(0.1, m*(q_emiss[i]-1)),nrow = m, byrow = TRUE)
      }
    } else if(is.numeric(var_emiss) & length(var_emiss) == n_dep){
      arg_var_emiss <- var_emiss
      var_emiss <- rep(list(NULL), n_dep)
      for(i in 1:n_dep){
        var_emiss[[i]] <- matrix(rep(arg_var_emiss[i], m*(q_emiss[i]-1)),nrow = m, byrow = TRUE)
      }
    } else if(is.list(var_emiss)){
      for(i in 1:n_dep){
        if(dim(var_emiss[[i]])[2] != q_emiss[i]-1){
          stop(paste("The number of columns of the between-subject variance for the emission distribution should be
                           equal to the number of observable categories minus one, which is", q_emiss[i], ". See emission distribution in element", i, "."))
        }
      }
    } else if(length(var_emiss) != n_dep){
      stop("The length of var_emiss specifying variance between subjects in each of the number of dependent variables should equal the number of dependent variables specified in n_dep. Note that var_emiss can either by a list of matrices, or a numeric vector.")
    }
  }
  
  #############
  # Simulating the data ---------------------
  #############
  
  states <- matrix(ncol = 2, nrow = n_t*n)
  states[,1] <- rep(1:n, each = n_t)
  obs <- matrix(ncol = 1 + n_dep, nrow = n_t*n)
  obs[,1] <- rep(1:n, each = n_t)
  sub_gamma <- rep(list(NULL), n)
  sub_emiss <- rep(list(vector("list", n_dep)), n)
  mnl_gamma <- prob_to_int(gamma)
  if(data_distr == "categorical"){
    mnl_emiss <- rep(list(NULL), n_dep)
    for(i in 1:n_dep){
      mnl_emiss[[i]] <- prob_to_int(emiss_distr[[i]])
    }
  }
  for(j in 1:n){
    sub_gamma[[j]] <- int_to_prob(mnl_gamma + xx_vec[[1]][j] * beta[[1]] +
                                    rnorm(n = m * (m-1), mean = 0, sd = sqrt(as.numeric(var_gamma))))
    for(i in 1:n_dep){
      if(data_distr == "categorical"){
        sub_emiss[[j]][[i]] <- int_to_prob(mnl_emiss[[i]] + xx_vec[[1+i]][j] * beta[[1+i]] +
                                             rnorm(n = m * (q_emiss[i]-1), mean = 0, sd = sqrt(as.numeric(var_emiss[[i]]))))
      } else if(data_distr == "continuous"){
        sub_emiss[[j]][[i]] <- emiss_distr[[i]]
        sub_emiss[[j]][[i]][,1] <- emiss_distr[[i]][,1] +  xx_vec[[1+i]][j] * beta[[1+i]] +
          rnorm(n = m, mean = 0, sd = sqrt(as.numeric(var_emiss[[i]])))
      } else if(data_distr == "tr_continuous"){
        sub_emiss[[j]][[i]] <- emiss_distr[[i]]
        sub_emiss[[j]][[i]][,1] <- emiss_distr[[i]][,1] +  xx_vec[[1+i]][j] * beta[[1+i]] +
          rnorm(n = m, mean = 0, sd = sqrt(as.numeric(var_emiss[[i]]))) # Add truncation here too?
      } else if(data_distr == "poisson-lognormal"){
        sub_emiss[[j]][[i]] <- emiss_distr[[i]]
        sub_emiss[[j]][[i]][,1] <- exp(emiss_distr[[i]][,1] +  xx_vec[[1+i]][j] * beta[[1+i]] +
                                         rnorm(n = m, mean = 0, sd = sqrt(as.numeric(var_emiss[[i]]))))
      }
    }
    
    if(n_t != 0){
      init <- solve(t(diag(m) - sub_gamma[[j]] + 1), rep(1, m))
      if (is.null(start_state)){
        states[((j-1) * n_t + 1), 2] <- sample(x = 1:m, size = 1, prob = init)
      } else {
        states[((j-1) * n_t + 1), 2] <- start_state
      }
      if(data_distr == "categorical"){
        for(i in 1:n_dep){
          obs[((j-1) * n_t + 1), (1+i)] <- sample(x = 1:q_emiss[i], size = 1, prob = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],])
        }
      } else if (data_distr == "continuous"){
        for(i in 1:n_dep){
          obs[((j-1) * n_t + 1), (1+i)] <- rnorm(1, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],1], sd = sqrt(sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],2]))
        }
      } else if (data_distr == "tr_continuous"){
        for(i in 1:n_dep){
          obs[((j-1) * n_t + 1), (1+i)] <- truncnorm::rtruncnorm(1, a = left_truncation[i], b = Inf, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],1], sd = sqrt(sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],2]))
        }
      } else if (data_distr == "poisson-lognormal"){
        for(i in 1:n_dep){
          obs[((j-1) * n_t + 1), (1+i)] <- rpois(1, lambda = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],1])
        }
      }
      for(t in 2:n_t){
        states[((j-1) * n_t + t), 2] <- sample(x = 1:m, size = 1, prob = sub_gamma[[j]][states[((j-1) * n_t + t - 1), 2],])
        if(data_distr == "categorical"){
          for(i in 1:n_dep){
            obs[((j-1) * n_t + t), (1+i)] <- sample(x = 1:q_emiss[i], size = 1, prob = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],])
          }
        } else if (data_distr == "continuous"){
          for(i in 1:n_dep){
            obs[((j-1) * n_t + t), (1+i)] <- rnorm(1, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],1], sd = sqrt(sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],2]))
          }
        } else if (data_distr == "tr_continuous"){
          for(i in 1:n_dep){
            obs[((j-1) * n_t + t), (1+i)] <- truncnorm::rtruncnorm(1, a = left_truncation[i], b = Inf, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],1], sd = sqrt(sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],2]))
          }
        } else if (data_distr == "poisson-lognormal"){
          for(i in 1:n_dep){
            obs[((j-1) * n_t + t), (1+i)] <- rpois(1, lambda = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],1])
          }
        }
      }
    }
  }
  
  #############
  # Returning output  ---------------------
  #############
  colnames(states) <- c("subj", "state")
  colnames(obs)    <- c("subj", paste("observation", 1:n_dep))
  if (return_ind_par == FALSE & n_t != 0){
    return(list(states = states, obs = obs))
  } else if (return_ind_par == TRUE & n_t != 0){
    return(list(states = states, obs = obs, subject_gamma = sub_gamma, subject_emiss = sub_emiss))
  } else if (n_t == 0){
    return(list(subject_gamma = sub_gamma, subject_emiss = sub_emiss))
  }
}
