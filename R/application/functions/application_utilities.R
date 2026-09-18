# Application utility functions for the public replication script.
# This file intentionally contains the implementation details used by R/application/run_application_analysis.R.


# ---- application_utilities.R ----
get_counts <- function(name, n_time = 499, winlist) {
  n_districts <- length(winlist)
  point_counts <- matrix(0, nrow = n_districts, ncol = n_time)
  for (t in 1:n_time) {
    p <- dat_hfr[[name]][[t]]
    for (d in 1:n_districts) {
      point_counts[d, t] <- sum(inside.owin(p, w = winlist[[d]]))
    }
  }
  point_counts
}

add_lagged_sum <- function(res_list, colname = "SAF", lag = 1, new_colname = NULL) {
  if (is.null(new_colname)) new_colname <- paste0("prev_", colname, "_", lag)
  c_array <- res_list[[colname]]
  new_array <- sapply(1:nrow(c_array), function(t) {
    if (t == 1) return(rep(0, ncol(c_array)))
    index_l <- max(1, t - lag)
    colSums(c_array[index_l:(t - 1), , drop = FALSE])
  })
  res_list[[new_colname]] <- t(new_array)
  return(res_list)
}


add_lagged_max <- function(res_list, colname = "SAF", lag = 1, new_colname = NULL) {
  if (is.null(new_colname)) new_colname <- paste0("prev_", colname, "_", lag)
  c_array <- res_list[[colname]]
  new_array <- sapply(1:nrow(c_array), function(t) {
    if (t == 1) return(rep(0, ncol(c_array)))
    index_l <- max(1, t - lag)
    apply(c_array[index_l:(t - 1), , drop = FALSE], 2, max)
  })
  res_list[[new_colname]] <- t(new_array)
  res_list
}


take_log <- function(res_list, colname_vec){
  for (name in colname_vec) {
    res_list[[name]] <- log(res_list[[name]] + 1)
  }
  return(res_list)
}




quantile_factor <- function(x,nlevel=5) {
  # Ensure x is numeric
  x <- as.numeric(x)
  
  # Initialize output as all "0"
  result <- rep("0", length(x))
  
  # Nonzero values
  x_nonzero <- x[x != 0]
  
  # Quantile cutoffs
  
  q <- quantile(x_nonzero, probs = 1:(nlevel-1)/nlevel, na.rm = TRUE)
  
  # Assign levels
  for (l in 1:nlevel) {
    if(l==1){
      result[x > 0 & x <= q[1]] <- paste(l)
    }
    
    if(l>1){
      result[x > q[l-1] & x <= q[l]] <- paste(l)
    }
    
    if(l==nlevel){
      result[x > q[l-1]] <- paste(l)
    }
  }
  
  return(as.factor(result))
}





#------------------------compute asd-----------------------------------------


compute_asd_cat <- function(var, treat, weight = NULL) {
  # var: categorical or binary variable (factor or character or integer)
  # treat: binary treatment variable (0/1)
  # weight: optional weights (default NULL means unweighted)
  
  # Convert to factor to handle categorical properly
  var <- as.factor(var)
  treat <- as.numeric(treat)
  
  levels_var <- levels(var)
  
  # Get proportions for each level by group (treated/control)
  prop_table <- function(v, tr, w) {
    tmp <- xtabs(w ~ v + tr)
    return(sapply(c(0,1), function(x) tmp[,x+1]/sum(w[tr==x])))
  }
  # Unweighted
  w_unweighted <- rep(1, length(var))
  p_unweighted <- prop_table(var, treat, w_unweighted)
  
  # Weighted
  if (is.null(weight)) {
    p_weighted <- p_unweighted
  } else {
    p_weighted <- prop_table(var, treat, weight)
  }
  
  # Compute ASD for each level
  asd <- function(p) {
    abs(p[, 2] - p[, 1]) / sqrt((p[, 2] * (1 - p[, 2]) + p[, 1] * (1 - p[, 1])) / 2)
  }
  
  asd_unweighted <- asd(p_unweighted)
  asd_weighted   <- asd(p_weighted)
  
  if(length(levels_var)>2){
    return( list(data.frame(
      Level = levels_var,
      ASD = asd_unweighted
    ),
    data.frame(
      Level = levels_var,
      ASD = asd_weighted
    )
    ))
  }else{
    return(list(asd_unweighted,asd_weighted))
  }
  

}


compute_asd_num <- function(var, treat, weight = NULL,remove_zero = FALSE){
  
  
  if(remove_zero){
    index <- var>0
  }else{
    index <- 1:length(var)
  }
  
  ASD <- abs(mean(var[treat==1 & index])-mean(var[treat==0 & index]))/sqrt(1/2*(var(var[treat==1 & index])+var(var[treat==0 & index])))
  
  if (is.null(weight)) {
    weight <- rep(1, length(var))
  }

  w_mean_treat <- weighted.mean(var[treat == 1 & index], weight[treat == 1 & index])
  w_mean_control <- weighted.mean(var[treat == 0 & index], weight[treat == 0 & index])
  
  w_var <- function(x, w) {
    m <- weighted.mean(x, w)
    sum(w * (x - m)^2) / sum(w)
  }
  
  # Weighted pooled SD
  w_sd_treat <- w_var(var[treat == 1 & index], weight[treat == 1 & index])
  w_sd_control <- w_var(var[treat == 0 & index], weight[treat == 0 & index])
  w_pooled_sd <- sqrt((w_sd_treat + w_sd_control) / 2)
  
  ASD_weighted <- abs(w_mean_treat - w_mean_control) / w_pooled_sd
  
  return(list(ASD,ASD_weighted))
  
}


compute_ASD <- function(df,name_vector,treatment_name, weight_name){
  treatment <- df[,treatment_name]
  weights <- df[,weight_name]
  ASD <- list()
  ASD_weighted <- list()
  for (i in 1:length(name_vector)) {
    predictor <- df[,name_vector[i]]
    
    if(is.factor(predictor)){
      tmp <- compute_asd_cat(predictor, treatment, weights)
    }else{
      remove_zero <- ifelse(grepl("costamount",name_vector[i]),1,0)
      tmp <- compute_asd_num(predictor, treatment, weights,remove_zero)
    }
    
    ASD[[i]] <- tmp[[1]]
    ASD_weighted[[i]] <- tmp[[2]]
    
  }
  names(ASD) <- name_vector
  names(ASD_weighted) <- name_vector
  return(list(ASD = ASD, ASD_weighted = ASD_weighted))
}





# Helper to tidy ASD list into long format with variable-level label
tidy_asd_list <- function(asd_list, type_label) {
  bind_rows(
    lapply(names(asd_list), function(var) {
      value <- asd_list[[var]]
      if (is.data.frame(value)) {
        # For categorical: paste level onto variable name
        data.frame(
          Variable = var,
          VariableLevel = paste0(var, value[[1]]),
          ASD = value[[2]],
          Type = type_label
        )
      } else {
        # For continuous: use "var" only
        data.frame(
          Variable = var,
          VariableLevel = var,
          ASD = value,
          Type = type_label
        )
      }
    }),
    .id = NULL
  )
}





add_std <- function(data, vars) {
  # data   : named list of 2D-arrays (matrices)
  # vars   : character vector of element names in data to transform
  # offset : constant to add before taking log (default = 1)
  
  for (nm in vars) {
    if (! nm %in% names(data)) {
      warning(sprintf("'%s' not found in data; skipping.", nm))
      next
    }
    
    mat <- data[[nm]]
    
    # standardize
    std_name <- paste0(nm, "_std")
    data[[std_name]] <- (mat-mean(mat))/sd(mat)
    
  }
  
  data
}



add_log_binary_std <- function(data, vars, offset = 1) {
  # data   : named list of 2D-arrays (matrices)
  # vars   : character vector of element names in data to transform
  # offset : constant to add before taking log (default = 1)
  
  for (nm in vars) {
    if (! nm %in% names(data)) {
      warning(sprintf("'%s' not found in data; skipping.", nm))
      next
    }
    
    mat <- data[[nm]]
    # compute log-scale (log(x + offset))
    log_mat <- log(mat + offset)
    
    # name for the log version
    log_name <- paste0(nm, "_log")
    data[[log_name]] <- log_mat
    
    # binary indicator: is log_mat == 0?  (i.e. original mat was 0)
    binary_ind_name <- paste0(nm, "_binary")
    data[[binary_ind_name]] <- (log_mat != 0)
    
    # standardize
    std_name <- paste0(nm, "_log_std")
    data[[std_name]] <- (log_mat-mean(log_mat))/sd(log_mat)
    
    std_name <- paste0(nm, "_binary_std")
    data[[std_name]] <- (data[[binary_ind_name]]-mean(data[[binary_ind_name]]))/sd(data[[binary_ind_name]])
  }
  
  data
}


# ---- application_data_functions.R ----
# Helpers for application policy histories and value-estimator recursion.
get_history_app <- function(data, t, X_cols = c("all_airstrike"), Z_cols = c("cities_dist"),W_name = "aid_array", Y_name = "all_outcome") {
  H <- list()
  
  # --- W and Y ---
  W_now <- data[[W_name]][t, ]
  Y_now <- data[[Y_name]][t, ]
  
  # --- Combine X matrices ---
  X_now <- do.call(cbind, lapply(X_cols, function(x) data[[x]][t, ]))
  
  # --- Combine Z matrices ---
  Z_now <- do.call(cbind, lapply(Z_cols, function(z) data[[z]]))

  # --- Current time slice ---
  H[[1]] <- rbind(W = W_now, t(X_now), t(Z_now), Y = Y_now)
  rownames(H[[1]]) <- c("W", X_cols, Z_cols, "Y")
  
  # --- Previous time slice ---
  if (t == 1) {
    total_rows <- length(rownames(H[[1]]))
    H[[2]] <- array(0, c(total_rows, ncol(data$aid_array)))
    rownames(H[[2]]) <- rownames(H[[1]])
  } else {
    W_prev <- data[[W_name]][t - 1, ]
    Y_prev <- data[[Y_name]][t - 1, ]
    X_prev <- do.call(cbind, lapply(X_cols, function(x) data[[x]][t - 1, ]))
    Z_prev <- do.call(cbind, lapply(Z_cols, function(z) data[[z]]))
    
    H[[2]] <- rbind(W = W_prev, t(X_prev), t(Z_prev), Y = Y_prev)
    rownames(H[[2]]) <- rownames(H[[1]])
  }
  
  return(H)
}



get_X_all_app <- function(data,t, X_cols = c("all_airstrike"), Z_cols = c("cities_dist"), policy_cols = c("W","all_airstrike","cities_dist"),W_name = "aid_array", Y_name = "all_outcome"){
  
  H <- get_history_app(data,t, X_cols , Z_cols, W_name = W_name, Y_name = Y_name)
  X_all <- cbind(1,t(H[[2]][policy_cols,,drop=FALSE]))  
  
  return(X_all)
}







# Stabilized addIPW value contributions used by the application analysis.
get_Y_est2 <- function(Y, weight_ls, M, method, prob_arr = NULL, prob_ls = NULL,
                       district = 1:9) {
  n <- ncol(Y)
  time_points <- nrow(Y)
  a <- NULL
  weight_t_arr1 <- NULL
  weight_t_arr2 <- NULL
  
  Y_est_arr <- array(NA,c(M,time_points-M+1))
  weight_t_arr <- array(NA, c(M,time_points-M+1))

  
  if(method == "addIPW"){
    
    for (m in 1:M) {
      weight_t_arr[m,] <- sapply(1:(time_points-M+1), function(t) (sum(weight_ls[[m]][t,]-1,na.rm = TRUE)+1))
       
    }
    
    if(M>1){
      for (m in (M-1):1) {
        weight_t_arr[m,] <- weight_t_arr[m,]*weight_t_arr[m+1,]
      }
    }
    
    
    a <- 1/(1+exp(-(apply(weight_t_arr,1,mean)^2-1)/0.1))
    weight_t_arr1 <- t(sapply(1:length(a), function(x) 
      a[x]*weight_t_arr[x,]/mean(weight_t_arr[x,])))
    weight_t_arr2 <- t(sapply(1:length(a), function(x) 
      (1-a[x])*(weight_t_arr[x,]+1-mean(weight_t_arr[x,]))))
    weight_t_arr <- weight_t_arr1+weight_t_arr2
    
    for (t in 1:(time_points-M+1)) {
      
      for(m in 1:M){
        Y_est_arr[m,t] <- weight_t_arr[m,t]*mean(Y[t+M-m,district])
        if(m==M){
          prob_arr[m,t] <- prob_arr[m,t]
        }else{
          prob_arr[m,t] <- weight_t_arr[m+1,t]*prob_arr[m,t]
          prob_ls[[m]][t,] <- weight_t_arr[m+1,t]*prob_ls[[m]][t,]
        }
        
      }
    }
    
  }
  
  if(method == "addIPW2"){
    
    subsets <- combn(1:n, 2, simplify = TRUE)
    
    for (m in 1:M) {
      
      weight_t_arr[m,] <- sapply(1:(time_points-M+1), function(t) (sum(weight_ls[[m]][t,]-1,na.rm = TRUE)+sum(apply(subsets, 2, function(indices) prod(weight_ls[[m]][t,indices] - 1,na.rm = TRUE)))+1))
      
    }
    
    if(M>1){
      for (m in (M-1):1) {
        weight_t_arr[m,] <- weight_t_arr[m,]*weight_t_arr[m+1,]
      }
    }
    
    
    a <- 1/(1+exp(-(apply(weight_t_arr,1,mean)^2-1)/0.1))
    weight_t_arr1 <- t(sapply(1:length(a), function(x) 
      a[x]*weight_t_arr[x,]/mean(weight_t_arr[x,])))
    weight_t_arr2 <- t(sapply(1:length(a), function(x) 
      (1-a[x])*(weight_t_arr[x,]+1-mean(weight_t_arr[x,]))))
    weight_t_arr <- weight_t_arr1+weight_t_arr2
    
    for (t in 1:(time_points-M+1)) {

        for(m in 1:M){
          Y_est_arr[m,t] <- weight_t_arr[m,t]*mean(Y[t+M-m,district])
          if(m==M){
            prob_arr[m,t] <- prob_arr[m,t]
          }else{
            prob_arr[m,t] <- weight_t_arr[m+1,t]*prob_arr[m,t]
            prob_ls[[m]][t,] <- weight_t_arr[m+1,t]*prob_ls[[m]][t,]
          }
          
        }
      }
      

  }
  
  res <- list(Y_est_arr = Y_est_arr,
              weight_t_arr = weight_t_arr,
              prob_arr=prob_arr,
              prob_ls=prob_ls,
              a=a,
              weight_t_arr1=weight_t_arr1,
              weight_t_arr2=weight_t_arr2)

  return(res)

}



# ---- application_policy_value.R ----
# The application workflow compares addIPW and addIPW2 and uses stabilized
# time-level weights from get_Y_est2().




get_estimated_policy_value_app <- function(data, p, ga, method = "addIPW",
                                           X_cols, Z_cols, policy_cols,
                                           save_Y_est = FALSE, save_weight_t = FALSE, save_prob = FALSE,
                                           no_aid = NULL, al_aid = NULL, district = NULL,
                                           W_name = "aid_array_week", Y_name = "all_outcome",
                                           message = FALSE) {

  W <- data[[W_name]]

  if(length(p)==1){
    p <- array(p,dim(W))
  }
  
  ga <- matrix(ga,ncol=(1+length(policy_cols)))
  M <- nrow(ga)
  time_points <- nrow(W)
  n <- ncol(W)
  
  if(is.null(district)){
    district <- 1:n
  }
  
  weight_ls <- list()
  prob_ls <- list()
  prob_arr <- array(NA,c(M,time_points-M+1))
  
  if(length(method)==1 & M>1 ){
    method <- rep(method,M)
  }
  
  # get propensity score
  e <- array(NA,c(time_points,n))
  e[W==1] <- p[W==1]
  e[W==0] <- 1-p[W==0]
  for (m in 1:M) {
    
    # get the treatment assignment
    W_counter <- array(NA,c(time_points,n))
    for (t in 1:(time_points-m+1)) {
      X_all <- get_X_all_app(data,t,X_cols,Z_cols,policy_cols,W_name = W_name,Y_name = Y_name)
      W_counter[t,] <- X_all%*%ga[m,]
      W_counter[t,no_aid] <- -Inf
      W_counter[t,al_aid] <- Inf
      
    }

    
    # IPW weights for each unit
    prob <- (1+exp(-W_counter))^(-1)
    weight <- prob^W*(1-prob)^(1-W)/e
  
    
    prob_arr[m,] <- apply(prob, 1, mean)[(M-m+1):(time_points-m+1)]
    
    prob_ls[[m]] <- prob[(M-m+1):(time_points-m+1),]
    
    weight_ls[[m]] <- weight[(M-m+1):(time_points-m+1),]

    
    
    
  }

  tmp <- get_Y_est2(data[[Y_name]], weight_ls, M, method[M],
                    prob_arr = prob_arr, prob_ls = prob_ls,
                    district = district)
  Y_est_arr <- tmp$Y_est_arr
  
  
  
  V_est <- rowMeans(Y_est_arr,na.rm = TRUE)
  if(message){
    print(V_est)
  }
  Y_arr <- rowMeans(data[[Y_name]])[M:time_points]
    
  if(M>1){
    for (m in (M-1):1) {
      Y_arr <- rbind(Y_arr,rowMeans(data[[Y_name]])[m:(time_points+m-M)])
    }
  }
    
  if(is.vector(Y_arr)){
    Y_arr <- matrix(Y_arr,nrow = 1)
  }

  value_center <- tmp$a * V_est + (1 - tmp$a) * rowMeans(Y_arr)
  value_if <- Y_est_arr - sweep(tmp$weight_t_arr, 1, value_center, `*`)
  V_var <- mean(colSums(value_if)^2)/(time_points-M+1)
  

  CI_lb <- sum(V_est) - 1.96*sqrt(V_var)
  CI_ub <- sum(V_est) + 1.96*sqrt(V_var)
  
  if(!save_prob){
    tmp$prob_arr <- NULL
    tmp$prob_ls <- NULL
  }
  
  if(save_Y_est){
    return(list(tmp,Y_arr,V_est))
  }
  
  res <- list(V_est = sum(V_est), CI_lb = CI_lb, CI_ub = CI_ub,
              V_var = V_var, V_est_vec = V_est,
              prob_arr = tmp$prob_arr, prob_ls = tmp$prob_ls)


  
  if(save_weight_t){
    res$weight_t_arr <- tmp$weight_t_arr
    res$weight_ls <- weight_ls
    res$a <- tmp$a
    res$weight_t_arr1 <- tmp$weight_t_arr1
    res$weight_t_arr2 <- tmp$weight_t_arr2
  }
  
  
  return(res)
  
  
  
  
}
  
  
    
    

# ---- application_optimize_policy.R ----
# method: addIPW or addIPW2
# M: which period we want to optimize
# If M>0, need to set the policy for previous treatment ga_prev


get_optimal_policy_app <- function(data, p, method = "addIPW", ub = 30, lb = -30,
                               ga_prev = NULL, X_cols, Z_cols, policy_cols,
                               initial_par = NULL, no_aid = NULL, al_aid = NULL,
                               W_name = "aid_array_week", Y_name = "all_outcome",
                               value_barrier = 1e5) {
  W <- data[[ W_name ]]

  if(length(p)==1){
    p <- array(p,dim(W))
  }
  
  time_points <- nrow(W)
  n <- ncol(W)
  if(is.null(ga_prev)){
    M <- 1
  }else{
    M <- nrow(ga_prev)+1
  }
  if(method%in%c("addIPW", "addIPW2")){
    
    X_all <- do.call(rbind,lapply(1:(time_points-M+1),function(t) get_X_all_app(data,t,X_cols,Z_cols,policy_cols, W_name = W_name, Y_name = Y_name)))
    
    
    # Initial guesses
    par_len <- ncol(X_all)
    if(is.null(initial_par)){
      par <- rep(0,par_len)
      par[1] <- -5
    }else{
      par <- initial_par
    }
    
    
    objective_function <- function(par,X_all) {

       tmp <- get_estimated_policy_value_app(
         data, p, rbind(ga_prev, par),
         method = method, X_cols = X_cols, Z_cols = Z_cols,
         policy_cols = policy_cols, no_aid = no_aid, al_aid = al_aid,
         save_prob = TRUE, W_name = W_name, Y_name = Y_name
       )

       base_obj <- tmp$V_est
       if (value_barrier > 0) {
         neg_value <- pmin(tmp$V_est_vec, 0)
         neg_value[is.na(neg_value)] <- 0
         base_obj <- base_obj + value_barrier * sum(neg_value^2)
       }
       return(base_obj)
    }
    

    if(length(lb)==1){
      rep(lb,par_len)
    }
    
    if(length(ub)==1){
      rep(ub,par_len)
    }
    
    # Solve the smooth policy-search objective.
    result <- optim(
      par = par,
      fn = objective_function,
      method = "L-BFGS-B",
      X_all = X_all,
      lower = lb,  # Lower bounds
      upper = ub    # Upper bounds 
    )
    print("finish")

    
    res_policy <- rbind(ga_prev,result$par)
    value_stochastic <- get_estimated_policy_value_app(
      data, p, res_policy, method = method,
      X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
      no_aid = no_aid, al_aid = al_aid,
      W_name = W_name, Y_name = Y_name, message = TRUE
    )
    value_stochastic2 <- get_estimated_policy_value_app(
      data, p, res_policy, method = method,
      X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
      no_aid = no_aid, al_aid = al_aid, save_prob = TRUE,
      W_name = W_name, Y_name = Y_name
    )

    colnames(res_policy) <- c("intercept", policy_cols)
    return(list(policy = res_policy,value_stochastic = value_stochastic,
                avg_treat_district = rowMeans(value_stochastic2$prob_arr))
           )
  }
  
  

  
  
  
  return(rbind(ga_prev,ga))
}

# ---- application_test.R ----
perform_test_app <- function(data, p = NULL, ga_prev = NULL, method_prev = NULL,
                             lb = -100, ub = 100, var_percentage = NULL,
                             alpha = 0.05,
                             total_points = 10, B = 1000,
                             X_cols, Z_cols, policy_cols,
                             Y_name, W_name = "aid_array_week",
                             seed = 1234, grid_seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  W <- data[[W_name]]
  time_points <- nrow(W)
  M <- if (is.null(ga_prev)) 1 else nrow(as.matrix(ga_prev)) + 1
  n_gamma <- length(policy_cols) + 1

  if (!is.null(grid_seed)) {
    rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    set.seed(grid_seed)
  }

  gamma_grid <- as.data.frame(
    replicate(n_gamma, runif(total_points, min = lb, max = ub))
  )
  colnames(gamma_grid) <- paste0("gamma", 0:(n_gamma - 1))

  if (!is.null(grid_seed)) {
    if (!is.null(rng_state)) {
      assign(".Random.seed", rng_state, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }

  N <- nrow(gamma_grid)
  GP <- rep(NA_real_, N)
  tmp_func <- array(NA_real_, c(2 * N, time_points))

  get_if_contribution <- function(fit) {
    tmp <- fit[[1]]
    Y_arr <- fit[[2]]
    V_est <- fit[[3]]

    out <- apply(
      tmp$Y_est_arr -
        tmp$weight_t_arr * (tmp$a * V_est + (1 - tmp$a) * rowMeans(Y_arr)),
      2,
      sum
    )

    if (M >= 2) {
      out <- c(out, rep(NA_real_, M - 1))
    }
    out
  }

  for (i in seq_len(N)) {
    ga <- rbind(ga_prev, unlist(gamma_grid[i, ]))

    value1 <- get_estimated_policy_value_app(
      data, p, ga,
      method = c(method_prev, "addIPW"),
      X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
      Y_name = Y_name, W_name = W_name
    )

    value2 <- get_estimated_policy_value_app(
      data, p, ga,
      method = c(method_prev, "addIPW2"),
      X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
      Y_name = Y_name, W_name = W_name
    )

    GP[i] <- value1$V_est - value2$V_est

    fit1 <- get_estimated_policy_value_app(
      data, p, ga,
      method = c(method_prev, "addIPW"),
      X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
      Y_name = Y_name, W_name = W_name,
      save_Y_est = TRUE
    )

    fit2 <- get_estimated_policy_value_app(
      data, p, ga,
      method = c(method_prev, "addIPW2"),
      X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
      Y_name = Y_name, W_name = W_name,
      save_Y_est = TRUE
    )

    tmp_func[i, ] <- get_if_contribution(fit1)
    tmp_func[i + N, ] <- get_if_contribution(fit2)
  }

  n_eval_times <- time_points - M + 1
  contrast_func <- tmp_func[1:N, 1:n_eval_times, drop = FALSE] -
    tmp_func[(N + 1):(2 * N), 1:n_eval_times, drop = FALSE]
  contrast_func[!is.finite(contrast_func)] <- 0

  Z <- matrix(rnorm(n_eval_times * B), nrow = n_eval_times, ncol = B)
  simulated_paths <- contrast_func %*% Z / sqrt(n_eval_times)
  observed_process <- sqrt(n_eval_times) * GP

  observed_l2 <- mean(observed_process^2, na.rm = TRUE)
  bootstrap_l2 <- colMeans(simulated_paths^2, na.rm = TRUE)
  critical_value <- as.numeric(quantile(bootstrap_l2, 1 - alpha, na.rm = TRUE))
  p_value <- mean(bootstrap_l2 >= observed_l2, na.rm = TRUE)
  test <- as.integer(observed_l2 > critical_value)

  list(
    test = test,
    p_value = p_value,
    observed_l2 = observed_l2,
    critical_value = critical_value
  )
}

# ---- application_analysis_helpers.R ----
# ================================
# Helpers for policy learning + tables + plots
# ================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(sf)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(forcats)
  library(tibble)
  library(tidytext)   # reorder_within / scale_y_reordered
})

# ---- fixed aid-assignment district sets from empirical treatment rates ----
select_fixed_aid_districts <- function(data, always_aid_threshold, no_aid_threshold) {
  always_aid <- which(colMeans(data$aid_array_week) > always_aid_threshold)
  no_aid <- which(colMeans(data$aid_array_week) < no_aid_threshold)
  list(always_aid = always_aid, no_aid = no_aid)
}

# ---- run one time step ----
run_one_time <- function(data, ga_prev_policy, method_prev, initial_par,
                         X_cols, Z_cols, policy_cols, Y_name,
                         always_aid_threshold, no_aid_threshold, iraq_district,
                         value_barrier = 1e5,
                         lb = -100, ub = 100, var_percentage = 0.99,
                         alpha = 0.05, total_points = 10, B = 1000,
                         W_name = "aid_array_week") {
  
  # 1) pick districts with empirically near-deterministic aid histories
  fixed_districts <- select_fixed_aid_districts(
    data,
    always_aid_threshold = always_aid_threshold,
    no_aid_threshold = no_aid_threshold
  )
  
  # 2) run test to choose method for this step
  test_res <- perform_test_app(
    data = data, p = data$ps, ga_prev = ga_prev_policy, method_prev = method_prev,
    lb = lb, ub = ub, var_percentage = var_percentage,
    alpha = alpha, total_points = total_points, B = B,
    X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
    Y_name = Y_name, W_name = W_name
  )
  
  method_this <- if (isTRUE(test_res$test == 1)) "addIPW2" else "addIPW"
  cat("P-value:", test_res$p_value, "\n")
  cat("Use method:", method_this, "\n")
  
  # 3) optimize policy under chosen method
  cat(
    "Optimize:",
    "method =", method_this,
    "| outcome =", Y_name,
    "\n"
  )
  res_policy <- get_optimal_policy_app(
    data, data$ps, method_this,
    ub = 100, lb = -100,
    ga_prev = ga_prev_policy,
    X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
    initial_par = initial_par,
    W_name = W_name, Y_name = Y_name,
    al_aid = fixed_districts$always_aid, no_aid = fixed_districts$no_aid,
    value_barrier = value_barrier
  )
  
  # 4) whole-region value under chosen method
  tmp <- get_estimated_policy_value_app(
    data, data$ps,
    res_policy$policy,
    method = method_this,
    X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
    save_Y_est = FALSE, save_weight_t = TRUE, save_prob = TRUE,
    W_name = W_name, Y_name = Y_name,
    al_aid = fixed_districts$always_aid, no_aid = fixed_districts$no_aid
  )
  
  # 5) district-level
  n_district <- ncol(data[[W_name]])
  tmp_district <- list()
  tmp_district$V_est <- sapply(
    1:n_district,
    function(x) {
      get_estimated_policy_value_app(
        data, data$ps,
        res_policy$policy,
        method = method_this,
        X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
        save_Y_est = FALSE, save_weight_t = TRUE, save_prob = TRUE,
        W_name = W_name, Y_name = Y_name,
        al_aid = fixed_districts$always_aid, no_aid = fixed_districts$no_aid, district = x
      )$V_est
    }
  )
  tmp_district$ub <- sapply(
    1:n_district,
    function(x) {
      get_estimated_policy_value_app(
        data, data$ps,
        res_policy$policy,
        method = method_this,
        X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
        save_Y_est = FALSE, save_weight_t = TRUE, save_prob = TRUE,
        W_name = W_name, Y_name = Y_name,
        al_aid = fixed_districts$always_aid, no_aid = fixed_districts$no_aid, district = x
      )$CI_ub
    }
  )
  tmp_district$lb <- sapply(
    1:n_district,
    function(x) {
      get_estimated_policy_value_app(
        data, data$ps,
        res_policy$policy,
        method = method_this,
        X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
        save_Y_est = FALSE, save_weight_t = TRUE, save_prob = TRUE,
        W_name = W_name, Y_name = Y_name,
        al_aid = fixed_districts$always_aid, no_aid = fixed_districts$no_aid, district = x
      )$CI_lb
    }
  )
  
  # 6) governorate-level (average over its districts)
  districts <- colnames(data[[W_name]])
  govs <- unique(trimws(iraq_district$ADM2NAME))
  n_gov <- length(govs)
  tmp_gov <- list()
  tmp_gov$V_est <- sapply(
    1:n_gov,
    function(x) {
      get_estimated_policy_value_app(
        data, data$ps,
        res_policy$policy,
        method = method_this,
        X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
        save_Y_est = FALSE, save_weight_t = TRUE, save_prob = TRUE,
        W_name = W_name, Y_name = Y_name,
        al_aid = fixed_districts$always_aid, no_aid = fixed_districts$no_aid,
        district = which(districts %in% iraq_district$ADM3NAME[trimws(iraq_district$ADM2NAME) == govs[x]])
      )$V_est
    }
  )
  tmp_gov$ub <- sapply(
    1:n_gov,
    function(x) {
      get_estimated_policy_value_app(
        data, data$ps,
        res_policy$policy,
        method = method_this,
        X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
        save_Y_est = FALSE, save_weight_t = TRUE, save_prob = TRUE,
        W_name = W_name, Y_name = Y_name,
        al_aid = fixed_districts$always_aid, no_aid = fixed_districts$no_aid,
        district = which(districts %in% iraq_district$ADM3NAME[trimws(iraq_district$ADM2NAME) == govs[x]])
      )$CI_ub
    }
  )
  tmp_gov$lb <- sapply(
    1:n_gov,
    function(x) {
      get_estimated_policy_value_app(
        data, data$ps,
        res_policy$policy,
        method = method_this,
        X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols,
        save_Y_est = FALSE, save_weight_t = TRUE, save_prob = TRUE,
        W_name = W_name, Y_name = Y_name,
        al_aid = fixed_districts$always_aid, no_aid = fixed_districts$no_aid,
        district = which(districts %in% iraq_district$ADM3NAME[trimws(iraq_district$ADM2NAME) == govs[x]])
      )$CI_lb
    }
  )
  tmp_gov$ADM2NAME <- govs
  
  list(
    res_policy = res_policy,
    tmp = tmp,
    tmp_district = tmp_district,
    tmp_gov = tmp_gov,
    method = method_this,
    test_p = test_res$p_value,
    test_flag = test_res$test
  )
}

# ---- long tables from 'out' (district-level; unchanged) ----
make_aid_policy_long <- function(tmp_list, district_names) {
  Ms <- seq_along(tmp_list)
  purrr::map_dfr(Ms, function(M) {
    mat <- do.call(cbind, lapply(seq_len(M), function(m) {
      colMeans(tmp_list[[M]]$prob_ls[[m]], na.rm = TRUE)
    }))
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
    avg <- rowMeans(mat, na.rm = TRUE)
    if (!is.null(names(avg))) avg <- avg[district_names]
    tibble(
      ADM3NAME = district_names,
      policy_length = factor(paste0("M", M), levels = paste0("M", Ms)),
      avg_aid = as.numeric(avg)
    )
  })
}


# ---- driver over T periods ----
run_all_times <- function(data, initial_par_list, X_cols, Z_cols, policy_cols, Y_name,
                          always_aid_threshold = 0.95,
                          no_aid_threshold = 0.05,
                          iraq_district, value_barrier = 1e5,
                          # optional knobs forwarded to perform_test_app:
                          lb = -100, ub = 100, var_percentage = 0.99,
                          alpha = 0.05, total_points = 10, B = 1000,
                          W_name = "aid_array_week") {
  
  T <- length(initial_par_list)
  if (length(always_aid_threshold) == 1) always_aid_threshold <- rep(always_aid_threshold, T)
  if (length(no_aid_threshold) == 1) no_aid_threshold <- rep(no_aid_threshold, T)
  
  res_policy   <- vector("list", T)
  tmp          <- vector("list", T)
  tmp_district <- vector("list", T)
  tmp_gov      <- vector("list", T)
  method_used  <- character(T)
  test_pvals   <- numeric(T)
  test_flags   <- integer(T)
  
  ga_prev_policy <- NULL            # will become a matrix with rows = previous policies
  method_prev    <- character(0)    # vector of methods chosen so far
  
  
  
  for (t in seq_len(T)) {
    out <- run_one_time(
      data = data,
      ga_prev_policy = ga_prev_policy,
      method_prev = method_prev,
      initial_par = initial_par_list[[t]],
      X_cols = X_cols, Z_cols = Z_cols, policy_cols = policy_cols, Y_name = Y_name,
      always_aid_threshold = always_aid_threshold[[t]],
      no_aid_threshold = no_aid_threshold[[t]],
      iraq_district = iraq_district,
      value_barrier = value_barrier,
      lb = lb, ub = ub, var_percentage = var_percentage,
      alpha = alpha, total_points = total_points, B = B, W_name = W_name
    )
    
    res_policy[[t]]   <- out$res_policy
    tmp[[t]]          <- out$tmp
    tmp_district[[t]] <- out$tmp_district
    tmp_gov[[t]]      <- out$tmp_gov
    
    method_used[t]    <- out$method
    test_pvals[t]     <- out$test_p
    test_flags[t]     <- out$test_flag
    
    # carry forward for next step
    # --- append method history ---
    method_prev <- c(method_prev, out$method)
    
    
    ga_prev_policy <- out$res_policy$policy

    cat(
      "Finished M =", t,
      "| outcome =", Y_name,
      "| selected methods =", paste(method_prev, collapse = ", "),
      "| % treated =", paste(format_treated_percent(out$res_policy$avg_treat_district), collapse = ", "),
      "\n"
    )
    cat("Policy coefficients:\n")
    print(ga_prev_policy)
    
  }
  
  
  list(
    res_policy = res_policy,
    tmp = tmp,
    tmp_district = tmp_district,
    tmp_gov = tmp_gov,
    method_used = method_used,
    test_flags = test_flags,
    test_pvals = test_pvals
  )
}

format_treated_percent <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }
  scales::percent(x, accuracy = 0.1)
}


make_outcome_policy_long <- function(tmp_district_list, district_names) {
  Ms <- seq_along(tmp_district_list)
  purrr::map_dfr(Ms, function(M) {
    tibble(
      ADM3NAME = district_names,
      policy_length = factor(paste0("M", M), levels = paste0("M", Ms)),
      avg_outcome = as.numeric(tmp_district_list[[M]]$V_est) / M,
      avg_outcome_ub = as.numeric(tmp_district_list[[M]]$ub) / M,
      avg_outcome_lb = as.numeric(tmp_district_list[[M]]$lb) / M
    )
  })
}

attach_to_sf <- function(sf_obj, long_tbl, key = "ADM3NAME") {
  sf_obj %>%
    mutate(!!key := trimws(.data[[key]])) %>%
    inner_join(long_tbl, by = key)
}

complete_diff_for_full_map <- function(diff_tbl, iraq_district, value_cols) {
  policy_levels <- levels(diff_tbl$policy_length)
  if (is.null(policy_levels)) {
    policy_levels <- unique(as.character(diff_tbl$policy_length))
  }

  full_grid <- tidyr::expand_grid(
    ADM3NAME = trimws(iraq_district$ADM3NAME),
    policy_length = factor(policy_levels, levels = policy_levels)
  )

  out <- full_grid %>%
    left_join(
      diff_tbl %>%
        mutate(
          ADM3NAME = trimws(ADM3NAME),
          policy_length = factor(as.character(policy_length), levels = policy_levels)
        ),
      by = c("ADM3NAME", "policy_length")
    )

  for (col in value_cols) {
    out[[col]][is.na(out[[col]])] <- 0
  }
  out
}

# ---- prepare ALL mapping data (district; unchanged) ----
prepare_all_data <- function(iraq_district, data, out, Y_mat_name) {
  district_names <- colnames(data$aid_array_week)
  
  aid_obs <- colMeans(data$aid_array_week, na.rm = TRUE)
  y_obs   <- colMeans(data[[Y_mat_name]],  na.rm = TRUE)
  
  aid_policy_long <- make_aid_policy_long(out$tmp, district_names)
  y_policy_long   <- make_outcome_policy_long(out$tmp_district, district_names)
  
  aid_diff_long <- aid_policy_long %>%
    left_join(tibble(ADM3NAME = district_names, aid_obs = as.numeric(aid_obs)), by = "ADM3NAME") %>%
    mutate(diff_aid = avg_aid - aid_obs)
  
  y_diff_long <- y_policy_long %>%
    left_join(tibble(ADM3NAME = district_names, y_obs = as.numeric(y_obs)), by = "ADM3NAME") %>%
    mutate(
      diff_outcome     = avg_outcome     - y_obs,
      diff_outcome_lb  = avg_outcome_lb  - y_obs,
      diff_outcome_ub  = avg_outcome_ub  - y_obs
    )

  aid_diff_map <- complete_diff_for_full_map(
    aid_diff_long, iraq_district,
    value_cols = c("avg_aid", "aid_obs", "diff_aid")
  )
  y_diff_map <- complete_diff_for_full_map(
    y_diff_long, iraq_district,
    value_cols = c("avg_outcome", "avg_outcome_ub", "avg_outcome_lb",
                   "y_obs", "diff_outcome", "diff_outcome_lb", "diff_outcome_ub")
  )
  
  list(
    aid_policy_long = aid_policy_long,
    outcome_policy_long = y_policy_long,
    aid_diff_long = aid_diff_long,
    outcome_diff_long = y_diff_long,
    a_obs_aid_sf  = attach_to_sf(iraq_district, tibble(ADM3NAME = district_names, aid_obs = as.numeric(aid_obs))),
    a_obs_out_sf  = attach_to_sf(iraq_district, tibble(ADM3NAME = district_names, outcome_obs = as.numeric(y_obs))),
    a_policy_aid_sf = attach_to_sf(iraq_district, aid_policy_long),
    a_policy_out_sf = attach_to_sf(iraq_district, y_policy_long),
    a_diff_aid_sf = attach_to_sf(iraq_district, aid_diff_map),
    a_diff_out_sf = attach_to_sf(iraq_district, y_diff_map)
  )
}

# ---- governorate edges (one per ADM2NAME) ----
make_governorate_edges <- function(iraq_district) {
  iraq_district %>%
    mutate(ADM2NAME = trimws(ADM2NAME)) %>%
    group_by(ADM2NAME) %>%
    summarise(.groups = "drop")
}

.norm <- function(x) trimws(tolower(x))

dissolve_governorates_boundary <- function(iraq_district, remove_holes = TRUE, planar_crs = 3857) {
  crs0 <- sf::st_crs(iraq_district)
  ap   <- sf::st_make_valid(iraq_district)
  ap_p <- sf::st_transform(ap, planar_crs)
  adm2_poly <- ap_p %>%
    mutate(ADM2NAME = trimws(ADM2NAME)) %>%
    group_by(ADM2NAME) %>%
    summarise(geometry = sf::st_union(geometry), .groups = "drop")
  if (remove_holes) {
    if (requireNamespace("sfheaders", quietly = TRUE)) {
      adm2_poly <- sfheaders::sf_remove_holes(adm2_poly)
    } else if (requireNamespace("nngeo", quietly = TRUE)) {
      adm2_poly <- nngeo::st_remove_holes(adm2_poly, units = "km^2", threshold = 0.1)
    }
  }
  bnd <- sf::st_boundary(adm2_poly)
  bnd <- sf::st_collection_extract(bnd, "LINESTRING", warn = FALSE)
  bnd <- sf::st_cast(bnd, "MULTILINESTRING", warn = FALSE)
  bnd <- sf::st_line_merge(bnd)
  sf::st_transform(bnd, crs0)
}

# --- From a governorate-level diff table, pick top/bottom N by horizon M ---
top_bottom_governorates <- function(gov_diff_tbl, N = 5) {
  gov_diff_tbl %>%
    group_by(policy_length) %>%
    { bind_rows(
      slice_min(., order_by = value, n = N, with_ties = FALSE),
      slice_max(., order_by = value, n = N, with_ties = FALSE)
    ) } %>%
    distinct(policy_length, ADM2NAME, .keep_all = TRUE) %>%
    ungroup() %>%
    transmute(
      policy_length = factor(policy_length, levels = levels(gov_diff_tbl$policy_length)),
      ADM2NAME = trimws(ADM2NAME)
    )
}

make_highlight_edges <- function(iraq_district, selected_govs_by_M) {
  gov_edges <- make_governorate_edges(iraq_district) %>% mutate(ADM2NAME = trimws(ADM2NAME))
  inner_join(gov_edges, selected_govs_by_M, by = "ADM2NAME")
}
make_highlight_edges_clean <- function(iraq_district, selected_govs_by_M) {
  gov_bnd <- dissolve_governorates_boundary(iraq_district) %>%
    mutate(ADM2NAME = trimws(ADM2NAME))
  inner_join(gov_bnd, selected_govs_by_M, by = "ADM2NAME")
}

# ---- difference maps (district borders only) ----
make_diff_maps <- function(a_diff_aid_sf, a_diff_out_sf, iraq_district, event_type, p_aid_obs_title = "",
                           highlight_aid_edges = NULL, highlight_out_edges = NULL) {
  district_edges <- iraq_district
  lim_aid     <- max(abs(a_diff_aid_sf$diff_aid), na.rm = TRUE);     if (!is.finite(lim_aid)) lim_aid <- 0
  lim_outcome <- max(abs(a_diff_out_sf$diff_outcome), na.rm = TRUE); if (!is.finite(lim_outcome)) lim_outcome <- 0
  
  p_aid_diff <- ggplot(a_diff_aid_sf) +
    geom_sf(aes(fill = diff_aid), color = "white", linewidth = 0.2) +
    geom_sf(data = district_edges, fill = NA, color = "grey40", linewidth = 0.25) +
    { if (!is.null(highlight_aid_edges) && nrow(highlight_aid_edges) > 0)
      geom_sf(data = highlight_aid_edges, fill = NA, color = "black", linewidth = 0.7) else NULL } +
    scale_fill_gradient2("", low = "#2c7bb6", mid = "white", high = "#d7191c",
                         midpoint = 0, limits = c(-lim_aid, lim_aid),
                         oob = scales::squish, na.value = "grey90") +
    facet_wrap(~ policy_length, nrow = 1) +
    coord_sf(datum = NA) +
    theme_minimal() +
    theme(
      axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
      legend.text = element_text(size = 8), legend.title = element_text(size = 9),
      legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")
    ) +
    labs(title = paste0("Difference in aid (", event_type, "policy - observed)"))
  
  p_out_diff <- ggplot(a_diff_out_sf) +
    geom_sf(aes(fill = diff_outcome), color = "white", linewidth = 0.2) +
    geom_sf(data = district_edges, fill = NA, color = "grey40", linewidth = 0.25) +
    { if (!is.null(highlight_out_edges) && nrow(highlight_out_edges) > 0)
      geom_sf(data = highlight_out_edges, fill = NA, color = "black", linewidth = 0.7) else NULL } +
    scale_fill_gradient2("", low = "#2c7bb6", mid = "white", high = "#d7191c",
                         midpoint = 0, limits = c(-lim_outcome, lim_outcome),
                         oob = scales::squish, na.value = "grey90") +
    facet_wrap(~ policy_length, nrow = 1) +
    coord_sf(datum = NA) +
    theme_minimal() +
    theme(
      axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
      legend.text = element_text(size = 8), legend.title = element_text(size = 9),
      legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")
    ) +
    labs(title = paste0("Difference in outcome (", event_type, "policy - observed)"))
  
  leg_aid <- cowplot::get_legend(p_aid_diff)
  leg_out <- cowplot::get_legend(p_out_diff)
  p_aid_nl <- p_aid_diff + theme(legend.position = "none")
  p_out_nl <- p_out_diff + theme(legend.position = "none")
  aligned <- cowplot::align_plots(p_aid_nl, p_out_nl, align = "h")
  stack   <- cowplot::plot_grid(aligned[[1]], aligned[[2]], ncol = 1, rel_heights = c(1, 1))
  legcol  <- cowplot::plot_grid(leg_aid, leg_out, ncol = 1, rel_heights = c(1, 1))
  combined <- cowplot::plot_grid(stack, legcol, ncol = 2, rel_widths = c(1, 0.18))
  
  list(p_aid_diff = p_aid_diff, p_out_diff = p_out_diff, combined_diff = combined)
}

# ---- observed & policy level maps (district; unchanged) ----
make_level_maps <- function(a_obs_aid_sf, a_policy_aid_sf,
                            a_obs_out_sf, a_policy_out_sf, iraq_district, event_type = "",
                            p_out_obs_title ="", p_aid_obs_title = "") {
  district_edges <- iraq_district
  governorate_edges <- make_governorate_edges(iraq_district)
  baghdad_outline <- dissolve_governorates_boundary(iraq_district, remove_holes = TRUE) %>%
    mutate(ADM2NAME = trimws(ADM2NAME)) %>%
    filter(tolower(ADM2NAME) == "baghdad")
  mosul_basrah_dist <- iraq_district %>%
    mutate(ADM3NAME = trimws(ADM3NAME)) %>%
    filter(.norm(ADM3NAME) %in% c("mosul", "basrah"))
  
  rng_aid <- range(a_policy_aid_sf$avg_aid, na.rm = TRUE)
  rng_out <- range(a_policy_out_sf$avg_outcome, na.rm = TRUE)
  
  p_aid_obs <- ggplot(a_obs_aid_sf) +
    geom_sf(aes(fill = aid_obs), color = "grey30", linewidth = 0.2) +
    geom_sf(data = district_edges, fill = NA, color = "grey40", linewidth = 0.25) +
    geom_sf(data = baghdad_outline, fill = NA, color = "black", linewidth = 0.7, linejoin = "round") +
    geom_sf(data = mosul_basrah_dist, fill = NA, color = "black", linewidth = 0.7, linejoin = "round") +
    scale_fill_viridis_c("") +
    coord_sf(datum = NA) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) +
    labs(title = p_aid_obs_title)
  
  p_aid_policy <- ggplot(a_policy_aid_sf) +
    geom_sf(aes(fill = avg_aid), color = "white", linewidth = 0.2) +
    geom_sf(data = district_edges, fill = NA, color = "grey40", linewidth = 0.25) +
    geom_sf(data = governorate_edges, fill = NA, color = "grey40",   linewidth = 0.25) +
    scale_fill_viridis_c("Average aid", limits = rng_aid) +
    facet_wrap(~ policy_length, nrow = 1) +
    coord_sf(datum = NA) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) +
    labs(title = paste0(event_type,"Aid under learned policy (by horizon M)"))
  
  p_out_obs <- ggplot(a_obs_out_sf) +
    geom_sf(aes(fill = outcome_obs), color = "grey30", linewidth = 0.2) +
    geom_sf(data = district_edges, fill = NA, color = "grey40", linewidth = 0.25) +
    geom_sf(data = baghdad_outline, fill = NA, color = "black", linewidth = 0.7, linejoin = "round") +
    geom_sf(data = mosul_basrah_dist, fill = NA, color = "black", linewidth = 0.7, linejoin = "round") +
    scale_fill_viridis_c("") +
    coord_sf(datum = NA) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) +
    labs(title = p_out_obs_title)
  
  p_out_policy <- ggplot(a_policy_out_sf) +
    geom_sf(aes(fill = avg_outcome), color = "white", linewidth = 0.2) +
    geom_sf(data = district_edges, fill = NA, color = "grey40", linewidth = 0.25) +
    geom_sf(data = governorate_edges, fill = NA, color = "grey40",   linewidth = 0.25) +
    scale_fill_viridis_c("Average outcome", limits = rng_out) +
    facet_wrap(~ policy_length, nrow = 1) +
    coord_sf(datum = NA) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) +
    labs(title = paste0(event_type,"Outcome under learned policy (by horizon M)"))
  
  list(p_aid_obs = p_aid_obs, p_aid_policy = p_aid_policy,
       p_out_obs = p_out_obs, p_out_policy = p_out_policy)
}

# ---- mapping District -> Governorate (no geometry) ----
district_to_governorate <- function(iraq_district) {
  iraq_district %>%
    sf::st_drop_geometry() %>%
    transmute(
      ADM3NAME = trimws(ADM3NAME),
      ADM2NAME = trimws(ADM2NAME)
    ) %>%
    distinct(ADM3NAME, .keep_all = TRUE)
}

# =========================
# GOVERNORATE BARS

top_change_bars <- function(diff_tbl, value_col, N = 5, title = "Largest changes") {
  v <- rlang::enquo(value_col)

  top_tbl <- diff_tbl %>%
    filter(!is.na(!!v)) %>%
    group_by(policy_length) %>%
    { bind_rows(
      slice_min(., order_by = !!v, n = N, with_ties = FALSE),
      slice_max(., order_by = !!v, n = N, with_ties = FALSE)
    ) } %>%
    distinct(policy_length, ADM3NAME, .keep_all = TRUE) %>%
    ungroup() %>%
    mutate(sign = case_when(
      !!v > 0 ~ "Increase",
      !!v < 0 ~ "Decrease",
      TRUE ~ "No change"
    )) %>%
    group_by(policy_length) %>%
    mutate(ADM3NAME_infacet = fct_rev(reorder_within(ADM3NAME, !!v, policy_length))) %>%
    ungroup()

  ggplot(top_tbl, aes(x = !!v, y = ADM3NAME_infacet, fill = sign)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
    scale_fill_manual(
      values = c("Decrease" = "#2c7bb6", "Increase" = "#d7191c", "No change" = "grey70"),
      drop = FALSE
    ) +
    facet_wrap(~ policy_length, scales = "free_y", ncol = 4) +
    scale_y_reordered() +
    labs(title = title, x = "Change", y = NULL, fill = NULL) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      axis.text.y = element_text(size = 9)
    )
}

aggregate_diff_to_governorate <- function(diff_tbl, iraq_district, value_col,
                                          weight_col = NULL,
                                          agg = c("sum", "mean")) {
  agg <- match.arg(agg)
  v <- rlang::ensym(value_col)
  wq <- if (is.null(weight_col)) NULL else rlang::ensym(weight_col)
  key <- district_to_governorate(iraq_district)

  df <- diff_tbl %>%
    mutate(ADM3NAME = trimws(ADM3NAME)) %>%
    left_join(key, by = "ADM3NAME") %>%
    filter(!is.na(ADM2NAME))

  summarise_value <- function(x, w) {
    if (agg == "sum") {
      return(sum(x, na.rm = TRUE))
    }
    if (is.null(w)) {
      return(mean(x, na.rm = TRUE))
    }
    if (all(!is.finite(w)) || sum(w, na.rm = TRUE) == 0) {
      return(mean(x, na.rm = TRUE))
    }
    sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
  }

  df %>%
    group_by(policy_length, ADM2NAME) %>%
    summarise(
      value = summarise_value(!!v, if (is.null(wq)) NULL else !!wq),
      .groups = "drop"
    ) %>%
    mutate(sign = case_when(
      value > 0 ~ "Increase",
      value < 0 ~ "Decrease",
      TRUE ~ "No change"
    ))
}

aggregate_diff_to_governorate_ci <- function(diff_tbl, iraq_district, value_col, lb_col, ub_col,
                                             weight_col = NULL,
                                             agg = c("sum", "mean")) {
  agg <- match.arg(agg)
  v <- rlang::ensym(value_col)
  vl <- rlang::ensym(lb_col)
  vu <- rlang::ensym(ub_col)
  wq <- if (is.null(weight_col)) NULL else rlang::ensym(weight_col)
  key <- district_to_governorate(iraq_district)

  df <- diff_tbl %>%
    mutate(ADM3NAME = trimws(ADM3NAME)) %>%
    left_join(key, by = "ADM3NAME") %>%
    filter(!is.na(ADM2NAME))

  summarise_value <- function(x, w) {
    if (agg == "sum") {
      return(sum(x, na.rm = TRUE))
    }
    if (is.null(w)) {
      return(mean(x, na.rm = TRUE))
    }
    if (all(!is.finite(w)) || sum(w, na.rm = TRUE) == 0) {
      return(mean(x, na.rm = TRUE))
    }
    sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
  }

  df %>%
    group_by(policy_length, ADM2NAME) %>%
    summarise(
      value = summarise_value(!!v, if (is.null(wq)) NULL else !!wq),
      lb = summarise_value(!!vl, if (is.null(wq)) NULL else !!wq),
      ub = summarise_value(!!vu, if (is.null(wq)) NULL else !!wq),
      .groups = "drop"
    ) %>%
    mutate(sign = case_when(
      value > 0 ~ "Increase",
      value < 0 ~ "Decrease",
      TRUE ~ "No change"
    ))
}

build_gov_out_diff_ci_from_tmpgov <- function(iraq_district, data, out, Y_mat_name,
                                               weight_col = NULL,
                                               agg = c("sum", "mean")) {
  agg <- match.arg(agg)
  district_names <- colnames(data$aid_array_week)
  key <- district_to_governorate(iraq_district)

  obs_tbl <- tibble(
    ADM3NAME = district_names,
    obs = as.numeric(colMeans(data[[Y_mat_name]], na.rm = TRUE))
  ) %>%
    left_join(key, by = "ADM3NAME")

  if (!is.null(weight_col)) {
    w <- data[[weight_col]]
    obs_tbl$weight <- if (is.matrix(w)) colMeans(w, na.rm = TRUE) else as.numeric(w)
  } else {
    obs_tbl$weight <- NA_real_
  }

  obs_gov <- obs_tbl %>%
    filter(!is.na(ADM2NAME)) %>%
    group_by(ADM2NAME) %>%
    summarise(
      n_district = dplyr::n(),
      observed = if (agg == "sum") {
        sum(obs, na.rm = TRUE)
      } else if (is.null(weight_col)) {
        mean(obs, na.rm = TRUE)
      } else {
        sum(obs * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE)
      },
      .groups = "drop"
    )

  map_dfr(seq_along(out$tmp_gov), function(M) {
    govs <- if (!is.null(out$tmp_gov[[M]]$ADM2NAME)) {
      trimws(out$tmp_gov[[M]]$ADM2NAME)
    } else {
      unique(trimws(iraq_district$ADM2NAME))
    }
    tibble(
      ADM2NAME = govs,
      policy_length = factor(paste0("M", M), levels = paste0("M", seq_along(out$tmp_gov))),
      policy_value_mean = as.numeric(out$tmp_gov[[M]]$V_est) / M,
      lb_mean = as.numeric(out$tmp_gov[[M]]$lb) / M,
      ub_mean = as.numeric(out$tmp_gov[[M]]$ub) / M
    )
  }) %>%
    left_join(obs_gov, by = "ADM2NAME") %>%
    mutate(
      scale = if (agg == "sum") n_district else 1,
      policy_value = policy_value_mean * scale,
      lb = lb_mean * scale,
      ub = ub_mean * scale,
      value = policy_value - observed,
      lb = lb - observed,
      ub = ub - observed,
      sign = case_when(
        value > 0 ~ "Increase",
        value < 0 ~ "Decrease",
        TRUE ~ "No change"
      )
    )
}

top_change_bars_governorate <- function(gov_diff_tbl, N = 5,
                                        title = "Largest changes by governorate") {
  top_tbl <- gov_diff_tbl %>%
    group_by(policy_length) %>%
    { bind_rows(
      slice_min(., order_by = value, n = N, with_ties = FALSE),
      slice_max(., order_by = value, n = N, with_ties = FALSE)
    ) } %>%
    distinct(policy_length, ADM2NAME, .keep_all = TRUE) %>%
    ungroup() %>%
    group_by(policy_length) %>%
    mutate(ADM2_infacet = fct_rev(reorder_within(ADM2NAME, value, policy_length))) %>%
    ungroup()

  ggplot(top_tbl, aes(x = value, y = ADM2_infacet, fill = sign)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
    scale_fill_manual(
      values = c("Decrease" = "#2c7bb6", "Increase" = "#d7191c", "No change" = "grey70"),
      drop = FALSE
    ) +
    annotate("segment", x = -Inf, xend = Inf, y = N + 0.5, yend = N + 0.5,
             linetype = "dashed", color = "black") +
    facet_wrap(~ policy_length, scales = "free_y", ncol = 4) +
    scale_y_reordered() +
    labs(title = title, x = "Change (policy - observed)", y = NULL, fill = NULL) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      axis.text.y = element_text(size = 9)
    )
}

top_change_bars_governorate_ci <- function(gov_diff_tbl_ci, N = 5,
                                           title = "Largest changes by governorate") {
  top_tbl <- gov_diff_tbl_ci %>%
    group_by(policy_length) %>%
    { bind_rows(
      slice_min(., order_by = value, n = N, with_ties = FALSE),
      slice_max(., order_by = value, n = N, with_ties = FALSE)
    ) } %>%
    distinct(policy_length, ADM2NAME, .keep_all = TRUE) %>%
    ungroup() %>%
    group_by(policy_length) %>%
    mutate(ADM2_infacet = fct_rev(reorder_within(ADM2NAME, value, policy_length))) %>%
    ungroup()

  ggplot(top_tbl, aes(x = value, y = ADM2_infacet, fill = sign)) +
    geom_col(width = 0.7) +
    geom_errorbarh(aes(xmin = lb, xmax = ub), height = 0.25, linewidth = 0.5, color = "black") +
    scale_fill_manual(
      values = c("Decrease" = "#2c7bb6", "Increase" = "#d7191c", "No change" = "grey70"),
      drop = FALSE
    ) +
    annotate("segment", x = -Inf, xend = Inf, y = N + 0.5, yend = N + 0.5,
             linetype = "dashed", color = "black") +
    facet_wrap(~ policy_length, scales = "free_y", ncol = 4) +
    scale_y_reordered() +
    labs(title = title, x = "Change (policy - observed)", y = NULL, fill = NULL) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      axis.text.y = element_text(size = 9)
    )
}

.pretty_coef_names <- function(nm) {
  nm <- sub("^tf_array_US_missing$", "tf array US (miss)", nm)
  nm <- sub("_week_log_pop$", " (log) x pop", nm)
  nm <- sub("_week_log_missing$", " (log) x tf (miss)", nm)
  nm <- sub("_week_log$", " (log)", nm)
  gsub("_", " ", nm, fixed = TRUE)
}

policy_coef_table <- function(out, label) {
  imap_dfr(out$res_policy, function(res, M) {
    as.data.frame(res$policy) %>%
      rownames_to_column("Time t") %>%
      mutate(
        M = M,
        `% treated` = if (!is.null(res$avg_treat_district)) {
          scales::percent(res$avg_treat_district, accuracy = 0.1)
        } else {
          NA_character_
        },
        Outcome = label,
        .before = 1
      )
  }) %>%
    relocate(Outcome, M, `Time t`)
}

value_table_from_out <- function(out, y_mat, label) {
  imap_dfr(out$res_policy, function(res, M) {
    tibble(
      M = M,
      Outcome = label,
      Vhat = res$value_stochastic$V_est,
      Lower = res$value_stochastic$CI_lb,
      Upper = res$value_stochastic$CI_ub,
      Observed = M * mean(y_mat, na.rm = TRUE)
    )
  })
}

coef_table_render_dual <- function(tbl_policy,
                                   caption,
                                   M_select = 4,
                                   drop_M_col = TRUE,
                                   print_latex = TRUE) {
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("Package `kableExtra` is required to render application tables.")
  }

  tbl <- tbl_policy %>%
    filter(M == M_select) %>%
    select(-Outcome) %>%
    arrange(`Time t`)

  if (drop_M_col) {
    tbl <- select(tbl, -M)
  }

  id_keep <- c("M", "Time t", "% treated")
  coef_cols <- setdiff(names(tbl), id_keep)
  tbl <- mutate(tbl, across(all_of(coef_cols), ~ round(., 3)))
  names(tbl) <- .pretty_coef_names(names(tbl))

  kb_html <- kableExtra::kbl(
    tbl, format = "html", booktabs = TRUE, escape = FALSE,
    caption = caption
  ) |>
    kableExtra::kable_styling(full_width = FALSE, font_size = 11)

  if (print_latex) {
    kb_tex <- kableExtra::kbl(
      tbl, format = "latex", booktabs = TRUE, escape = FALSE,
      caption = caption
    ) |>
      kableExtra::kable_styling(
        full_width = FALSE, font_size = 11,
        latex_options = "hold_position"
      )
    cat(as.character(kb_tex), "\n")
  }

  invisible(kb_html)
}


# ---- Public application workflow utilities ----

read_dataverse_rds <- function(fileid, server = "https://dataverse.harvard.edu") {
  tf <- tempfile(fileext = ".rds")
  httr::GET(
    paste0(server, "/api/access/datafile/", fileid),
    query = list(format = "original"),
    httr::write_disk(tf, overwrite = TRUE)
  ) |>
    httr::stop_for_status()
  readRDS(tf)
}

count_points_by_window <- function(hfr, name, winlist, n_time = 499) {
  point_patterns <- hfr[, name, drop = TRUE]
  n_districts <- length(winlist)
  point_counts <- matrix(0, nrow = n_districts, ncol = n_time)
  for (t in seq_len(n_time)) {
    p <- point_patterns[[t]]
    for (d in seq_len(n_districts)) {
      point_counts[d, t] <- sum(spatstat.geom::inside.owin(p, w = winlist[[d]]))
    }
  }
  point_counts
}

prepare_application_data <- function(aid_csv, troop_density_csv) {
  data("airstrikes", package = "geocausal")
  data("insurgencies", package = "geocausal")
  data("airstrikes_base", package = "geocausal")
  data("iraq_window", package = "geocausal")

  iraq_district <- read_dataverse_rds(8138903)
  iraq_district_sf <- sf::st_as_sf(iraq_district)
  winlist <- lapply(seq_len(nrow(iraq_district_sf)), function(z) spatstat.geom::as.owin(iraq_district_sf$geometry[[z]]))
  names(winlist) <- iraq_district_sf$ADM3NAME
  district_names <- names(winlist)

  aid <- data.table::fread(aid_csv)
  aid <- subset(aid, Actual_Start != "")
  aid[, USE_DATE := lubridate::dmy(Actual_Start)]
  data.table::setorder(aid, USE_DATE)

  aid_district <- aid[
    USE_DATE >= as.Date("2007-01-01") & USE_DATE <= as.Date("2008-07-30"),
    .(Total_cost = sum(Construction_Cost, na.rm = TRUE), Count = .N),
    by = .(District, USE_DATE)
  ]
  aid_district <- subset(aid_district, Total_cost > 0)

  all_dates <- seq(min(aid_district$USE_DATE), max(aid_district$USE_DATE), by = 1)
  all_aid <- data.table::data.table(expand.grid(District = district_names, USE_DATE = all_dates))
  all_aid <- merge(all_aid, aid_district, by = c("District", "USE_DATE"), all.x = TRUE)
  aid_district <- subset(all_aid, USE_DATE >= as.Date("2007-02-23") & USE_DATE <= as.Date("2008-07-05"))
  aid_district[, New_aid := ifelse(is.na(Total_cost), 0, 1)]
  aid_district[, Total_cost := ifelse(is.na(Total_cost), 0, Total_cost)]
  aid_district[, Count := ifelse(is.na(Count), 0, Count)]
  aid_district[, District := factor(District, levels = district_names)]

  aid_matrix <- data.table::dcast(aid_district, District ~ USE_DATE, value.var = "New_aid")
  aid_array <- t(as.matrix(aid_matrix[, -1]))
  colnames(aid_array) <- aid_matrix$District

  airstr <- airstrikes
  insurg <- insurgencies
  airstr$time <- as.numeric(airstr$date - min(airstr$date) + 1)
  insurg$time <- as.numeric(insurg$date - min(insurg$date) + 1)

  treatment_hfr <- geocausal::get_hfr(
    data = airstr, col = "type", window = iraq_window,
    time_col = "time", time_range = c(1, 499),
    coordinates = c("longitude", "latitude"), combine = TRUE
  )
  outcome_hfr <- geocausal::get_hfr(
    data = insurg, col = "type", window = iraq_window,
    time_col = "time", time_range = c(1, 499),
    coordinates = c("longitude", "latitude"), combine = TRUE
  )

  dat_hfr <- spatstat.geom::cbind.hyperframe(treatment_hfr, outcome_hfr[, -1])
  names(dat_hfr)[names(dat_hfr) == "all_combined"] <- "all_airstrike"
  names(dat_hfr)[names(dat_hfr) == "all_combined.1"] <- "all_outcome"

  res_list <- list()
  time_splines <- splines::ns(seq_len(499), df = 3)
  res_list$time_1 <- time_splines[, 1]
  res_list$time_2 <- time_splines[, 2]
  res_list$time_3 <- time_splines[, 3]

  min_s <- as.numeric(as.Date("2007-03-25") - min(airstr$date) + 1)
  max_s <- as.numeric(as.Date("2008-01-01") - min(airstr$date) + 1)
  res_list$surge <- ifelse(dat_hfr$time >= min_s & dat_hfr$time <= max_s, 1, 0)

  for (nm in c("SAF", "all_outcome", "SOF", "all_airstrike", "IED")) {
    res_list[[nm]] <- t(count_points_by_window(dat_hfr, nm, winlist = winlist))
    for (lag in c(1, 7, 14, 28, 30)) {
      res_list <- add_lagged_sum(res_list, colname = nm, lag = lag)
    }
  }

  cities_dist <- read_dataverse_rds(8079278)
  cities_tmp <- lapply(seq_len(5), function(x) {
    tmp <- sapply(seq_along(winlist), function(d) mean(cities_dist$distance[[x]][winlist[[d]]]))
    exp(-(2 * x) * tmp)
  })
  res_list$cities_dist <- Reduce("+", cities_tmp)

  rivers_dist <- read_dataverse_rds(8079275)
  routes_dist <- read_dataverse_rds(8079274)
  ethnicity <- read_dataverse_rds(8079276)

  get_mode <- function(x) {
    ux <- unique(x[!is.na(x)])
    ux[which.max(tabulate(match(x, ux)))]
  }
  res_list$rivers_dist <- sapply(seq_along(winlist), function(d) exp(-3 * mean(rivers_dist[winlist[[d]]])))
  res_list$routes_dist <- sapply(seq_along(winlist), function(d) exp(-3 * mean(routes_dist[winlist[[d]]])))
  res_list$pop <- sapply(seq_along(winlist), function(d) log(get_mode(ethnicity$total_pop[winlist[[d]]])))

  aid_matrix <- data.table::dcast(aid_district, District ~ USE_DATE, value.var = "Total_cost")
  aid_array_cost <- t(as.matrix(aid_matrix[, -1]))
  colnames(aid_array_cost) <- aid_matrix$District
  aid_array_costamount <- aid_array_cost
  aid_array_cost1 <- aid_array_cost <= 1e5 & aid_array_cost > 0
  aid_array_cost2 <- aid_array_cost > 1e5
  aid_array_cost[aid_array_cost1 == 1] <- 1
  aid_array_cost[aid_array_cost2 == 1] <- 2

  aid_matrix <- data.table::dcast(aid_district, District ~ USE_DATE, value.var = "Count")
  aid_array_count <- t(as.matrix(aid_matrix[, -1]))
  colnames(aid_array_count) <- aid_matrix$District

  res_list$aid_array <- aid_array
  res_list$aid_costamount <- aid_array_costamount
  res_list$aid_array_cost <- aid_array_cost
  res_list$aid_array_cost1 <- aid_array_cost1
  res_list$aid_array_cost2 <- aid_array_cost2
  res_list$aid_array_count <- aid_array_count

  for (lag in c(1, 7, 15, 14, 28, 30)) {
    res_list <- add_lagged_max(res_list, colname = "aid_array", lag = lag)
    res_list <- add_lagged_sum(res_list, colname = "aid_array_count", lag = lag)
    res_list <- add_lagged_sum(res_list, colname = "aid_costamount", lag = lag)
  }
  res_list <- take_log(
    res_list,
    c("aid_costamount", "prev_aid_costamount_1", "prev_aid_costamount_7",
      "prev_aid_costamount_15", "prev_aid_costamount_30")
  )

  tf_raw <- data.table::fread(troop_density_csv, colClasses = list(integer = "week"))
  required_tf_columns <- c("District", "week", "TForce")
  missing_tf_columns <- setdiff(required_tf_columns, names(tf_raw))
  if (length(missing_tf_columns) > 0L) {
    stop("Troop-density input is missing required columns: ",
         paste(missing_tf_columns, collapse = ", "))
  }
  origin_date <- lubridate::ymd("2005-05-02") - lubridate::weeks(2358 - 1)
  tf_raw[, week_start := origin_date + lubridate::weeks(week - 1)]

  all_dates <- seq(as.Date("2007-02-23"), as.Date("2008-07-05"), by = "day")
  tf_daily <- data.table::CJ(District = district_names, USE_DATE = all_dates)
  tf_daily[, week_start := lubridate::floor_date(USE_DATE, unit = "week", week_start = 1)]

  make_tf_array <- function(raw_sub) {
    dt <- merge(
      tf_daily,
      raw_sub[, .(District, week_start, TForce)],
      by = c("District", "week_start"),
      all.x = TRUE
    )
    data.table::setorder(dt, District, USE_DATE)
    dt[, TForce := zoo::na.locf(TForce, na.rm = FALSE), by = District]
    dt[, TForce := zoo::na.locf(TForce, fromLast = TRUE), by = District]
    wide <- data.table::dcast(dt, District ~ USE_DATE, value.var = "TForce")
    mat <- t(as.matrix(wide[, -1, with = FALSE]))
    mat2 <- mat[, match(district_names, wide$District)]
    colnames(mat2) <- district_names
    mat2
  }

  # CSM_filter.csv is already restricted to observations with british == 0.
  tf_array_US <- make_tf_array(tf_raw)
  res_list$tf_array_US <- tf_array_US
  res_list$tf_array_US[is.na(tf_array_US)] <- 0
  res_list$tf_array_US_missing <- is.na(tf_array_US)

  res_list$aid_array_week <- apply(res_list$aid_array, 2, function(x) {
    zoo::rollapply(x, width = 7, FUN = function(v) as.integer(any(v == 1)), align = "left", fill = NA)
  })
  res_list$SAF_week <- apply(res_list$SAF, 2, function(x) {
    zoo::rollapply(x, width = 7, FUN = sum, align = "left", fill = NA)
  })
  res_list$IED_week <- apply(res_list$IED, 2, function(x) {
    zoo::rollapply(x, width = 7, FUN = sum, align = "left", fill = NA)
  })
  res_list$all_airstrike_week <- apply(res_list$all_airstrike, 2, function(x) {
    zoo::rollapply(x, width = 7, FUN = sum, align = "left", fill = NA)
  })

  ps_vars <- c(
    "prev_SAF_7", "prev_SAF_14", "prev_SAF_28",
    "prev_all_outcome_7", "prev_all_outcome_14", "prev_all_outcome_28",
    "prev_SOF_7", "prev_SOF_14", "prev_SOF_28",
    "prev_all_airstrike_7", "prev_all_airstrike_14", "prev_all_airstrike_28",
    "pop", "time_1", "time_2", "time_3", "surge", "tf_array_US", "tf_array_US_missing",
    "rivers_dist", "routes_dist", "cities_dist",
    "prev_aid_costamount_7", "prev_aid_costamount_14", "prev_aid_costamount_28",
    "prev_aid_array_7", "prev_aid_array_14", "prev_aid_array_28"
  )

  keep_time <- seq(1, 493, by = 7)
  keep_index <- which(!district_names %in% names(which(colSums(res_list$aid_array_week[keep_time, ], na.rm = TRUE) <= 0)))
  n_time <- length(keep_time)
  n_district <- length(keep_index)

  ps_df <- data.table::data.table(
    aid = as.vector(t(res_list$aid_array[keep_time, keep_index])),
    aid_cost = as.vector(t(res_list$aid_array_cost[keep_time, keep_index])),
    aid_cost1 = as.vector(t(res_list$aid_array_cost1[keep_time, keep_index])),
    aid_cost2 = as.vector(t(res_list$aid_array_cost2[keep_time, keep_index])),
    time = rep(seq_len(n_time), each = n_district),
    district = rep(keep_index, times = n_time)
  )

  for (v in ps_vars) {
    val <- res_list[[v]]
    if (is.matrix(val) && all(dim(val) == c(499, ncol(aid_array)))) {
      ps_df[[v]] <- as.vector(t(val[keep_time, keep_index]))
    } else if (length(val) == ncol(aid_array)) {
      ps_df[[v]] <- rep(val[keep_index], times = n_time)
    } else if (length(val) == 499) {
      ps_df[[v]] <- rep(val[keep_time], each = n_district)
    } else {
      stop(sprintf("Dimension mismatch in propensity-score variable: %s", v))
    }
  }

  policy_vars <- c(
    "aid_array_week", "all_airstrike_week", "SAF_week", "IED_week",
    "tf_array_US", "tf_array_US_missing", "pop", "cities_dist"
  )
  data <- list()
  for (v in policy_vars) {
    val <- res_list[[v]]
    if (is.matrix(val) && all(dim(val) == c(499, ncol(aid_array)))) {
      data[[v]] <- val[keep_time, keep_index]
    } else if (length(val) == ncol(aid_array)) {
      data[[v]] <- val[keep_index]
    } else if (length(val) == 499) {
      data[[v]] <- val[keep_time]
    } else {
      stop(sprintf("Dimension mismatch in policy variable: %s", v))
    }
  }
  data <- add_log_binary_std(data, names(data)[2:4], offset = 1)

  list(
    data = data,
    iraq_district = iraq_district_sf,
    iraq_district_analysis = iraq_district_sf[iraq_district_sf$ADM3NAME %in% colnames(data$aid_array_week), ],
    ps_df = ps_df,
    ps_vars = ps_vars,
    balance_vars = ps_vars,
    keep_time = keep_time,
    keep_index = keep_index,
    n_time = n_time,
    n_district = n_district
  )
}

fit_application_propensity <- function(application, ps_vars = NULL) {
  ps_df <- application$ps_df
  if (is.null(ps_vars)) {
    ps_vars <- application$ps_vars
  } else {
    ps_vars <- unique(ps_vars)
  }
  missing_vars <- setdiff(ps_vars, names(ps_df))
  if (length(missing_vars) > 0) {
    stop(
      "The requested propensity-score variables are not in ps_df: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  application$ps_vars <- ps_vars

  ps_formula <- paste("aid ~", paste(ps_vars, collapse = " + "))
  ps_fit <- geepack::geeglm(
    as.formula(ps_formula),
    data = as.data.frame(ps_df),
    id = time,
    family = binomial(link = "logit"),
    corstr = "independence"
  )

  ps_df$ps_hat <- predict(ps_fit, type = "response")
  ps_df$weights <- 1 / (ps_df$ps_hat * ps_df$aid + (1 - ps_df$ps_hat) * (1 - ps_df$aid))

  data <- application$data
  data$ps <- t(matrix(ps_df$ps_hat, nrow = application$n_district, ncol = application$n_time))
  colnames(data$ps) <- colnames(data$aid_array_week)
  rownames(data$ps) <- rownames(data$aid_array_week)

  application$data <- data
  application$ps_df <- ps_df
  application$ps_formula <- ps_formula
  application$ps_fit <- ps_fit
  application
}

add_application_policy_terms <- function(data) {
  data$SAF_week_log_pop <- data$SAF_week_log * data$pop
  data$SAF_week_log_missing <- data$SAF_week_log * data$tf_array_US_missing
  data$IED_week_log_pop <- data$IED_week_log * data$pop
  data$IED_week_log_missing <- data$IED_week_log * data$tf_array_US_missing
  data
}

application_policy_spec <- function(include_population = FALSE) {
  X_cols <- c(
    "all_airstrike_week_log", "all_airstrike_week_binary",
    "SAF_week_log", "SAF_week_binary",
    "tf_array_US", "tf_array_US_missing",
    "IED_week_log", "IED_week_binary",
    "SAF_week_log_missing", "SAF_week_log_pop",
    "IED_week_log_missing", "IED_week_log_pop"
  )
  Z_cols <- c("pop", "cities_dist")
  policy_cols_SAF <- c("W", "tf_array_US", "all_airstrike_week_log", "tf_array_US_missing", "SAF_week_log")
  policy_cols_IED <- c("W", "tf_array_US", "all_airstrike_week_log", "tf_array_US_missing", "IED_week_log")
  if (isTRUE(include_population)) {
    policy_cols_SAF <- c(policy_cols_SAF, "pop", "SAF_week_log_pop")
    policy_cols_IED <- c(policy_cols_IED, "pop", "IED_week_log_pop")
  }

  missing_policy_cols <- setdiff(c(policy_cols_SAF, policy_cols_IED), c("W", X_cols, Z_cols))
  if (length(missing_policy_cols) > 0) {
    stop(
      "The requested policy variables are not available in X_cols/Z_cols: ",
      paste(missing_policy_cols, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    X_cols = X_cols,
    Z_cols = Z_cols,
    policy_cols_SAF = policy_cols_SAF,
    policy_cols_IED = policy_cols_IED
  )
}

make_initial_par_list <- function(policy_cols, n_periods = 4) {
  rep(list(c(-1, rep(0, length(policy_cols)))), n_periods)
}

make_warm_start_initial_par_list <- function(previous_out, policy_cols, n_periods = 4) {
  initial_par_list <- make_initial_par_list(policy_cols, n_periods)
  target_names <- c("intercept", policy_cols)

  if (is.null(previous_out) || is.null(previous_out$res_policy)) {
    return(initial_par_list)
  }

  n_available <- min(n_periods, length(previous_out$res_policy))
  for (m in seq_len(n_available)) {
    policy <- previous_out$res_policy[[m]]$policy
    if (is.null(policy)) next

    policy <- as.matrix(policy)
    old_par <- policy[nrow(policy), ]
    new_par <- initial_par_list[[m]]
    names(new_par) <- target_names

    common_names <- intersect(names(old_par), target_names)
    new_par[common_names] <- old_par[common_names]
    initial_par_list[[m]] <- unname(new_par)
  }

  initial_par_list
}

make_ps_overlap_plot <- function(ps_df) {
  ggplot2::ggplot(ps_df, ggplot2::aes(x = ps_hat, fill = factor(aid))) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(x = "Estimated Propensity Score", fill = "Aid", title = "") +
    ggplot2::theme_minimal()
}

make_application_asd <- function(ps_df, balance_vars, ps_model_vars = NULL) {
  missing_vars <- setdiff(balance_vars, names(ps_df))
  if (length(missing_vars) > 0) {
    stop(
      "The requested balance-check variables are not in ps_df: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  asd <- compute_ASD(as.data.frame(ps_df), balance_vars, "aid", "weights")
  asd_df <- dplyr::bind_rows(
    tidy_asd_list(asd$ASD, "Before"),
    tidy_asd_list(asd$ASD_weighted, "After")
  )
  if (!is.null(ps_model_vars)) {
    asd_df$InPropensityModel <- asd_df$Variable %in% ps_model_vars
  }
  asd_df
}

make_asd_plot <- function(asd_df) {
  asd_labels <- c(
    prev_SAF_7 = "prev SAF 1 week",
    prev_SAF_14 = "prev SAF 2 weeks",
    prev_SAF_28 = "prev SAF 4 weeks",
    prev_all_outcome_7 = "prev outcome 1 week",
    prev_all_outcome_14 = "prev outcome 2 weeks",
    prev_all_outcome_28 = "prev outcome 4 weeks",
    prev_SOF_7 = "prev SOF 1 week",
    prev_SOF_14 = "prev SOF 2 weeks",
    prev_SOF_28 = "prev SOF 4 weeks",
    prev_all_airstrike_7 = "prev airstrike 1 week",
    prev_all_airstrike_14 = "prev airstrike 2 weeks",
    prev_all_airstrike_28 = "prev airstrike 4 weeks",
    pop = "population density",
    time_1 = "time spline 1",
    time_2 = "time spline 2",
    time_3 = "time spline 3",
    surge = "surge",
    tf_array_US = "troop",
    tf_array_US_missing = "troop NA",
    rivers_dist = "rivers dist",
    routes_dist = "routes dist",
    cities_dist = "cities dist",
    prev_aid_costamount_7 = "prev aid amount 1 week",
    prev_aid_costamount_14 = "prev aid amount 2 weeks",
    prev_aid_costamount_28 = "prev aid amount 4 weeks",
    prev_aid_array_7 = "prev aid 1 week",
    prev_aid_array_14 = "prev aid 2 weeks",
    prev_aid_array_28 = "prev aid 4 weeks"
  )
  ggplot2::ggplot(asd_df, ggplot2::aes(x = ASD, y = reorder(VariableLevel, ASD), color = Type)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(title = "", x = "Absolute Standardized Difference", y = "", color = "Weighting") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed", color = "gray") +
    ggplot2::scale_y_discrete(labels = function(x) {
      dplyr::recode(x, !!!asd_labels, .default = x)
    })
}

save_application_diagnostics <- function(ps_df, balance_vars, plot_dir,
                                         ps_model_vars = NULL) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  p_ps_overlap <- make_ps_overlap_plot(ps_df)
  ggplot2::ggsave(file.path(plot_dir, "p_ps_overlap.pdf"), p_ps_overlap,
                  width = 4.7, height = 3, units = "in", dpi = 600)
  asd_df <- make_application_asd(ps_df, balance_vars, ps_model_vars = ps_model_vars)

  p_asd <- make_asd_plot(asd_df)
  ggplot2::ggsave(file.path(plot_dir, "p_ASD.pdf"), p_asd,
                  width = 4.7, height = 5, units = "in", dpi = 600)

  print(p_ps_overlap)
  print(p_asd)

  invisible(asd_df)
}

make_application_column_header <- function(label) {
  ggplot() +
    annotate("text", x = 0.62, y = 0.72, label = label, size = 3.4) +
    annotate("segment", x = 0.28, xend = 0.98, y = 0.43, yend = 0.43,
             linewidth = 0.35) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_application_group_label <- function(label) {
  ggplot() +
    annotate("segment", x = 0.68, xend = 0.68, y = 0.08, yend = 0.92, linewidth = 0.35) +
    annotate("segment", x = 0.68, xend = 0.88, y = 0.92, yend = 0.92, linewidth = 0.35) +
    annotate("segment", x = 0.68, xend = 0.88, y = 0.08, yend = 0.08, linewidth = 0.35) +
    annotate("text", x = 0.22, y = 0.5, label = label, angle = 90, size = 3.8) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_application_blank_plot <- function() {
  ggplot() + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}

application_governorate_label_specs <- function() {
  tibble::tribble(
    ~ADM2NAME,        ~label_x, ~label_y, ~hjust, ~vjust, ~draw_segment,
    "Anbar",             41.15,   32.25,   0.5,    0.5, FALSE,
    "Babylon",           50.05,   32.95,   1.0,    0.5, TRUE,
    "Baghdad",           50.05,   33.70,   1.0,    0.5, TRUE,
    "Basrah",            50.05,   29.20,   1.0,    0.5, FALSE,
    "Dahuk",             44.70,   38.70,   0.5,    0.5, FALSE,
    "Diyala",            50.05,   34.25,   1.0,    0.5, TRUE,
    "Erbil",             45.25,   36.75,   0.0,    0.5, TRUE,
    "Kerbala",           43.15,   32.20,   1.0,    0.5, TRUE,
    "Missan",            50.05,   32.35,   1.0,    0.5, FALSE,
    "Muthanna",          46.00,   29.65,   0.0,    0.5, TRUE,
    "Najaf",             42.55,   29.90,   1.0,    0.5, TRUE,
    "Ninewa",            40.30,   38.65,   0.5,    0.5, FALSE,
    "Qadissiya",         44.4,    31.4,    1.0,    0.5, TRUE,
    "Salah al-Din",      39.05,   37.35,   0.0,    0.5, TRUE,
    "Sulaymaniyah",      46.35,   35.80,   0.0,    0.5, TRUE,
    "Tameem",            50.05,   36.20,   1.0,    0.5, TRUE,
    "Thi-Qar",           49.80,   30.80,   1.0,    0.5, TRUE,
    "Wassit",            50.05,   32.75,   1.0,    0.5, TRUE
  )
}

application_governorate_label_overrides <- function() {
  tibble::tribble(
    ~plot_policy, ~plot_outcome, ~plot_map_type, ~plot_level, ~ADM2NAME,       ~label_x, ~label_y, ~hjust, ~vjust, ~draw_segment, ~segment_x, ~segment_y, ~segment_xend, ~segment_yend,
    NA,           "SAF",         "aid",          "M4",        "Salah al-Din",    38.92,    37.40,   0.0,    0.5, TRUE,              43.8,       34.9,       42.00,         37.32,
    NA,           "SAF",         "aid",          "M4",        "Ninewa",          40.20,    38.75,   0.5,    0.5, FALSE,             NA,         NA,         NA,            NA,
    NA,           "SAF",         "aid",          "M4",        "Dahuk",           45.05,    38.95,   0.5,    0.5, FALSE,             NA,         NA,         NA,            NA,
    NA,           "SAF",         "outcome",      "M4",        "Thi-Qar",         47.85,    30.55,   0.0,    0.5, FALSE,             NA,         NA,         NA,            NA,
    NA,           "SAF",         "outcome",      "M4",        "Missan",          50.10,    32.25,   1.0,    0.5, FALSE,             NA,         NA,         NA,            NA,
    NA,           "SAF",         "outcome",      "M4",        "Basrah",          50.10,    28.85,   1.0,    0.5, FALSE,             NA,         NA,         NA,            NA,
    "wpop",       "SAF",         "outcome",      "M1",        "Salah al-Din",    38.95,    37.55,   0.0,    0.5, TRUE,              43.2,       34.9,       42.00,         37.42,
    "wpop",       "SAF",         "outcome",      "M1",        "Erbil",           45.55,    37.10,   0.0,    0.5, TRUE,              NA,         NA,         45.25,         36.95,
    "wpop",       "IED",         "outcome",      "M4",        "Najaf",           42.10,    29.72,   1.0,    0.5, TRUE,              43.9,       31.1,       42.50,         29.88,
    "wopop",      "IED",         "outcome",      "M4",        "Najaf",           42.10,    29.72,   1.0,    0.5, TRUE,              43.9,       31.1,       42.50,         29.88,
    NA,           NA,            NA,             NA,          "Baghdad",         50.05,    33.70,   1.0,    0.5, TRUE,              NA,         NA,         47.55,         33.65,
    NA,           NA,            NA,             NA,          "Babylon",         50.05,    32.95,   1.0,    0.5, TRUE,              NA,         NA,         47.85,         32.95,
    NA,           NA,            NA,             NA,          "Tameem",          50.05,    36.20,   1.0,    0.5, TRUE,              NA,         NA,         47.70,         36.05,
    NA,           NA,            NA,             NA,          "Diyala",          50.05,    34.35,   1.0,    0.5, TRUE,              NA,         NA,         47.80,         34.30,
    NA,           NA,            NA,             NA,          "Wassit",          50.05,    32.10,   1.0,    0.5, TRUE,              NA,         NA,         47.75,         32.10,
    NA,           NA,            NA,             NA,          "Thi-Qar",         49.80,    30.80,   1.0,    0.5, TRUE,              NA,         NA,         47.95,         30.80,
    NA,           NA,            NA,             NA,          "Najaf",           42.55,    29.90,   1.0,    0.5, TRUE,              NA,         NA,         42.90,         30.05
  )
}

make_application_map_label_sets <- function(selected_tbl, policy_levels) {
  seen <- character()
  label_sets <- stats::setNames(vector("list", length(policy_levels)), policy_levels)
  for (level in rev(policy_levels)) {
    names_level <- selected_tbl %>%
      filter(as.character(policy_length) == level) %>%
      pull(ADM2NAME) %>%
      unique()
    label_sets[[level]] <- setdiff(names_level, seen)
    seen <- union(seen, names_level)
  }
  label_sets
}

make_application_map_labels <- function(iraq_district, label_names,
                                        policy_name = NA_character_,
                                        event_label = NA_character_,
                                        map_type = NA_character_,
                                        policy_level = NA_character_) {
  label_names <- unique(label_names)
  if (length(label_names) == 0) {
    return(tibble::tibble())
  }

  gov_sf <- iraq_district %>%
    filter(ADM2NAME %in% label_names) %>%
    group_by(ADM2NAME) %>%
    summarise(geometry = sf::st_union(geometry), .groups = "drop")

  if (nrow(gov_sf) == 0) {
    return(tibble::tibble())
  }

  anchor <- sf::st_coordinates(sf::st_point_on_surface(gov_sf))
  gov_centers <- sf::st_drop_geometry(gov_sf) %>%
    mutate(anchor_x = anchor[, 1], anchor_y = anchor[, 2])

  specs <- application_governorate_label_specs()
  overrides <- application_governorate_label_overrides() %>%
    mutate(
      specificity =
        as.integer(!is.na(plot_policy)) +
        as.integer(!is.na(plot_outcome)) +
        as.integer(!is.na(plot_map_type)) +
        as.integer(!is.na(plot_level))
    ) %>%
    filter(
      is.na(plot_policy) | plot_policy == policy_name,
      is.na(plot_outcome) | plot_outcome == event_label,
      is.na(plot_map_type) | plot_map_type == map_type,
      is.na(plot_level) | plot_level == policy_level
    ) %>%
    group_by(ADM2NAME) %>%
    slice_max(specificity, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-plot_policy, -plot_outcome, -plot_map_type, -plot_level, -specificity)

  gov_centers %>%
    left_join(specs, by = "ADM2NAME") %>%
    left_join(overrides, by = "ADM2NAME", suffix = c("", "_override")) %>%
    mutate(
      label_x = dplyr::coalesce(label_x_override, label_x),
      label_y = dplyr::coalesce(label_y_override, label_y),
      hjust = dplyr::coalesce(hjust_override, hjust),
      vjust = dplyr::coalesce(vjust_override, vjust),
      draw_segment = dplyr::coalesce(draw_segment_override, draw_segment),
      label_x = ifelse(is.na(label_x), anchor_x + 0.35, label_x),
      label_y = ifelse(is.na(label_y), anchor_y, label_y),
      hjust = ifelse(is.na(hjust), 0, hjust),
      vjust = ifelse(is.na(vjust), 0.5, vjust),
      draw_segment = ifelse(is.na(draw_segment), TRUE, draw_segment),
      text_width = 0.38 * nchar(ADM2NAME),
      segment_x = dplyr::coalesce(segment_x, anchor_x),
      segment_y = dplyr::coalesce(segment_y, anchor_y),
      segment_xend = dplyr::case_when(
        !is.na(segment_xend) ~ segment_xend,
        hjust == 1 & label_x > segment_x ~ label_x - text_width - 0.28,
        hjust == 0 & label_x < segment_x ~ label_x + text_width + 0.28,
        hjust == 0 ~ label_x - 0.28,
        hjust == 1 ~ label_x + 0.28,
        TRUE ~ label_x
      ),
      segment_yend = dplyr::coalesce(segment_yend, label_y)
    )
}

make_application_governorate_label_palette <- function(governorates = NULL) {
  if (is.null(governorates)) {
    governorates <- sort(application_governorate_label_specs()$ADM2NAME)
  }
  n_per_col <- ceiling(length(governorates) / 2)
  label_tbl <- tibble::tibble(
    ADM2NAME = governorates,
    idx = seq_along(governorates),
    col = (idx - 1) %/% n_per_col,
    row = n_per_col - ((idx - 1) %% n_per_col),
    x0 = col * 4.35,
    x1 = x0 + 0.75,
    x_text = x0 + 0.95
  )

  ggplot(label_tbl) +
    geom_segment(
      aes(x = x0, xend = x1, y = row, yend = row),
      linewidth = 0.28,
      color = "grey25"
    ) +
    geom_text(
      aes(x = x_text, y = row, label = ADM2NAME),
      hjust = 0,
      vjust = 0.5,
      size = 3.15,
      color = "grey25"
    ) +
    coord_cartesian(
      xlim = c(-0.05, max(label_tbl$x_text) + 2.8),
      ylim = c(0.45, n_per_col + 0.55),
      expand = FALSE
    ) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

save_application_governorate_label_palette <- function(output_dir,
                                                       filename = "governorate_label_palette.pdf") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  p <- make_application_governorate_label_palette()
  ggplot2::ggsave(
    file.path(output_dir, filename),
    p,
    width = 5.6,
    height = 2.6,
    units = "in",
    dpi = 600
  )
  invisible(file.path(output_dir, filename))
}

make_governorate_change_panel <- function(gov_diff_tbl, policy_level,
                                          x_limits, include_ci = FALSE,
                                          show_legend = FALSE, N = 3) {
  panel_tbl <- gov_diff_tbl %>%
    filter(as.character(policy_length) == policy_level) %>%
    { bind_rows(
      slice_min(., order_by = value, n = N, with_ties = FALSE),
      slice_max(., order_by = value, n = N, with_ties = FALSE)
    ) } %>%
    distinct(ADM2NAME, .keep_all = TRUE) %>%
    mutate(ADM2_infacet = fct_rev(reorder(ADM2NAME, value)))

  p <- ggplot(panel_tbl, aes(x = value, y = ADM2_infacet, fill = sign)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey65") +
    annotate("segment", x = -Inf, xend = Inf, y = N + 0.5, yend = N + 0.5,
             linetype = "dashed", color = "black") +
    scale_fill_manual(
      values = c("Decrease" = "#2c7bb6", "Increase" = "#d7191c", "No change" = "grey70"),
      breaks = c("Decrease", "Increase"),
      drop = FALSE
    ) +
    scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.35, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9.6),
      axis.ticks.length = unit(0, "pt"),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )

  if (include_ci) {
    p <- p + geom_errorbarh(aes(xmin = lb, xmax = ub),
                            height = 0.22, linewidth = 0.42, color = "black")
  }

  p
}

make_district_diff_map_panel <- function(diff_sf, value_col, policy_level,
                                         fill_title, fill_limit,
                                         iraq_district, highlight_edges = NULL,
                                         label_tbl = NULL, show_legend = FALSE) {
  panel_sf <- diff_sf %>% filter(as.character(policy_length) == policy_level)
  map_bbox <- sf::st_bbox(iraq_district)

  ggplot(panel_sf) +
    geom_sf(aes(fill = .data[[value_col]]), color = "white", linewidth = 0.16) +
    geom_sf(data = iraq_district, fill = NA, color = "grey45", linewidth = 0.22) +
    { if (!is.null(highlight_edges) && nrow(highlight_edges) > 0)
      geom_sf(data = highlight_edges, fill = NA, color = "black", linewidth = 0.65) else NULL } +
    { if (!is.null(label_tbl) && nrow(label_tbl) > 0)
      geom_segment(
        data = filter(label_tbl, draw_segment),
        aes(x = segment_x, y = segment_y, xend = segment_xend, yend = segment_yend),
        inherit.aes = FALSE, linewidth = 0.28, color = "grey25"
      ) else NULL } +
    { if (!is.null(label_tbl) && nrow(label_tbl) > 0)
      geom_text(
        data = label_tbl,
        aes(x = label_x, y = label_y, label = ADM2NAME, hjust = hjust, vjust = vjust),
        inherit.aes = FALSE, size = 3.15, color = "grey25"
      ) else NULL } +
    scale_fill_gradient2(
      fill_title,
      low = "#2c7bb6", mid = "white", high = "#d7191c",
      midpoint = 0, limits = c(-fill_limit, fill_limit),
      oob = scales::squish, na.value = "grey90"
    ) +
    coord_sf(
      xlim = c(map_bbox[["xmin"]] - 0.05, map_bbox[["xmax"]] + 0.25),
      ylim = c(map_bbox[["ymin"]] - 0.05, map_bbox[["ymax"]] + 0.15),
      datum = NA,
      expand = FALSE
    ) +
    theme_void(base_size = 9) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.46, "cm"),
      legend.key.width = unit(0.35, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    )
}

make_application_composite_plot <- function(data, iraq_district, out, Y_mat_name,
                                            event_label, N = 3,
                                            policy_name = NA_character_,
                                            layout = c("maps_together", "paired")) {
  layout <- match.arg(layout)
  dat <- prepare_all_data(
    iraq_district = iraq_district,
    data = data,
    out = out,
    Y_mat_name = Y_mat_name
  )

  gov_aid_diff <- aggregate_diff_to_governorate(
    diff_tbl = dat$aid_diff_long,
    iraq_district = iraq_district,
    value_col = diff_aid,
    weight_col = NULL,
    agg = "sum"
  )

  gov_out_diff_ci <- if (!is.null(out$tmp_gov) && length(out$tmp_gov) > 0) {
    build_gov_out_diff_ci_from_tmpgov(
      iraq_district = iraq_district,
      data = data,
      out = out,
      Y_mat_name = Y_mat_name,
      weight_col = NULL,
      agg = "sum"
    )
  } else {
    aggregate_diff_to_governorate_ci(
      diff_tbl = dat$outcome_diff_long,
      iraq_district = iraq_district,
      value_col = diff_outcome,
      lb_col = diff_outcome_lb,
      ub_col = diff_outcome_ub,
      weight_col = NULL,
      agg = "sum"
    )
  }

  selected_aid <- top_bottom_governorates(gov_aid_diff, N = N)
  selected_out <- top_bottom_governorates(gov_out_diff_ci, N = N)
  policy_levels <- paste0("M", seq_along(out$res_policy))

  aid_bar_limit <- max(abs(gov_aid_diff$value), na.rm = TRUE)
  outcome_bar_limit <- max(abs(c(gov_out_diff_ci$lb, gov_out_diff_ci$ub)), na.rm = TRUE)
  aid_map_limit <- max(abs(dat$a_diff_aid_sf$diff_aid), na.rm = TRUE)
  outcome_map_limit <- max(abs(dat$a_diff_out_sf$diff_outcome), na.rm = TRUE)
  for (nm in c("aid_bar_limit", "outcome_bar_limit", "aid_map_limit", "outcome_map_limit")) {
    if (!is.finite(get(nm)) || get(nm) == 0) assign(nm, 1)
  }

  legend_width <- 0.46
  panel_widths <- c(rep(1, length(policy_levels)), legend_width)
  group_label_width <- 0.048

  header_row <- cowplot::plot_grid(
    plotlist = c(
      list(make_application_blank_plot()),
      lapply(seq_along(policy_levels), function(M) {
        make_application_column_header(paste0("M = ", M))
      }),
      list(make_application_blank_plot())
    ),
    nrow = 1,
    rel_widths = c(group_label_width, panel_widths)
  )

  aid_bar_plots <- lapply(policy_levels, function(level) {
    make_governorate_change_panel(
      gov_aid_diff, level,
      x_limits = c(-aid_bar_limit, aid_bar_limit),
      include_ci = FALSE,
      show_legend = FALSE,
      N = N
    )
  })
  aid_bar_legend <- cowplot::get_legend(
    make_governorate_change_panel(
      gov_aid_diff, policy_levels[[1]],
      x_limits = c(-aid_bar_limit, aid_bar_limit),
      include_ci = FALSE,
      show_legend = TRUE,
      N = N
    )
  )

  aid_map_plots <- lapply(policy_levels, function(level) {
    selected_aid_level <- selected_aid %>% filter(as.character(policy_length) == level)
    make_district_diff_map_panel(
      dat$a_diff_aid_sf, "diff_aid", level,
      fill_title = "Aid",
      fill_limit = aid_map_limit,
      iraq_district = iraq_district,
      highlight_edges = make_highlight_edges_clean(iraq_district, selected_aid_level),
      show_legend = FALSE
    )
  })
  aid_map_legend <- cowplot::get_legend(
    make_district_diff_map_panel(
      dat$a_diff_aid_sf, "diff_aid", policy_levels[[1]],
      fill_title = "Aid",
      fill_limit = aid_map_limit,
      iraq_district = iraq_district,
      highlight_edges = NULL,
      show_legend = TRUE
    )
  )

  outcome_map_plots <- lapply(policy_levels, function(level) {
    selected_out_level <- selected_out %>% filter(as.character(policy_length) == level)
    make_district_diff_map_panel(
      dat$a_diff_out_sf, "diff_outcome", level,
      fill_title = "Outcome",
      fill_limit = outcome_map_limit,
      iraq_district = iraq_district,
      highlight_edges = make_highlight_edges_clean(iraq_district, selected_out_level),
      show_legend = FALSE
    )
  })
  outcome_map_legend <- cowplot::get_legend(
    make_district_diff_map_panel(
      dat$a_diff_out_sf, "diff_outcome", policy_levels[[1]],
      fill_title = "Outcome",
      fill_limit = outcome_map_limit,
      iraq_district = iraq_district,
      highlight_edges = NULL,
      show_legend = TRUE
    )
  )

  outcome_bar_plots <- lapply(policy_levels, function(level) {
    make_governorate_change_panel(
      gov_out_diff_ci, level,
      x_limits = c(-outcome_bar_limit, outcome_bar_limit),
      include_ci = TRUE,
      show_legend = FALSE,
      N = N
    )
  })
  outcome_bar_legend <- cowplot::get_legend(
    make_governorate_change_panel(
      gov_out_diff_ci, policy_levels[[1]],
      x_limits = c(-outcome_bar_limit, outcome_bar_limit),
      include_ci = TRUE,
      show_legend = TRUE,
      N = N
    )
  )

  make_panel_row <- function(plot_list, legend) {
    cowplot::plot_grid(
      plotlist = c(plot_list, list(legend)),
      nrow = 1,
      rel_widths = panel_widths,
      align = "none"
    )
  }

  row_gap_height <- 0.055
  row_gap <- make_application_blank_plot()

  make_group <- function(label, first_row, second_row, rel_heights) {
    cowplot::plot_grid(
      make_application_group_label(label),
      cowplot::plot_grid(
        first_row,
        row_gap,
        second_row,
        ncol = 1,
        rel_heights = c(rel_heights[1], row_gap_height, rel_heights[2])
      ),
      nrow = 1,
      rel_widths = c(group_label_width, 1),
      align = "none"
    )
  }

  aid_group <- make_group(
    "Aid",
    make_panel_row(aid_bar_plots, aid_bar_legend),
    make_panel_row(aid_map_plots, aid_map_legend),
    rel_heights = c(0.42, 0.58)
  )

  if (layout == "maps_together") {
    outcome_group <- make_group(
      "Outcome",
      make_panel_row(outcome_map_plots, outcome_map_legend),
      make_panel_row(outcome_bar_plots, outcome_bar_legend),
      rel_heights = c(0.56, 0.44)
    )
  } else {
    outcome_group <- make_group(
      "Outcome",
      make_panel_row(outcome_bar_plots, outcome_bar_legend),
      make_panel_row(outcome_map_plots, outcome_map_legend),
      rel_heights = c(0.44, 0.56)
    )
  }

  cowplot::plot_grid(
    header_row,
    aid_group,
    row_gap,
    outcome_group,
    ncol = 1,
    rel_heights = c(0.055, 0.46, 0.025, 0.46),
    align = "none"
  )
}

save_application_plots <- function(data, iraq_district, out_SAF = NULL, out_IED = NULL, policy_name, output_dir) {
  out_dir <- file.path(output_dir, policy_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plot_one_outcome <- function(out, Y_mat_name, event_label) {
    dat <- prepare_all_data(iraq_district = iraq_district, data = data, out = out, Y_mat_name = Y_mat_name)
    levels <- make_level_maps(dat$a_obs_aid_sf, dat$a_policy_aid_sf, dat$a_obs_out_sf, dat$a_policy_out_sf, iraq_district, event_type = paste0(event_label, ": "))
    gov_aid_diff <- aggregate_diff_to_governorate(diff_tbl = dat$aid_diff_long, iraq_district = iraq_district, value_col = diff_aid, weight_col = NULL, agg = "sum")
    p_gov_aid <- top_change_bars_governorate(gov_aid_diff, N = 3, title = "")
    gov_out_diff_ci <- if (!is.null(out$tmp_gov) && length(out$tmp_gov) > 0) {
      build_gov_out_diff_ci_from_tmpgov(iraq_district = iraq_district, data = data, out = out, Y_mat_name = Y_mat_name, weight_col = NULL, agg = "sum")
    } else {
      aggregate_diff_to_governorate_ci(diff_tbl = dat$outcome_diff_long, iraq_district = iraq_district, value_col = diff_outcome, lb_col = diff_outcome_lb, ub_col = diff_outcome_ub, weight_col = NULL, agg = "sum")
    }
    p_gov_out_ci <- top_change_bars_governorate_ci(gov_out_diff_ci, N = 3, title = "")
    selected_aid <- top_bottom_governorates(gov_aid_diff, N = 3)
    selected_out <- top_bottom_governorates(gov_out_diff_ci, N = 3)
    maps <- make_diff_maps(
      a_diff_aid_sf = dat$a_diff_aid_sf,
      a_diff_out_sf = dat$a_diff_out_sf,
      iraq_district = iraq_district,
      event_type = paste0(event_label, ": "),
      highlight_aid_edges = make_highlight_edges_clean(iraq_district, selected_aid),
      highlight_out_edges = make_highlight_edges_clean(iraq_district, selected_out)
    )
    ggplot2::ggsave(file.path(out_dir, paste0("p_gov_aid_", event_label, "_", policy_name, ".pdf")), p_gov_aid, width = 10, height = 2.8, units = "in", dpi = 600)
    ggplot2::ggsave(file.path(out_dir, paste0("p_gov_out_ci_", event_label, "_", policy_name, ".pdf")), p_gov_out_ci, width = 10, height = 2.8, units = "in", dpi = 600)
    ggplot2::ggsave(file.path(out_dir, paste0("combined_diff_", event_label, "_", policy_name, ".pdf")), maps$combined_diff, width = 10, height = 5, units = "in", dpi = 600)
    composite <- make_application_composite_plot(
      data = data,
      iraq_district = iraq_district,
      out = out,
      Y_mat_name = Y_mat_name,
      event_label = event_label,
      N = 3,
      policy_name = policy_name
    )
    ggplot2::ggsave(file.path(out_dir, paste0("plot_diff_", tolower(event_label), "_", policy_name, "_grid.pdf")), composite, width = 8.6, height = 7.6, units = "in", dpi = 600)
    composite_grid_2 <- make_application_composite_plot(
      data = data,
      iraq_district = iraq_district,
      out = out,
      Y_mat_name = Y_mat_name,
      event_label = event_label,
      N = 3,
      policy_name = policy_name,
      layout = "paired"
    )
    ggplot2::ggsave(file.path(out_dir, paste0("plot_diff_", tolower(event_label), "_", policy_name, "_grid_2.pdf")), composite_grid_2, width = 8.6, height = 7.6, units = "in", dpi = 600)
    if (event_label == "SAF" && startsWith(policy_name, "wopop")) {
      ggplot2::ggsave(file.path(out_dir, "p_aid_obs_SAF.pdf"), levels$p_aid_obs, width = 4.7, height = 3, units = "in", dpi = 600)
      ggplot2::ggsave(file.path(out_dir, "p_out_obs_SAF.pdf"), levels$p_out_obs, width = 4.7, height = 3, units = "in", dpi = 600)
    }
  }
  if (!is.null(out_SAF)) {
    plot_one_outcome(out_SAF, "SAF_week", "SAF")
  }
  if (!is.null(out_IED)) {
    plot_one_outcome(out_IED, "IED_week", "IED")
  }
  invisible(out_dir)
}

print_application_tables <- function(data, out_SAF = NULL, out_IED = NULL, case_name) {
  if (!is.null(out_SAF)) {
    tbl_policy_SAF <- policy_coef_table(out_SAF, "SAF")
    tbl_val_SAF <- value_table_from_out(out_SAF, data$SAF_week, "SAF")
  }
  if (!is.null(out_IED)) {
    tbl_policy_IED <- policy_coef_table(out_IED, "IED")
    tbl_val_IED <- value_table_from_out(out_IED, data$IED_week, "IED")
  }
  if (!is.null(out_SAF) && !is.null(out_IED)) {
    wide_values <- values_twoheader_wide(tbl_val_IED, tbl_val_SAF, M_select = seq_len(length(out_SAF$res_policy)), drop_M_cols = FALSE)
  }
  if (requireNamespace("kableExtra", quietly = TRUE)) {
    if (!is.null(out_SAF)) {
      coef_table_render_dual(tbl_policy_SAF, caption = paste0("Learned policy coefficients for SAF (", case_name, ")"), M_select = length(out_SAF$res_policy), drop_M_col = TRUE, print_latex = TRUE)
    }
    if (!is.null(out_IED)) {
      coef_table_render_dual(tbl_policy_IED, caption = paste0("Learned policy coefficients for IED (", case_name, ")"), M_select = length(out_IED$res_policy), drop_M_col = TRUE, print_latex = TRUE)
    }
    if (!is.null(out_SAF) && !is.null(out_IED)) {
      values_twoheader_render_dual(wide_values, caption = paste0("Estimated stochastic policy values (", case_name, ")"), print_latex = TRUE)
    }
  } else {
    if (!is.null(out_SAF)) {
      print(tbl_policy_SAF)
      print(tbl_val_SAF)
    }
    if (!is.null(out_IED)) {
      print(tbl_policy_IED)
      print(tbl_val_IED)
    }
    if (!is.null(out_SAF) && !is.null(out_IED)) {
      print(wide_values)
    }
  }
  invisible(list(
    SAF_policy = if (exists("tbl_policy_SAF")) tbl_policy_SAF else NULL,
    SAF_value = if (exists("tbl_val_SAF")) tbl_val_SAF else NULL,
    IED_policy = if (exists("tbl_policy_IED")) tbl_policy_IED else NULL,
    IED_value = if (exists("tbl_val_IED")) tbl_val_IED else NULL,
    combined_values = if (exists("wide_values")) wide_values else NULL
  ))
}

values_twoheader_wide <- function(tbl_val_IED,
                                  tbl_val_SAF,
                                  M_select = 4,
                                  drop_M_cols = TRUE) {
  clean_vals <- function(df) {
    df %>%
      transmute(
        M,
        Vhat = sprintf("%.2f", Vhat),
        Lower = sprintf("%.2f", Lower),
        Upper = sprintf("%.2f", Upper),
        Observed = sprintf("%.3f", Observed)
      )
  }

  ied <- clean_vals(tbl_val_IED) %>% filter(M %in% M_select)
  saf <- clean_vals(tbl_val_SAF) %>% filter(M %in% M_select)
  stopifnot(nrow(ied) == nrow(saf), all(ied$M == saf$M))

  wide <- tibble(
    M = ied$M,
    IED_Vhat = ied$Vhat,
    IED_Lower = ied$Lower,
    IED_Upper = ied$Upper,
    IED_Observed = ied$Observed,
    SAF_Vhat = saf$Vhat,
    SAF_Lower = saf$Lower,
    SAF_Upper = saf$Upper,
    SAF_Observed = saf$Observed
  )

  if (drop_M_cols) {
    wide <- select(wide, -M)
    colnames(wide) <- c(
      "$\\hat V$", "Lower", "Upper", "Observed",
      "$\\hat V$", "Lower", "Upper", "Observed"
    )
  } else {
    colnames(wide) <- c(
      "$M$",
      "$\\hat V$", "Lower", "Upper", "Observed",
      "$\\hat V$", "Lower", "Upper", "Observed"
    )
  }

  wide
}

values_twoheader_render_dual <- function(wide_tbl, caption, print_latex = TRUE) {
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("Package `kableExtra` is required to render application tables.")
  }

  has_m_col <- ncol(wide_tbl) == 9
  header <- if (has_m_col) {
    c(" " = 1, "IED" = 4, "SAF" = 4)
  } else {
    c("IED" = 4, "SAF" = 4)
  }
  align <- if (has_m_col) c("c", rep("r", 8)) else rep("r", 8)

  kb_html <- kableExtra::kbl(
    wide_tbl, format = "html", booktabs = TRUE, escape = FALSE,
    align = align,
    caption = caption
  ) |>
    kableExtra::kable_styling(full_width = FALSE, font_size = 11) |>
    kableExtra::add_header_above(header, escape = FALSE)

  if (print_latex) {
    kb_tex <- kableExtra::kbl(
      wide_tbl, format = "latex", booktabs = TRUE, escape = FALSE,
      align = align,
      caption = caption
    ) |>
      kableExtra::kable_styling(
        full_width = FALSE, font_size = 11,
        latex_options = "hold_position"
      ) |>
      kableExtra::add_header_above(header, escape = FALSE)
    cat(as.character(kb_tex), "\n")
  }

  invisible(kb_html)
}
