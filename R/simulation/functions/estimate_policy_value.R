# Estimate the value of a fixed policy using IPW/addIPW variants.
# Public models 1, 2, and 3 all use the four-parameter policy
#   gamma0 + gamma1 W_{t-1} + gamma2 X_{t-1} + gamma3 Z_{t-1}.




get_estimated_policy_value <- function(data,p,ga, method = "addIPW",
                                       model,save_Y_est=FALSE,
                                       stochastic = FALSE, 
                                       stabilize = FALSE, get_var = FALSE, stabilize2 = FALSE,
                                       shift = NULL,bound_weight = FALSE,save_weight=FALSE,stabilize3=FALSE,
                                       mix_par1=1.5,mix_par2=0.1){
  W <- data$W

  if(length(p)==1){
    p <- array(p,dim(W))
  }
  
  if(model %in% c(1, 2, 3)){
    ga <- matrix(ga,ncol=4)
  }else{
    stop("Unsupported model. Use model = 1, 2, or 3.", call. = FALSE)
  }
  

  time_points <- nrow(data$W)
  n <- ncol(data$W)
  if(length(method)==1 & nrow(ga)>1 ){
    method <- rep(method,nrow(ga))
  }
  X_all_matrix <- get_X_all_matrix(data, model, time_points)
  
  # get propensity score
  e <- array(NA,c(time_points,n))
  e[data$W==1] <- p[data$W==1]
  e[data$W==0] <- 1-p[data$W==0]
  
  Y_est <- data$Y
  
  if(!is.null(shift)){
    Y_est <- data$Y+shift
  }
  

  
  for (M in 1:nrow(ga)) {
    
    # get the treatment assignment
    W_counter <- array(NA,c(time_points,n))
    idx <- 1:(time_points-M+1)
    row_idx <- 1:(length(idx) * n)
    W_counter[idx,] <- matrix(X_all_matrix[row_idx, , drop = FALSE] %*% ga[M,],
                              nrow = length(idx), ncol = n, byrow = TRUE)
    if(stochastic==FALSE){
      W_counter[idx,] <- W_counter[idx,] >= 0
    }
    

    # IPW weights for each unit
    if(stochastic==FALSE){
      weight <- as.numeric(data$W==W_counter)/e
    }else{
      prob <- plogis(W_counter)
      weight <- array(NA_real_, c(time_points,n))
      weight[data$W == 1] <- prob[data$W == 1] / e[data$W == 1]
      weight[data$W == 0] <- (1 - prob[data$W == 0]) / e[data$W == 0]
    }
    
    if(stabilize){
      weight <- weight/mean(weight,na.rm = TRUE)
    }
    if(stabilize2|stabilize3| ! method[M]%in%c("addIPW","addIPW2","IPW")){
      assign(paste0("weight",M),weight)
    }else{
      Y_est <- get_Y_est(Y_est, weight, M, method[M])
    }
    
    
  }
  
  if(stabilize2|stabilize3|! method[M]%in%c("addIPW","addIPW2","IPW")){
    Y_est <- get_Y_est2(data$Y, weight1, weight2, M, method[M],shift = shift,stabilize2 = stabilize2,bound_weight = bound_weight,save_weight=save_weight,stabilize3 = stabilize3,save_Y_est=save_Y_est,mix_par1 = mix_par1,mix_par2 = mix_par2)
  }
  
  if(save_Y_est){
    return(Y_est)
  }
  
  V_est <- mean(Y_est,na.rm = TRUE)

  if(get_var){
    
    V_var <- mean(apply(Y_est,1,mean)^2,na.rm = TRUE)/(time_points-nrow(ga)+1)
    CI_lb <- V_est - 1.96*sqrt(V_var)
    CI_ub <- V_est + 1.96*sqrt(V_var)
    
  }
  
  
  if(save_weight){
    return(Y_est)
  }
  


  
  if(save_Y_est==FALSE & get_var==FALSE){
    return(V_est)
  }
  
  if(save_Y_est==FALSE & get_var==TRUE){
    return(list(V_est = V_est, V_var = V_var,CI_lb = CI_lb, CI_ub = CI_ub, sta_coeff = mean(weight)))
  }
  
  
  
  
}
  
  
    
    
