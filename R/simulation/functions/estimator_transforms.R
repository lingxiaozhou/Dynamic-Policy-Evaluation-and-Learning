# Estimator pseudo-outcome transformations used by the simulation code.
# Extracted from the research code and kept separate from the data-generating process.

pair_weight_sum <- function(weight){
  z <- weight - 1
  missing_n <- rowSums(is.na(z))
  z[is.na(z)] <- 0
  z_sum <- rowSums(z)
  z_pair_sum <- 0.5 * (z_sum^2 - rowSums(z^2))
  z_sum + z_pair_sum + missing_n * z_sum + missing_n * (missing_n - 1) / 2 + 1
}

first_weight_sum <- function(weight){
  rowSums(weight - 1, na.rm = TRUE) + 1
}

f_smooth <- function(x, eps = 0.01) {
  y <- numeric(length(x))
  y[x <= 0] <- 0

  idx1 <- x > 0 & x < eps
  y[idx1] <- (x[idx1]^2) / (2 * eps)

  idx2 <- x >= eps & x <= 1 - eps
  y[idx2] <- x[idx2]

  idx3 <- x > 1 - eps & x < 1
  y[idx3] <- 1 - (1 - x[idx3])^2 / (2 * eps)

  y[x >= 1] <- 1
  y
}

get_Y_est <- function(Y, weight, M, method, truncation = 0, stabilize2 = FALSE){
  
  n <- ncol(Y)
  time_points <- nrow(Y)
  Y_est <- array(NA,dim(Y))
  weight_t <- NULL
  
  if(method == "addIPW"){
    
    weight_t <- first_weight_sum(weight)
    
    if(stabilize2){
      weight_t <- weight_t/mean(weight_t)
    }
    
    if(truncation>0){
      ub <- quantile(weight_t,1-truncation,na.rm = TRUE)
      lb <- quantile(weight_t,truncation,na.rm = TRUE)
      weight_t <- ifelse(weight_t > lb, weight_t, lb)
      weight_t <- ifelse(weight_t < ub, weight_t, ub)
    }
    
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * weight_t[idx]
  }
  
  if(method == "addIPW2"){
    weight_t <- pair_weight_sum(weight)
    if(stabilize2){
      weight_t <- weight_t/mean(weight_t)
    }
    
    if(truncation>0){
      ub <- quantile(weight_t,1-truncation)
      lb <- quantile(weight_t,truncation)
      weight_t <- ifelse(weight_t > lb, weight_t, lb)
      weight_t <- ifelse(weight_t < ub, weight_t, ub)
    }
    
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * weight_t[idx]
  }
  
  
  if(method == "IPWnoinf"){
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- weight[idx,] * Y[idx + M - 1,]
  }
  
  
  if(method == "IPW"){
    
    weight_t <- sapply(1:(time_points), function(t) prod(weight[t,],na.rm = TRUE))
    
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * weight_t[idx]
  }
  
  if(M>1){
    Y_est[(time_points-M+2):time_points,] <- NA
  }
  
  return(Y_est)
  
}




# For stabilization2
get_Y_est2 <- function(Y, weight1, weight2=NULL, M, method, shift = NULL, stabilize2 = FALSE,bound_weight = FALSE,save_weight=FALSE,stabilize3=FALSE, save_Y_est=FALSE,mix_par1=1.5,mix_par2=0.1){

  if(!is.null(shift)){
    Y <- Y+shift
  }
  n <- ncol(Y)
  time_points <- nrow(Y)
  Y_est <- array(NA,dim(Y))
  weight_t <- NULL
  mean_weight <- NULL
  
  if(method == "addIPW"){
    if(M==2){
      weight_t1 <- first_weight_sum(weight1)
      weight_t2 <- first_weight_sum(weight2)
      weight_t <- weight_t2[-time_points]*weight_t1[-1]
    }else{
      weight_t1 <- first_weight_sum(weight1)
      weight_t <- weight_t1
    }
    mean_weight <- mean(weight_t)
    
    
    if(stabilize2){

      if(!bound_weight){
        weight_t <- weight_t/mean(weight_t)
      }else{
        weight_t <- weight_t/(mean(weight_t1)*mean(weight_t2))
      }
    }
    
    if(stabilize3){
      
      a <- 1/(1+exp(-(mean(weight_t)^2-mix_par1)/mix_par2))
      # a <- 1/(1+exp(-(mean(weight_t)^2-1)/0.1))

      
      if(!bound_weight){
        weight_t <- a*weight_t/mean(weight_t)+(1-a)*(weight_t+1-mean(weight_t))
      }else{
        weight_t <- weight_t/(mean(weight_t1)*mean(weight_t2))
      }
      
      
      if(save_Y_est){
        return(list(a,weight_t))
      }
    }
    
    # if(stabilize2){
    #   
    #   if(!bound_weight){
    #     weight_t <- weight_t+1-mean(weight_t)
    #   }else{
    #     weight_t <- weight_t/(mean(weight_t1)*mean(weight_t2))
    #   }
    # }
    
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * weight_t[idx]
  }
  
  
  if(method == "addIPW2"){
    
    if(M==2){

      weight_t1 <- pair_weight_sum(weight1)
      weight_t2 <- pair_weight_sum(weight2)
      weight_t <- weight_t2[-time_points]*weight_t1[-1]
    }else{
      weight_t1 <- pair_weight_sum(weight1)
      weight_t <- weight_t1
    }
    
    mean_weight <- mean(weight_t)
    if(stabilize2){
      # if(!bound_weight){
      #   weight_t <- weight_t+1-mean(weight_t)
      # }else{
      #   weight_t <- weight_t/(mean(weight_t1)*mean(weight_t2))
      # }

      if(!bound_weight){
        weight_t <- weight_t/mean(weight_t)
      }else{
        weight_t <- weight_t/(mean(weight_t1)*mean(weight_t2))
      }
      
    }
    
    
    if(stabilize3){
      a <- 1/(1+exp(-(mean(weight_t)^2-mix_par1)/mix_par2))
      # a <- 1/(1+exp(-(mean(weight_t)^2-1)/0.1))
      raw_weight_t <- weight_t
      if(!bound_weight){
        weight_t <- a*weight_t/mean(weight_t)+(1-a)*(weight_t+1-mean(weight_t))
      }else{
        weight_t <- weight_t/(mean(weight_t1)*mean(weight_t2))
      }
      
      
      if(save_Y_est){
        return(list(a,weight_t,raw_weight_t))
      }

    }
    
    
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * weight_t[idx]
  }
  

  
  if(method == "IPW"){
    
   
    
    if(M==2){
      weight_t1 <- sapply(1:(time_points), function(t) prod(weight1[t,],na.rm = TRUE))
      weight_t2 <- sapply(1:(time_points), function(t) prod(weight2[t,],na.rm = TRUE))
      weight_t <- weight_t2[-time_points]*weight_t1[-1]
    }else{
      weight_t1 <- sapply(1:(time_points), function(t) prod(weight1[t,],na.rm = TRUE))
      weight_t <- weight_t1
    }
    
    mean_weight <- mean(weight_t)
    if(stabilize2|stabilize3){
      
      if(!bound_weight){
        weight_t <- weight_t/mean(weight_t)
      }else{
        weight_t <- weight_t/(mean(weight_t1)*mean(weight_t2))
      }
    }
    
    
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * weight_t[idx]
  }
  
  if(method == "adaptive"){
    
    
    if(M==2){
      weight_t1 <- first_weight_sum(weight1)
      weight_t2 <- first_weight_sum(weight2)
      weight_t <- weight_t2[-time_points]*weight_t1[-1]
    }else{
      weight_t1 <- first_weight_sum(weight1)
      weight_t <- weight_t1
    }
    
    A <- sapply(1:(time_points-M+1), function(t) weight_t[t]^2*mean(Y[t+M-1,]))
    B <- sapply(1:(time_points-M+1), function(t) weight_t[t]*mean(Y[t+M-1,]))
    C <- mean(weight_t)
    D <- sapply(1:(time_points-M+1), function(t) weight_t[t]^2)
    
    # lambda <- (mean(A)-mean(B)*C)/(mean(D)-C^2)
    lambda <- (mean(A)-mean(B)*1)/(mean(D)-1)
    
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * weight_t[idx] + lambda * (1 - C)
    
    if(save_Y_est){
      return(list(A=A,B=B,C=C,D=D,lambda=lambda,weight_t=weight_t))
    }
  }
  
  if(method == "adaptive2"){
    
    if(M==2){
      
      weight_t1 <- pair_weight_sum(weight1)
      weight_t2 <- pair_weight_sum(weight2)
      weight_t <- weight_t2[-time_points]*weight_t1[-1]
    }else{
      weight_t1 <- pair_weight_sum(weight1)
      weight_t <- weight_t1
    }
    
    A <- sapply(1:(time_points-M+1), function(t) weight_t[t]^2*mean(Y[t+M-1,]))
    B <- sapply(1:(time_points-M+1), function(t) weight_t[t]*mean(Y[t+M-1,]))
    C <- mean(weight_t)
    D <- sapply(1:(time_points-M+1), function(t) weight_t[t]^2)
    
    lambda <- (mean(A)-mean(B)*C)/(mean(D)-C^2)
    
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * weight_t[idx] + lambda * (1 - C)
    
    if(save_Y_est){
      return(list(A=A,B=B,C=C,D=D,lambda=lambda,weight_t=weight_t))
    }
  }
  
  
  if(method == "shift_adaptive"){
    
    
    if(M==2){
      weight_t1 <- first_weight_sum(weight1)
      weight_t2 <- first_weight_sum(weight2)
      weight_t <- weight_t2[-time_points]*weight_t1[-1]
    }else{
      weight_t1 <- first_weight_sum(weight1)
      weight_t <- weight_t1
    }
    
    mu_hat <- mean(sapply(1:(time_points-M+1), function(t) weight_t[t]*mean(Y[t+M-1,])/mean(weight_t)))
    Y_bar <- mean(Y)
    
    A <- sapply(1:(time_points-M+1), function(t) weight_t[t]*(mean(Y[t+M-1,])-mu_hat))
    B <- sapply(1:(time_points-M+1), function(t) weight_t[t]*(mu_hat-Y_bar))
    C <- mean(weight_t)

    
    # lambda <- (mean(A)-mean(B)*C)/(mean(D)-C^2)
    lambda <- max(0,-mean(A*(B-(mu_hat-Y_bar)))/mean((B-(mu_hat-Y_bar))^2))
    lambda <- min(1,lambda)
    shift_weight_t <- sapply(1:(time_points-M+1), function(t) weight_t[t]+lambda*(1-C))
    shift_weight_t <- shift_weight_t/mean(shift_weight_t)
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * shift_weight_t[idx]
    
    mu_hat <- mean(Y_est,na.rm = TRUE)
    
    if(save_Y_est){
      return(list(lambda= lambda, shift_weight_t=shift_weight_t,mu_hat=mu_hat))
    }
  }
  
  if(method == "shift_adaptive2"){
    
    
    if(M==2){
      
      weight_t1 <- pair_weight_sum(weight1)
      weight_t2 <- pair_weight_sum(weight2)
      weight_t <- weight_t2[-time_points]*weight_t1[-1]
    }else{
      weight_t1 <- pair_weight_sum(weight1)
      weight_t <- weight_t1
    }
    
    mu_hat <- mean(sapply(1:(time_points-M+1), function(t) weight_t[t]*mean(Y[t+M-1,])/mean(weight_t)))
    
    Y_bar <- mean(Y)
    
    A <- sapply(1:(time_points-M+1), function(t) weight_t[t]*(mean(Y[t+M-1,])-mu_hat))
    B <- sapply(1:(time_points-M+1), function(t) weight_t[t]*(mu_hat-Y_bar))
    C <- mean(weight_t)
    
    
    # lambda <- (mean(A)-mean(B)*C)/(mean(D)-C^2)
    lambda <- max(0,-mean(A*(B-(mu_hat-Y_bar)))/mean((B-(mu_hat-Y_bar))^2))
    lambda <- min(1,lambda)
    
    
    shift_weight_t <- sapply(1:(time_points-M+1), function(t) weight_t[t]+lambda*(1-C))
    shift_weight_t <- shift_weight_t/mean(shift_weight_t)
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * shift_weight_t[idx]
    
    mu_hat <- mean(Y_est,na.rm = TRUE)
    
    if(save_Y_est){
      return(list( lambda = lambda, shift_weight_t=shift_weight_t, mu_hat = mu_hat, mean_weight = C, A=A, B=B,weight_t=weight_t))
    }
  }
  
  
  
  
  if(method %in% c("mix_adaptive","mix_adaptive2")){
    
    if(method=="mix_adaptive"){
      if(M==2){
        weight_t1 <- first_weight_sum(weight1)
        weight_t2 <- first_weight_sum(weight2)
        weight_t <- weight_t2[-time_points]*weight_t1[-1]
      }else{
        weight_t1 <- first_weight_sum(weight1)
        weight_t <- weight_t1
      }
    }else{
      if(M==2){
        
        weight_t1 <- pair_weight_sum(weight1)
        weight_t2 <- pair_weight_sum(weight2)
        weight_t <- weight_t2[-time_points]*weight_t1[-1]
      }else{
        weight_t1 <- pair_weight_sum(weight1)
        weight_t <- weight_t1
      }
    }
    
    
    
    mu_hat <- mean(sapply(1:(time_points-M+1), function(t) weight_t[t]*mean(Y[t+M-1,])/mean(weight_t)))
    Y_bar <- mean(Y)
    
    A <- sapply(1:(time_points-M+1), function(t) weight_t[t]*(mean(Y[t+M-1,])-Y_bar))
    B <- sapply(1:(time_points-M+1), function(t) weight_t[t]*(mu_hat-Y_bar))
    C <- mean(weight_t)
    
    
    lambda <- mean(A*(B-(mu_hat-Y_bar)))/mean((B-(mu_hat-Y_bar))^2)
    # lambda <- min(max(0,lambda),1)
    lambda <- f_smooth(lambda)
    shift_weight_t <- sapply(1:(time_points-M+1), function(t) lambda*weight_t[t]/mean(weight_t)+(1-lambda)*(weight_t[t]+1-C))
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * shift_weight_t[idx]
    
    mu_hat <- mean(Y_est,na.rm = TRUE)
    
    if(save_Y_est){
      return(list(lambda= lambda, shift_weight_t=shift_weight_t,mu_hat=mu_hat,C=C,weight_t=weight_t))
    }
  }
  
  
  if(method %in% c("shift","shift2")){
    
    if(method=="shift"){
      if(M==2){
        weight_t1 <- first_weight_sum(weight1)
        weight_t2 <- first_weight_sum(weight2)
        weight_t <- weight_t2[-time_points]*weight_t1[-1]
      }else{
        weight_t1 <- first_weight_sum(weight1)
        weight_t <- weight_t1
      }
    }else{
      if(M==2){
        
        weight_t1 <- pair_weight_sum(weight1)
        weight_t2 <- pair_weight_sum(weight2)
        weight_t <- weight_t2[-time_points]*weight_t1[-1]
      }else{
        weight_t1 <- pair_weight_sum(weight1)
        weight_t <- weight_t1
      }
    }

    mean_weight <- mean(weight_t)
    
    

    weight_t <- weight_t+(1-mean_weight)
    
    idx <- 1:(time_points-M+1)
    Y_est[idx,] <- Y[idx + M - 1,] * weight_t[idx]
    
    
    if(save_Y_est){
      bound <- Y_est-(mean_weight)*mean(data$Y)
      return(bound)
    }
  }
  
  
  if(M>1){
    Y_est[(time_points-M+2):time_points,] <- NA
  }
  if(save_weight){

    return(mean_weight)
  }
  return(Y_est)
  
}
