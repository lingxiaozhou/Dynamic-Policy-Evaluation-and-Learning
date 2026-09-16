# If stochastic is TRUE, use stochastic policy

get_true_policy_value <- function(data,model,ga,stochastic = FALSE,MC = FALSE,nrep=100,save_array = FALSE,shift = NULL){
  
  W <- data$W
  if (is.vector(ga)) {
    ga <- matrix(ga, nrow = 1)
  }
  time_points <- nrow(W)
  n <- ncol(W)
  M <- nrow(ga)
  X_all_list <- get_X_all_cached(data, model, time_points)
  

  
  Q <- rep(NA,time_points)
  
  if(!stochastic){
    
    for (t in 1:(time_points-M+1)) {
      X_all <- X_all_list[[t]]
      W_counter <- c(X_all%*%ga[M,]>=0)
      Q[t] <- mean(compute_true_Q(data,t,model,M,W_counter,ga = ga))
    }

  }else{
    
    
    if(model==1){
      if(MC==FALSE){
        
        
        W_counter_all <- expand.grid(rep(list(c(0, 1)), n))
        
        for (t in 1:(time_points-M+1)) {
          
          Qt <- rep(NA,nrow(W_counter_all))
          
          X_all <- X_all_list[[t]]
          W_log_odds <- X_all%*%ga[M,]
          prob <- (1+exp(-W_log_odds))^(-1)
          
          
          for (index in 1:nrow(W_counter_all)){
            
            Qt[index] <- mean(compute_true_Q(data,t,model,M,unlist(W_counter_all[index,]),ga))*
              prod(prob^W_counter_all[index,])*prod((1-prob)^(1-W_counter_all[index,]))
            
          }
          
          
          Q[t] <- sum(Qt)
          
        }
        
      }else{
        
        for (t in 1:(time_points-M+1)) {
          
          Qt <- rep(NA,nrep)
          
          X_all <- X_all_list[[t]]
          W_log_odds <- X_all%*%ga[M,]
          prob <- (1+exp(-W_log_odds))^(-1)
          W_counter_all <- t(sapply(1:nrep, function(x) rbinom(n, size = 1, prob = prob)))
          
          for (index in 1:nrep){
            
            Qt[index] <- mean(compute_true_Q(data,t,model,M,unlist(W_counter_all[index,]),ga))
            
          }
          
          
          Q[t] <- mean(Qt)
          
        }
        
        
      }
    }
    
    if(model!=1){
      for (t in 1:(time_points-M+1)) {
        X_all <- X_all_list[[t]]
        W_counter <- (1+exp(-X_all%*%ga[M,]))^(-1)
        Q[t] <- mean(compute_true_Q(data,t,model,M,W_counter,ga = ga,stochastic = stochastic))
      }
    }

    
    

  }
  
  if(!is.null(shift)){
    Q <- Q+shift
  }
  
  if(save_array){
    return(Q)
  }
  
  return(mean(Q,na.rm = TRUE))
}
