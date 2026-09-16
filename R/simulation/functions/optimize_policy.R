# Optimize a stochastic policy under an estimated value criterion.
# Public models 1, 2, and 3 all use four policy parameters.


get_optimal_policy <- function(data,p, method = "addIPW",stochastic = FALSE,ub=30,lb=-30,model,
                               ga_prev = NULL,stabilize = FALSE,stabilize2 = FALSE,shift=NULL,stabilize3 = FALSE,mix_par1=1.5,mix_par2=0.1){
  W <- data$W
  
  if(length(p)==1){
    p <- array(p,dim(W))
  }
  
  time_points <- nrow(W)
  n <- ncol(W)
  X_all_matrix <- get_X_all_matrix(data, model, time_points)
  if(is.null(ga_prev)){
    M <- 1
  }else{
    M <- nrow(ga_prev)+1
  }
  
  finite_objective <- function(value){
    if(length(value) != 1 || !is.finite(value)){
      return(1e100)
    }
    -value
  }
  
  if(method%in%c("addIPW", "IPWnoinf") & stochastic==FALSE){
    weight <- array(NA,c(time_points,n))
    weight[W==1] <- 1/p[W==1]
    weight[W==0] <- -1/(1-p[W==0])
    
    if(M==1){
      Y_est <- data$Y
    }else{
      Y_est <- get_estimated_policy_value(data,p,ga_prev, method = method,model=model,save_Y_est=TRUE,stabilize2 = stabilize2,shift=shift,stabilize3 = stabilize3,mix_par1 = mix_par1,mix_par2 = mix_par2)
    }

    
    get_pseudo_out <- function(weight, Y, t, method="addIPW",M=NULL){
      
      if(method == "addIPW"){
        
        return(matrix(weight[t,]*mean(Y[t+M-1,]),ncol=1))
        
      }else{
        
        return(matrix(weight[t,]*Y[t,],ncol=1))
        
      }
    }
        
    
    
    
    # optimization
    X_all <- X_all_matrix[1:((time_points - M + 1) * n), , drop = FALSE]
    
    nn_all <- nrow(X_all)
    d <- ncol(X_all)
    f <- c(rep(0,d), 
           do.call(rbind,lapply(1:(time_points-M+1), function(t) get_pseudo_out(weight, Y_est, t, method,M))))  # objective function coefficients
    
    B <- 100  # bounds on coefficients
    C <- B * apply(abs(X_all), 1, sum)  # maximum values of x'beta
    # minmargin <- max(1, C) * (1e-8)  # prevent non-integer numbers the integrality constraint of integers from being counted as integers
    # minmargin <-  max(1, C) * (1e-2)/2   # handle the strict inequality
    minmargin <- 0.1
    
    
    # if(model==2){
    #   B <- 10  # bounds on coefficients
    #   C <- B * apply(abs(X_all), 1, sum)  # maximum values of x'beta
    #   # minmargin <- max(1, C) * (1e-8)  # prevent non-integer numbers the integrality constraint of integers from being counted as integers
    #   # minmargin <-  max(1, C) * (1e-2)/2   # handle the strict inequality
    #   minmargin <- max(1, C) * (1e-4)
    # }

    
    Aineq <- rbind(cbind(-X_all, diag(C)), cbind(X_all, -diag(C)) )
    dir <-  c(rep("<=", 0.5*nrow(Aineq)),rep("<=", 0.5*nrow(Aineq)))
    bineq <- c(C, 0*C-minmargin)

    lb <- c(-B * rep(1,d), rep(0, nn_all))
    ub <- c(B * rep(1, d), rep(1, nn_all))

    # Variable type string
    ctype <- c(rep("C", d), rep("B", nn_all))


    bounds <- list(lower = list(ind = 1:length(lb),val = lb),
                 upper = list(ind =  1:length(lb),val = ub))

    
    sol_pd <- Rglpk_solve_LP(obj = f, mat = Aineq, dir = dir, rhs = bineq, bounds = bounds, types = ctype,max = TRUE)
    ga <- sol_pd$solution[c(1:d)]
    
    return(rbind(ga_prev,ga))

  }
  
  
  if(method=="IPW"){
    
    if(model %in% c(1, 2, 3)){
      # Initial guesses
      par <- c(0, 0, 0, 0)
    }else{
      stop("Unsupported model. Use model = 1, 2, or 3.", call. = FALSE)
    }
    
    d <- length(par)

    # Objective function to maximize (negate for optimization)
    objective_function <- function(par) {
      value <- get_estimated_policy_value(data, p, rbind(ga_prev,par),method,model=model,stochastic=TRUE,stabilize = stabilize,stabilize2 = stabilize2,shift=shift,stabilize3 = stabilize3,mix_par1 = mix_par1,mix_par2 = mix_par2)
      finite_objective(value)
    }
    
    # Gradient function (optional, for efficiency)
    gradient_function <- function(par) {
      -1 * get_estimated_policy_value_gradient(data, p, rbind(ga_prev,par),method,model=model,stabilize2 = stabilize2,shift=shift,stabilize3 = stabilize3)
    }
    
    
    # Solve with constraints
    # result <- optim(
    #   par = par,
    #   fn = objective_function,
    #   gr = gradient_function,
    #   method = "L-BFGS-B",
    #   lower = rep(lb,3),  # Lower bounds
    #   upper = rep(ub,3)    # Upper bounds 
    # )
    
    result <- optim(
      par = par,
      fn = objective_function,
      method = "L-BFGS-B",
      lower = rep(lb,d),  # Lower bounds
      upper = rep(ub,d)    # Upper bounds 
    )
    
    
    return(rbind(ga_prev,result$par))
    
    
    
  }
  
  if(method!="IPW" & stochastic==TRUE){


    if(model %in% c(1, 2, 3)){
      # Initial guesses
      par <- c(-1, 0, 0, 0)
    }else{
      stop("Unsupported model. Use model = 1, 2, or 3.", call. = FALSE)
    }
    d <- length(par)
    
    # Objective function to maximize (negate for optimization)
    
    if(stabilize3){
      objective_function <- function(par) {
        value <- get_estimated_policy_value(data, p, rbind(ga_prev,par),method = method,model=model,stochastic=TRUE,stabilize = stabilize,stabilize2 = stabilize2,shift=shift,stabilize3 = stabilize3,mix_par1 = mix_par1,mix_par2 = mix_par2)
        finite_objective(value)
      }
    }else{
      objective_function <- function(par) {
        value <- get_estimated_policy_value(data, p, rbind(ga_prev,par),method = method,model=model,stochastic=TRUE,stabilize = stabilize,stabilize2 = stabilize2,shift=shift,stabilize3 = stabilize3,mix_par1 = mix_par1,mix_par2 = mix_par2)
        finite_objective(value)
      }
    }
 
    
    


    
    
    # Solve with constraints
    result <- optim(
      par = par,
      fn = objective_function,
      method = "L-BFGS-B",
      lower = rep(lb,d),  # Lower bounds
      upper = rep(ub,d)    # Upper bounds 
    )
    
    
    return(rbind(ga_prev,result$par))
  }
  
  

  
  
  
  return(rbind(ga_prev,ga))
}
