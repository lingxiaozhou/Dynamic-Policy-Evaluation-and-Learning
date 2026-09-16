# Optimize the true policy value for the public simulation models.
get_true_optimal_policy <- function(data,model,ga_prev = NULL,lb=-30,ub=30,MC=TRUE,shift = NULL){
  
  W <- data$W
  if(model %in% c(1, 2, 3)){
    # Initial guesses
    par <- c(0, 0, 0, 0)
  }else{
    stop("Unsupported model. Use model = 1, 2, or 3.", call. = FALSE)
  }
  
  # Objective function to maximize (negate for optimization)
  objective_function <- function(par) {
    -1 * get_true_policy_value(data,model = model,rbind(ga_prev,par),stochastic = TRUE,MC=MC,shift=shift)
  }
  
  # Gradient function (optional, for efficiency)
  # gradient_function <- function(par) {
  #   -1 * get_true_policy_value(W,Y,X=X,model = model,rbind(ga_prev,par),stochastic = TRUE,gradient=TRUE)
  # }
  
  
  # Solve with constraints
  result <- optim(
    par = par,
    fn = objective_function,
    #gr = gradient_function,
    method = "L-BFGS-B",
    lower = rep(lb,length(par)),  # Lower bounds
    upper = rep(ub,length(par))    # Upper bounds 
  )
  
 
  return(list(value = -result$value, ga = rbind(ga_prev,result$par)))
  
  
}
