get_estimated_ps <- function(data, model, save_coeff = FALSE){
  validate_public_model(model)
  
  W <- data$W
  
  n <- ncol(W)
  time_points <- nrow(W)
  W_neigh <- t(apply(W,1, function(x){sapply(1:n, function(y) mean(x[-y]))}))
  
  df <- data.frame(
    W1 = c(rbind(0,data$W[-time_points,])),
    Y = c(rbind(0,data$Y[-time_points,])),
    X = c(rbind(0,data$X[-time_points,])),
    Z = c(rbind(0,data$Z[-time_points,])),
    W_neigh = c(rbind(0,W_neigh[-time_points,])),
    W = c(W)
  )
  
  fit <- glm(W ~.,family=binomial(link='logit'),data=df)
  p <- array(fit$fitted.values,dim(W))
  if(save_coeff==FALSE){
    return(p)
  }else{
    return(list(p=p,coeff=fit$coefficients))
  }
  
}
