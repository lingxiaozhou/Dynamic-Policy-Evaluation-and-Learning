# Data-generating process and policy-evaluation helpers for the public simulations.
#
# Public model codes:
#   1: no interaction term, displayed as C = 0
#   2: interaction coefficient C = -2
#   3: interaction coefficient C = -10

public_model_labels <- c("1" = "C = 0", "2" = "C = -2", "3" = "C = -10")

interaction_coefficient <- function(model){
  coeff <- c("1" = 0, "2" = -2, "3" = -10)[as.character(model)]
  if(is.na(coeff)){
    stop("Unsupported model. Use model = 1, 2, or 3.", call. = FALSE)
  }
  unname(coeff)
}

validate_public_model <- function(model){
  if(!model %in% c(1, 2, 3)){
    stop("Unsupported model. Use model = 1, 2, or 3.", call. = FALSE)
  }
  invisible(model)
}

logistic <- function(x){
  1 / (1 + exp(-x))
}

get_history <- function(data, t, model){
  validate_public_model(model)
  n <- ncol(data$W)
  current <- rbind(W = data$W[t,], X = data$X[t,], Z = data$Z[t,], Y = data$Y[t,])
  if(t == 1){
    previous <- array(0, c(4, n))
    rownames(previous) <- c("W", "X", "Z", "Y")
  }else{
    previous <- rbind(
      W = data$W[t - 1,],
      X = data$X[t - 1,],
      Z = data$Z[t - 1,],
      Y = data$Y[t - 1,]
    )
  }
  list(current, previous)
}

get_X_all <- function(data, t, model){
  validate_public_model(model)
  H <- get_history(data, t, model)
  cbind(1, H[[2]]["W",], H[[2]]["X",], H[[2]]["Z",])
}

get_X_all_list <- function(data, model, max_t = nrow(data$W)){
  validate_public_model(model)
  W_prev <- rbind(0, data$W[-nrow(data$W), , drop = FALSE])
  X_prev <- rbind(0, data$X[-nrow(data$X), , drop = FALSE])
  Z_prev <- rbind(0, data$Z[-nrow(data$Z), , drop = FALSE])
  lapply(1:max_t, function(t) cbind(1, W_prev[t,], X_prev[t,], Z_prev[t,]))
}

prepare_policy_cache <- function(data, model){
  time_points <- nrow(data$W)
  data$.X_all_list <- get_X_all_list(data, model, time_points)
  data$.X_all_matrix <- do.call(rbind, data$.X_all_list)
  data
}

get_X_all_cached <- function(data, model, max_t = nrow(data$W)){
  if(!is.null(data$.X_all_list) && length(data$.X_all_list) >= max_t){
    return(data$.X_all_list[1:max_t])
  }
  get_X_all_list(data, model, max_t)
}

get_X_all_matrix <- function(data, model, max_t = nrow(data$W)){
  n <- ncol(data$W)
  n_rows <- max_t * n
  if(!is.null(data$.X_all_matrix) && nrow(data$.X_all_matrix) >= n_rows){
    return(data$.X_all_matrix[1:n_rows, , drop = FALSE])
  }
  do.call(rbind, get_X_all_cached(data, model, max_t))
}

compute_ps <- function(data, t, model){
  validate_public_model(model)
  n <- ncol(data$W)
  H <- get_history(data, t, model)
  W_prev <- H[[2]]["W",]
  X_prev <- H[[2]]["X",]
  Z_prev <- H[[2]]["Z",]
  Y_prev <- H[[2]]["Y",]
  W_prev_excl_mean <- (sum(W_prev) - W_prev) / (n - 1)

  logits <- (-0.8 + 0.06 * X_prev - 0.1 * Z_prev) -
    0.2 * W_prev_excl_mean + 0.2 * W_prev + 0.08 * Y_prev
  logistic(logits)
}

compute_mu <- function(data, t, model, W_counter = NULL, W_counter2 = NULL){
  validate_public_model(model)
  n <- ncol(data$W)
  H <- get_history(data, t, model)

  if(!is.null(W_counter)){
    H[[1]]["W",] <- W_counter
  }
  if(!is.null(W_counter2)){
    H[[2]]["W",] <- W_counter2
  }

  W_cur <- H[[1]]["W",]
  W_prev <- H[[2]]["W",]
  X_prev <- H[[2]]["X",]
  Z_prev <- H[[2]]["Z",]
  Y_prev <- H[[2]]["Y",]

  mean_W_cur_excl <- (sum(W_cur) - W_cur) / (n - 1)
  mean_W_prev_excl <- (sum(W_prev) - W_prev) / (n - 1)
  sum_XW_excl <- sum(X_prev * W_cur) - X_prev * W_cur

  base_mu <- (-0.5 * X_prev + Z_prev) * W_cur +
    sum_XW_excl / (n - 1) +
    0.2 * Y_prev - 0.3 * mean_W_prev_excl + Z_prev -
    0.3 * W_prev * W_cur

  base_mu + interaction_coefficient(model) * Z_prev * W_cur * mean_W_cur_excl
}

compute_W_value <- function(ga, H, ind = NULL, stochastic = FALSE){
  if(is.null(ind)){
    X_all <- cbind(1, H[[2]]["W",], H[[2]]["X",], H[[2]]["Z",])
  }else{
    X_all <- cbind(1, H[[2]]["W", ind], H[[2]]["X", ind], H[[2]]["Z", ind])
  }
  W_value <- X_all %*% ga[nrow(ga),]
  if(stochastic){
    return(logistic(W_value))
  }
  as.numeric(W_value >= 0)
}

compute_utilities <- function(ga, W_counter, interaction){
  get_ga <- function(k) if(is.matrix(ga)) ga[1, k] else ga[k]
  A <- get_ga(1)
  B <- get_ga(2)
  C <- get_ga(3)
  D <- get_ga(4)

  n <- length(W_counter)
  sdT <- sqrt(C^2 + D^2)

  if(sdT == 0){
    return(list(
      E_linear_indicator = (-0.3 * W_counter) * as.numeric(A + B * W_counter >= 0),
      E_cross = rep(0, n)
    ))
  }

  alpha <- (-(A + B * W_counter)) / sdT
  term1 <- ((0.5 * C + D) / sdT) * dnorm(alpha)
  term2 <- (-0.3 * W_counter) * (1 - pnorm(alpha))
  E_linear_indicator <- term1 + term2

  EI <- 1 - pnorm(alpha)
  EZI <- (D / sdT) * dnorm(alpha)
  sum_EI_excl_i <- sum(EI) - EI
  E_cross <- (interaction / (n - 1)) * EZI * sum_EI_excl_i

  list(E_linear_indicator = E_linear_indicator, E_cross = E_cross)
}

compute_true_Q <- function(data, t, model, M = 1, W_counter, ga, stochastic = FALSE){
  validate_public_model(model)
  n <- ncol(data$W)
  H <- get_history(data, t, model)
  interaction <- interaction_coefficient(model)

  if(M == 1){
    H[[1]]["W",] <- W_counter
    return(compute_mu(data, t, model, W_counter = W_counter))
  }

  if(M != 2){
    stop("This public simulation code supports M = 1 or M = 2.", call. = FALSE)
  }

  tmp <- compute_utilities(ga, W_counter, interaction = interaction)
  mean_counter_excl <- (sum(W_counter) - W_counter) / (n - 1)

  part1 <- tmp$E_linear_indicator
  part2 <- (0.1 * H[[2]]["X",] + 0.2 * H[[2]]["Z",] - 0.3) * W_counter
  part3 <- -0.06 * H[[2]]["W",] + 0.02 * H[[2]]["Z",] +
    0.2 * interaction * H[[2]]["Z",] * W_counter * mean_counter_excl
  part4 <- tmp$E_cross
  part1 + part2 + part3 + part4
}

generate_data <- function(p_w = 0.3, n, time_points, model){
  validate_public_model(model)
  W <- array(NA, c(time_points, n))
  X <- array(rnorm(time_points * n), c(time_points, n))
  Z <- array(rnorm(time_points * n), c(time_points, n))
  Y <- array(NA, c(time_points, n))
  p <- array(NA, c(time_points, n))
  data <- list(W = W, X = X, Z = Z, Y = Y, p = p)

  for(t in 1:time_points){
    data$p[t,] <- compute_ps(data, t, model = model)
    data$W[t,] <- rbinom(n, 1, data$p[t,])
    mu <- compute_mu(data, t, model = model)
    data$Y[t,] <- rnorm(n, mu, 1)
  }

  data
}

# Backward-compatible alias for readers comparing with the research scripts.
generateData <- generate_data
