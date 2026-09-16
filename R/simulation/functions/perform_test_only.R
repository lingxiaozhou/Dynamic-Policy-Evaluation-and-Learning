# Recompute only the L2 bootstrap estimator-selection test for an existing
# policy-learning fit. This is the test used in the paper.

perform_test_only <- function(data, model, p = 0.3, fit, M = 2,
                              lb = -5, ub = 5, n_points = 10,
                              var_percentage = 0.99,
                              alpha = c(0.05, 0.1), B = 1000,
                              stabilize = FALSE, stabilize2 = FALSE,
                              shift = NULL,
                              methods = c("addIPW", "addIPW2"),
                              total_points = 1000,
                              stabilize3 = FALSE, bound2 = FALSE,
                              mix_par1 = 1.5, mix_par2 = 0.1,
                              seed = NULL, grid_seed = NULL){
  validate_public_model(model)
  if(length(stabilize2) == 1){
    stabilize2 <- rep(stabilize2, 2)
  }
  if(length(stabilize3) == 1){
    stabilize3 <- rep(stabilize3, 2)
  }
  if(is.null(shift)){
    shift <- 0
  }
  if(!is.null(seed)){
    set.seed(seed)
  }

  time_points <- nrow(data$W)
  data <- prepare_policy_cache(data, model)
  policy_list <- fit$policy_list
  test <- array(NA, c(length(alpha), M), dimnames = list(alpha, 1:M))

  if(!is.null(grid_seed)){
    rng_state <- if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
      get(".Random.seed", envir = .GlobalEnv)
    }else{
      NULL
    }
    set.seed(grid_seed)
  }
  gamma_grid <- data.frame(
    gamma0 = runif(total_points, min = lb, max = ub),
    gamma1 = runif(total_points, min = lb, max = ub),
    gamma2 = runif(total_points, min = lb, max = ub),
    gamma3 = runif(total_points, min = lb, max = ub)
  )
  if(!is.null(grid_seed)){
    if(!is.null(rng_state)){
      assign(".Random.seed", rng_state, envir = .GlobalEnv)
    }else if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }

  N <- nrow(gamma_grid)
  tmp_func <- array(NA, c(2 * N, time_points))

  get_contribution <- function(tmp, method, m, Y_value){
    method_family <- sub("2$", "", method)

    if(method_family == "adaptive"){
      out <- tmp$weight_t[1:(time_points - m + 1)] * (Y_value - tmp$lambda)
    }else if(method_family == "shift_adaptive"){
      out <- tmp$shift_weight_t[1:(time_points - m + 1)] *
        (Y_value - tmp$mu_hat + tmp$lambda * (tmp$mu_hat - mean(data$Y)))
    }else if(method_family == "mix_adaptive"){
      out <- tmp$shift_weight_t[1:(time_points - m + 1)] *
        (Y_value - mean(data$Y) + tmp$lambda * (mean(data$Y) - tmp$mu_hat))
    }else if(method_family == "shift"){
      out <- rowMeans(tmp)
    }else{
      out <- tmp[[2]] *
        (Y_value - (tmp[[1]] * mean(tmp[[2]] * Y_value) +
                      (1 - tmp[[1]]) * mean(data$Y + shift)))
    }

    if(m == 2 && method_family != "shift"){
      out <- c(out, NA)
    }
    out
  }

  for(m in 1:M){
    if(m == 1){
      method_prev_combinations <- array(1, c(1, 1))
    }else{
      method_prev_combinations <- unique(test)[, 1:(m - 1), drop = FALSE] + 1
      combo_name <- apply(method_prev_combinations, 1, function(x) paste(x, collapse = ""))
      test_tmp <- array(
        NA,
        c(length(alpha), nrow(method_prev_combinations)),
        dimnames = list(alpha, combo_name)
      )
    }

    for(vv in 1:nrow(method_prev_combinations)){
      if(m == 1){
        ga_prev <- NULL
        method_prev <- NULL
        prev_name <- ""
      }else{
        prev_name <- paste(method_prev_combinations[vv, ], collapse = "")
        ga_prev <- policy_list[[m - 1]][[prev_name]]
        method_prev <- methods[method_prev_combinations[vv, ]]
      }

      GP <- rep(NA_real_, N)
      for(i in 1:N){
        ga <- rbind(ga_prev, unlist(gamma_grid[i, ]))

        value1 <- get_estimated_policy_value(
          data, p, ga, method = c(method_prev, methods[1]), model = model,
          stochastic = TRUE, save_Y_est = FALSE, stabilize = stabilize,
          stabilize2 = stabilize2[1], shift = shift, stabilize3 = stabilize3[1],
          mix_par1 = mix_par1, mix_par2 = mix_par2
        )
        value2 <- get_estimated_policy_value(
          data, p, ga, method = c(method_prev, methods[2]), model = model,
          stochastic = TRUE, save_Y_est = FALSE, stabilize = stabilize,
          stabilize2 = stabilize2[2], shift = shift, stabilize3 = stabilize3[2],
          mix_par1 = mix_par1, mix_par2 = mix_par2
        )
        GP[i] <- value1 - value2

        if(stabilize2[1]){
          if(bound2){
            bound_shift <- get_estimated_policy_value(
              data, p, ga, method = c(method_prev, methods[1]), model = model,
              stochastic = TRUE, save_Y_est = FALSE, stabilize = stabilize,
              stabilize2 = stabilize2[1], shift = shift, stabilize3 = stabilize3[1],
              bound_weight = TRUE, mix_par1 = mix_par1, mix_par2 = mix_par2
            )
            tmp1 <- rowMeans(get_estimated_policy_value(
              data, p, ga, method = c(method_prev, methods[1]), model = model,
              stochastic = TRUE, save_Y_est = TRUE, stabilize = stabilize,
              stabilize2 = stabilize2[1], shift = shift - bound_shift,
              mix_par1 = mix_par1, mix_par2 = mix_par2
            ))
          }else{
            tmp1 <- rowMeans(get_estimated_policy_value(
              data, p, ga, method = c(method_prev, methods[1]), model = model,
              stochastic = TRUE, save_Y_est = TRUE, stabilize = stabilize,
              stabilize2 = stabilize2[1], shift = shift - value1,
              mix_par1 = mix_par1, mix_par2 = mix_par2
            ))
          }
        }else if(!stabilize3[1] && !sub("2$", "", methods[1]) %in% c("adaptive", "shift_adaptive", "mix_adaptive", "shift")){
          tmp1 <- rowMeans(get_estimated_policy_value(
            data, p, ga, method = c(method_prev, methods[1]), model = model,
            stochastic = TRUE, save_Y_est = TRUE, stabilize = stabilize,
            stabilize2 = stabilize2[1], shift = shift,
            mix_par1 = mix_par1, mix_par2 = mix_par2
          ))
        }else{
          raw1 <- get_estimated_policy_value(
            data, p, ga, method = c(method_prev, methods[1]), model = model,
            stochastic = TRUE, save_Y_est = TRUE, stabilize = stabilize,
            stabilize2 = stabilize2[1], shift = shift, stabilize3 = stabilize3[1],
            mix_par1 = mix_par1, mix_par2 = mix_par2
          )
          Y_value <- if(m == 1) rowMeans(data$Y) + shift else rowMeans(data$Y)[-1] + shift
          tmp1 <- get_contribution(raw1, methods[1], m, Y_value)
        }

        if(stabilize2[2]){
          if(bound2){
            bound_shift <- get_estimated_policy_value(
              data, p, ga, method = c(method_prev, methods[2]), model = model,
              stochastic = TRUE, save_Y_est = FALSE, stabilize = stabilize,
              stabilize2 = stabilize2[2], shift = shift, stabilize3 = stabilize3[2],
              bound_weight = TRUE, mix_par1 = mix_par1, mix_par2 = mix_par2
            )
            tmp2 <- rowMeans(get_estimated_policy_value(
              data, p, ga, method = c(method_prev, methods[2]), model = model,
              stochastic = TRUE, save_Y_est = TRUE, stabilize = stabilize,
              stabilize2 = stabilize2[2], shift = shift - bound_shift,
              mix_par1 = mix_par1, mix_par2 = mix_par2
            ))
          }else{
            tmp2 <- rowMeans(get_estimated_policy_value(
              data, p, ga, method = c(method_prev, methods[2]), model = model,
              stochastic = TRUE, save_Y_est = TRUE, stabilize = stabilize,
              stabilize2 = stabilize2[2], shift = shift - value2,
              mix_par1 = mix_par1, mix_par2 = mix_par2
            ))
          }
        }else if(!stabilize3[2] && !sub("2$", "", methods[2]) %in% c("adaptive", "shift_adaptive", "mix_adaptive", "shift")){
          tmp2 <- rowMeans(get_estimated_policy_value(
            data, p, ga, method = c(method_prev, methods[2]), model = model,
            stochastic = TRUE, save_Y_est = TRUE, stabilize = stabilize,
            stabilize2 = stabilize2[2], shift = shift,
            mix_par1 = mix_par1, mix_par2 = mix_par2
          ))
        }else{
          raw2 <- get_estimated_policy_value(
            data, p, ga, method = c(method_prev, methods[2]), model = model,
            stochastic = TRUE, save_Y_est = TRUE, stabilize = stabilize,
            stabilize2 = stabilize2[2], shift = shift, stabilize3 = stabilize3[2],
            mix_par1 = mix_par1, mix_par2 = mix_par2
          )
          Y_value <- if(m == 1) rowMeans(data$Y) + shift else rowMeans(data$Y)[-1] + shift
          tmp2 <- get_contribution(raw2, methods[2], m, Y_value)
        }

        tmp_func[i, ] <- tmp1
        tmp_func[i + N, ] <- tmp2
      }

      n_eval_times <- time_points - m + 1
      contrast_func <- tmp_func[1:N, 1:n_eval_times, drop = FALSE] -
        tmp_func[(N + 1):(2 * N), 1:n_eval_times, drop = FALSE]
      contrast_func[!is.finite(contrast_func)] <- 0

      Z <- matrix(rnorm(n_eval_times * B), nrow = n_eval_times, ncol = B)
      simulated_paths <- contrast_func %*% Z / sqrt(n_eval_times)
      observed_process <- sqrt(n_eval_times) * GP

      observed_l2 <- mean(observed_process^2, na.rm = TRUE)
      bootstrap_l2 <- colMeans(simulated_paths^2, na.rm = TRUE)
      critical_value <- as.numeric(quantile(bootstrap_l2, 1 - alpha / M, na.rm = TRUE))
      reject <- as.numeric(observed_l2 > critical_value)

      if(m == 1){
        test[, m] <- reject
      }else{
        test_tmp[, prev_name] <- reject
      }
    }

    if(m > 1){
      prev_test <- apply(test, 1, function(x) paste(x[1:(m - 1)] + 1, collapse = ""))
      test[, m] <- sapply(seq_along(alpha), function(x) test_tmp[x, prev_test[x]])
    }
  }

  list(
    test = test,
    test_1 = test
  )
}
