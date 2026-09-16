# Policy-learning simulation core.
#
# `perform_test()` estimates the candidate policies and then calls
# `perform_test_only()` for the final L2 bootstrap estimator-selection test.

perform_test <- function(data, model, p = 0.3, M = 1,
                         lb = -5, ub = 5, n_points = 10,
                         var_percentage = 0.99,
                         alpha = c(0.05, 0.1), B = 1000,
                         stabilize = FALSE, stabilize2 = FALSE,
                         shift = NULL,
                         methods = c("addIPW", "addIPW2"),
                         test_sta1 = FALSE, test_sta2 = FALSE,
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

  data <- prepare_policy_cache(data, model)
  policy_list <- vector("list", M)
  policy_values <- vector("list", M)
  true_policy_list <- vector("list", M)
  true_policy_values <- vector("list", M)

  for(m in 1:M){
    message("Start time period: ", m)
    method_combinations <- expand.grid(rep(list(c(1, 2)), m))

    policy_list[[m]] <- list()
    policy_values[[m]] <- list()

    true_ga_prev <- if(m == 1) NULL else true_policy_list[[m - 1]]
    true_policy_tmp <- get_true_optimal_policy(
      data = data,
      model = model,
      ga_prev = true_ga_prev,
      lb = lb,
      ub = ub,
      shift = shift
    )
    true_policy_list[[m]] <- true_policy_tmp$ga
    true_policy_values[[m]] <- true_policy_tmp$value

    for(vv in seq_len(nrow(method_combinations))){
      combo <- unlist(method_combinations[vv,])
      combo_name <- paste(combo, collapse = "")
      ga_prev <- if(m == 1){
        NULL
      }else{
        policy_list[[m - 1]][[paste(combo[1:(m - 1)], collapse = "")]]
      }

      method_index <- combo[m]
      policy_list[[m]][[combo_name]] <- get_optimal_policy(
        data = data,
        p = p,
        method = methods[method_index],
        stochastic = TRUE,
        ub = ub,
        lb = lb,
        model = model,
        ga_prev = ga_prev,
        stabilize = stabilize,
        stabilize2 = stabilize2[method_index],
        shift = shift,
        stabilize3 = stabilize3[method_index],
        mix_par1 = mix_par1,
        mix_par2 = mix_par2
      )

      policy_values[[m]][[combo_name]] <- get_true_policy_value(
        data = data,
        model = model,
        ga = policy_list[[m]][[combo_name]],
        stochastic = TRUE,
        MC = FALSE,
        shift = shift
      )
    }
  }

  test_result <- perform_test_only(
    data = data,
    model = model,
    p = p,
    fit = list(policy_list = policy_list),
    M = M,
    lb = lb,
    ub = ub,
    n_points = n_points,
    var_percentage = var_percentage,
    alpha = alpha,
    B = B,
    stabilize = stabilize,
    stabilize2 = stabilize2,
    shift = shift,
    methods = methods,
    total_points = total_points,
    stabilize3 = stabilize3,
    bound2 = bound2,
    mix_par1 = mix_par1,
    mix_par2 = mix_par2,
    grid_seed = grid_seed
  )

  list(
    test = test_result$test,
    test_1 = test_result$test_1,
    policy_list = policy_list,
    policy_values = policy_values,
    true_policy_list = true_policy_list,
    true_policy_values = true_policy_values
  )
}
