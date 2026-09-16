# Policy evaluation simulation.
#
# Usage:
#   Rscript R/simulation/policy_evaluation_simulation.R <sim> <model> <time_points> <n> [output_dir]
#
# Public model codes:
#   1: C = 0
#   2: C = -2
#   3: C = -10

args <- commandArgs(TRUE)
if(length(args) < 4){
  stop("Usage: Rscript R/simulation/policy_evaluation_simulation.R <sim> <model> <time_points> <n> [output_dir]",
       call. = FALSE)
}

sim <- as.integer(args[1])
model <- as.integer(args[2])
time_points <- as.integer(args[3])
n <- as.integer(args[4])
output_dir <- if(length(args) >= 5) args[5] else file.path("outputs", "results", "simulation", "policy_evaluation")

script_path <- function(){
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if(length(file_arg)){
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  normalizePath(file.path("R", "simulation", "policy_evaluation_simulation.R"), mustWork = FALSE)
}

source_all_functions <- function(){
  function_dir <- file.path(dirname(script_path()), "functions")
  if(!dir.exists(function_dir)){
    function_dir <- file.path("R", "simulation", "functions")
  }
  invisible(lapply(sort(list.files(function_dir, pattern = "[.]R$", full.names = TRUE)), source))
}
source_all_functions()

suppressPackageStartupMessages({
  library(mvtnorm)
})

validate_public_model(model)
set.seed(sim)

p <- 0.3
shift <- 0

compute_v_adaptive <- function(ga, V, method, time_points, V_true, data, model, save_test = FALSE){
  tmp <- get_estimated_policy_value(
    data, p = data$p, ga = ga, method = method, model = model,
    shift = 0, save_Y_est = TRUE
  )
  m <- ifelse(is.null(nrow(ga)), 1, 2)
  Y_value <- if(m == 1) rowMeans(data$Y) else rowMeans(data$Y)[-1]

  if(method %in% c("adaptive", "adaptive2")){
    bound <- tmp$weight_t[1:(time_points - m + 1)] * (Y_value - tmp$lambda)
  }else if(method %in% c("shift_adaptive", "shift_adaptive2")){
    bound <- tmp$shift_weight_t[1:(time_points - m + 1)] *
      (Y_value - tmp$mu_hat + tmp$lambda * (tmp$mu_hat - mean(data$Y)))
  }else{
    bound <- tmp$shift_weight_t[1:(time_points - m + 1)] *
      (Y_value - mean(data$Y) + tmp$lambda * (mean(data$Y) - tmp$mu_hat))
  }

  tmp_bound <- bound
  bound <- mean(bound^2, na.rm = TRUE) / time_points
  CI <- V + qnorm(c(0.025, 0.975)) * sqrt(bound)
  coverage <- V_true >= CI[1] & V_true <= CI[2]

  if(save_test){
    return(list(bound = bound, CI = CI, coverage = coverage, tmp_bound = tmp_bound, V = V))
  }
  list(bound = bound, CI = CI, coverage = coverage)
}

compute_v_ipw <- function(ga, V, method, time_points, V_true, data, model){
  bound <- rowMeans(
    get_estimated_policy_value(
      data, p = data$p, ga = ga, method = method, model = model,
      shift = 0, save_Y_est = TRUE
    ),
    na.rm = TRUE
  )
  bound <- mean(bound^2, na.rm = TRUE) / time_points
  CI <- V + qnorm(c(0.025, 0.975)) * sqrt(bound)
  coverage <- V_true >= CI[1] & V_true <= CI[2]
  list(bound = bound, CI = CI, coverage = coverage)
}

compute_v_hajek <- function(ga, V, method, time_points, V_true, data, model){
  bound <- rowMeans(
    get_estimated_policy_value(
      data, p = data$p, ga = ga, method = method, model = model,
      stabilize2 = TRUE, shift = (-V), save_Y_est = TRUE
    ),
    na.rm = TRUE
  )
  bound <- mean(bound^2, na.rm = TRUE) / time_points
  CI <- V + qnorm(c(0.025, 0.975)) * sqrt(bound)
  coverage <- V_true >= CI[1] & V_true <= CI[2]
  list(bound = bound, CI = CI, coverage = coverage)
}

perform_value_test <- function(res_ls1, res_ls2, time_points){
  bound <- mean((res_ls1$tmp_bound - res_ls2$tmp_bound)^2, na.rm = TRUE) / time_points
  CI <- (res_ls1$V - res_ls2$V) + qnorm(c(0.025, 0.975)) * sqrt(bound)
  CI[1] <= 0 && 0 <= CI[2]
}

data <- generate_data(p_w = p, n = n, time_points = time_points, model = model)
if(is.null(data$p)){
  data$p <- array(p, c(time_points, n))
}
ps_fit <- get_estimated_ps(data, model, save_coeff = TRUE)
data$p <- ps_fit$p

ga1 <- c(0, -0.3, 0.5, 1)
ga11 <- rbind(ga1, c(-0.5, -0.06, 0.1, 0.2))
ga2 <- c(-1, 0, 0, 0)
ga22 <- rbind(ga2, ga2)
ga_list <- list(ga1 = ga1, ga11 = ga11, ga2 = ga2, ga22 = ga22)

estimator_vec <- c(
  "truth", "shift", "shift2", "addIPW(H)", "addIPW2(H)",
  "mix_adaptive", "mix_adaptive2", "test-based(adaptive)"
)
method_vec <- c("truth", "shift", "shift2", "addIPW", "addIPW2", "mix_adaptive", "mix_adaptive2", "test")
stabilize2_vec <- c(NA, rep(FALSE, 2), rep(TRUE, 2), rep(FALSE, 3))

Value_arr <- array(NA, c(length(ga_list), length(method_vec)))
mean_weight_arr <- array(NA, c(length(ga_list), length(method_vec)))
colnames(Value_arr) <- estimator_vec
rownames(Value_arr) <- names(ga_list)
colnames(mean_weight_arr) <- estimator_vec
rownames(mean_weight_arr) <- names(ga_list)
ci_res_ls <- list()

for(i in seq_along(ga_list)){
  ci_res_tmp <- list()
  saved_adaptive <- list()

  for(j in seq_along(method_vec)){
    estimator <- estimator_vec[j]
    method <- method_vec[j]
    ga <- ga_list[[i]]

    if(method == "truth"){
      Value_arr[i, j] <- get_true_policy_value(data, model = model, ga = ga)
      ci_res_tmp[[j]] <- NA
      next
    }

    if(estimator == "test-based(adaptive)"){
      test_res <- perform_value_test(saved_adaptive[["mix_adaptive"]], saved_adaptive[["mix_adaptive2"]], time_points)
      selected <- ifelse(test_res, which(estimator_vec == "mix_adaptive"), which(estimator_vec == "mix_adaptive2"))
      Value_arr[i, j] <- Value_arr[i, selected]
      ci_res_tmp[[j]] <- ci_res_tmp[[selected]]
      next
    }

    Value_arr[i, j] <- get_estimated_policy_value(
      data, p = data$p, ga = ga, method = method, model = model,
      stabilize2 = stabilize2_vec[j]
    )

    if(method %in% c("mix_adaptive", "mix_adaptive2")){
      ci_res_tmp[[j]] <- compute_v_adaptive(
        ga = ga, V = Value_arr[i, j], method = method, time_points = time_points,
        V_true = Value_arr[i, 1], data = data, model = model
      )
      saved_adaptive[[method]] <- compute_v_adaptive(
        ga = ga, V = Value_arr[i, j], method = method, time_points = time_points,
        V_true = Value_arr[i, 1], data = data, model = model, save_test = TRUE
      )
    }else if(isTRUE(stabilize2_vec[j])){
      ci_res_tmp[[j]] <- compute_v_hajek(
        ga = ga, V = Value_arr[i, j], method = method, time_points = time_points,
        V_true = Value_arr[i, 1], data = data, model = model
      )
    }else{
      ci_res_tmp[[j]] <- compute_v_ipw(
        ga = ga, V = Value_arr[i, j], method = method, time_points = time_points,
        V_true = Value_arr[i, 1], data = data, model = model
      )
    }
  }

  names(ci_res_tmp) <- estimator_vec
  ci_res_ls[[i]] <- ci_res_tmp
}
names(ci_res_ls) <- names(ga_list)

result <- list(
  sim = sim,
  model = model,
  time_points = time_points,
  n = n,
  propensity_coefficients = ps_fit$coeff,
  mean_weight = mean_weight_arr,
  value = Value_arr,
  ci = ci_res_ls
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  output_dir,
  sprintf("policy_evaluation_sim_%s_model%s_T%s_n%s.rds", sim, model, time_points, n)
)
saveRDS(result, output_file)
message("Saved policy evaluation simulation result to: ", output_file)
