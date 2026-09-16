# Policy learning simulation.
#
# Usage:
#   Rscript R/simulation/policy_learning_simulation.R <sim> <model> <time_points> <n> [output_dir]
#
# Public model codes:
#   1: C = 0
#   2: C = -2
#   3: C = -10

args <- commandArgs(TRUE)
if(length(args) < 4){
  stop("Usage: Rscript R/simulation/policy_learning_simulation.R <sim> <model> <time_points> <n> [output_dir]",
       call. = FALSE)
}

sim <- as.integer(args[1])
model <- as.integer(args[2])
time_points <- as.integer(args[3])
n <- as.integer(args[4])
output_dir <- if(length(args) >= 5) args[5] else file.path("outputs", "results", "simulation", "policy_learning")

script_path <- function(){
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if(length(file_arg)){
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  normalizePath(file.path("R", "simulation", "policy_learning_simulation.R"), mustWork = FALSE)
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
lb <- -5
ub <- 5
shift <- 0

data <- generate_data(p_w = p, n = n, time_points = time_points, model = model)
if(is.null(data$p)){
  data$p <- array(p, c(time_points, n))
}

methods_hajek <- c("addIPW", "addIPW2")
methods_mix <- c("mix_adaptive", "mix_adaptive2")
methods_ipw <- c("addIPW", "IPW")

res_true_ps_hajek <- perform_test(
  data, model, p = data$p, M = 2,
  stabilize2 = c(TRUE, TRUE), stabilize3 = c(FALSE, FALSE),
  shift = shift, methods = methods_hajek, ub = ub, lb = lb,
  mix_par1 = 1, mix_par2 = 0.05
)
res_true_ps_mix <- perform_test(
  data, model, p = data$p, M = 2,
  shift = shift, methods = methods_mix, ub = ub, lb = lb
)
res_true_ps_ipw <- perform_test(
  data, model, p = data$p, M = 2,
  shift = shift, methods = methods_ipw, ub = ub, lb = lb
)

ps_fit <- get_estimated_ps(data, model, save_coeff = TRUE)
ps <- ps_fit$p
coeff <- ps_fit$coeff

res_estimated_ps_hajek <- perform_test(
  data, model, p = ps, M = 2,
  stabilize2 = c(TRUE, TRUE), stabilize3 = c(FALSE, FALSE),
  shift = shift, methods = methods_hajek, ub = ub, lb = lb,
  mix_par1 = 1, mix_par2 = 0.05
)
res_estimated_ps_mix <- perform_test(
  data, model, p = ps, M = 2,
  shift = shift, methods = methods_mix, ub = ub, lb = lb
)
res_estimated_ps_ipw <- perform_test(
  data, model, p = ps, M = 2,
  shift = shift, methods = methods_ipw, ub = ub, lb = lb
)

result <- list(
  sim = sim,
  model = model,
  time_points = time_points,
  n = n,
  propensity_coefficients = coeff,
  true_ps = list(
    hajek = res_true_ps_hajek,
    mix = res_true_ps_mix,
    ipw = res_true_ps_ipw
  ),
  estimated_ps = list(
    hajek = res_estimated_ps_hajek,
    mix = res_estimated_ps_mix,
    ipw = res_estimated_ps_ipw
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  output_dir,
  sprintf("policy_learning_sim_%s_model%s_T%s_n%s.rds", sim, model, time_points, n)
)
saveRDS(result, output_file)
message("Saved policy learning simulation result to: ", output_file)
