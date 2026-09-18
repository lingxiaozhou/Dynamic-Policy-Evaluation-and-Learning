suppressPackageStartupMessages({
  library(cowplot)
  library(data.table)
  library(dplyr)
  library(forcats)
  library(geepack)
  library(geocausal)
  library(ggplot2)
  library(httr)
  library(lubridate)
  library(purrr)
  library(readr)
  library(sf)
  library(spatstat)
  library(tibble)
  library(tidytext)
  library(tidyr)
  library(viridis)
  library(zoo)
})

get_script_path <- function() {
  file_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]])))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    path <- rstudioapi::getSourceEditorContext()$path
    if (!is.null(path) && nzchar(path)) {
      return(normalizePath(path))
    }
  }
  NA_character_
}

script_path <- get_script_path()
repo_dir <- if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."))
} else {
  normalizePath(".")
}

source(file.path(repo_dir, "R/application/functions/application_utilities.R"))

data_dir <- file.path(repo_dir, "data/application")
plot_dir <- file.path(repo_dir, "outputs/plots/application")
result_dir <- file.path(repo_dir, "outputs/results")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

aid_csv <- file.path(data_dir, "Maaws_PublicData.csv")
troop_density_csv <- file.path(data_dir, "CSM_filter.csv")
prepared_data_file <- file.path(data_dir, "prepared_application_data.rds")
analysis_data_file <- file.path(data_dir, "application_data.rds")
iraq_sf_file <- file.path(data_dir, "iraq_district_sf.rds")
ps_diagnostic_file <- file.path(data_dir, "propensity_score_diagnostics.rds")

if (!file.exists(aid_csv)) {
  stop("Missing required local input: ", aid_csv)
}
if (!file.exists(troop_density_csv)) {
  stop("Missing required local input: ", troop_density_csv)
}

application_result_file <- function(outcome_name, policy_name) {
  file.path(result_dir, paste0("res_", outcome_name, "_", policy_name, ".dat"))
}

load_application_data <- function() {
  saved <- readRDS(analysis_data_file)
  list(
    data = saved$data,
    iraq_district = readRDS(iraq_sf_file)
  )
}

load_policy_result <- function(outcome_name, policy_name) {
  result_path <- application_result_file(outcome_name, policy_name)
  if (!file.exists(result_path)) {
    stop("Cannot find saved result: ", result_path, call. = FALSE)
  }

  result_env <- new.env(parent = emptyenv())
  load(result_path, envir = result_env)

  object_name <- paste0("out_", outcome_name)
  if (!exists(object_name, envir = result_env, inherits = FALSE)) {
    stop("Saved result does not contain ", object_name, ": ", result_path, call. = FALSE)
  }
  get(object_name, envir = result_env)
}

load_warm_start <- function(outcome_name, policy_cols, n_periods = 4) {
  result_path <- application_result_file(outcome_name, "wopop")
  if (!file.exists(result_path)) {
    return(NULL)
  }
  previous_out <- load_policy_result(outcome_name, "wopop")
  make_warm_start_initial_par_list(previous_out, policy_cols, n_periods = n_periods)
}


# Part 1: Prepare Data -------------------------------------------------------

message("Part 1: preparing district-week application data")

application <- prepare_application_data(
  aid_csv = aid_csv,
  troop_density_csv = troop_density_csv
)

saveRDS(application, prepared_data_file)


# Part 2: Fit Propensity Model And Diagnostics -------------------------------

message("Part 2: fitting propensity model and saving diagnostics")

if (!exists("application")) {
  application <- readRDS(prepared_data_file)
}

# Edit this vector when checking alternative propensity-score specifications.
ps_model_vars <- c(
  "prev_SAF_7", "prev_SAF_14", "prev_SAF_28",
  "prev_all_outcome_7","prev_all_outcome_14","prev_all_outcome_28",
  "prev_all_airstrike_7", "prev_all_airstrike_14", "prev_all_airstrike_28",
  "pop", "time_1", "time_2", "time_3", "surge", "tf_array_US", "tf_array_US_missing",
  "rivers_dist", "routes_dist", "cities_dist",
  "prev_aid_costamount_7", "prev_aid_costamount_14", "prev_aid_costamount_28",
  "prev_aid_array_7", "prev_aid_array_14", "prev_aid_array_28"
)

application <- fit_application_propensity(
  application,
  ps_vars = ps_model_vars
)

data <- application$data
iraq_district <- application$iraq_district

saveRDS(
  list(
    data = data,
    ps_formula = application$ps_formula,
    ps_vars = application$ps_vars,
    balance_vars = application$balance_vars,
    source_files = list(
      aid_csv = normalizePath(aid_csv),
      troop_density_csv = normalizePath(troop_density_csv)
    ),
    dataverse_file_ids = list(
      iraq_district = 8138903,
      cities_dist = 8079278,
      rivers_dist = 8079275,
      routes_dist = 8079274,
      ethnicity = 8079276
    )
  ),
  analysis_data_file
)
saveRDS(iraq_district, iraq_sf_file)

saveRDS(
  list(
    df = application$ps_df,
    ps_vars = application$ps_vars,
    balance_vars = application$balance_vars,
    ps_formula = application$ps_formula,
    diagnostic_weight = "raw inverse-probability weight"
  ),
  ps_diagnostic_file
)

save_application_diagnostics(
  ps_df = application$ps_df,
  balance_vars = application$balance_vars,
  plot_dir = file.path(plot_dir, "diagnostics"),
  ps_model_vars = application$ps_vars
)


# Part 3: Perform Policy Learning -------------------------------------------

message("Part 3: performing policy learning")

if (!exists("data") || !exists("iraq_district")) {
  loaded_application <- load_application_data()
  data <- loaded_application$data
  iraq_district <- loaded_application$iraq_district
}

policy_search_points <- 1000
bootstrap_reps <- 1000
policy_seed <- 1234

# This is on by default to discourage policies with negative horizon-specific
# estimated values. Set to 0 here for debugging if needed.
value_barrier <- 1e5

wopop_spec <- application_policy_spec(include_population = FALSE)
wpop_spec <- application_policy_spec(include_population = TRUE)
data <- add_application_policy_terms(data)

message("  Running SAF without-population policy")
set.seed(policy_seed)
out_SAF_wopop <- run_all_times(
  data = data,
  initial_par_list = make_initial_par_list(wopop_spec$policy_cols_SAF, 4),
  X_cols = wopop_spec$X_cols,
  Z_cols = wopop_spec$Z_cols,
  policy_cols = wopop_spec$policy_cols_SAF,
  Y_name = "SAF_week",
  always_aid_threshold = 0.95,
  no_aid_threshold = 0.05,
  value_barrier = value_barrier,
  iraq_district = iraq_district,
  alpha = 0.05 / 4,
  total_points = policy_search_points,
  B = bootstrap_reps
)
out_SAF <- out_SAF_wopop
save(out_SAF, file = application_result_file("SAF", "wopop"))

message("  Running IED without-population policy")
set.seed(policy_seed)
out_IED_wopop <- run_all_times(
  data = data,
  initial_par_list = make_initial_par_list(wopop_spec$policy_cols_IED, 4),
  X_cols = wopop_spec$X_cols,
  Z_cols = wopop_spec$Z_cols,
  policy_cols = wopop_spec$policy_cols_IED,
  Y_name = "IED_week",
  always_aid_threshold = 0.95,
  no_aid_threshold = 0.05,
  value_barrier = value_barrier,
  iraq_district = iraq_district,
  alpha = 0.05 / 4,
  total_points = policy_search_points,
  B = bootstrap_reps
)
out_IED <- out_IED_wopop
save(out_IED, file = application_result_file("IED", "wopop"))

message("  Running SAF with-population policy")
initial_par_list_SAF_wpop <- load_warm_start("SAF", wpop_spec$policy_cols_SAF)
if (is.null(initial_par_list_SAF_wpop)) {
  initial_par_list_SAF_wpop <- make_initial_par_list(wpop_spec$policy_cols_SAF, 4)
}
set.seed(policy_seed)
out_SAF_wpop <- run_all_times(
  data = data,
  initial_par_list = initial_par_list_SAF_wpop,
  X_cols = wpop_spec$X_cols,
  Z_cols = wpop_spec$Z_cols,
  policy_cols = wpop_spec$policy_cols_SAF,
  Y_name = "SAF_week",
  always_aid_threshold = 0.95,
  no_aid_threshold = 0.05,
  value_barrier = value_barrier,
  iraq_district = iraq_district,
  alpha = 0.05 / 4,
  total_points = policy_search_points,
  B = bootstrap_reps
)
out_SAF <- out_SAF_wpop
save(out_SAF, file = application_result_file("SAF", "wpop"))

message("  Running IED with-population policy")
initial_par_list_IED_wpop <- load_warm_start("IED", wpop_spec$policy_cols_IED)
if (is.null(initial_par_list_IED_wpop)) {
  initial_par_list_IED_wpop <- make_initial_par_list(wpop_spec$policy_cols_IED, 4)
}
set.seed(policy_seed)
out_IED_wpop <- run_all_times(
  data = data,
  initial_par_list = initial_par_list_IED_wpop,
  X_cols = wpop_spec$X_cols,
  Z_cols = wpop_spec$Z_cols,
  policy_cols = wpop_spec$policy_cols_IED,
  Y_name = "IED_week",
  always_aid_threshold = 0.95,
  no_aid_threshold = 0.05,
  value_barrier = value_barrier,
  iraq_district = iraq_district,
  alpha = 0.05 / 4,
  total_points = policy_search_points,
  B = bootstrap_reps
)
out_IED <- out_IED_wpop
save(out_IED, file = application_result_file("IED", "wpop"))


# Part 4: Plot And Print Tables ---------------------------------------------

message("Part 4: saving application plots and printing tables")

if (!exists("data") || !exists("iraq_district")) {
  loaded_application <- load_application_data()
  data <- loaded_application$data
  iraq_district <- loaded_application$iraq_district
}

if (!exists("out_SAF_wopop")) out_SAF_wopop <- load_policy_result("SAF", "wopop")
if (!exists("out_IED_wopop")) out_IED_wopop <- load_policy_result("IED", "wopop")
if (!exists("out_SAF_wpop")) out_SAF_wpop <- load_policy_result("SAF", "wpop")
if (!exists("out_IED_wpop")) out_IED_wpop <- load_policy_result("IED", "wpop")

save_application_plots(
  data = data,
  iraq_district = iraq_district,
  out_SAF = out_SAF_wopop,
  out_IED = out_IED_wopop,
  policy_name = "wopop",
  output_dir = plot_dir
)

print_application_tables(
  data = data,
  out_SAF = out_SAF_wopop,
  out_IED = out_IED_wopop,
  case_name = "wopop"
)

save_application_plots(
  data = data,
  iraq_district = iraq_district,
  out_SAF = out_SAF_wpop,
  out_IED = out_IED_wpop,
  policy_name = "wpop",
  output_dir = plot_dir
)

save_application_governorate_label_palette(plot_dir)

print_application_tables(
  data = data,
  out_SAF = out_SAF_wpop,
  out_IED = out_IED_wpop,
  case_name = "wpop"
)

message("Application workflow complete.")
message("Generated data:  ", analysis_data_file)
message("Generated plots: ", plot_dir)
