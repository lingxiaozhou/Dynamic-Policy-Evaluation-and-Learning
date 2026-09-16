# Generate policy learning plots from prepared public plotting data.
#
# Usage:
#   Rscript R/simulation/plots/policy_learning_plot.R [data_file] [output_dir]

args <- commandArgs(TRUE)
data_file <- if(length(args) >= 1) args[1] else file.path("data", "simulation", "policy_learning_plot_data.rds")
output_dir <- if(length(args) >= 2) args[2] else file.path("outputs", "plots", "simulation", "policy_learning")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

plot_data <- readRDS(data_file)
full_plot_df <- plot_data$full_plot_df
plot_df_hajek <- plot_data$plot_df_hajek
keep_time <- plot_data$keep_time
keep_n <- plot_data$keep_n
model_labels <- plot_data$model_labels
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

keep_estimator_full <- c("addIPW", "addIPW2", "addIPW12", "addIPW21", "test_1_5", "oracle")
keep_estimator_hajek <- c("IPW", "addIPW", "addIPW2", "oracle")

pal <- c(
  "IPW" = "#E6C200",
  "addIPW" = "#1b9e77",
  "addIPW2" = "#729ECE",
  "addIPW12" = "#7570b3",
  "addIPW21" = "#C080A1",
  "test_1_5" = "#e41a1c",
  "test_1_10" = "#e41a1c",
  "oracle" = "#000000"
)

label_map_m1 <- list(
  IPW = "Hajek",
  addIPW = expression(hat(V)^"(1)"),
  addIPW2 = expression(hat(V)^"(2)"),
  addIPW12 = expression(hat(V)^"(1,2)"),
  addIPW21 = expression(hat(V)^"(2,1)"),
  test_1_5 = "Test-based",
  test_1_10 = "Test-based",
  oracle = "Oracle"
)

label_map_m2 <- list(
  IPW = "Hajek",
  addIPW = expression(hat(V)^"(1,1)"),
  addIPW2 = expression(hat(V)^"(2,2)"),
  addIPW12 = expression(hat(V)^"(1,2)"),
  addIPW21 = expression(hat(V)^"(2,1)"),
  test_1_5 = "Test-based",
  test_1_10 = "Test-based",
  oracle = "Oracle"
)

label_map_combined_m <- list(
  IPW = "Hajek",
  addIPW = expression(hat(V)^"(1)"*","~hat(V)^"(1,1)"),
  addIPW2 = expression(hat(V)^"(2)"*","~hat(V)^"(2,2)"),
  addIPW12 = expression(hat(V)^"(1,2)"),
  addIPW21 = expression(hat(V)^"(2,1)"),
  test_1_5 = "Test-based",
  test_1_10 = "Test-based",
  oracle = "Oracle"
)

get_labels_m1 <- function(keys) unname(label_map_m1[keys])
get_labels_m2 <- function(keys) unname(label_map_m2[keys])
get_labels_combined_m <- function(keys) unname(label_map_combined_m[keys])

paper_theme <- function(){
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(size = 9),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold")
    )
}

plot_policy <- function(df, yvar, legend_breaks, labels_fun){
  ggplot(df, aes(x = time_point, y = .data[[yvar]], fill = estimator)) +
    geom_boxplot(width = 0.7, outlier.alpha = 0.35, linewidth = 0.3) +
    facet_grid(
      model ~ ps,
      scales = "free_y",
      labeller = labeller(model = as_labeller(model_labels))
    ) +
    scale_fill_manual(
      values = pal,
      breaks = legend_breaks,
      labels = labels_fun(legend_breaks),
      name = "Estimator"
    ) +
    labs(title = "", x = "Time points", y = "Policy value") +
    paper_theme() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
}

save_pdf <- function(filename, plot, width = 8.2, height = 4.7){
  ggsave(filename, plot = plot, width = width, height = height, units = "in", dpi = 600)
}

for(n_value in keep_n){
  df_n <- full_plot_df %>% filter(n == n_value)
  if(!nrow(df_n)) next

  df_m1 <- df_n %>% filter(!is.na(policy_value1))
  breaks_m1 <- intersect(keep_estimator_full, as.character(unique(df_m1$estimator)))
  p1 <- plot_policy(df_m1, yvar = "policy_value1", legend_breaks = breaks_m1, labels_fun = get_labels_m1)

  df_m2 <- df_n %>% filter(!is.na(policy_value2))
  breaks_m2 <- intersect(keep_estimator_full, as.character(unique(df_m2$estimator)))
  p2 <- plot_policy(df_m2, yvar = "policy_value2", legend_breaks = breaks_m2, labels_fun = get_labels_m2)

  save_pdf(file.path(output_dir, paste0("sim_res_ps1_M1_n", n_value, ".pdf")), p1)
  save_pdf(file.path(output_dir, paste0("sim_res_ps1_M2_n", n_value, ".pdf")), p2)

  df_mix_est <- df_n %>%
    filter(ps == "Estimated PS") %>%
    select(-ps) %>%
    pivot_longer(cols = starts_with("policy_value"), names_to = "M", values_to = "policy_value") %>%
    filter(!is.na(policy_value)) %>%
    mutate(M = factor(M, levels = c("policy_value1", "policy_value2"), labels = c("M = 1", "M = 2")))

  breaks_mix <- intersect(keep_estimator_full, as.character(unique(df_mix_est$estimator)))
  p_mix_combined <- ggplot(df_mix_est, aes(x = time_point, y = policy_value, fill = estimator)) +
    geom_boxplot(width = 0.7, outlier.alpha = 0.35, linewidth = 0.3) +
    facet_grid(
      model ~ M,
      scales = "free_y",
      labeller = labeller(model = as_labeller(model_labels))
    ) +
    scale_fill_manual(
      values = pal,
      breaks = breaks_mix,
      labels = get_labels_combined_m(breaks_mix),
      name = "Estimator"
    ) +
    labs(title = "", x = "Time points", y = "Policy value") +
    paper_theme() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

  save_pdf(file.path(output_dir, paste0("sim_res_M1M2_mix_estPS_n", n_value, ".pdf")), p_mix_combined)

  df_est_long <- plot_df_hajek %>%
    filter(n == n_value, ps == "Estimated PS") %>%
    select(-ps) %>%
    pivot_longer(cols = starts_with("policy_value"), names_to = "M", values_to = "policy_value") %>%
    filter(!is.na(policy_value)) %>%
    mutate(M = factor(M, levels = c("policy_value1", "policy_value2"), labels = c("M = 1", "M = 2")))

  breaks_hajek <- intersect(keep_estimator_hajek, as.character(unique(df_est_long$estimator)))
  p_est_combined <- ggplot(df_est_long, aes(x = time_point, y = policy_value, fill = estimator)) +
    geom_boxplot(width = 0.7, outlier.alpha = 0.35, linewidth = 0.3) +
    facet_grid(
      model ~ M,
      scales = "free_y",
      labeller = labeller(model = as_labeller(model_labels))
    ) +
    scale_fill_manual(
      values = pal,
      breaks = breaks_hajek,
      labels = get_labels_combined_m(breaks_hajek),
      name = "Estimator"
    ) +
    labs(title = "", x = "Time points", y = "Policy value") +
    paper_theme() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

  save_pdf(file.path(output_dir, paste0("sim_res_M1M2_Hajek_compare_estPS_n", n_value, ".pdf")), p_est_combined)
}

df_T1000_est_n <- full_plot_df %>%
  filter(time_point == "1000", ps == "Estimated PS") %>%
  select(-ps, -time_point) %>%
  pivot_longer(cols = starts_with("policy_value"), names_to = "M", values_to = "policy_value") %>%
  filter(!is.na(policy_value)) %>%
  mutate(M = factor(M, levels = c("policy_value1", "policy_value2"), labels = c("M = 1", "M = 2")))

if(nrow(df_T1000_est_n)){
  breaks_T1000 <- intersect(keep_estimator_full, as.character(unique(df_T1000_est_n$estimator)))
  p_T1000_est_n <- ggplot(df_T1000_est_n, aes(x = n, y = policy_value, fill = estimator)) +
    geom_boxplot(width = 0.7, outlier.alpha = 0.35, linewidth = 0.3) +
    facet_grid(
      model ~ M,
      scales = "free_y",
      labeller = labeller(model = as_labeller(model_labels))
    ) +
    scale_fill_manual(
      values = pal,
      breaks = breaks_T1000,
      labels = get_labels_combined_m(breaks_T1000),
      name = "Estimator"
    ) +
    labs(title = "", x = "n", y = "Policy value") +
    paper_theme() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

  save_pdf(file.path(output_dir, "sim_res_M1M2_mix_estPS_T1000_by_n.pdf"), p_T1000_est_n)
}

message("Saved policy learning plots to: ", output_dir)
