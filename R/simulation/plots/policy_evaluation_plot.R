# Generate policy evaluation plots from prepared public plotting data.
#
# Usage:
#   Rscript R/simulation/plots/policy_evaluation_plot.R [data_file] [output_dir]

args <- commandArgs(TRUE)
data_file <- if(length(args) >= 1) args[1] else file.path("data", "simulation", "policy_evaluation_plot_data.rds")
output_dir <- if(length(args) >= 2) args[2] else file.path("outputs", "plots", "simulation", "policy_evaluation")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

plot_data <- readRDS(data_file)
policy_summary <- plot_data$policy_summary
truth_summary <- plot_data$truth_summary
model_labels <- plot_data$model_labels
T_levels <- plot_data$keep_time
n_levels <- plot_data$keep_n
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

estimator_labels <- c(
  "mix_adaptive" = "addIPW(mix)",
  "mix_adaptive2" = "addIPW2(mix)",
  "shift" = "addIPW(S)",
  "shift2" = "addIPW2(S)",
  "addIPW(H)" = "Hajek",
  "addIPW2(H)" = "Hajek2"
)

estimator_shapes <- c(
  "mix_adaptive" = 16,
  "mix_adaptive2" = 17,
  "shift" = 15,
  "shift2" = 18,
  "addIPW(H)" = 15,
  "addIPW2(H)" = 18
)

plot_theme <- function(){
  theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold")
    )
}

prep_plot_data <- function(n_value, ga_keep, estimator_keep, ga_labels){
  policy_df <- policy_summary %>%
    filter(n == n_value, ga %in% ga_keep, estimator %in% estimator_keep) %>%
    mutate(
      T = factor(T, levels = T_levels),
      model = factor(model, levels = c(1, 2, 3)),
      ga = factor(ga, levels = ga_keep),
      estimator = factor(estimator, levels = estimator_keep)
    )

  truth_df <- truth_summary %>%
    filter(n == n_value, ga %in% ga_keep) %>%
    mutate(
      T = factor(T, levels = T_levels),
      model = factor(model, levels = c(1, 2, 3)),
      ga = factor(ga, levels = ga_keep)
    )

  list(policy = policy_df, truth = truth_df, ga_labels = ga_labels)
}

old_default_colors <- function(estimator_order){
  stats::setNames(scales::hue_pal()(length(estimator_order)), estimator_order)
}

make_policy_plot <- function(plot_data){
  dodge_w <- 0.6
  estimator_order <- levels(plot_data$policy$estimator)
  ggplot(plot_data$policy, aes(x = T, y = mean_value, color = estimator)) +
    geom_point(position = position_dodge(width = dodge_w), size = 1.9) +
    geom_errorbar(
      aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
      position = position_dodge(width = dodge_w),
      width = 0.22
    ) +
    geom_segment(
      data = plot_data$truth,
      aes(
        x = as.numeric(T) - (dodge_w / 2 - 0.01),
        xend = as.numeric(T) + (dodge_w / 2 - 0.01),
        y = truth_mean,
        yend = truth_mean
      ),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      linewidth = 0.7
    ) +
    facet_grid(
      rows = vars(model),
      cols = vars(ga),
      scales = "free_y",
      labeller = labeller(model = as_labeller(model_labels), ga = as_labeller(plot_data$ga_labels))
    ) +
    scale_color_manual(
      values = old_default_colors(estimator_order),
      breaks = estimator_order,
      labels = estimator_labels[estimator_order],
      drop = FALSE
    ) +
    labs(x = "Time points", y = "Estimated value") +
    plot_theme()
}

make_rmse_plot <- function(plot_data){
  estimator_order <- levels(plot_data$policy$estimator)
  ggplot(plot_data$policy, aes(x = T, y = rmse, color = estimator, shape = estimator)) +
    geom_line(aes(group = estimator), linetype = "dashed", linewidth = 0.5, alpha = 0.75) +
    geom_point(size = 2.6) +
    facet_grid(
      rows = vars(model),
      cols = vars(ga),
      scales = "free_y",
      labeller = labeller(model = as_labeller(model_labels), ga = as_labeller(plot_data$ga_labels))
    ) +
    scale_color_manual(
      values = old_default_colors(estimator_order),
      breaks = estimator_order,
      labels = estimator_labels[estimator_order],
      drop = FALSE
    ) +
    scale_shape_manual(
      values = estimator_shapes,
      breaks = estimator_order,
      labels = estimator_labels[estimator_order],
      drop = FALSE
    ) +
    labs(x = "Time points", y = "RMSE") +
    plot_theme()
}

plot_specs <- list(
  mix_shift = list(
    estimators = c("mix_adaptive", "shift", "mix_adaptive2", "shift2"),
    prefix = "mix_shift"
  ),
  mix_hajek = list(
    estimators = c("mix_adaptive", "addIPW(H)", "mix_adaptive2", "addIPW2(H)"),
    prefix = "mix_hajek"
  )
)

ga_specs <- list(
  M1 = list(ga = c("ga1", "ga2"), labels = c("ga1" = "Policy1", "ga2" = "Policy2")),
  M2 = list(ga = c("ga11", "ga22"), labels = c("ga11" = "Policy1", "ga22" = "Policy2"))
)

for(n_value in n_levels){
  n_dir <- file.path(output_dir, paste0("n", n_value))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)

  for(spec in plot_specs){
    for(m_name in names(ga_specs)){
      ga_spec <- ga_specs[[m_name]]
      plot_input <- prep_plot_data(
        n_value = n_value,
        ga_keep = ga_spec$ga,
        estimator_keep = spec$estimators,
        ga_labels = ga_spec$labels
      )
      if(!nrow(plot_input$policy)) next

      ggsave(
        file.path(n_dir, paste0("policy_", spec$prefix, "_", m_name, "_n", n_value, ".pdf")),
        make_policy_plot(plot_input),
        width = 10,
        height = 5,
        dpi = 300
      )
      ggsave(
        file.path(n_dir, paste0("RMSE_", spec$prefix, "_", m_name, "_n", n_value, ".pdf")),
        make_rmse_plot(plot_input),
        width = 10,
        height = 5,
        dpi = 300
      )
    }
  }
}

message("Saved policy evaluation plots to: ", output_dir)
