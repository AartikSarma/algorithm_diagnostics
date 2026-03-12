#' Create Conditional Bias Diagnostic Plots
#'
#' Generates diagnostic plots that evaluate clinical algorithms for conditional
#' biases across subgroups. For each combination of dependent variable,
#' independent variable, and grouping variable, the function produces a panel
#' showing percentile-level trends, ntile-aggregated point estimates with
#' confidence intervals, LOESS smoothers, and a beeswarm distribution strip.
#' Panels are assembled into a single figure using \pkg{patchwork}.
#'
#' @param data A data frame containing all variables referenced by the other
#'   arguments.
#' @param dependent_vars Character vector of column names for the dependent
#'   (outcome) variables to plot on the y-axis.
#' @param independent_vars Character vector of column names for the independent
#'   (predictor) variables whose percentiles define the x-axis.
#' @param grouping_vars Character vector of column names used to stratify the
#'   data. Numeric variables are automatically split into terciles; categorical
#'   variables use their existing factor levels.
#' @param color_palettes Optional list of named color vectors, one per grouping
#'   variable. Names must match the factor levels (or tercile labels for numeric
#'   variables). When \code{NULL} (default), the Okabe-Ito colorblind-friendly
#'   palette is used.
#' @param x_size Numeric size of the percentile-level scatter points
#'   (default \code{0.5}).
#' @param n_tiles Integer number of quantile bins for the aggregated point
#'   estimates and confidence intervals (default \code{10}).
#' @param conf_level Numeric confidence level for the interval estimates
#'   (default \code{0.95}). Wilson intervals are used for binary outcomes;
#'   t-based intervals for continuous outcomes.
#' @param y_limits Optional list of length-2 numeric vectors giving
#'   \code{c(lower, upper)} y-axis limits for each dependent variable. When
#'   \code{NULL}, limits are computed automatically.
#' @param min_n Optional minimum number of observations required in an ntile
#'   bin; bins with fewer observations are plotted as \code{NA}.
#' @param dep_var_labels Optional character vector of display labels for the
#'   dependent variables (must match length of \code{dependent_vars}).
#' @param strat_var_labels Optional character vector of display labels for the
#'   grouping variables (must match length of \code{grouping_vars}).
#' @param x_labels Optional character vector of x-axis labels (must match
#'   length of \code{independent_vars}). Defaults to
#'   \code{"Percentile of <var>"}.
#' @param title_size Numeric text size for axis titles and legend titles
#'   (default \code{11}).
#' @param axis_text_size Numeric text size for axis tick labels
#'   (default \code{10}).
#' @param legend_text_size Numeric text size for legend item labels
#'   (default \code{10}).
#'
#' @return A \code{\link[patchwork]{wrap_plots}} object combining all panels
#'   across dependent variables, independent variables, and grouping variables.
#'
#' @examples
#' \dontrun{
#' library(algorithmDiagnostics)
#'
#' result <- conditional_bias_plot(
#'   data = my_data,
#'   dependent_vars = c("outcome1", "outcome2"),
#'   independent_vars = c("predictor1"),
#'   grouping_vars = c("race", "sex"),
#'   n_tiles = 10,
#'   conf_level = 0.95
#' )
#' print(result)
#' }
#'
#' @importFrom dplyr filter group_by summarise mutate n ntile case_when `%>%`
#' @importFrom tidyr complete unnest_wider
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_smooth
#'   coord_cartesian scale_color_manual scale_shape_manual scale_x_continuous
#'   guides labs theme_minimal theme element_text element_blank
#'   position_dodge
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom scales squish
#' @importFrom rlang sym `:=`
#' @importFrom stats sd qt qnorm quantile
#' @export
conditional_bias_plot <- function(data,
                                  dependent_vars,
                                  independent_vars,
                                  grouping_vars,
                                  color_palettes = NULL,
                                  x_size = 0.5,
                                  n_tiles = 10,
                                  conf_level = 0.95,
                                  y_limits = NULL,
                                  min_n = NULL,
                                  dep_var_labels = NULL,
                                  strat_var_labels = NULL,
                                  x_labels = NULL,
                                  title_size = 11,
                                  axis_text_size = 10,
                                  legend_text_size = 10) {

  # Input validation
  if (!is.null(x_labels) && length(x_labels) != length(independent_vars)) {
    stop("Length of x_labels must match number of independent variables")
  }
  if (!is.null(dep_var_labels) && length(dep_var_labels) != length(dependent_vars)) {
    stop("Length of dep_var_labels must match number of dependent variables")
  }
  if (!is.null(strat_var_labels) && length(strat_var_labels) != length(grouping_vars)) {
    stop("Length of strat_var_labels must match number of grouping variables")
  }

  # Define theme for consistent text sizes across plots
  size_theme <- theme(
    axis.title = element_text(size = title_size),
    axis.text = element_text(size = axis_text_size),
    legend.title = element_text(size = title_size),
    legend.text = element_text(size = legend_text_size)
  )

  # Okabe-Ito color palette (colorblind-friendly default)
  okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#000000")

  # Validate color palettes if provided
  if (!is.null(color_palettes)) {
    if (length(color_palettes) != length(grouping_vars)) {
      stop("Length of color_palettes must match number of grouping variables")
    }

    for (g in seq_along(grouping_vars)) {
      if (!is.null(color_palettes[[g]])) {
        if (is.numeric(data[[grouping_vars[g]]])) {
          expected_levels <- c("Tercile 1", "Tercile 2", "Tercile 3")
        } else {
          expected_levels <- levels(factor(data[[grouping_vars[g]]]))
        }

        if (length(color_palettes[[g]]) != length(expected_levels) ||
            !all(names(color_palettes[[g]]) %in% expected_levels)) {
          stop(sprintf(
            "Color palette for stratification variable %s must be a named vector with names matching levels: %s",
            grouping_vars[g], paste(expected_levels, collapse = ", ")
          ))
        }
      }
    }
  }

  # Helper function to calculate confidence intervals
  calc_ci <- function(x, binary = FALSE) {
    if (binary) {
      n <- length(x)
      p <- mean(x, na.rm = TRUE)
      z <- qnorm((1 + conf_level) / 2)

      center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
      spread <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / (1 + z^2 / n)

      return(list(lower = center - spread, upper = center + spread))
    } else {
      n <- length(x)
      se <- sd(x, na.rm = TRUE) / sqrt(n)
      t_val <- qt((1 + conf_level) / 2, df = n - 1)

      return(list(
        lower = mean(x, na.rm = TRUE) - t_val * se,
        upper = mean(x, na.rm = TRUE) + t_val * se
      ))
    }
  }

  # Initialize list to store all plots
  all_strat_plots <- list()

  # Process each stratification variable
  for (g in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[g]
    group_label <- if (!is.null(strat_var_labels)) strat_var_labels[g] else grouping_var

    # Create group factor for current stratification variable
    if (is.numeric(data[[grouping_var]])) {
      data$group_factor <- cut(
        data[[grouping_var]],
        breaks = quantile(data[[grouping_var]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
        labels = c("Tercile 1", "Tercile 2", "Tercile 3"),
        include.lowest = TRUE
      )
    } else {
      data$group_factor <- factor(data[[grouping_var]])
    }

    # Get color scale for this stratification variable
    if (!is.null(color_palettes) && !is.null(color_palettes[[g]])) {
      color_scale <- scale_color_manual(
        name = group_label,
        values = color_palettes[[g]],
        limits = names(color_palettes[[g]])
      )
    } else {
      n_levels <- length(levels(data$group_factor))
      color_scale <- scale_color_manual(
        name = group_label,
        values = okabe_ito[seq_len(n_levels)],
        limits = levels(data$group_factor)
      )
    }

    # Initialize lists for summaries
    all_perc_summaries <- list()
    all_ntile_summaries <- list()
    all_y_ranges <- list()

    # Process each dependent variable
    for (d in seq_along(dependent_vars)) {
      dep_var <- dependent_vars[d]

      perc_summaries <- list()
      ntile_summaries <- list()

      # Process each independent variable
      for (i in seq_along(independent_vars)) {
        var_name <- independent_vars[i]
        perc_col <- paste0("perc_var", i)
        ntile_col <- paste0("ntile_var", i)

        # Create percentile and ntile versions
        if (!(perc_col %in% names(data))) {
          data[[perc_col]] <- ntile(data[[var_name]], 100)
          data[[ntile_col]] <- ntile(data[[var_name]], n_tiles)
        }

        # Calculate percentile summaries
        perc_summaries[[i]] <- data %>%
          filter(!is.na(!!sym(dep_var))) %>%
          group_by(group_factor, !!sym(perc_col)) %>%
          summarise(avg_dep = mean(!!sym(dep_var), na.rm = TRUE), .groups = "drop") %>%
          complete(group_factor, !!sym(perc_col) := 1:100, fill = list(avg_dep = NA))

        # Calculate ntile summaries with confidence intervals
        ntile_summaries[[i]] <- data %>%
          filter(!is.na(!!sym(dep_var))) %>%
          group_by(group_factor, !!sym(ntile_col)) %>%
          summarise(
            n = n(),
            avg_dep = if (is.null(min_n) || n >= min_n) mean(!!sym(dep_var), na.rm = TRUE) else NA,
            ci = if (is.null(min_n) || n >= min_n)
              list(calc_ci(!!sym(dep_var), binary = all(!!sym(dep_var) %in% c(0, 1))))
            else list(list(lower = NA, upper = NA)),
            .groups = "drop"
          ) %>%
          unnest_wider(ci) %>%
          mutate(
            x_pos = !!sym(ntile_col) * (100 / n_tiles) - (50 / n_tiles)
          )
      }

      all_perc_summaries[[d]] <- perc_summaries
      all_ntile_summaries[[d]] <- ntile_summaries

      # Calculate y-axis limits for this dependent variable
      if (is.null(y_limits) || is.null(y_limits[[d]])) {
        y_values <- c(
          sapply(perc_summaries, function(x) x$avg_dep),
          sapply(ntile_summaries, function(x) c(x$lower, x$upper))
        )
        y_values <- unlist(y_values)
        y_range <- range(y_values, na.rm = TRUE)
        y_sd <- sd(y_values, na.rm = TRUE)
        y_mean <- mean(y_values, na.rm = TRUE)
        all_y_ranges[[d]] <- c(
          max(y_mean - (3 * y_sd), min(y_range)) * 0.95,
          min(y_mean + (3 * y_sd), max(y_range)) * 1.05
        )
      } else {
        all_y_ranges[[d]] <- y_limits[[d]]
      }
    }

    # Plot generation parameters
    dodge_width <- 5
    point_size <- 2

    # Function to create plot for each dependent variable
    create_dep_var_panel <- function(d, i, show_y_axis = TRUE) {
      dep_var <- dependent_vars[d]
      y_range <- all_y_ranges[[d]]
      perc_col <- paste0("perc_var", i)
      ntile_col <- paste0("ntile_var", i)

      p <- ggplot() +
        geom_point(
          data = all_perc_summaries[[d]][[i]],
          aes(x = !!sym(perc_col),
              y = squish(avg_dep, y_range),
              color = group_factor,
              shape = case_when(
                squish(avg_dep, y_range) > avg_dep ~ "OOBUp",
                squish(avg_dep, y_range) < avg_dep ~ "OOBDown",
                TRUE ~ "x"
              )),
          size = x_size, alpha = 0.5
        ) +
        geom_point(
          data = all_ntile_summaries[[d]][[i]],
          aes(x = x_pos, y = avg_dep, color = group_factor),
          position = position_dodge(width = dodge_width)
        ) +
        geom_errorbar(
          data = all_ntile_summaries[[d]][[i]],
          aes(x = x_pos, ymin = lower, ymax = upper, color = group_factor),
          position = position_dodge(width = dodge_width),
          width = 2
        ) +
        geom_smooth(
          data = filter(data, !is.na(!!sym(dep_var))),
          aes(x = !!sym(perc_col), y = !!sym(dep_var), color = group_factor),
          method = "loess", se = FALSE, linewidth = 0.5
        ) +
        coord_cartesian(xlim = c(0, 100), ylim = y_range) +
        color_scale +
        scale_shape_manual(values = c("OOBUp" = 6, "x" = 4, "OOBDown" = 2)) +
        guides(shape = "none") +
        labs(
          x = NULL,
          y = if (show_y_axis)
            if (!is.null(dep_var_labels)) dep_var_labels[d] else dep_var
          else NULL
        ) +
        theme_minimal() +
        size_theme

      if (!show_y_axis) {
        p <- p + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())
      }
      return(p)
    }

    # Function to create bottom panel plot (bee swarm distribution)
    create_bottom_panel <- function(i, show_y_axis = TRUE) {
      perc_col <- paste0("perc_var", i)

      p <- ggplot(data, aes(x = !!sym(perc_col), y = group_factor)) +
        geom_quasirandom(aes(color = group_factor), size = 0.1, alpha = 0.5) +
        scale_x_continuous(limits = c(0, 100)) +
        color_scale +
        labs(
          x = if (is.null(x_labels)) paste("Percentile of", independent_vars[i]) else x_labels[i],
          y = if (show_y_axis) group_label else NULL
        ) +
        theme_minimal() +
        theme(legend.position = "none") +
        guides(color = "none") +
        size_theme

      if (!show_y_axis) {
        p <- p + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())
      }
      return(p)
    }

    # Create plots for this stratification variable
    strat_plots <- list()

    for (d in seq_along(dependent_vars)) {
      row_plots <- list()
      for (i in seq_along(independent_vars)) {
        row_plots[[i]] <- create_dep_var_panel(d, i, show_y_axis = (i == 1))
      }
      strat_plots[[d]] <- Reduce(`+`, row_plots)
    }

    # Create bottom row (bee swarm plots)
    bottom_plots <- list()
    for (i in seq_along(independent_vars)) {
      bottom_plots[[i]] <- create_bottom_panel(i, show_y_axis = (i == 1))
    }
    strat_plots[[length(strat_plots) + 1]] <- Reduce(`+`, bottom_plots)

    # Combine plots for this stratification variable with proper height ratios
    all_strat_plots[[g]] <- wrap_plots(
      strat_plots,
      ncol = 1,
      heights = c(rep(4, length(dependent_vars)), 1)
    ) +
      plot_layout(guides = "collect")
  }

  # Create final combined plot with all stratification variables
  final_plot <- wrap_plots(
    all_strat_plots,
    ncol = 1,
    heights = rep(1, length(grouping_vars))
  ) +
    plot_layout(guides = "collect")

  return(final_plot)
}
