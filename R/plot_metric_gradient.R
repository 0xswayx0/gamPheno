#' Plot Phenological Metric Across a Continuous Predictor Gradient
#'
#' Evaluates and plots a chosen phenological metric (with 95% Credible Intervals)
#' across the observed range of a continuous environmental predictor, extracting
#' the training dataset directly from `model$model` while holding other predictors
#' at their mean (or user-specified presets via `...`).
#'
#' @param model A fitted `mgcv` GAM or BAM object.
#' @param metric Character, the metric to evaluate: `"centroid"` (default), `"duration"`, `"modality"`, `"peak"`, `"divergence"`, `"onset"`, or `"offset"`.
#' @param x_axis Character, name of the continuous predictor in `model$model` to vary along the x-axis (e.g., `"Latitude"` or `"GDD5_YTD"`). Matched case-insensitively.
#' @param doy Character, name of the Day-of-Year column (default = `"DOY"`).
#' @param n_steps Integer, number of evaluation points along the predictor gradient (default = 30).
#' @param n_sims Integer, number of posterior draws per evaluation point (default = 1000).
#' @param doby_origin Character (`"auto"`, `"scenario"`, `"joint"`, or `"ref"`), integer/numeric calendar `DOY` in `(0, 366)`, or reference simulation matrix passed to `DOBY`-dependent metric estimators (default = `"auto"`).
#' @param add Logical, whether to return `ggplot2` layers to add to an existing plot (default = FALSE).
#' @param color Character, hex color for the line and credible interval ribbon (default = "#009E73").
#' @param ... Optional named arguments passed to [create_newdata()] to hold other covariates at specific values or quantiles (default = `"mean"`).
#' @return A `ggplot` object (if `add = FALSE`) or a list of `ggplot2` layers (if `add = TRUE`).
#' @examples
#' data(sim_carabids)
#' mod <- mgcv::gam(
#'   Count ~ s(DOY, bs = "cc", k = 10) + ti(DOY, Latitude, bs = c("cc", "cr"), k = c(8, 4)) +
#'     s(Latitude, k = 4) + offset(LogEffort),
#'   family = poisson(), data = sim_carabids,
#'   knots = list(DOY = c(0.5, 365.5))
#' )
#' p <- plot_metric_gradient(
#'   mod, metric = "centroid", x_axis = "Latitude", n_steps = 8, n_sims = 50
#' )
#' @export
plot_metric_gradient <- function(model, metric = "centroid", x_axis = "Latitude",
                                 doy = "DOY", n_steps = 30, n_sims = 1000,
                                 doby_origin = "auto", add = FALSE,
                                 color = "#009E73", ...) {
  if (is.null(model$model)) {
    stop("The supplied model object does not contain training data in 'model$model'.")
  }

  train_df <- model$model
  col_idx <- which(tolower(names(train_df)) == tolower(x_axis))
  if (length(col_idx) == 0) {
    stop(sprintf(
      "Predictor '%s' not found in model$model. Available predictors: %s",
      x_axis, paste(names(train_df)[-1], collapse = ", ")
    ))
  }
  x_col <- names(train_df)[col_idx[1]]
  x_seq <- seq(min(train_df[[x_col]], na.rm = TRUE),
               max(train_df[[x_col]], na.rm = TRUE),
               length.out = n_steps)

  # Baseline reference newdata (all predictors at mean or '...' presets)
  base_nd <- create_newdata(model, doy = doy, ...)
  metric_lc <- tolower(metric)
  need_ref <- metric_lc %in% c("divergence", "noa", "onset", "offset")
  ref_sim <- if (need_ref) {
    simulate_posterior_density(model, base_nd, doy = doy, n_sims = n_sims, warnings = FALSE)
  } else {
    NULL
  }

  metric_res <- data.frame(X_Val = x_seq, Mean = 0, Lower = 0, Upper = 0)

  for (i in seq_len(n_steps)) {
    step_nd <- base_nd
    step_nd[[x_col]] <- x_seq[i]
    attr(step_nd, "mean_baseline") <- attr(base_nd, "mean_baseline")
    sim_i <- simulate_posterior_density(model, step_nd, doy = doy, n_sims = n_sims, warnings = FALSE)

    out <- switch(
      metric_lc,
      "centroid"   = estimate_centroid(sim_i, doby_origin = doby_origin),
      "duration"   = estimate_duration(sim_i),
      "modality"   = estimate_modality(sim_i, doby_origin = doby_origin),
      "peak"       = estimate_peak(sim_i, doby_origin = doby_origin),
      "divergence" = estimate_divergence(sim_i, ref_sim),
      "noa"        = estimate_divergence(sim_i, ref_sim),
      "onset"      = estimate_onset(sim_i, ref_sim, doby_origin = doby_origin),
      "offset"     = estimate_offset(sim_i, ref_sim, doby_origin = doby_origin),
      stop(sprintf("Unknown metric '%s'.", metric))
    )

    metric_res[["Mean"]][i]  <- out$Mean
    metric_res[["Lower"]][i] <- out$Lower_CI
    metric_res[["Upper"]][i] <- out$Upper_CI
  }

  layer_ribbon <- ggplot2::geom_ribbon(
    data = metric_res,
    ggplot2::aes(x = X_Val, ymin = Lower, ymax = Upper),
    fill = color, alpha = 0.2, inherit.aes = FALSE
  )
  layer_line <- ggplot2::geom_line(
    data = metric_res,
    ggplot2::aes(x = X_Val, y = Mean),
    color = color, linewidth = 1, inherit.aes = FALSE
  )

  if (add) {
    return(list(layer_ribbon, layer_line))
  } else {
    y_label <- switch(
      metric_lc,
      "centroid"   = "Centroid (DOY)",
      "duration"   = "Active Duration (Days)",
      "modality"   = "Modality (Peaks)",
      "peak"       = "Peak Timing (DOY)",
      "divergence" = "Divergence (NOA)",
      "noa"        = "Divergence (NOA)",
      "onset"      = "Onset Shift Volume",
      "offset"     = "Offset Shift Volume",
      tools::toTitleCase(metric)
    )
    p <- ggplot2::ggplot() +
      layer_ribbon +
      layer_line +
      ggplot2::theme_classic() +
      ggplot2::labs(x = x_col, y = y_label)
    return(p)
  }
}
