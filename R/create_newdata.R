#' Create a 365-Day Prediction Grid (newdata) for a Fitted GAM/BAM
#'
#' Constructs a complete 365-day prediction `data.frame` compatible with a fitted
#' `mgcv` GAM or BAM object. Automatically populates factor levels for random
#' effects, zeroes out sampling effort offsets, and sets continuous environmental
#' predictors to `"mean"`, hull-bounded `"min"` or `"max"`, empirical quantiles
#' (`"5%"`, `"10%"`, `"25%"`, `"75%"`, `"90%"`, `"95%"`), or custom numeric
#' values (bounded to within `hull_limit` of the empirical training envelope).
#' When `interactive = TRUE`, launches an interactive visual scenario builder
#' showing the live phenological curve (contrasted against a light grey mean
#' baseline curve) and an interactive predictor barplot.
#'
#' @param model A fitted `mgcv` GAM or BAM object.
#' @param doy Character, name of the Day-of-Year column (default = `"DOY"`). Automatically detects `"Yday"`, `"doy"`, or cyclic smooth terms if `"DOY"` is not found in `model`.
#' @param ... Named arguments specifying values for continuous predictors: `"mean"` (default), `"min"`, `"5%"`, `"10%"`, `"25%"`, `"75%"`, `"90%"`, `"95%"`, `"max"`, a single numeric value, or a numeric vector of length 365.
#' @param interactive Logical, whether to open the interactive Shiny scenario builder window (default = FALSE).
#' @param hull_limit Numeric, maximum allowed extrapolation beyond the empirical data envelope as a fraction of the training range (default = 0.10, i.e., 10%).
#' @param assign_name Optional character string; if provided, assigns the created `newdata` object to this variable name in `envir`.
#' @param envir Environment in which to assign `assign_name` (or interactive `"Create"` exports); defaults to `parent.frame()` (i.e., the calling workspace).
#' @return A 365-row `data.frame` ready to pass to [simulate_posterior_density()].
#' @examples
#' data(sim_carabids)
#' mod <- mgcv::gam(
#'   Count ~ s(DOY, bs = "cc", k = 10) + s(Latitude, k = 4) +
#'     s(GDD5_YTD, k = 4) + offset(LogEffort),
#'   family = poisson(), data = sim_carabids,
#'   knots = list(DOY = c(0.5, 365.5))
#' )
#' nd_warm <- create_newdata(mod, doy = "DOY", Latitude = "5%", GDD5_YTD = "95%")
#' head(nd_warm)
#' @export
create_newdata <- function(model, doy = "DOY", ..., interactive = FALSE, hull_limit = 0.10, assign_name = NULL, envir = parent.frame()) {
  if (is.null(model$model)) {
    stop("The supplied model object does not contain training data in 'model$model'.")
  }

  train_df <- model$model
  resp_var <- names(train_df)[1]
  # Exclude response column and special model.frame columns like '(weights)' or '(offset)'
  pred_vars <- setdiff(names(train_df), resp_var)
  pred_vars <- pred_vars[!grepl("^\\(.*\\)$", pred_vars)]

  doy_col <- resolve_doy_column(model, train_df, doy)
  if (is.null(doy_col)) {
    stop("Could not identify the Day-of-Year column in 'model'. Please specify its column name via the 'doy' argument.")
  }

  # Classify predictor columns
  other_vars <- setdiff(pred_vars, doy_col)
  offset_vars <- other_vars[grepl("^offset\\(|Effort|Offset", other_vars, ignore.case = TRUE)]
  factor_vars <- other_vars[vapply(other_vars, function(v) is.factor(train_df[[v]]) || is.character(train_df[[v]]), logical(1))]
  env_vars <- setdiff(other_vars, c(offset_vars, factor_vars))
  env_vars <- env_vars[vapply(env_vars, function(v) is.numeric(train_df[[v]]), logical(1))]

  # Determine which environmental predictors vary across DOY within grouping units
  is_doy_varying <- stats::setNames(rep(FALSE, length(env_vars)), env_vars)
  if (length(env_vars) > 0) {
    grp_key <- if (length(factor_vars) > 0) {
      interaction(train_df[, factor_vars, drop = FALSE], drop = TRUE)
    } else {
      rep(1L, nrow(train_df))
    }
    for (v in env_vars) {
      sds <- tapply(train_df[[v]], grp_key, function(x) if (length(x) > 1) stats::sd(x, na.rm = TRUE) else 0)
      is_doy_varying[v] <- any(sds > 1e-6, na.rm = TRUE)
    }
  }

  # Initialize baseline 365-day data.frame
  nd <- data.frame(row.names = seq_len(365))
  nd[[doy_col]] <- seq_len(365)
  if (doy != doy_col && !(doy %in% names(nd))) {
    nd[[doy]] <- seq_len(365)
  }
  for (v in offset_vars) {
    nd[[v]] <- 0
    inner_v <- sub("^offset\\((.*)\\)$", "\\1", v)
    if (inner_v != v) nd[[inner_v]] <- 0
  }
  for (v in factor_vars) nd[[v]] <- train_df[[v]][1]

  # Helper: normalize preset string (e.g., "5%", "q05", "0.05" -> "5%")
  normalize_preset <- function(val_str) {
    s <- tolower(trimws(as.character(val_str)))
    map <- c(
      "min" = "min", "0%" = "min",
      "5%" = "5%", "q5" = "5%", "q05" = "5%", "0.05" = "5%",
      "10%" = "10%", "q10" = "10%", "0.10" = "10%", "0.1" = "10%",
      "25%" = "25%", "q25" = "25%", "0.25" = "25%",
      "mean" = "mean",
      "75%" = "75%", "q75" = "75%", "0.75" = "75%",
      "90%" = "90%", "q90" = "90%", "0.90" = "90%", "0.9" = "90%",
      "95%" = "95%", "q95" = "95%", "0.95" = "95%",
      "max" = "max", "100%" = "max"
    )
    if (s %in% names(map)) unname(map[[s]]) else s
  }

  # Helper: compute exact 2D convex-hull slice bounds, quantiles, and 10% extrapolation bounds for predictor v
  get_predictor_bounds <- function(v, current_nd) {
    obs_v <- train_df[[v]]
    full_rng <- range(obs_v, na.rm = TRUE)
    span <- diff(full_rng)

    hull_min <- full_rng[1]
    hull_max <- full_rng[2]

    # Intersect with the 2D convex hull of any other scalar environmental predictors
    other_env <- setdiff(env_vars, v)
    for (ov in other_env) {
      if (!is_doy_varying[ov] && (ov %in% names(current_nd))) {
        x0 <- mean(current_nd[[ov]], na.rm = TRUE)
        ok <- stats::complete.cases(train_df[[ov]], obs_v)
        pts <- unique(cbind(train_df[[ov]][ok], obs_v[ok]))
        if (nrow(pts) >= 3) {
          h_idx <- grDevices::chull(pts[, 1], pts[, 2])
          poly <- pts[c(h_idx, h_idx[1]), , drop = FALSE]
          y_cross <- numeric(0)
          for (k in seq_len(nrow(poly) - 1)) {
            x1 <- poly[k, 1]; y1 <- poly[k, 2]
            x2 <- poly[k + 1, 1]; y2 <- poly[k + 1, 2]
            if ((x1 <= x0 && x0 <= x2) || (x2 <= x0 && x0 <= x1)) {
              if (abs(x2 - x1) < 1e-12) {
                y_cross <- c(y_cross, y1, y2)
              } else {
                y_cross <- c(y_cross, y1 + (x0 - x1) * (y2 - y1) / (x2 - x1))
              }
            }
          }
          if (length(y_cross) >= 2) {
            hull_min <- max(hull_min, min(y_cross))
            hull_max <- min(hull_max, max(y_cross))
          }
        }
      }
    }
    if (hull_min > hull_max) {
      hull_min <- full_rng[1]
      hull_max <- full_rng[2]
    }

    # Subset training observations within the hull bounds to compute empirical quantiles
    in_hull_obs <- obs_v[!is.na(obs_v) & obs_v >= hull_min & obs_v <= hull_max]
    if (length(in_hull_obs) < 10) {
      in_hull_obs <- seq(hull_min, hull_max, length.out = 100)
    }
    qs <- as.numeric(stats::quantile(
      in_hull_obs,
      probs = c(0.05, 0.10, 0.25, 0.75, 0.90, 0.95),
      na.rm = TRUE
    ))

    presets <- c(
      "min"  = hull_min,
      "5%"   = qs[1],
      "10%"  = qs[2],
      "25%"  = qs[3],
      "mean" = mean(in_hull_obs, na.rm = TRUE),
      "75%"  = qs[4],
      "90%"  = qs[5],
      "95%"  = qs[6],
      "max"  = hull_max
    )

    list(
      hull_min = hull_min,
      hull_max = hull_max,
      hull_mean = presets[["mean"]],
      presets = presets,
      ext_min = max(full_rng[1] - hull_limit * span, hull_min - hull_limit * span),
      ext_max = min(full_rng[2] + hull_limit * span, hull_max + hull_limit * span),
      full_min = full_rng[1],
      full_max = full_rng[2]
    )
  }

  # Helper: compute daily baseline profile for DOY-varying predictors
  get_doy_profile <- function(v, mode = "mean") {
    m <- normalize_preset(mode)
    d_obs <- train_df[[doy_col]]
    v_obs <- train_df[[v]]
    ok <- stats::complete.cases(d_obs, v_obs)
    d_obs <- d_obs[ok]; v_obs <- v_obs[ok]

    prob_map <- c(
      "min" = 0.0, "5%" = 0.05, "10%" = 0.10, "25%" = 0.25,
      "75%" = 0.75, "90%" = 0.90, "95%" = 0.95, "max" = 1.0
    )

    raw_prof <- numeric(365)
    for (d in 1:365) {
      cdist <- pmin(abs(d_obs - d), 365 - abs(d_obs - d))
      sub_v <- v_obs[cdist <= 15]
      if (length(sub_v) < 3) sub_v <- v_obs
      if (m %in% names(prob_map)) {
        raw_prof[d] <- as.numeric(stats::quantile(sub_v, prob_map[[m]], na.rm = TRUE))
      } else {
        raw_prof[d] <- mean(sub_v, na.rm = TRUE)
      }
    }
    # Apply light circular Gaussian smoothing (sigma = 7 days) to avoid discrete boxcar steps
    days <- 1:365
    prof <- numeric(365)
    for (d in days) {
      kdist <- pmin(abs(days - d), 365 - abs(days - d))
      w <- exp(-0.5 * (kdist / 7)^2)
      prof[d] <- sum(w * raw_prof) / sum(w)
    }
    pmin(pmax(prof, min(v_obs, na.rm = TRUE)), max(v_obs, na.rm = TRUE))
  }

  # Set initial mean values for all environmental predictors
  for (v in env_vars) {
    if (is_doy_varying[v]) {
      nd[[v]] <- get_doy_profile(v, "mean")
    } else {
      nd[[v]] <- mean(train_df[[v]], na.rm = TRUE)
    }
  }

  # Store a copy of the all-mean baseline newdata for reference comparison in interactive mode
  nd_mean_baseline <- nd

  # Apply user-supplied arguments in '...' (also support passing a named list via 'values')
  user_args <- list(...)
  if ("values" %in% names(user_args) && is.list(user_args$values)) {
    user_args <- utils::modifyList(user_args$values, user_args[setdiff(names(user_args), "values")])
  }
  valid_presets <- c("min", "5%", "10%", "25%", "mean", "75%", "90%", "95%", "max")

  for (v in intersect(names(user_args), env_vars)) {
    val <- user_args[[v]]
    bnds <- get_predictor_bounds(v, nd)
    if (is.character(val) && length(val) == 1) {
      m <- normalize_preset(val)
      if (!(m %in% valid_presets)) {
        stop(sprintf(
          "Unknown preset '%s' for predictor '%s'. Choose one of: %s",
          val, v, paste(valid_presets, collapse = ", ")
        ))
      }
      if (is_doy_varying[v]) {
        nd[[v]] <- get_doy_profile(v, m)
      } else if (m == "mean") {
        nd[[v]] <- mean(train_df[[v]], na.rm = TRUE)
      } else {
        nd[[v]] <- unname(bnds$presets[[m]])
      }
    } else if (is.numeric(val)) {
      clamped <- pmin(pmax(val, bnds$ext_min), bnds$ext_max)
      if (any(clamped != val, na.rm = TRUE)) {
        warning(sprintf(
          "Predictor '%s' was clamped to [%.2f, %.2f] to stay within %.0f%% of the empirical data envelope.",
          v, bnds$ext_min, bnds$ext_max, hull_limit * 100
        ), call. = FALSE)
      }
      nd[[v]] <- clamped
    }
  }

  attr(nd, "doy_col") <- doy_col
  attr(nd, "mean_baseline") <- nd_mean_baseline

  # Launch interactive Shiny scenario builder if requested
  if (isTRUE(interactive)) {
    if (!requireNamespace("shiny", quietly = TRUE)) {
      stop("Package 'shiny' is required for interactive = TRUE. Please install it via install.packages('shiny').")
    }
    nd <- launch_newdata_builder(
      model = model, nd = nd, nd_mean_baseline = nd_mean_baseline,
      train_df = train_df, doy_col = doy_col,
      env_vars = env_vars, is_doy_varying = is_doy_varying,
      get_predictor_bounds = get_predictor_bounds,
      get_doy_profile = get_doy_profile,
      hull_limit = hull_limit, assign_name = assign_name,
      target_env = envir
    )
  } else if (!is.null(assign_name) && nzchar(assign_name)) {
    assign(assign_name, nd, envir = envir)
  }

  return(nd)
}

# Internal helper to launch the Shiny interactive newdata builder
launch_newdata_builder <- function(model, nd, nd_mean_baseline, train_df, doy_col, env_vars,
                                   is_doy_varying, get_predictor_bounds, get_doy_profile,
                                   hull_limit, assign_name, target_env) {
  default_name <- if (!is.null(assign_name) && nzchar(assign_name)) assign_name else "newdata_scenario"
  valid_presets <- c("min", "5%", "10%", "25%", "mean", "75%", "90%", "95%", "max")

  # Fast mean curve evaluator for live preview
  compute_mean_curve <- function(cur_nd) {
    X <- stats::predict(model, newdata = cur_nd, type = "lpmatrix")
    re_cols <- integer(0)
    if (!is.null(model$smooth)) {
      for (sm in model$smooth) {
        if (inherits(sm, "random.effect")) re_cols <- c(re_cols, sm$first.para:sm$last.para)
      }
    }
    re_cols <- unique(c(re_cols, grep("s\\(fYear\\)|s\\(SubSite\\)", colnames(X))))
    if (length(re_cols) > 0) X[, re_cols] <- 0
    linkinv <- if (!is.null(model$family$linkinv)) model$family$linkinv else exp
    raw_dens <- as.numeric(linkinv(X %*% stats::coef(model)))
    s <- sum(raw_dens)
    if (s > 0) raw_dens / s else rep(1 / 365, 365)
  }

  # Pre-compute the static all-mean reference curve (light grey background curve)
  baseline_mean_dens <- compute_mean_curve(nd_mean_baseline)

  ui <- shiny::fluidPage(
    shiny::titlePanel("gamPheno: Interactive Scenario Builder (create_newdata)"),
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::plotOutput("pheno_plot", height = "260px")
      )
    ),
    shiny::hr(),
    shiny::fluidRow(
      shiny::column(
        width = 4,
        shiny::selectInput("pred_var", "Select Predictor:", choices = env_vars),
        shiny::radioButtons(
          "preset_mode", "Envelope Preset / Quantile (Single Choice):",
          choices = c(
            "Min (inside hull)"  = "min",
            "5% Quantile"        = "5%",
            "10% Quantile"       = "10%",
            "25% Quantile"       = "25%",
            "Mean (inside hull)" = "mean",
            "75% Quantile"       = "75%",
            "90% Quantile"       = "90%",
            "95% Quantile"       = "95%",
            "Max (inside hull)"  = "max",
            "Custom"             = "custom"
          ),
          selected = "mean"
        ),
        shiny::uiOutput("dynamic_control"),
        shiny::hr(),
        shiny::textInput("obj_name", "Newdata Object Name:", value = default_name),
        shiny::actionButton("create_btn", "Create & Save to Environment", class = "btn-primary")
      ),
      shiny::column(
        width = 8,
        shiny::plotOutput(
          "bar_plot", height = "300px",
          click = "bar_click",
          brush = shiny::brushOpts(id = "bar_brush", direction = "xy", resetOnNew = TRUE)
        ),
        shiny::helpText(
          "For DOY-varying predictors: click or drag a box across any DOY range on the barplot to pull bars up or down (bounded to max 10% extrapolation beyond the empirical hull)."
        )
      )
    )
  )

  server <- function(input, output, session) {
    rv <- shiny::reactiveValues(
      nd = nd,
      modes = stats::setNames(rep("mean", length(env_vars)), env_vars)
    )

    shiny::observeEvent(input$pred_var, {
      v <- input$pred_var
      shiny::updateRadioButtons(session, "preset_mode", selected = rv$modes[[v]])
    })

    shiny::observeEvent(input$preset_mode, {
      v <- input$pred_var
      shiny::req(v)
      m <- input$preset_mode
      rv$modes[[v]] <- m
      if (m %in% valid_presets) {
        bnds <- get_predictor_bounds(v, rv$nd)
        if (is_doy_varying[v]) {
          rv$nd[[v]] <- get_doy_profile(v, m)
        } else if (m == "mean") {
          rv$nd[[v]] <- mean(train_df[[v]], na.rm = TRUE)
        } else {
          rv$nd[[v]] <- unname(bnds$presets[[m]])
        }
      }
    })

    output$dynamic_control <- shiny::renderUI({
      v <- input$pred_var
      shiny::req(v)
      if (!is_doy_varying[v]) {
        bnds <- get_predictor_bounds(v, rv$nd)
        cur_val <- round(mean(rv$nd[[v]]), 3)
        shiny::sliderInput(
          "scalar_val",
          label = sprintf("%s Value (Hull + %.0f%% max):", v, hull_limit * 100),
          min = round(bnds$ext_min, 3),
          max = round(bnds$ext_max, 3),
          value = cur_val,
          step = signif((bnds$ext_max - bnds$ext_min) / 100, 2)
        )
      }
    })

    shiny::observeEvent(input$scalar_val, {
      v <- input$pred_var
      shiny::req(v, !is_doy_varying[v])
      if (abs(mean(rv$nd[[v]]) - input$scalar_val) > 1e-4) {
        rv$nd[[v]] <- input$scalar_val
        rv$modes[[v]] <- "custom"
        shiny::updateRadioButtons(session, "preset_mode", selected = "custom")
      }
    })

    shiny::observeEvent(input$bar_brush, {
      v <- input$pred_var
      shiny::req(v, is_doy_varying[v])
      br <- input$bar_brush
      bnds <- get_predictor_bounds(v, rv$nd)
      d_idx <- which(seq_len(365) >= max(1, round(br$xmin)) & seq_len(365) <= min(365, round(br$xmax)))
      if (length(d_idx) > 0) {
        target_y <- min(max(br$ymax, bnds$ext_min), bnds$ext_max)
        cur_vec <- rv$nd[[v]]
        cur_vec[d_idx] <- target_y
        k <- stats::dnorm(-4:4, 0, 2); k <- k / sum(k)
        pad <- c(cur_vec[362:365], cur_vec, cur_vec[1:4])
        sm_vec <- as.numeric(stats::filter(pad, k, sides = 2))[5:369]
        rv$nd[[v]] <- pmin(pmax(sm_vec, bnds$ext_min), bnds$ext_max)
        rv$modes[[v]] <- "custom"
        shiny::updateRadioButtons(session, "preset_mode", selected = "custom")
      }
    })

    shiny::observeEvent(input$bar_click, {
      v <- input$pred_var
      shiny::req(v)
      clk <- input$bar_click
      bnds <- get_predictor_bounds(v, rv$nd)
      target_y <- min(max(clk$y, bnds$ext_min), bnds$ext_max)
      if (!is_doy_varying[v]) {
        rv$nd[[v]] <- target_y
        rv$modes[[v]] <- "custom"
        shiny::updateRadioButtons(session, "preset_mode", selected = "custom")
      } else {
        d0 <- min(max(1, round(clk$x)), 365)
        cdist <- pmin(abs((1:365) - d0), 365 - abs((1:365) - d0))
        w <- exp(-(cdist^2) / (2 * 15^2))
        new_vec <- rv$nd[[v]] * (1 - w) + target_y * w
        rv$nd[[v]] <- pmin(pmax(new_vec, bnds$ext_min), bnds$ext_max)
        rv$modes[[v]] <- "custom"
        shiny::updateRadioButtons(session, "preset_mode", selected = "custom")
      }
    })

    # Top plot: Target Phenological Curve (green) overlaid on All-Mean Baseline Curve (light grey)
    output$pheno_plot <- shiny::renderPlot({
      dens <- compute_mean_curve(rv$nd)
      df_p <- data.frame(
        DOY = 1:365,
        Target = dens,
        Baseline = baseline_mean_dens
      )
      ggplot2::ggplot(df_p, ggplot2::aes(x = DOY)) +
        ggplot2::geom_area(ggplot2::aes(y = Baseline), fill = "grey80", alpha = 0.45) +
        ggplot2::geom_line(ggplot2::aes(y = Baseline), color = "grey55", linewidth = 0.9, linetype = "dashed") +
        ggplot2::geom_area(ggplot2::aes(y = Target), fill = "#009E73", alpha = 0.25) +
        ggplot2::geom_line(ggplot2::aes(y = Target), color = "#009E73", linewidth = 1.15) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::labs(
          title = "Predicted Mean Phenological Activity Curve",
          subtitle = "Green = Target Scenario  |  Light Grey (dashed) = All-Mean Baseline Prediction",
          x = "Day of Year (DOY)",
          y = "Probability Density"
        )
    })

    # Bottom plot: Predictor Barplot (1 bar if constant, 365 bars if DOY-varying)
    output$bar_plot <- shiny::renderPlot({
      v <- input$pred_var
      shiny::req(v)
      bnds <- get_predictor_bounds(v, rv$nd)
      if (!is_doy_varying[v]) {
        df_b <- data.frame(Predictor = v, Value = mean(rv$nd[[v]]))
        ggplot2::ggplot(df_b, ggplot2::aes(x = Predictor, y = Value)) +
          ggplot2::geom_col(fill = "#CC79A7", width = 0.35) +
          ggplot2::geom_hline(yintercept = c(bnds$hull_min, bnds$hull_max), linetype = "dashed", color = "grey30") +
          ggplot2::geom_hline(yintercept = c(bnds$presets[["5%"]], bnds$presets[["95%"]]), linetype = "dotdash", color = "grey60") +
          ggplot2::geom_hline(yintercept = c(bnds$ext_min, bnds$ext_max), linetype = "dotted", color = "#D55E00") +
          ggplot2::coord_cartesian(ylim = c(min(0, bnds$ext_min), bnds$ext_max * 1.05)) +
          ggplot2::theme_classic(base_size = 13) +
          ggplot2::labs(
            title = sprintf("Predictor: %s (Site/Annual Value -- Click Bar or Use Slider)", v),
            subtitle = "Dashed = Hull [Min, Max] | Dot-Dash = [5%, 95%] | Dotted = 10% Extrapolation Limit",
            x = NULL, y = v
          )
      } else {
        df_b <- data.frame(DOY = 1:365, Value = rv$nd[[v]])
        ggplot2::ggplot(df_b, ggplot2::aes(x = DOY, y = Value)) +
          ggplot2::geom_col(fill = "#56B4E9", width = 1.0) +
          ggplot2::geom_hline(yintercept = c(bnds$hull_min, bnds$hull_max), linetype = "dashed", color = "grey30") +
          ggplot2::geom_hline(yintercept = c(bnds$presets[["5%"]], bnds$presets[["95%"]]), linetype = "dotdash", color = "grey60") +
          ggplot2::geom_hline(yintercept = c(bnds$ext_min, bnds$ext_max), linetype = "dotted", color = "#D55E00") +
          ggplot2::coord_cartesian(ylim = c(min(0, bnds$ext_min), bnds$ext_max * 1.05)) +
          ggplot2::theme_classic(base_size = 13) +
          ggplot2::labs(
            title = sprintf("Daily Predictor Profile: %s (Click or Drag Across Bars to Adjust)", v),
            subtitle = "Dashed = Hull [Min, Max] | Dot-Dash = [5%, 95%] | Dotted = 10% Extrapolation Limit",
            x = "Day of Year (DOY)", y = v
          )
      }
    })

    shiny::observeEvent(input$create_btn, {
      nm <- trimws(input$obj_name)
      if (!nzchar(nm)) nm <- "newdata_scenario"
      assign(nm, rv$nd, envir = target_env)
      arg_strs <- vapply(env_vars, function(v) {
        if (rv$modes[[v]] %in% valid_presets) {
          sprintf('%s = "%s"', v, rv$modes[[v]])
        } else if (!is_doy_varying[v]) {
          sprintf("%s = %.3f", v, mean(rv$nd[[v]]))
        } else {
          sprintf("%s = <custom_365_vector>", v)
        }
      }, character(1))
      message(sprintf("\nSaved '%s' (365 rows) to calling environment.", nm))
      message(sprintf("Reproducible R command:\n%s <- create_newdata(model, doy = \"%s\", %s)\n",
                      nm, doy_col, paste(arg_strs, collapse = ", ")))
      shiny::stopApp(rv$nd)
    })
  }

  shiny::runApp(shiny::shinyApp(ui, server))
}
