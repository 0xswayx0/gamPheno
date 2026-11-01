#' Simulate Posterior Probability Density Curves from a Fitted GAM/BAM
#'
#' Draws simulated parameter vectors from the posterior covariance matrix of a
#' fitted `mgcv` GAM or BAM object, transforms predictions to the response scale
#' via the model's inverse link function, and normalizes each simulated trajectory
#' into a 365-day probability density function summing to 1. Standardizes the
#' temporal axis to `DOY` (1 to 365) for downstream metric extraction functions
#' and optionally checks whether `newdata` extrapolates beyond the empirical
#' predictor ranges or 2D interaction envelopes of the training data.
#'
#' @param model A fitted `mgcv` GAM or BAM object.
#' @param newdata A `data.frame` containing a 365-day prediction grid for the scenario of interest (such as one produced by [create_newdata()]).
#' @param doy Character, name of the Day-of-Year column in `newdata` and `model` (default = `"DOY"`). If the specified name is not found, the function automatically checks for `"Yday"`, `"doy"`, `"yday"`, or cyclic smooth terms in `model`.
#' @param n_sims Integer, number of posterior simulation draws (default = `1000`).
#' @param seed Optional integer random seed for reproducibility (default = `42`). Using a shared integer seed across scenarios evaluated on the same fitted model pairs identical posterior coefficient draws across scenarios (preserving within-model covariance). Set `seed = NULL` to draw independent random samples without resetting the RNG state.
#' @param warnings Logical, whether to check for univariate and 2D interaction envelope extrapolation in `newdata` relative to the model training data (default = `TRUE`). Set to `FALSE` to skip extrapolation checks during large batch or raster runs.
#' @return A numeric matrix of dimensions `nrow(newdata)` by `n_sims` with row names corresponding to `DOY` and attribute `"DOY"`, where each column is a normalized probability density curve summing to `1`.
#' @examples
#' data(sim_carabids)
#' mod <- mgcv::gam(
#'   Count ~ s(DOY, bs = "cc", k = 10) + s(Latitude, k = 4) + offset(LogEffort),
#'   family = poisson(), data = sim_carabids,
#'   knots = list(DOY = c(0.5, 365.5))
#' )
#' nd <- create_newdata(mod, doy = "DOY", Latitude = "5%")
#' sim_mat <- simulate_posterior_density(mod, nd, n_sims = 100)
#' dim(sim_mat)
#' @export
simulate_posterior_density <- function(model, newdata, doy = "DOY", n_sims = 1000, seed = 42, warnings = TRUE) {
  doy_col <- resolve_doy_column(model, newdata, doy)

  # If the user supplied a column named 'DOY' in newdata when the model expects e.g. 'Yday' (or vice versa),
  # ensure the model's expected column name is present for predict.gam()
  model_doy_col <- resolve_doy_column(model, model$model, doy)
  if (!is.null(doy_col) && !is.null(model_doy_col) && !(model_doy_col %in% names(newdata)) && (doy_col %in% names(newdata))) {
    newdata[[model_doy_col]] <- newdata[[doy_col]]
    doy_col <- model_doy_col
  }

  # Order newdata rows by DOY (1 to 365) so all downstream functions receive a standardized DOY matrix
  if (!is.null(doy_col) && (doy_col %in% names(newdata))) {
    ord <- order(newdata[[doy_col]])
    newdata <- newdata[ord, , drop = FALSE]
  }

  if (isTRUE(warnings) && !is.null(model$model)) {
    check_extrapolation_envelope(model, newdata, doy_col = doy_col)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }
  beta <- stats::coef(model)
  V <- stats::vcov(model, unconditional = TRUE)
  br <- mgcv::rmvn(n_sims, beta, V)

  X <- stats::predict(model, newdata = newdata, type = "lpmatrix")

  # Zero out random effect smooths (bs = "re") to obtain marginal population predictions
  re_cols <- integer(0)
  if (!is.null(model$smooth)) {
    for (sm in model$smooth) {
      if (inherits(sm, "random.effect")) {
        re_cols <- c(re_cols, sm$first.para:sm$last.para)
      }
    }
  }
  re_cols <- unique(c(re_cols, grep("s\\(fYear\\)|s\\(SubSite\\)", colnames(X))))
  if (length(re_cols) > 0) X[, re_cols] <- 0

  sim_link <- X %*% t(br)
  linkinv <- if (!is.null(model$family$linkinv)) model$family$linkinv else exp
  sim_val <- linkinv(sim_link)

  # Normalize each simulated curve to integrate to 1
  sim_mat <- apply(sim_val, 2, function(x) {
    s <- sum(x)
    if (s > 0) x / s else rep(1 / nrow(newdata), nrow(newdata))
  })

  doy_vec <- if (!is.null(doy_col) && (doy_col %in% names(newdata))) newdata[[doy_col]] else seq_len(nrow(newdata))
  rownames(sim_mat) <- as.character(doy_vec)
  attr(sim_mat, "DOY") <- doy_vec
  return(sim_mat)
}

# Internal helper to resolve the Day-of-Year column name flexibly
resolve_doy_column <- function(model, df, doy = "DOY") {
  if (!is.null(df) && !is.null(attr(df, "doy_col")) && (attr(df, "doy_col") %in% names(df))) {
    return(attr(df, "doy_col"))
  }
  if (!is.null(df) && (doy %in% names(df))) {
    return(doy)
  }
  if (!is.null(model$smooth)) {
    for (sm in model$smooth) {
      if (inherits(sm, "cyclic.smooth") && length(sm$term) == 1) {
        if (is.null(df) || (sm$term %in% names(df))) return(sm$term)
      }
    }
  }
  if (!is.null(df)) {
    candidates <- intersect(c("DOY", "Yday", "doy", "yday"), names(df))
    if (length(candidates) > 0) return(candidates[1])
  }
  return(NULL)
}

# Internal helper to detect univariate and 2D convex-hull envelope extrapolation
check_extrapolation_envelope <- function(model, newdata, doy_col = NULL) {
  train_df <- model$model
  common_vars <- intersect(names(newdata), names(train_df))
  if (length(common_vars) == 0) return(invisible(NULL))

  resp_var <- names(train_df)[1]
  all_num <- common_vars[vapply(common_vars, function(v) {
    v != resp_var && !grepl("Effort|Offset", v, ignore.case = TRUE) &&
      is.numeric(train_df[[v]]) && is.numeric(newdata[[v]])
  }, logical(1))]

  if (length(all_num) == 0) return(invisible(NULL))

  # Identify cyclic / Day-of-Year variables (spanning 1:365 in newdata)
  doy_vars <- character(0)
  if (!is.null(doy_col)) doy_vars <- c(doy_vars, doy_col)
  if (!is.null(model$smooth)) {
    for (sm in model$smooth) {
      if (inherits(sm, "cyclic.smooth")) {
        doy_vars <- unique(c(doy_vars, sm$term))
      }
    }
  }
  doy_vars <- unique(c(doy_vars, intersect(all_num, c("Yday", "DOY", "doy", "yday"))))
  env_vars <- setdiff(all_num, doy_vars)

  warn_msgs <- character(0)

  # 1. Univariate range checks on continuous environmental predictors
  for (v in env_vars) {
    obs_rng <- range(train_df[[v]], na.rm = TRUE)
    nd_rng <- range(newdata[[v]], na.rm = TRUE)
    if (nd_rng[1] < obs_rng[1] || nd_rng[2] > obs_rng[2]) {
      warn_msgs <- c(warn_msgs, sprintf(
        "Univariate extrapolation in '%s': newdata range [%.2f, %.2f] exceeds observed training range [%.2f, %.2f].",
        v, nd_rng[1], nd_rng[2], obs_rng[1], obs_rng[2]
      ))
    }
  }

  # Helper to test whether 2D points fall inside the convex hull of training points
  check_hull_2d <- function(obs_x, obs_y, nd_x, nd_y, is_x_cyclic = FALSE, is_y_cyclic = FALSE) {
    ok <- stats::complete.cases(obs_x, obs_y)
    obs_pts <- unique(cbind(obs_x[ok], obs_y[ok]))
    if (nrow(obs_pts) < 3) return(rep(TRUE, length(nd_x)))

    if (is_x_cyclic) {
      obs_pts <- rbind(
        cbind(obs_pts[, 1] - 365, obs_pts[, 2]),
        obs_pts,
        cbind(obs_pts[, 1] + 365, obs_pts[, 2])
      )
    } else if (is_y_cyclic) {
      obs_pts <- rbind(
        cbind(obs_pts[, 1], obs_pts[, 2] - 365),
        obs_pts,
        cbind(obs_pts[, 1], obs_pts[, 2] + 365)
      )
    }

    h_idx <- grDevices::chull(obs_pts[, 1], obs_pts[, 2])
    hull <- obs_pts[c(h_idx, h_idx[1]), , drop = FALSE]
    center <- colMeans(obs_pts[h_idx, , drop = FALSE])
    hull[, 1] <- center[1] + (hull[, 1] - center[1]) * (1 + 1e-5)
    hull[, 2] <- center[2] + (hull[, 2] - center[2]) * (1 + 1e-5)

    mgcv::in.out(hull, cbind(nd_x, nd_y))
  }

  # 2. 2D Convex Hull check for pairs of continuous environmental covariates
  checked_pairs <- character(0)
  if (length(env_vars) >= 2) {
    env_pairs <- utils::combn(env_vars, 2, simplify = FALSE)
    for (pr in env_pairs) {
      v1 <- pr[1]
      v2 <- pr[2]
      key <- paste(sort(c(v1, v2)), collapse = ":")
      checked_pairs <- c(checked_pairs, key)
      inside <- check_hull_2d(train_df[[v1]], train_df[[v2]], newdata[[v1]], newdata[[v2]])
      if (any(!inside)) {
        pct_out <- round(100 * mean(!inside), 1)
        warn_msgs <- c(warn_msgs, sprintf(
          "Covariate envelope extrapolation in (%s, %s): %.1f%% of newdata rows fall outside the 2D convex hull of the training data.",
          v1, v2, pct_out
        ))
      }
    }
  }

  # 3. 2D Smooth / Tensor interaction envelope check (including DOY x EnvDriver)
  if (!is.null(model$smooth)) {
    for (sm in model$smooth) {
      if (length(sm$term) >= 2 && all(sm$term %in% all_num)) {
        t_pairs <- utils::combn(sm$term, 2, simplify = FALSE)
        for (pr in t_pairs) {
          v1 <- pr[1]
          v2 <- pr[2]
          key <- paste(sort(c(v1, v2)), collapse = ":")
          if (key %in% checked_pairs) next
          checked_pairs <- c(checked_pairs, key)

          inside <- check_hull_2d(
            train_df[[v1]], train_df[[v2]],
            newdata[[v1]], newdata[[v2]],
            is_x_cyclic = (v1 %in% doy_vars),
            is_y_cyclic = (v2 %in% doy_vars)
          )
          if (any(!inside)) {
            pct_out <- round(100 * mean(!inside), 1)
            warn_msgs <- c(warn_msgs, sprintf(
              "Interaction envelope extrapolation in smooth(%s, %s): %.1f%% of newdata rows fall outside the observed 2D interaction envelope.",
              v1, v2, pct_out
            ))
          }
        }
      }
    }
  }

  if (length(warn_msgs) > 0) {
    warning(
      paste(c(
        "gamPheno extrapolation warning (set warnings = FALSE to suppress):",
        paste0("  - ", warn_msgs)
      ), collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(NULL)
}
