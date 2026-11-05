# Internal shared utilities for circular Day-of-Biological-Year (DOBY) transformations,
# path-tracked circular steepest descent, and temporal grid resolution handling.

utils::globalVariables(c(
  "Yday", "Lower", "Upper", "Density",
  "X_Val", "Mean",
  "DOY", "Baseline", "Target", "Predictor", "Value"
))

# Extract standardized Day-of-Year vector and step width (in days) from sim_mat
get_doy_grid <- function(sim_mat) {
  if (!is.matrix(sim_mat) || nrow(sim_mat) < 3L) {
    stop("'sim_mat' must be a numeric matrix with at least 3 temporal rows.")
  }
  n_rows <- nrow(sim_mat)
  doy_vec <- attr(sim_mat, "DOY")
  if (is.null(doy_vec) || length(doy_vec) != n_rows) {
    rn <- suppressWarnings(as.numeric(rownames(sim_mat)))
    if (!any(is.na(rn)) && length(rn) == n_rows) {
      doy_vec <- rn
    } else {
      doy_vec <- seq_len(n_rows) * (365 / n_rows)
    }
  }
  dt <- 365 / n_rows
  list(doy = as.numeric(doy_vec), n_rows = n_rows, dt = dt)
}

# Circular global minimum index (averaging ties on the circle for flat zero troughs)
circular_global_min_idx <- function(vec) {
  n_rows <- length(vec)
  min_val <- min(vec)
  min_idx <- which(vec <= min_val + 1e-9)
  if (length(min_idx) == 1L) return(min_idx)
  theta <- min_idx * 2 * pi / n_rows
  mean_theta <- atan2(mean(sin(theta)), mean(cos(theta)))
  if (mean_theta < 0) mean_theta <- mean_theta + 2 * pi
  idx_min <- as.integer(round(mean_theta * n_rows / (2 * pi)))
  if (idx_min <= 0L || idx_min > n_rows) idx_min <- n_rows
  idx_min
}

# Sub-day 3-point circular parabolic vertex refinement around discrete trough index
parabolic_trough_doy <- function(vec, idx_min, doy_vec) {
  n_rows <- length(vec)
  dt <- 365 / n_rows
  t1 <- doy_vec[1]
  prv <- if (idx_min == 1L) n_rows else idx_min - 1L
  nxt <- if (idx_min == n_rows) 1L else idx_min + 1L
  denom <- vec[prv] - 2 * vec[idx_min] + vec[nxt]
  if (is.finite(denom) && denom > 1e-12 &&
      vec[idx_min] <= vec[prv] + 1e-12 &&
      vec[idx_min] <= vec[nxt] + 1e-12) {
    delta <- max(-0.5, min(0.5, 0.5 * (vec[prv] - vec[nxt]) / denom))
  } else {
    delta <- 0
  }
  ((doy_vec[idx_min] + delta * dt - t1) %% 365) + t1
}

# Circular steepest-descent walk to local minimum on 1:n_rows
circular_steepest_descent <- function(vec, start_idx) {
  n_rows <- length(vec)
  cur <- as.integer(start_idx)
  for (step in seq_len(n_rows)) {
    prv <- if (cur == 1L) n_rows else cur - 1L
    nxt <- if (cur == n_rows) 1L else cur + 1L
    vp <- vec[prv]
    vc <- vec[cur]
    vn <- vec[nxt]
    if (vc <= vp && vc <= vn) break
    cur <- if (vp < vn) prv else nxt
  }
  cur
}

# Homotopy path-tracked circular steepest descent from all-mean baseline eta_0 to eta_target
track_path_doby_origin <- function(eta_base, eta_target, doy_vec, n_path = 10L) {
  if (max(eta_base) - min(eta_base) < 1e-12 && max(eta_target) - min(eta_target) < 1e-12) {
    return(list(d_star = doy_vec[1], d_base = doy_vec[1]))
  }
  idx0 <- circular_global_min_idx(eta_base)
  d_base <- parabolic_trough_doy(eta_base, idx0, doy_vec)
  cur_idx <- idx0
  for (k in seq_len(n_path)) {
    s <- k / n_path
    eta_s <- (1 - s) * eta_base + s * eta_target
    cur_idx <- circular_steepest_descent(eta_s, cur_idx)
  }
  d_star <- parabolic_trough_doy(eta_target, cur_idx, doy_vec)
  list(d_star = d_star, d_base = d_base)
}

# Convert continuous DOY cut-point d_star into integer row m, fraction alpha, and rotation orders
dstar_to_grid_parts <- function(d_star, grid) {
  n_rows <- grid$n_rows
  dt <- grid$dt
  t1 <- grid$doy[1]
  d_norm <- ((as.numeric(d_star) - t1) %% 365) + t1
  u_star <- 1 + (d_norm - t1) / dt
  m <- as.integer(floor(u_star + 1e-12))
  if (m < 1L) m <- 1L
  if (m > n_rows) {
    m <- 1L
    u_star <- 1
  }
  alpha <- max(0, min(1 - 1e-12, u_star - m))
  idx_nearest <- as.integer(((round(u_star) - 1L) %% n_rows) + 1L)
  shift_order_m <- order(((seq_len(n_rows) - m) %% n_rows) + 1L)
  shift_order_nearest <- order(((seq_len(n_rows) - idx_nearest) %% n_rows) + 1L)
  list(
    D_min = d_norm,
    m = m,
    alpha = alpha,
    idx_min = idx_nearest,
    shift_order = shift_order_m,
    shift_order_nearest = shift_order_nearest
  )
}

# Identify biological antimode (D_min) from a density vector and return grid decomposition
compute_doby_origin <- function(base_dens, doy_vec) {
  grid <- list(doy = as.numeric(doy_vec), n_rows = length(base_dens), dt = 365 / length(base_dens))
  if (max(base_dens) - min(base_dens) < 1e-12) {
    return(dstar_to_grid_parts(doy_vec[1], grid))
  }
  idx_min <- circular_global_min_idx(base_dens)
  if (sum(base_dens <= min(base_dens) + 1e-9) > 1L) {
    d_star <- doy_vec[idx_min]
  } else {
    d_star <- parabolic_trough_doy(base_dens, idx_min, doy_vec)
  }
  dstar_to_grid_parts(d_star, grid)
}

# Order-invariant shortest circular arc of grid indices connecting d1 and d2
shortest_circular_arc_indices <- function(d1, d2, grid) {
  n_rows <- grid$n_rows
  dt <- grid$dt
  t1 <- grid$doy[1]
  u1 <- 1 + (((d1 - t1) %% 365)) / dt
  u2 <- 1 + (((d2 - t1) %% 365)) / dt
  i1 <- as.integer(((round(u1) - 1L) %% n_rows) + 1L)
  i2 <- as.integer(((round(u2) - 1L) %% n_rows) + 1L)
  fwd <- (i2 - i1) %% n_rows
  bwd <- (i1 - i2) %% n_rows
  if (fwd < bwd) {
    ((i1 - 1L + 0:fwd) %% n_rows) + 1L
  } else if (bwd < fwd) {
    ((i2 - 1L + 0:bwd) %% n_rows) + 1L
  } else {
    i_start <- min(i1, i2)
    ((i_start - 1L + 0:fwd) %% n_rows) + 1L
  }
}

# Unified DOBY origin resolver for single-curve and two-curve estimators
resolve_doby_origin <- function(doby_origin, sim_mat, ref_sim_mat = NULL, two_curve = FALSE) {
  grid <- get_doy_grid(sim_mat)

  # Explicit numeric DOY cut-point
  if (is.numeric(doby_origin) && !is.matrix(doby_origin) && length(doby_origin) == 1L) {
    if (!is.finite(doby_origin) || doby_origin <= 0 || doby_origin >= 366) {
      stop("Numeric 'doby_origin' must be a finite calendar Day of Year in (0, 366).")
    }
    return(dstar_to_grid_parts(doby_origin, grid))
  }

  # Matrix passed directly to doby_origin (convenience reference anchoring)
  if (is.matrix(doby_origin)) {
    if (nrow(doby_origin) != grid$n_rows) {
      stop("Matrix supplied to 'doby_origin' must have the same number of temporal rows as 'sim_mat'.")
    }
    d_attr <- attr(doby_origin, "doby_origin")
    if (!is.null(d_attr) && is.finite(d_attr)) {
      return(dstar_to_grid_parts(d_attr, grid))
    }
    return(compute_doby_origin(rowMeans(doby_origin), grid$doy))
  }

  mode_str <- if (is.null(doby_origin)) "auto" else tolower(trimws(as.character(doby_origin[1])))

  if (!two_curve) {
    if (!(mode_str %in% c("auto", "scenario"))) {
      stop("For single-curve estimators, 'doby_origin' must be \"auto\", \"scenario\", a numeric/integer DOY in (0, 366), or a reference simulation matrix.")
    }
    if (mode_str == "scenario") {
      return(compute_doby_origin(rowMeans(sim_mat), grid$doy))
    }
    d_attr <- attr(sim_mat, "doby_origin")
    if (!is.null(d_attr) && is.finite(d_attr)) {
      return(dstar_to_grid_parts(d_attr, grid))
    }
    return(compute_doby_origin(rowMeans(sim_mat), grid$doy))
  }

  # Two-curve resolution (estimate_onset / estimate_offset)
  if (!(mode_str %in% c("auto", "scenario", "joint", "ref"))) {
    stop("For two-curve estimators, 'doby_origin' must be \"auto\", \"scenario\", \"joint\", \"ref\", a numeric/integer DOY in (0, 366), or a reference simulation matrix.")
  }

  d_sim <- attr(sim_mat, "doby_origin")
  d_ref <- attr(ref_sim_mat, "doby_origin")
  id_sim <- attr(sim_mat, "model_id")
  id_ref <- attr(ref_sim_mat, "model_id")

  if (mode_str == "ref") {
    if (!is.null(d_ref) && is.finite(d_ref)) {
      return(dstar_to_grid_parts(d_ref, grid))
    }
    return(compute_doby_origin(rowMeans(ref_sim_mat), grid$doy))
  }

  # Unconstrained scenario-only joint mean curve (ignores model steepest-descent metadata)
  if (mode_str == "scenario") {
    p_joint <- 0.5 * (rowMeans(sim_mat) + rowMeans(ref_sim_mat))
    if (min(p_joint) > 0.05 * max(p_joint) && (max(p_joint) - min(p_joint) > 1e-9)) {
      warning(
        "Joint mean density lacks a deep shared inactivity trough (min > 5% of max); onset/offset direction depends on the shallowest inter-cohort saddle.",
        call. = FALSE
      )
    }
    return(compute_doby_origin(p_joint, grid$doy))
  }

  diff_models <- !is.null(id_sim) && !is.null(id_ref) && !identical(id_sim, id_ref)

  if (mode_str == "auto" && diff_models) {
    d_use <- if (!is.null(d_ref) && is.finite(d_ref)) d_ref else compute_doby_origin(rowMeans(ref_sim_mat), grid$doy)$D_min
    d_other <- if (!is.null(d_sim) && is.finite(d_sim)) d_sim else compute_doby_origin(rowMeans(sim_mat), grid$doy)$D_min
    cdist <- min(abs(d_other - d_use), 365 - abs(d_other - d_use))
    if (cdist > 45) {
      warning(sprintf(
        paste0(
          "Comparing posterior matrices from two distinct models whose tracked dormancy troughs differ by %.1f days ",
          "(DOY %.1f vs. reference DOY %.1f). Anchoring DOBY origin to 'ref_sim_mat' (DOY %.1f); ",
          "override via doby_origin = \"joint\" or a numeric DOY if needed."
        ),
        cdist, d_other, d_use, d_use
      ), call. = FALSE)
    }
    return(dstar_to_grid_parts(d_use, grid))
  }

  # Joint mean density curve (used by "joint" and intraspecific/raw-matrix "auto")
  p_joint <- 0.5 * (rowMeans(sim_mat) + rowMeans(ref_sim_mat))

  if (!diff_models && !is.null(d_sim) && !is.null(d_ref) && is.finite(d_sim) && is.finite(d_ref)) {
    cdist <- min(abs(d_sim - d_ref), 365 - abs(d_sim - d_ref))
    if (cdist < 1e-9) {
      return(dstar_to_grid_parts(0.5 * (d_sim + d_ref), grid))
    }
    arc_idx <- shortest_circular_arc_indices(d_sim, d_ref, grid)
    sub_vals <- p_joint[arc_idx]
    idx_min <- arc_idx[which.min(sub_vals)]
    d_star <- parabolic_trough_doy(p_joint, idx_min, grid$doy)
    return(dstar_to_grid_parts(d_star, grid))
  }

  if (min(p_joint) > 0.05 * max(p_joint) && (max(p_joint) - min(p_joint) > 1e-9)) {
    warning(
      "Joint mean density lacks a deep shared inactivity trough (min > 5% of max); onset/offset direction depends on the shallowest inter-cohort saddle.",
      call. = FALSE
    )
  }
  compute_doby_origin(p_joint, grid$doy)
}

# Map unwrapped calendar coordinate (or 1-based DOBY grid coordinate) back onto calendar DOY in [1, 366)
wrap_unwrapped_doy <- function(unwrapped_doy) {
  ((unwrapped_doy - 1) %% 365) + 1
}

wrap_doby_to_doy <- function(doby_idx, D_min, dt = 1) {
  ((D_min + (doby_idx - 1) * dt - 1) %% 365) + 1
}
