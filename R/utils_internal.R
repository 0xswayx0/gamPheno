# Internal shared utilities for circular Day-of-Biological-Year (DOBY) transformations
# and temporal grid resolution handling.

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

# Identify biological antimode (D_min) and circular rotation order
compute_doby_origin <- function(base_dens, doy_vec) {
  n_rows <- length(base_dens)
  min_val <- min(base_dens)
  min_idx <- which(base_dens <= min_val + 1e-9)
  theta <- min_idx * 2 * pi / n_rows
  mean_theta <- atan2(mean(sin(theta)), mean(cos(theta)))
  if (mean_theta < 0) mean_theta <- mean_theta + 2 * pi
  idx_min <- round(mean_theta * n_rows / (2 * pi))
  if (idx_min == 0L) idx_min <- n_rows
  D_min <- doy_vec[idx_min]
  shift_order <- order(((seq_len(n_rows) - idx_min) %% n_rows) + 1L)
  list(idx_min = idx_min, D_min = D_min, shift_order = shift_order)
}

# Map unwrapped 1-based DOBY grid coordinate back onto calendar DOY in (0, 365]
wrap_doby_to_doy <- function(doby_idx, D_min, dt = 1) {
  ((D_min + (doby_idx - 1) * dt - 1) %% 365) + 1
}
