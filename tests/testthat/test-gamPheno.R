test_that("Circular DOBY protects winter-active centroids, peaks, and modality across Jan 1", {
  # Construct a synthetic 365 x 100 posterior matrix peaking across Dec 31 - Jan 2
  # Half the draws peak at DOY 365 (Dec 31), half peak at DOY 3 (Jan 3); summer dormant (D_min ~ 183)
  n_days <- 365
  n_sims <- 100
  doy <- seq_len(n_days)
  sim_winter <- matrix(0, nrow = n_days, ncol = n_sims)

  for (j in seq_len(n_sims)) {
    mu <- if (j <= 50) 365 else 3
    cdist <- pmin(abs(doy - mu), 365 - abs(doy - mu))
    dens <- exp(-0.5 * (cdist / 12)^2)
    # Zero out mid-summer dormancy around DOY 183
    dens[150:215] <- 0
    sim_winter[, j] <- dens / sum(dens)
  }
  attr(sim_winter, "DOY") <- doy

  cent <- estimate_centroid(sim_winter)
  pk <- estimate_peak(sim_winter)
  mod_cnt <- estimate_modality(sim_winter)

  # True circular centroid and peak lie at ~DOY 1.0 (NOT artificial mid-summer DOY ~184!)
  expect_true(cent$Mean <= 4 || cent$Mean >= 364)
  expect_true(pk$Mean <= 4 || pk$Mean >= 364)
  expect_equal(as.numeric(cent$Lower_CI), 365, tolerance = 1.5)
  expect_equal(as.numeric(cent$Upper_CI), 3, tolerance = 1.5)

  # Single winter peak straddling Dec 31 / Jan 1 should be detected as unimodal (1 peak)
  expect_equal(mod_cnt$Mean, 1)
  expect_equal(mod_cnt$Prob_Multimodal, 0)
})

test_that("Bimodal aestivation duration excludes mid-summer dormancy and detects 2 peaks", {
  doy <- 1:365
  spring <- exp(-0.5 * ((doy - 120) / 14)^2)
  autumn <- exp(-0.5 * ((doy - 280) / 14)^2)
  dens <- spring + autumn
  dens <- dens / sum(dens)
  sim_bi <- matrix(rep(dens, 40), nrow = 365, ncol = 40)
  attr(sim_bi, "DOY") <- doy

  mod_bi <- estimate_modality(sim_bi)
  dur_bi <- estimate_duration(sim_bi, threshold_fraction = 0.05)

  expect_equal(mod_bi$Mean, 2)
  expect_equal(mod_bi$Prob_Multimodal, 1)
  # Each Gaussian window above 5% max (~2.45 * sigma * 2 ~ 68 days) -> ~136 active days,
  # strictly less than naive (Offset - Onset) span (> 220 days)
  expect_lt(dur_bi$Mean, 160)
  expect_gt(dur_bi$Mean, 110)
})

test_that("Sub-daily decimal DOY grids yield consistent calendar-day metrics", {
  doy_daily <- 1:365
  doy_half  <- seq(0.5, 365.0, by = 0.5) # 730 half-day steps

  f_dens <- function(d) {
    v <- exp(-0.5 * ((d - 160) / 20)^2)
    v / sum(v)
  }

  sim_365 <- matrix(rep(f_dens(doy_daily), 20), nrow = 365, ncol = 20)
  attr(sim_365, "DOY") <- doy_daily

  sim_730 <- matrix(rep(f_dens(doy_half), 20), nrow = 730, ncol = 20)
  attr(sim_730, "DOY") <- doy_half

  expect_equal(estimate_duration(sim_365)$Mean, estimate_duration(sim_730)$Mean, tolerance = 1.0)
  expect_equal(estimate_centroid(sim_365)$Mean, estimate_centroid(sim_730)$Mean, tolerance = 0.6)
  expect_equal(estimate_peak(sim_365)$Mean, estimate_peak(sim_730)$Mean, tolerance = 0.6)
})

test_that("Comparative metrics satisfy zero-identity and directional signs", {
  doy <- 1:365
  early <- exp(-0.5 * ((doy - 130) / 18)^2); early <- early / sum(early)
  base  <- exp(-0.5 * ((doy - 155) / 18)^2); base  <- base / sum(base)
  late  <- exp(-0.5 * ((doy - 180) / 18)^2); late  <- late / sum(late)

  sim_early <- matrix(rep(early, 30), nrow = 365, ncol = 30)
  sim_base  <- matrix(rep(base,  30), nrow = 365, ncol = 30)
  sim_late  <- matrix(rep(late,  30), nrow = 365, ncol = 30)

  # Comparing identical matrices yields exact 0 divergence, 0 onset, 0 offset
  expect_equal(estimate_divergence(sim_base, sim_base)$Mean, 0)
  expect_equal(estimate_onset(sim_base, sim_base)$Mean, 0)
  expect_equal(estimate_offset(sim_base, sim_base)$Mean, 0)

  # Advanced curve -> positive onset, negative offset; delayed curve -> negative onset, positive offset
  expect_gt(estimate_onset(sim_early, sim_base)$Mean, 0.2)
  expect_lt(estimate_offset(sim_early, sim_base)$Mean, -0.2)
  expect_lt(estimate_onset(sim_late, sim_base)$Mean, -0.2)
  expect_gt(estimate_offset(sim_late, sim_base)$Mean, 0.2)
  expect_true(estimate_divergence(sim_early, sim_late)$Mean <= 1 &&
              estimate_divergence(sim_early, sim_late)$Mean > 0.4)
})

test_that("Full GAM workflow, extrapolation warning, and plotting helpers succeed on sim_carabids", {
  data(sim_carabids)
  mod <- mgcv::gam(
    Count ~ s(DOY, bs = "cc", k = 10) +
      s(Latitude, k = 4) +
      s(GDD5_YTD, k = 4) +
      ti(DOY, Latitude, bs = c("cc", "cr"), k = c(6, 4)) +
      s(SubSite, bs = "re") +
      offset(LogEffort),
    family = poisson(),
    data = sim_carabids,
    knots = list(DOY = c(0.5, 365.5))
  )

  nd_mean <- create_newdata(mod, doy = "DOY", Latitude = "mean", GDD5_YTD = "mean")
  nd_warm <- create_newdata(mod, doy = "DOY", Latitude = "5%",  GDD5_YTD = "95%")
  expect_equal(nrow(nd_mean), 365)

  # Paired default seed yields exact zero NOA against itself
  s1 <- simulate_posterior_density(mod, nd_mean, n_sims = 50, warnings = TRUE)
  s2 <- simulate_posterior_density(mod, nd_mean, n_sims = 50, warnings = FALSE)
  expect_equal(estimate_divergence(s1, s2)$Mean, 0)
  expect_equal(colSums(s1), rep(1, 50), tolerance = 1e-10)

  # Out-of-range newdata triggers extrapolation warning when warnings = TRUE
  nd_extrap <- nd_mean
  nd_extrap$Latitude <- 65.0
  expect_warning(
    simulate_posterior_density(mod, nd_extrap, n_sims = 20, warnings = TRUE),
    "extrapolation"
  )

  # Plotting functions return valid ggplot objects
  p1 <- plot_phenology_curve(s1, log = TRUE)
  p2 <- plot_metric_gradient(mod, metric = "centroid", x_axis = "Latitude", n_steps = 5, n_sims = 30)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("ref_sim_mat anchors D_min consistently when bimodal aestivation trough deepens", {
  doy <- 1:365
  spring <- exp(-0.5 * ((doy - 140) / 18)^2)
  autumn <- exp(-0.5 * ((doy - 280) / 18)^2)

  # Reference baseline: winter trough (DOY ~40) lower than mid-summer trough (DOY ~210)
  base_d <- spring + autumn + 0.02 * exp(-0.5 * ((doy - 210) / 25)^2) + 0.005
  base_d[35:45] <- 0.001
  base_d <- base_d / sum(base_d)
  sim_ref <- matrix(rep(base_d, 20), nrow = 365, ncol = 20)

  # Perturbed scenario: mid-summer trough dips slightly below winter trough
  pert_d <- spring + autumn + 0.005
  pert_d[205:215] <- 0.0005
  pert_d <- pert_d / sum(pert_d)
  sim_pert <- matrix(rep(pert_d, 20), nrow = 365, ncol = 20)

  cen_anchored <- estimate_centroid(sim_pert, ref_sim_mat = sim_ref)
  pek_anchored <- estimate_peak(sim_pert, ref_sim_mat = sim_ref)
  expect_equal(cen_anchored$D_min, estimate_centroid(sim_ref)$D_min)
  expect_equal(cen_anchored$Mean, 210, tolerance = 5)
  expect_true(!is.null(pek_anchored$Median))
})
