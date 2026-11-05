test_that("Circular DOBY protects winter-active centroids, peaks, and modality across Jan 1", {
  # Construct a synthetic 365 x 100 posterior matrix peaking across Dec 31 - Jan 3
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

test_that("Comparative metrics satisfy zero-identity, anti-symmetry, and disjoint edge cases", {
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

  # Exact anti-symmetry under default joint/auto DOBY origin
  expect_equal(
    estimate_onset(sim_early, sim_late)$Mean,
    -estimate_onset(sim_late, sim_early)$Mean,
    tolerance = 1e-12
  )
  expect_equal(
    estimate_offset(sim_early, sim_late)$Mean,
    -estimate_offset(sim_late, sim_early)$Mean,
    tolerance = 1e-12
  )

  # Advanced curve -> positive onset, negative offset; delayed curve -> negative onset, positive offset
  expect_gt(estimate_onset(sim_early, sim_base)$Mean, 0.2)
  expect_lt(estimate_offset(sim_early, sim_base)$Mean, -0.2)
  expect_lt(estimate_onset(sim_late, sim_base)$Mean, -0.2)
  expect_gt(estimate_offset(sim_late, sim_base)$Mean, 0.2)
  expect_true(estimate_divergence(sim_early, sim_late)$Mean <= 1 &&
              estimate_divergence(sim_early, sim_late)$Mean > 0.4)

  # Completely disjoint boxcar cohorts yield Onset = +1.0 and Offset = -1.0 (and anti-symmetric swap)
  box_spring <- numeric(365); box_spring[70:100] <- 1 / 31
  box_autumn <- numeric(365); box_autumn[230:260] <- 1 / 31
  sim_sp <- matrix(rep(box_spring, 10), nrow = 365, ncol = 10)
  sim_au <- matrix(rep(box_autumn, 10), nrow = 365, ncol = 10)
  expect_equal(estimate_onset(sim_sp, sim_au, doby_origin = 1)$Mean, 1.0, tolerance = 1e-12)
  expect_equal(estimate_offset(sim_sp, sim_au, doby_origin = 1)$Mean, -1.0, tolerance = 1e-12)
  expect_equal(estimate_onset(sim_au, sim_sp, doby_origin = 1)$Mean, -1.0, tolerance = 1e-12)
  expect_equal(estimate_offset(sim_au, sim_sp, doby_origin = 1)$Mean, 1.0, tolerance = 1e-12)
})

test_that("Fractional doby_origin yields C0-continuous centroid, onset, and offset across integer and New-Year boundaries", {
  doy <- 1:365
  # Non-zero winter background activity (p(D_min) > 0) where integer rotation would jump by 365 * p(D_min)
  p1 <- exp(-0.5 * ((doy - 170) / 35)^2) + 0.08
  p1 <- p1 / sum(p1)
  p2 <- exp(-0.5 * ((doy - 195) / 35)^2) + 0.05
  p2 <- p2 / sum(p2)
  s1 <- matrix(rep(p1, 10), nrow = 365, ncol = 10)
  s2 <- matrix(rep(p2, 10), nrow = 365, ncol = 10)

  # Across integer day boundary 54.999 -> 55.001
  c_left  <- estimate_centroid(s1, doby_origin = 54.999)$Mean
  c_right <- estimate_centroid(s1, doby_origin = 55.001)$Mean
  expect_lt(abs(c_right - c_left), 0.01)

  on_left  <- estimate_onset(s1, s2, doby_origin = 54.999)$Mean
  on_right <- estimate_onset(s1, s2, doby_origin = 55.001)$Mean
  expect_lt(abs(on_right - on_left), 1e-4)

  of_left  <- estimate_offset(s1, s2, doby_origin = 54.999)$Mean
  of_right <- estimate_offset(s1, s2, doby_origin = 55.001)$Mean
  expect_lt(abs(of_right - of_left), 1e-4)

  # Across New-Year boundary 365.999 -> 1.001
  c_ny1 <- estimate_centroid(s1, doby_origin = 365.999)$Mean
  c_ny2 <- estimate_centroid(s1, doby_origin = 1.001)$Mean
  cdist_ny <- min(abs(c_ny1 - c_ny2), 365 - abs(c_ny1 - c_ny2))
  expect_lt(cdist_ny, 0.01)
})

test_that("Full GAM workflow, gampheno_sim metadata ($), extrapolation warning, and plotting helpers succeed", {
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
  nd_cold <- create_newdata(mod, doy = "DOY", Latitude = "95%", GDD5_YTD = "5%")
  nd_warm <- create_newdata(mod, doy = "DOY", Latitude = "5%",  GDD5_YTD = "95%")
  expect_equal(nrow(nd_mean), 365)

  s1 <- simulate_posterior_density(mod, nd_mean, n_sims = 50, warnings = TRUE)
  s2 <- simulate_posterior_density(mod, nd_mean, n_sims = 50, warnings = FALSE)
  s_cold <- simulate_posterior_density(mod, nd_cold, n_sims = 50, warnings = FALSE)
  s_warm <- simulate_posterior_density(mod, nd_warm, n_sims = 50, warnings = FALSE)

  # gampheno_sim preserves base matrix behavior while exposing $ metadata
  expect_true(is.matrix(s1))
  expect_s3_class(s1, "gampheno_sim")
  expect_true(is.numeric(s1$doby_origin) && s1$doby_origin > 0 && s1$doby_origin < 366)
  expect_equal(length(s1$doy), 365)
  expect_true(is.data.frame(s1$scenario) && nrow(s1$scenario) == 1L)
  # Annual-mean scenario summary distinguishes DOY-varying GDD5_YTD ("5%" vs "95%")
  expect_lt(s_cold$scenario$GDD5_YTD, s_warm$scenario$GDD5_YTD)
  expect_identical(s1$model_id, s_warm$model_id)

  # Paired default seed yields exact zero NOA against itself
  expect_equal(estimate_divergence(s1, s2)$Mean, 0)
  expect_equal(colSums(s1), rep(1, 50), tolerance = 1e-10)

  # Intraspecific onset/offset anti-symmetry across tracked scenarios
  expect_equal(
    estimate_onset(s_warm, s_cold)$Mean,
    -estimate_onset(s_cold, s_warm)$Mean,
    tolerance = 1e-12
  )

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

test_that("doby_origin anchors consistently across bimodal aestivation deepening and interspecific comparisons", {
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

  cen_anchored <- estimate_centroid(sim_pert, doby_origin = sim_ref)
  pek_anchored <- estimate_peak(sim_pert, doby_origin = sim_ref)
  cen_numeric  <- estimate_centroid(sim_pert, doby_origin = 40L)
  expect_equal(cen_anchored$D_min, estimate_centroid(sim_ref)$D_min)
  expect_equal(cen_anchored$Mean, 210, tolerance = 5)
  expect_equal(cen_numeric$D_min, 40)
  expect_true(!is.null(pek_anchored$Median))

  # "auto" uses steepest-descent doby_origin attribute, whereas "scenario" uses unconstrained scenario D_min
  attr(sim_pert, "doby_origin") <- 40
  attr(sim_pert, "model_id") <- "same_species"
  attr(sim_ref,  "doby_origin") <- 42
  attr(sim_ref,  "model_id") <- "same_species"
  expect_equal(estimate_centroid(sim_pert, doby_origin = "auto")$D_min, 40)
  expect_equal(estimate_centroid(sim_pert, doby_origin = "scenario")$D_min, 210, tolerance = 2)

  # Two-curve "auto" searches within shortest arc [40, 42], whereas "scenario" uses unconstrained joint mean D_min
  ons_auto <- estimate_onset(sim_pert, sim_ref, doby_origin = "auto")
  expect_true(ons_auto$D_min >= 39.5 && ons_auto$D_min <= 42.5)
  # Create a second summer-trough scenario so unconstrained joint mean trough lies at ~210 while tracked origins are ~40
  sim_pert2 <- sim_pert
  attr(sim_pert2, "doby_origin") <- 42
  attr(sim_pert2, "model_id") <- "same_species"
  expect_true(estimate_onset(sim_pert, sim_pert2, doby_origin = "auto")$D_min <= 43)
  expect_equal(estimate_onset(sim_pert, sim_pert2, doby_origin = "scenario")$D_min, 210, tolerance = 2)
  expect_equal(estimate_offset(sim_pert, sim_pert2, doby_origin = 55L)$D_min, 55)

  # Two distinct models with > 45-day trough difference anchor to ref_sim_mat and warn
  attr(sim_pert, "doby_origin") <- 210
  attr(sim_pert, "model_id") <- "species_A"
  attr(sim_ref,  "doby_origin") <- 40
  attr(sim_ref,  "model_id") <- "species_B"
  expect_warning(
    ons_inter <- estimate_onset(sim_pert, sim_ref, doby_origin = "auto"),
    "distinct models"
  )
  expect_equal(ons_inter$D_min, 40)
})
