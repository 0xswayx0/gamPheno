# gamPheno: Posterior Extraction of Phenological Metrics from GAMs

`gamPheno` is an R package for extracting uncertainty-aware phenological metrics from pre-fitted Generalized Additive Models (`gam` and `bam` objects from the `mgcv` package). 

Rather than extracting deterministic point estimates from a single mean curve, `gamPheno` uses a two-step workflow:
1. **Posterior Density Simulation (`simulate_posterior_density`):** Draws simulated coefficient vectors from the model's posterior covariance matrix, applies the inverse link function, and normalizes each draw into a 365-day probability density function (PDF).
2. **Geometric Metric Extraction (`estimate_*`):** Extracts single-curve metrics (Day-of-Biological-Year centroid, active duration, modality, peak timing) and two-curve comparative metrics (Non-Overlap Area [NOA] divergence, onset shift volume, and offset shift volume) across all simulated curves, propagating spline variance into 95% Credible Intervals.

Because metric extraction operates directly on the simulated posterior density matrices, `gamPheno` supports both **within-species comparisons** (e.g., comparing one species across different climate scenarios, latitudes, or years) and **between-species comparisons** (e.g., quantifying trophic overlap or divergence between two interacting species under the same environmental regime).

---

## Installation

You can install the development version of `gamPheno` from GitHub using `devtools` or `remotes`:

```r
# install.packages("devtools")
devtools::install_github("0xswayx0/gamPheno")
```

---

## Quick-Start (Runnable Out-of-the-Box)

`gamPheno` (v0.7.0) ships with a lightweight simulated biweekly trapping dataset (`sim_carabids`) so you can test the full workflow immediately:

```r
library(gamPheno)
library(mgcv)

data(sim_carabids)
mod_sim <- gam(
  Count ~ s(DOY, bs = "cc", k = 12) + s(Latitude, k = 4) + s(GDD5_YTD, k = 4) +
    ti(DOY, Latitude, bs = c("cc", "cr"), k = c(8, 4)) + offset(LogEffort),
  family = poisson(), data = sim_carabids, knots = list(DOY = c(0.5, 365.5))
)

sim_warm <- simulate_posterior_density(mod_sim, create_newdata(mod_sim, Latitude = "5%", GDD5_YTD = "95%"))
sim_ref  <- simulate_posterior_density(mod_sim, create_newdata(mod_sim, Latitude = "mean", GDD5_YTD = "mean"))

estimate_centroid(sim_warm)
estimate_duration(sim_warm)
estimate_divergence(sim_warm, sim_ref)
plot_phenology_curve(sim_warm, log = TRUE, color = "#D55E00")
```

---

## Empirical Workflow Example (*Nebria brevicollis*)

```r
library(gamPheno)
library(mgcv)

# 1. Load a pre-fitted mgcv GAM/BAM object
mod <- readRDS("Nebria_brevicollis_Champion_GAMM.rds")

# 2. Create 365-day prediction grids using create_newdata()
# Option A: Programmatic creation using hull presets ("min", "5%", "10%", "25%", "mean", "75%", "90%", "95%", "max")
newdata_warm <- create_newdata(mod, doy = "DOY", Latitude = "5%",  GDD5_YTD = "95%")
newdata_cold <- create_newdata(mod, doy = "DOY", Latitude = "95%", GDD5_YTD = "5%")
newdata_ref  <- create_newdata(mod, doy = "DOY", Latitude = "mean", GDD5_YTD = "mean")

# Option B: Interactive visual scenario builder (launches Shiny app with live curve comparison)
# newdata_custom <- create_newdata(mod, doy = "DOY", interactive = TRUE)

# 3. Step 1: Generate normalized posterior probability density matrices (365 x n_sims)
sim_warm <- simulate_posterior_density(mod, newdata_warm, doy = "DOY", n_sims = 1000)
sim_cold <- simulate_posterior_density(mod, newdata_cold, doy = "DOY", n_sims = 1000)
sim_ref  <- simulate_posterior_density(mod, newdata_ref,  doy = "DOY", n_sims = 1000)

# 4. Step 2a: Extract single-curve phenological metrics (with 95% Credible Intervals)
estimate_centroid(sim_warm)
estimate_duration(sim_warm, threshold_fraction = 0.05)
estimate_modality(sim_warm, rel_prominence = 0.50, abs_threshold = 0.25)
estimate_peak(sim_warm)

# 5. Step 2b: Extract comparative metrics between two curves
# (Works for two scenarios of one species OR two different species under the same scenario)
estimate_divergence(sim_warm, sim_ref)
estimate_onset(sim_warm, sim_ref, threshold_fraction = 0.05)
estimate_offset(sim_warm, sim_ref, threshold_fraction = 0.05)

# 6. Visualize posterior density curves and continuous reaction norms
plot_phenology_curve(sim_warm, log = TRUE, color = "#D55E00")
plot_metric_gradient(mod, metric = "centroid", x_axis = "Latitude")
plot_metric_gradient(mod, metric = "onset", x_axis = "GDD5_YTD", color = "#CC79A7")
```

---

## Data Acknowledgment

The empirical case study and example workflows accompanying this package use long-term carabid beetle pitfall trap data from the **UK Environmental Change Network (ECN)**, made openly available by the UK Centre for Ecology & Hydrology (UKCEH):

* Rennie S, Adamson J, Anderson R, Andrews C, Bater J, Bayfield N, Beaton K, Beaumont D, Benham S, Bowmaker V, Britt C, Brooker R, Brooks D, Brunt J, Common G, Cooper R, Corbett S, Critchley N, Dennis P, Dick J, Dodd B, Dodd N, Donovan N, Easter J, Eaton E, Flexen M, Gardiner A, Hamilton D, Hargreaves P, Hatton-Ellis M, Howe M, Kahl J, Lane M, Langan S, Lloyd D, McElarney Y, McKenna C, McMillan S, Milne F, Milne L, Morecroft M, Murphy M, Nelson A, Nicholson H, Pallett D, Parry D, Pearce I, Pozsgai G, Rose R, Schafer S, Scott T, Sherrin L, Shortall C, Smith R, Smith P, Tait R, Taylor C, Taylor M, Thurlow M, Turner A, Tyson K, Watson H, Whittaker M, Wilkinson M (2015) *UK Environmental Change Network (ECN) ground predator (carabids and spiders) data: 1992–2012*. NERC Environmental Information Data Centre. https://doi.org/10.5285/8385f864-dd41-410f-b248-028f923cb281

---

## AI Assistance Acknowledgment

The development, code structuring, and documentation of the `gamPheno` package were assisted by **Gemini 4 Argon** (Google DeepMind). All statistical formulations, biological thresholds, and analytical outputs were reviewed and validated by the package authors.

---

## License

MIT License.
