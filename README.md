# gamPheno: Posterior Extraction of Phenological Metrics from GAMs

`gamPheno` is an R package for extracting uncertainty-aware phenological metrics from pre-fitted Generalized Additive Models (`gam` and `bam` objects from the `mgcv` package). 

Rather than extracting deterministic point estimates from a single mean curve, `gamPheno` uses a two-step workflow:
1. **Posterior Density Simulation (`simulate_posterior_density`):** Draws simulated coefficient vectors from the model's posterior covariance matrix, applies the inverse link function, and normalizes each draw into a 365-day probability density function (PDF).
2. **Geometric Metric Extraction (`estimate_*`):** Extracts single-curve metrics (Day-of-Biological-Year centroid, active duration, modality, peak timing) and two-curve comparative metrics (Non-Overlap Area [NOA] divergence, onset shift volume, and offset shift volume) across all simulated curves, propagating spline variance into 95% Credible Intervals.

Because metric extraction operates directly on the simulated posterior density matrices, `gamPheno` supports both **within-species comparisons** (e.g., comparing one species across different climate scenarios, latitudes, or years) and **between-species comparisons** (e.g., quantifying trophic overlap or divergence between two interacting species under the same environmental regime).

---

## Installation

You can install `gamPheno` from GitHub using `remotes`, `pak`, or `devtools`:

```r
# install.packages("remotes")
remotes::install_github("0xswayx0/gamPheno")

# Or via pak / devtools:
# pak::pak("0xswayx0/gamPheno")
# devtools::install_github("0xswayx0/gamPheno")
```

---

## Complete Runnable Workflow (`sim_carabids`)

`gamPheno` (v0.8.0) includes a simulated biweekly pitfall trapping dataset (`sim_carabids`) so you can run the entire workflow and every package function out of the box without external model files:

```r
library(gamPheno)
library(mgcv)

# 1. Load the built-in simulated dataset and fit a cyclic spatiotemporal GAMM
data(sim_carabids)
mod_sim <- gam(
  Count ~ s(DOY, bs = "cc", k = 12) +
    s(Latitude, bs = "cr", k = 4) +
    s(GDD5_YTD, bs = "cr", k = 4) +
    ti(DOY, Latitude, bs = c("cc", "cr"), k = c(8, 4)) +
    ti(DOY, GDD5_YTD, bs = c("cc", "cr"), k = c(8, 4)) +
    s(SubSite, bs = "re") +
    s(fYear, bs = "re") +
    offset(LogEffort),
  family = poisson(),
  data = sim_carabids,
  knots = list(DOY = c(0.5, 365.5))
)

# 2. Create 365-day prediction grids using create_newdata()
# Option A: Programmatic creation using hull presets ("min", "5%", "10%", "25%", "mean", "75%", "90%", "95%", "max")
newdata_warm <- create_newdata(mod_sim, doy = "DOY", Latitude = "5%",  GDD5_YTD = "95%")
newdata_cold <- create_newdata(mod_sim, doy = "DOY", Latitude = "95%", GDD5_YTD = "5%")
newdata_ref  <- create_newdata(mod_sim, doy = "DOY", Latitude = "mean", GDD5_YTD = "mean")

# Option B: Interactive visual scenario builder (launches Shiny app with live curve comparison)
# newdata_custom <- create_newdata(mod_sim, doy = "DOY", interactive = TRUE)

# 3. Step 1: Generate normalized posterior probability density matrices (365 x n_sims)
sim_warm <- simulate_posterior_density(mod_sim, newdata_warm, doy = "DOY", n_sims = 500)
sim_cold <- simulate_posterior_density(mod_sim, newdata_cold, doy = "DOY", n_sims = 500)
sim_ref  <- simulate_posterior_density(mod_sim, newdata_ref,  doy = "DOY", n_sims = 500)

# Inspect attached scenario metadata and path-tracked biological-year origin (DOBY)
sim_warm$doby_origin
sim_warm$scenario

# 4. Step 2a: Extract single-curve phenological metrics (with 95% Credible Intervals)
# doby_origin supports "auto" (default path-tracked steepest descent), "scenario", or a custom DOY (1:365)
estimate_centroid(sim_warm, doby_origin = "auto")
estimate_duration(sim_warm, threshold_fraction = 0.05)
estimate_modality(sim_warm, rel_prominence = 0.50, abs_threshold = 0.25)
estimate_peak(sim_warm)

# 5. Step 2b: Extract comparative metrics between two curves
# (Works for two scenarios of one species OR two different species under the same scenario)
estimate_divergence(sim_warm, sim_ref)
estimate_onset(sim_warm, sim_ref, threshold_fraction = 0.05, doby_origin = "auto")
estimate_offset(sim_warm, sim_ref, threshold_fraction = 0.05, doby_origin = "auto")

# 6. Visualize posterior density curves and continuous reaction norms
plot_phenology_curve(sim_warm, log = TRUE, color = "#D55E00")
plot_metric_gradient(mod_sim, metric = "centroid", x_axis = "Latitude", n_steps = 15, n_sims = 200)
plot_metric_gradient(mod_sim, metric = "onset", x_axis = "GDD5_YTD", n_steps = 15, n_sims = 200, color = "#CC79A7")
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
