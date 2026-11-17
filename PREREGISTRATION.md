# Pre-registration: Z boson mass from CMS open dimuon data

Written before the data file was downloaded.

## Questions

1. Using public CMS collision events with two muons, what Z boson mass do I measure with a proper likelihood fit, how big are the statistical and the systematic uncertainties, and does the result agree with the world average?
2. When I scan a mass region with no known particle for a bump, how large is the most significant fluctuation, and what is its significance after accounting for the number of places I looked?

## Data

CERN Open Data record 5201, "Events with two muons from 2011 (Primary dataset DoubleMu 2011A)", CMS experiment, licence CC0 1.0. Record page: https://opendata.cern.ch/record/5201. File: `Dimuon_DoubleMu.csv`. No other file, record or mirror is used. The record describes 100000 events with two muons, `|eta| < 2.4`, at least one global muon, opposite charge, invariant mass between 0.3 and 300 GeV.

Required columns: `Run`, `Event`, `type1`, `E1`, `px1`, `py1`, `pz1`, `pt1`, `eta1`, `phi1`, `Q1`, `type2`, `E2`, `px2`, `py2`, `pz2`, `pt2`, `eta2`, `phi2`, `Q2`, `M`.

## Reference values (fixed, not looked up again)

- Z mass 91.1879 +/- 0.0020 GeV, Z full width 2.4955 +/- 0.0023 GeV (Particle Data Group, Review of Particle Physics, 2026 edition, Int. J. Mod. Phys. A 41, 2630011 (2026)).
- J/psi mass 3.096900 +/- 0.000006 GeV (Particle Data Group).

## Constants (the same as in `config.py`)

```
ETA_MAX = 2.4
Z_PT_MIN = 20.0
Z_RANGE = (70.0, 110.0)
Z_BIN_WIDTH = 0.5
Z_RANGE_VARIANTS = [(75.0, 105.0), (66.0, 116.0)]
Z_PT_VARIANTS = [25.0, 30.0]
NORM_GRID_POINTS = 4001
N_TOYS_Z = 200
JPSI_RANGE = (2.9, 3.3)
JPSI_BIN_WIDTH = 0.005
CB_N = 3.0
SCAN_RANGE = (15.0, 60.0)
SCAN_BIN_WIDTH = 0.5
SCAN_MASS_FIRST = 16.0
SCAN_MASS_LAST = 59.0
SCAN_MASS_STEP = 0.5
N_TOYS_SCAN = 500
RANDOM_SEED = 42
GOF_ALPHA = 0.01
```

Other fixed numbers from the design: detector-region variant `|eta| < 1.2`; Pearson chi-square bins need an expected count of at least 5; model integrated per bin with the trapezoid rule on 21 points; start values and limits as listed in the fit definitions below (they are also in `config.py`).

## Selections

The invariant mass is always recomputed: `M = sqrt(max(0, (E1+E2)^2 - (px1+px2)^2 - (py1+py2)^2 - (pz1+pz2)^2))`. The share of rows where it differs from the file's `M` by more than 0.01 GeV is logged. The recomputed value is used everywhere.

- Base: `type1 == "G"` and `type2 == "G"`, `Q1 * Q2 == -1`, `|eta1| < 2.4`, `|eta2| < 2.4`.
- Z selection: base, `pt1 > 20` and `pt2 > 20`, `70 <= M < 110`.
- J/psi selection: base, `2.9 <= M < 3.3`.
- Scan selection: base, `15 <= M < 60`.

## Probability densities

All normalized to 1 on the fit range. Voigt: `scipy.special.voigt_profile(m - mass, sigma, gamma)` with `gamma = width / 2`, normalized by the trapezoid rule on 4001 points. Exponential background `exp(lam * m)`, analytic normalization (flat for `lam` near 0). Linear background `1 + c * (m - mid) / half`, `-1 < c < 1`. Second-order background `1 + c1 x + c2 (2 x^2 - 1)` with `x = (m - mid)/half`, `-0.9 < c1 < 0.9`, `-0.4 < c2 < 0.4`. Crystal Ball `scipy.stats.crystalball(beta, 3.0, loc=mean, scale=sigma)` normalized with its cdf, tail on the low-mass side. Gaussian `scipy.stats.norm` normalized with its cdf.

## Fits

Unbinned extended maximum likelihood with iminuit: `NLL = (Ns + Nb) - sum(log(Ns f_s(m_i) + Nb f_b(m_i)))`, `errordef = 0.5`. migrad, then hesse, then minos for the mass. I record validity, limit hits and the covariance accuracy flag. An invalid fit is rerun once from starting values moved by 10 percent. If it is still invalid it is kept, marked invalid, and excluded from the systematic totals.

Nominal Z fit: Z selection. Voigt signal with free `mass` (start 91.0, limits 88 to 94), free `sigma` (start 1.5, limits 0.2 to 5.0), width fixed to 2.4955. Exponential background with free `lam` (start -0.02, limits -0.5 to 0.5). Yields `Ns` (start 0.9 n) and `Nb` (start 0.1 n), both limited to between 0 and 2 n.

Goodness of fit: histogram with bin width 0.5 GeV, model integrated per bin, Pearson chi-square over bins with expected count at least 5, degrees of freedom = bins used minus free parameters, p-value from the chi-square distribution.

## Systematic variants (each changes one thing)

1. Background shape: linear; second-order.
2. Fit range: (75, 105) and (66, 116) (the selection range changes with it).
3. Muon momentum threshold: 25 and 30 GeV.
4. Signal model: width free (start 2.4955, limits 0.5 to 6.0).
5. Likelihood type: binned Poisson likelihood, 0.5 GeV bins, model integrated per bin.
6. Detector region: both muons with `|eta| < 1.2`. Reported in the table but NOT added to the systematic total, because it is a different sample and its difference includes statistical fluctuation.

For groups 1 to 5 the group systematic is the largest `|mass_variant - mass_nominal|` among its valid variants. The total systematic is the square root of the sum of squares of the five group values. I call it an indicative systematic estimate. It covers only my modelling choices, not detector calibration.

## Toy study for the nominal fit

200 toys from the fitted nominal model. Signal and background counts are Poisson with means `Ns` and `Nb`. Voigt sample: `normal(mass, sigma) + gamma * standard_cauchy()`, redrawn until inside the range. Exponential sample by inverse transform on the range. Each toy is fitted with the nominal model. I report the mean and standard deviation of the pull `(mass_fit - mass_true) / hesse_error` and the number of invalid toy fits. Generator: `default_rng(SeedSequence([42, 1, toy_index]))`.

## J/psi scale check

J/psi selection. Crystal Ball signal with free `mean` (start 3.097, limits 3.0 to 3.2), `sigma` (start 0.03, limits 0.005 to 0.1), `beta` (start 1.5, limits 0.3 to 5.0), `n` fixed to 3.0. Exponential background. One variant with a Gaussian signal. Goodness of fit with 0.005 GeV bins. I report the relative scale `(mass_fit - reference) / reference` at the J/psi and at the Z side by side. This is a check of the momentum scale at two masses, not a precision measurement.

## Bump scan with the look-elsewhere effect

Scan selection, 90 bins of 0.5 GeV. Background `B(m) = exp(p0 + p1 x + p2 x^2 + p3 x^3)` with `x = (m - 37.5) / 22.5`, binned Poisson likelihood, iminuit. Background-only Pearson chi-square with 86 degrees of freedom. Signal for a test mass `mt`: Gaussian with mean `mt` and standard deviation `r * mt`, where `r = sigma_Z / mass_Z` from the nominal Z fit. Expected signal per bin = `S` times the Gaussian probability in the bin, `S` between 0 and the total number of events. For test masses 16.0 to 59.0 in steps of 0.5 (87 masses) fit background plus signal with all five parameters free, starting from the background-only result and `S = 1`. `q0 = 2 (NLL_bkg - NLL_s+b)`, set to 0 if negative. `Z_local = sqrt(q0)`. Local p-value `0.5 (1 - erf(Z_local / sqrt 2))`.

Global significance: 500 pseudo-experiments drawn bin by bin from Poisson distributions with means from the background-only fit to the data (`default_rng(SeedSequence([42, 2, toy_index]))`). Each toy gets its own background-only fit and its own full scan. `p_global = (1 + number of toys with largest Z_local >= observed) / (1 + 500)`. The trials factor is `p_global / p_local` of the observed maximum. If the background-only goodness-of-fit p-value is below 0.01, the scan is still reported and I say that the smooth background model does not describe the region well. The scan region contains no known resonance. I do not describe any bump as a particle, a hint or a discovery.

## Primary result

The Z mass from the nominal fit with statistical and systematic uncertainty, and its difference from the reference value.

This dataset has no detector simulation and no momentum-scale calibration applied by me, and the fit model ignores final-state radiation and the Drell-Yan line shape, so a shift from the reference value larger than the statistical error is possible. Whatever shift appears is reported, not tuned away.

These choices were fixed before the data was downloaded and were not changed afterward.
