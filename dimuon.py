"""Library for the dimuon analysis: selections, densities, fits, toys, bump scan."""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from iminuit import Minuit
from scipy import stats
from scipy.special import erf, voigt_profile

import config

# ----------------------------------------------------------------------------
# Kinematics and selections
# ----------------------------------------------------------------------------


def invariant_mass(E1, px1, py1, pz1, E2, px2, py2, pz2) -> np.ndarray:
    """Invariant mass of two particles from their four-vectors (GeV)."""
    m2 = ((np.asarray(E1, dtype=float) + E2) ** 2 - (np.asarray(px1, dtype=float) + px2) ** 2
          - (np.asarray(py1, dtype=float) + py2) ** 2 - (np.asarray(pz1, dtype=float) + pz2) ** 2)
    return np.sqrt(np.maximum(0.0, m2))


def add_recomputed_mass(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of df with the recomputed invariant mass in column 'mass'."""
    out = df.copy()
    out["mass"] = invariant_mass(out["E1"], out["px1"], out["py1"], out["pz1"],
                                 out["E2"], out["px2"], out["py2"], out["pz2"])
    return out


def base_mask(df: pd.DataFrame, eta_max: float = config.ETA_MAX) -> np.ndarray:
    """Both muons global, opposite charge, both inside the pseudorapidity limit."""
    return ((df["type1"] == "G") & (df["type2"] == "G") & (df["Q1"] * df["Q2"] == -1)
            & (df["eta1"].abs() < eta_max) & (df["eta2"].abs() < eta_max)).to_numpy()


def range_mask(df: pd.DataFrame, rng: Tuple[float, float]) -> np.ndarray:
    """Recomputed mass in [low, high)."""
    return ((df["mass"] >= rng[0]) & (df["mass"] < rng[1])).to_numpy()


def z_selection(df: pd.DataFrame, pt_min: float = config.Z_PT_MIN,
                mrange: Tuple[float, float] = config.Z_RANGE,
                eta_max: float = config.ETA_MAX) -> np.ndarray:
    """Z selection: base, both pt above pt_min, mass in range."""
    return (base_mask(df, eta_max) & (df["pt1"] > pt_min).to_numpy()
            & (df["pt2"] > pt_min).to_numpy() & range_mask(df, mrange))


def jpsi_selection(df: pd.DataFrame) -> np.ndarray:
    """J/psi selection: base and mass in the J/psi range."""
    return base_mask(df) & range_mask(df, config.JPSI_RANGE)


def scan_selection(df: pd.DataFrame) -> np.ndarray:
    """Scan selection: base and mass in the scan range."""
    return base_mask(df) & range_mask(df, config.SCAN_RANGE)


def cutflow(df: pd.DataFrame) -> pd.DataFrame:
    """Number of events after every cut for the Z, J/psi and scan selections."""
    rows = []

    def add(sel, step, desc, mask):
        rows.append({"selection": sel, "step": step, "cut": desc, "events": int(np.sum(mask))})

    n = len(df)
    allm = np.ones(n, dtype=bool)
    glob = ((df["type1"] == "G") & (df["type2"] == "G")).to_numpy()
    opp = glob & (df["Q1"] * df["Q2"] == -1).to_numpy()
    eta = opp & (df["eta1"].abs() < config.ETA_MAX).to_numpy() & (df["eta2"].abs() < config.ETA_MAX).to_numpy()
    pt = eta & (df["pt1"] > config.Z_PT_MIN).to_numpy() & (df["pt2"] > config.Z_PT_MIN).to_numpy()
    for sel in ["Z", "Jpsi", "scan"]:
        add(sel, 0, "all events in file", allm)
        add(sel, 1, "both muons global", glob)
        add(sel, 2, "opposite charge", opp)
        add(sel, 3, f"|eta| < {config.ETA_MAX} for both", eta)
    add("Z", 4, f"pt > {config.Z_PT_MIN} GeV for both", pt)
    add("Z", 5, f"mass in [{config.Z_RANGE[0]}, {config.Z_RANGE[1]})", z_selection(df))
    add("Jpsi", 4, f"mass in [{config.JPSI_RANGE[0]}, {config.JPSI_RANGE[1]})", jpsi_selection(df))
    add("scan", 4, f"mass in [{config.SCAN_RANGE[0]}, {config.SCAN_RANGE[1]})", scan_selection(df))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Probability densities (each normalized to 1 on [a, b))
# ----------------------------------------------------------------------------


def voigt_pdf(m, mass: float, sigma: float, gamma: float, a: float, b: float,
              n_grid: int = config.NORM_GRID_POINTS) -> np.ndarray:
    """Voigt density (Gaussian sigma convolved with Lorentzian half width gamma), normalized on [a, b)."""
    grid = np.linspace(a, b, n_grid)
    norm = np.trapezoid(voigt_profile(grid - mass, sigma, gamma), grid)
    return voigt_profile(np.asarray(m, dtype=float) - mass, sigma, gamma) / norm


def exponential_pdf(m, lam: float, a: float, b: float) -> np.ndarray:
    """Density proportional to exp(lam * m) on [a, b), flat if lam * (b - a) is tiny."""
    m = np.asarray(m, dtype=float)
    span = b - a
    if abs(lam * span) < 1e-9:
        return np.full_like(m, 1.0 / span)
    return lam * np.exp(lam * (m - a)) / np.expm1(lam * span)


def linear_pdf(m, c: float, a: float, b: float) -> np.ndarray:
    """Density proportional to 1 + c * (m - mid) / half on [a, b)."""
    mid, half = 0.5 * (a + b), 0.5 * (b - a)
    x = (np.asarray(m, dtype=float) - mid) / half
    return (1.0 + c * x) / (b - a)


def second_order_pdf(m, c1: float, c2: float, a: float, b: float) -> np.ndarray:
    """Density proportional to 1 + c1 x + c2 (2 x^2 - 1), x = (m - mid) / half."""
    mid, half = 0.5 * (a + b), 0.5 * (b - a)
    x = (np.asarray(m, dtype=float) - mid) / half
    return (1.0 + c1 * x + c2 * (2.0 * x ** 2 - 1.0)) / (half * (2.0 - 2.0 * c2 / 3.0))


def crystalball_pdf(m, mean: float, sigma: float, beta: float, a: float, b: float,
                    n: float = config.CB_N) -> np.ndarray:
    """Crystal Ball density (tail on the low-mass side), normalized on [a, b)."""
    d = stats.crystalball(beta, n, loc=mean, scale=sigma)
    return d.pdf(np.asarray(m, dtype=float)) / (d.cdf(b) - d.cdf(a))


def gaussian_pdf(m, mean: float, sigma: float, a: float, b: float) -> np.ndarray:
    """Gaussian density normalized on [a, b)."""
    d = stats.norm(loc=mean, scale=sigma)
    return d.pdf(np.asarray(m, dtype=float)) / (d.cdf(b) - d.cdf(a))


# ----------------------------------------------------------------------------
# Models
# ----------------------------------------------------------------------------


def build_model(kind: str, a: float, b: float, n: int, bkg: str = "exp") -> Dict:
    """Describe a signal plus background model: parameter names, starts, limits, densities.

    kind: 'z_fixed_width', 'z_free_width', 'jpsi_cb', 'jpsi_gauss'.
    bkg: 'exp', 'lin', 'quad' (Z models only; J/psi models use 'exp').
    """
    names: List[str] = []
    start: Dict[str, float] = {}
    limits: Dict[str, Tuple[float, float]] = {}
    fixed: Dict[str, bool] = {}
    if kind in ("z_fixed_width", "z_free_width"):
        mass_name = "mass"
        names += ["mass", "sigma"]
        start.update({"mass": config.Z_MASS_START, "sigma": config.Z_SIGMA_START})
        limits.update({"mass": config.Z_MASS_LIMITS, "sigma": config.Z_SIGMA_LIMITS})
        if kind == "z_free_width":
            names.append("width")
            start["width"] = config.PDG_GAMMA_Z
            limits["width"] = config.Z_WIDTH_FREE_LIMITS

        def pdf_s(m, p):
            width = p["width"] if "width" in p else config.PDG_GAMMA_Z
            return voigt_pdf(m, p["mass"], p["sigma"], 0.5 * width, a, b)
        ns_frac, nb_frac = config.Z_SIG_FRAC_START, config.Z_BKG_FRAC_START
    elif kind in ("jpsi_cb", "jpsi_gauss"):
        mass_name = "mean"
        names += ["mean", "sigma"]
        start.update({"mean": config.JPSI_MEAN_START, "sigma": config.JPSI_SIGMA_START})
        limits.update({"mean": config.JPSI_MEAN_LIMITS, "sigma": config.JPSI_SIGMA_LIMITS})
        if kind == "jpsi_cb":
            names.append("beta")
            start["beta"] = config.JPSI_BETA_START
            limits["beta"] = config.JPSI_BETA_LIMITS

            def pdf_s(m, p):
                return crystalball_pdf(m, p["mean"], p["sigma"], p["beta"], a, b)
        else:
            def pdf_s(m, p):
                return gaussian_pdf(m, p["mean"], p["sigma"], a, b)
        ns_frac, nb_frac = config.JPSI_SIG_FRAC_START, config.JPSI_BKG_FRAC_START
    else:
        raise ValueError(kind)
    if kind.startswith("jpsi") or bkg == "exp":
        names.append("lam")
        if kind.startswith("jpsi"):
            start["lam"], limits["lam"] = config.JPSI_LAM_START, config.JPSI_LAM_LIMITS
        else:
            start["lam"], limits["lam"] = config.Z_LAM_START, config.Z_LAM_LIMITS

        def pdf_b(m, p):
            return exponential_pdf(m, p["lam"], a, b)
        bkg_names = ["lam"]
    elif bkg == "lin":
        names.append("c")
        start["c"], limits["c"] = 0.0, config.LINEAR_C_LIMITS

        def pdf_b(m, p):
            return linear_pdf(m, p["c"], a, b)
        bkg_names = ["c"]
    elif bkg == "quad":
        names += ["c1", "c2"]
        start.update({"c1": 0.0, "c2": 0.0})
        limits.update({"c1": config.SECOND_C1_LIMITS, "c2": config.SECOND_C2_LIMITS})

        def pdf_b(m, p):
            return second_order_pdf(m, p["c1"], p["c2"], a, b)
        bkg_names = ["c1", "c2"]
    else:
        raise ValueError(bkg)
    names += ["Ns", "Nb"]
    start.update({"Ns": ns_frac * n, "Nb": nb_frac * n})
    limits.update({"Ns": (0.0, config.Z_YIELD_MAX_FACTOR * n), "Nb": (0.0, config.Z_YIELD_MAX_FACTOR * n)})
    return {"kind": kind, "bkg": bkg, "a": a, "b": b, "n": n, "names": names, "start": start,
            "limits": limits, "fixed": fixed, "pdf_s": pdf_s, "pdf_b": pdf_b,
            "mass_name": mass_name, "bkg_names": bkg_names}


def model_density_components(model: Dict, m, p: Dict[str, float]) -> Tuple[np.ndarray, np.ndarray]:
    """Expected events per GeV from signal and from background at m."""
    return p["Ns"] * model["pdf_s"](m, p), p["Nb"] * model["pdf_b"](m, p)


def bin_expected(model: Dict, p: Dict[str, float], edges: np.ndarray,
                 npts: int = config.GOF_POINTS_PER_BIN) -> Tuple[np.ndarray, np.ndarray]:
    """Expected signal and background counts per bin, trapezoid rule with npts points per bin."""
    nb = len(edges) - 1
    t = np.linspace(0.0, 1.0, npts)
    pts = edges[:-1, None] + (edges[1:] - edges[:-1])[:, None] * t[None, :]
    flat = pts.ravel()
    fs = model["pdf_s"](flat, p).reshape(nb, npts)
    fb = model["pdf_b"](flat, p).reshape(nb, npts)
    w = (edges[1:] - edges[:-1])

    def trap(f):
        return w * (f.sum(axis=1) - 0.5 * (f[:, 0] + f[:, -1])) / (npts - 1)
    return p["Ns"] * trap(fs), p["Nb"] * trap(fb)


# ----------------------------------------------------------------------------
# Fitting
# ----------------------------------------------------------------------------


def _make_minuit(nll: Callable, model: Dict, scale: float) -> Minuit:
    names = model["names"]
    starts = []
    for k in names:
        v = model["start"][k]
        if scale != 1.0:
            lo, hi = model["limits"][k]
            v = v * scale
            span = hi - lo
            v = min(max(v, lo + 1e-3 * span), hi - 1e-3 * span)
        starts.append(v)
    m = Minuit(nll, *starts, name=names)
    m.errordef = 0.5
    for k in names:
        m.limits[k] = model["limits"][k]
        if model["fixed"].get(k, False):
            m.fixed[k] = True
    m.print_level = 0
    return m


def _unbinned_nll(model: Dict, data: np.ndarray) -> Callable:
    names = model["names"]

    def nll(*args):
        p = dict(zip(names, args))
        dens = p["Ns"] * model["pdf_s"](data, p) + p["Nb"] * model["pdf_b"](data, p)
        return (p["Ns"] + p["Nb"]) - float(np.sum(np.log(np.maximum(dens, 1e-300))))
    return nll


def _binned_nll(model: Dict, counts: np.ndarray, edges: np.ndarray) -> Callable:
    names = model["names"]

    def nll(*args):
        p = dict(zip(names, args))
        es, eb = bin_expected(model, p, edges)
        mu = np.maximum(es + eb, 1e-300)
        return float(np.sum(mu - counts * np.log(mu)))
    return nll


def _summarize_minuit(m: Minuit, model: Dict, do_minos: bool) -> Dict:
    names = model["names"]
    fm = m.fmin
    at_limit = bool(fm.has_parameters_at_limit)
    out = {
        "valid": bool(fm.is_valid), "at_limit": at_limit,
        "accurate_covariance": bool(fm.has_accurate_covar), "has_covariance": bool(fm.has_covariance),
        "nll": float(fm.fval), "edm": float(fm.edm), "nfcn": int(fm.nfcn),
        "values": {k: float(m.values[k]) for k in names},
        "errors": {k: float(m.errors[k]) for k in names},
        "fixed": {k: bool(m.fixed[k]) for k in names},
        "at_limit_params": [k for k in names if (not m.fixed[k]) and (
            min(m.values[k] - model["limits"][k][0], model["limits"][k][1] - m.values[k])
            < config.NEAR_LIMIT_FRACTION * (model["limits"][k][1] - model["limits"][k][0]))],
        "minos": None,
    }
    if do_minos:
        mn = model["mass_name"]
        try:
            m.minos(mn)
            me = m.merrors[mn]
            out["minos"] = {"lower": float(me.lower), "upper": float(me.upper), "valid": bool(me.is_valid)}
        except Exception:  # noqa: BLE001
            out["minos"] = {"lower": float("nan"), "upper": float("nan"), "valid": False}
    return out


def _fit(nll: Callable, model: Dict, do_minos: bool) -> Dict:
    rerun = False
    m = _make_minuit(nll, model, 1.0)
    m.migrad()
    m.hesse()
    if not m.fmin.is_valid:
        rerun = True
        m = _make_minuit(nll, model, 1.1)
        m.migrad()
        m.hesse()
    res = _summarize_minuit(m, model, do_minos)
    res["rerun"] = rerun
    res["n_free"] = int(sum(1 for k in model["names"] if not m.fixed[k]))
    return res


def fit_unbinned(data: np.ndarray, model: Dict, do_minos: bool = True) -> Dict:
    """Unbinned extended maximum likelihood fit with migrad, hesse and (optionally) minos."""
    res = _fit(_unbinned_nll(model, np.asarray(data, dtype=float)), model, do_minos)
    res["mass"] = res["values"][model["mass_name"]]
    res["mass_err"] = res["errors"][model["mass_name"]]
    return res


def fit_binned(data: np.ndarray, model: Dict, bin_width: float, do_minos: bool = True) -> Dict:
    """Binned Poisson extended likelihood fit, model integrated per bin."""
    edges = np.arange(model["a"], model["b"] + 0.5 * bin_width, bin_width)
    counts, _ = np.histogram(data, bins=edges)
    res = _fit(_binned_nll(model, counts.astype(float), edges), model, do_minos)
    res["mass"] = res["values"][model["mass_name"]]
    res["mass_err"] = res["errors"][model["mass_name"]]
    return res


def goodness_of_fit(data: np.ndarray, model: Dict, res: Dict, bin_width: float) -> Dict:
    """Pearson chi-square over bins with expected count >= 5, p-value from the chi-square distribution."""
    edges = np.arange(model["a"], model["b"] + 0.5 * bin_width, bin_width)
    counts, _ = np.histogram(data, bins=edges)
    es, eb = bin_expected(model, res["values"], edges)
    exp = es + eb
    use = exp >= config.MIN_EXPECTED_PER_BIN
    chi2 = float(np.sum((counts[use] - exp[use]) ** 2 / exp[use]))
    ndof = int(use.sum() - res["n_free"])
    p = float(stats.chi2.sf(chi2, ndof)) if ndof > 0 else float("nan")
    return {"chi2": chi2, "ndof": ndof, "p_value": p, "n_bins_used": int(use.sum()), "n_bins": int(len(counts))}


# ----------------------------------------------------------------------------
# Toy generation
# ----------------------------------------------------------------------------


def sample_voigt(rng: np.random.Generator, size: int, mass: float, sigma: float, gamma: float,
                 a: float, b: float) -> np.ndarray:
    """Voigt sample: normal(mass, sigma) + gamma * standard_cauchy, redrawn until inside [a, b)."""
    out = np.empty(0)
    while len(out) < size:
        need = size - len(out)
        draw = max(int(need * 1.5) + 10, 10)
        x = rng.normal(mass, sigma, draw) + gamma * rng.standard_cauchy(draw)
        out = np.concatenate([out, x[(x >= a) & (x < b)]])
    return out[:size]


def sample_exponential(rng: np.random.Generator, size: int, lam: float, a: float, b: float) -> np.ndarray:
    """Sample from exp(lam * m) on [a, b) by inverse transform."""
    u = rng.random(size)
    span = b - a
    if abs(lam * span) < 1e-9:
        return a + u * span
    return a + np.log1p(u * np.expm1(lam * span)) / lam


def sample_background(rng: np.random.Generator, size: int, model: Dict, p: Dict[str, float]) -> np.ndarray:
    """Background sample for the model (exponential, or by rejection for linear and second-order)."""
    a, b = model["a"], model["b"]
    if model["bkg"] == "exp":
        return sample_exponential(rng, size, p["lam"], a, b)
    out = np.empty(0)
    fmax = float(np.max(model["pdf_b"](np.linspace(a, b, 2001), p))) * 1.05
    while len(out) < size:
        x = rng.uniform(a, b, max(2 * (size - len(out)), 10))
        u = rng.uniform(0, fmax, len(x))
        out = np.concatenate([out, x[u < model["pdf_b"](x, p)]])
    return out[:size]


def generate_toy_z(model: Dict, p: Dict[str, float], toy_index: int,
                   fixed_counts: Optional[Tuple[int, int]] = None) -> np.ndarray:
    """One toy data set from the fitted Z model, from generator SeedSequence([seed, 1, toy_index]).

    Signal and background counts are Poisson with means Ns and Nb unless fixed_counts is given.
    """
    rng = np.random.default_rng(np.random.SeedSequence([config.RANDOM_SEED, config.TOY_STREAM_Z, toy_index]))
    if fixed_counts is None:
        ns, nb = int(rng.poisson(p["Ns"])), int(rng.poisson(p["Nb"]))
    else:
        ns, nb = fixed_counts
    width = p["width"] if "width" in p else config.PDG_GAMMA_Z
    sig = sample_voigt(rng, ns, p["mass"], p["sigma"], 0.5 * width, model["a"], model["b"])
    bkg = sample_background(rng, nb, model, p)
    return np.concatenate([sig, bkg])


# ----------------------------------------------------------------------------
# Systematics
# ----------------------------------------------------------------------------


def group_systematic(shifts: Sequence[float], valid: Sequence[bool]) -> float:
    """Largest absolute shift among valid variants of one group (0 if none is valid)."""
    vals = [abs(s) for s, v in zip(shifts, valid) if v and np.isfinite(s)]
    return float(max(vals)) if vals else 0.0


def total_systematic(group_values: Sequence[float]) -> float:
    """Square root of the sum of squares of the group values."""
    return float(np.sqrt(np.sum(np.square(np.asarray(group_values, dtype=float)))))


# ----------------------------------------------------------------------------
# Bump scan
# ----------------------------------------------------------------------------


def scan_edges() -> np.ndarray:
    """Bin edges of the scan histogram."""
    lo, hi = config.SCAN_RANGE
    n = int(round((hi - lo) / config.SCAN_BIN_WIDTH))
    return np.linspace(lo, hi, n + 1)


def scan_counts(masses: np.ndarray) -> np.ndarray:
    """Histogram of masses in the scan bins."""
    return np.histogram(masses, bins=scan_edges())[0].astype(float)


def _xcoord(centres: np.ndarray) -> np.ndarray:
    return (centres - config.SCAN_X_CENTRE) / config.SCAN_X_HALF


def background_curve(p: Sequence[float], centres: np.ndarray) -> np.ndarray:
    """Expected count in a bin with the given centre: exp(p0 + p1 x + p2 x^2 + p3 x^3)."""
    x = _xcoord(centres)
    return np.exp(p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3)


def _poisson_nll(mu: np.ndarray, counts: np.ndarray) -> float:
    mu = np.maximum(mu, 1e-300)
    return float(np.sum(mu - counts * np.log(mu)))


def fit_scan_background(counts: np.ndarray) -> Dict:
    """Binned Poisson fit of the four-parameter smooth background."""
    centres = 0.5 * (scan_edges()[1:] + scan_edges()[:-1])
    x = _xcoord(centres)
    start = np.polyfit(x, np.log(counts + 0.5), 3, w=np.sqrt(counts + 1.0))[::-1]

    def nll(p0, p1, p2, p3):
        return _poisson_nll(background_curve([p0, p1, p2, p3], centres), counts)

    m = Minuit(nll, *[float(v) for v in start], name=["p0", "p1", "p2", "p3"])
    m.errordef = 0.5
    m.migrad()
    if not m.fmin.is_valid:
        m.migrad()
    m.hesse()
    p = [float(m.values[k]) for k in ["p0", "p1", "p2", "p3"]]
    mu = background_curve(p, centres)
    use = mu >= config.MIN_EXPECTED_PER_BIN
    chi2 = float(np.sum((counts[use] - mu[use]) ** 2 / mu[use]))
    ndof = int(len(counts) - 4)
    return {"p": p, "errors": [float(m.errors[k]) for k in ["p0", "p1", "p2", "p3"]],
            "nll": float(m.fmin.fval), "valid": bool(m.fmin.is_valid),
            "chi2": float(np.sum((counts - mu) ** 2 / np.maximum(mu, 1e-300))),
            "ndof": ndof,
            "p_value": float(stats.chi2.sf(np.sum((counts - mu) ** 2 / np.maximum(mu, 1e-300)), ndof)),
            "chi2_bins_ge5": chi2, "n_events": float(counts.sum())}


def signal_template(mt: float, r: float) -> np.ndarray:
    """Probability inside each scan bin for a Gaussian with mean mt and standard deviation r * mt."""
    cdf = stats.norm.cdf(scan_edges(), loc=mt, scale=r * mt)
    return np.diff(cdf)


def local_p_value(z: float) -> float:
    """One-sided local p-value from a local significance."""
    return float(0.5 * (1.0 - erf(z / np.sqrt(2.0))))


def run_scan(counts: np.ndarray, r: float, bkg: Optional[Dict] = None,
             masses: Optional[np.ndarray] = None) -> Dict:
    """Scan the test masses: fit background plus signal, return S, q0 and the local significance."""
    if masses is None:
        masses = np.arange(config.SCAN_MASS_FIRST, config.SCAN_MASS_LAST + 0.5 * config.SCAN_MASS_STEP,
                           config.SCAN_MASS_STEP)
    if bkg is None:
        bkg = fit_scan_background(counts)
    centres = 0.5 * (scan_edges()[1:] + scan_edges()[:-1])
    ntot = float(counts.sum())
    S_fit, q0s, valid = [], [], []
    for mt in masses:
        tpl = signal_template(float(mt), r)

        def nll(p0, p1, p2, p3, S):
            mu = background_curve([p0, p1, p2, p3], centres) + S * tpl
            return _poisson_nll(mu, counts)

        m = Minuit(nll, *bkg["p"], 1.0, name=["p0", "p1", "p2", "p3", "S"])
        m.errordef = 0.5
        m.limits["S"] = (0.0, ntot)
        m.migrad()
        if not m.fmin.is_valid:
            m.migrad()
        q0 = max(0.0, 2.0 * (bkg["nll"] - float(m.fmin.fval)))
        S_fit.append(float(m.values["S"]))
        valid.append(bool(m.fmin.is_valid))
        q0s.append(q0)
    q0s = np.array(q0s)
    z = np.sqrt(q0s)
    return {"mass": np.asarray(masses, dtype=float), "S": np.array(S_fit), "valid": np.array(valid), "q0": q0s, "z_local": z,
            "p_local": np.array([local_p_value(v) for v in z]), "bkg": bkg}


def scan_toy_counts(bkg_p: Sequence[float], toy_index: int) -> np.ndarray:
    """One background-only pseudo-experiment, Poisson bin by bin, generator SeedSequence([seed, 2, toy_index])."""
    rng = np.random.default_rng(np.random.SeedSequence([config.RANDOM_SEED, config.TOY_STREAM_SCAN, toy_index]))
    centres = 0.5 * (scan_edges()[1:] + scan_edges()[:-1])
    return rng.poisson(background_curve(bkg_p, centres)).astype(float)


def global_p_value(observed: float, toy_maxima: Sequence[float]) -> float:
    """(1 + number of toys with maximum >= observed) / (1 + number of toys)."""
    toys = np.asarray(toy_maxima, dtype=float)
    return float((1 + np.sum(toys >= observed)) / (1 + len(toys)))
