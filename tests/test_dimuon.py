"""Unit tests for dimuon.py. Fixed seeds, no network."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402
import dimuon  # noqa: E402


def muon(pt, eta, phi):
    """Massless muon four-vector."""
    return pt * np.cosh(eta), pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)


def make_row(m1, m2, t1="G", t2="G", q1=-1, q2=1):
    """Table row from two (pt, eta, phi) tuples."""
    e1, x1, y1, z1 = muon(*m1)
    e2, x2, y2, z2 = muon(*m2)
    return {"Run": 1, "Event": 1, "type1": t1, "E1": e1, "px1": x1, "py1": y1, "pz1": z1, "pt1": m1[0],
            "eta1": m1[1], "phi1": m1[2], "Q1": q1, "type2": t2, "E2": e2, "px2": x2, "py2": y2,
            "pz2": z2, "pt2": m2[0], "eta2": m2[1], "phi2": m2[2], "Q2": q2, "M": 0.0}


def test_invariant_mass():
    assert dimuon.invariant_mass(45.594, 45.594, 0, 0, 45.594, -45.594, 0, 0) == pytest.approx(91.188)
    # E=10, p=(6,0,8) massless; E=5, p=(0,3,-4) massless; M^2 = 15^2 - 36 - 9 - 16 = 164
    m = dimuon.invariant_mass(10.0, 6.0, 0.0, 8.0, 5.0, 0.0, 3.0, -4.0)
    assert abs(float(m) - np.sqrt(164.0)) < 1e-9


def test_selections_and_cutflow():
    ref1, ref2 = (45.0, 0.3, 0.0), (45.0, -0.3, np.pi)  # roughly 90 GeV
    rows = [
        make_row(ref1, ref2),                                  # 0 passes everything
        make_row(ref1, ref2, t2="T"),                          # 1 tracker muon
        make_row(ref1, ref2, q1=1, q2=1),                      # 2 same sign
        make_row((10.0, 0.3, 0.0), (10.0, -0.3, np.pi)),       # 3 low pt, mass about 20
        make_row((45.0, 2.6, 0.0), (45.0, -0.3, np.pi)),       # 4 high eta
        make_row((70.0, 0.3, 0.0), (70.0, -0.3, np.pi)),       # 5 mass about 140, out of range
    ]
    df = dimuon.add_recomputed_mass(pd.DataFrame(rows))
    sel = dimuon.z_selection(df)
    assert list(sel) == [True, False, False, False, False, False]
    cf = dimuon.cutflow(df)
    z = cf[cf.selection == "Z"].set_index("step")["events"]
    assert list(z) == [6, 5, 4, 3, 2, 1]
    assert bool(dimuon.scan_selection(df)[3]) is True    # mass about 20 is inside (15, 60)
    assert bool(dimuon.jpsi_selection(df)[3]) is False


DENSITY_CASES = [
    ("voigt", lambda m: dimuon.voigt_pdf(m, 91.0, 1.5, 1.25, 70.0, 110.0), (70.0, 110.0)),
    ("voigt2", lambda m: dimuon.voigt_pdf(m, 92.5, 3.0, 1.0, 66.0, 116.0), (66.0, 116.0)),
    ("exp", lambda m: dimuon.exponential_pdf(m, -0.05, 70.0, 110.0), (70.0, 110.0)),
    ("exp2", lambda m: dimuon.exponential_pdf(m, 0.3, 15.0, 60.0), (15.0, 60.0)),
    ("lin", lambda m: dimuon.linear_pdf(m, 0.4, 70.0, 110.0), (70.0, 110.0)),
    ("lin2", lambda m: dimuon.linear_pdf(m, -0.9, 66.0, 116.0), (66.0, 116.0)),
    ("quad", lambda m: dimuon.second_order_pdf(m, 0.3, 0.2, 70.0, 110.0), (70.0, 110.0)),
    ("quad2", lambda m: dimuon.second_order_pdf(m, -0.7, -0.3, 75.0, 105.0), (75.0, 105.0)),
    ("cb", lambda m: dimuon.crystalball_pdf(m, 3.097, 0.03, 1.5, 2.9, 3.3), (2.9, 3.3)),
    ("cb2", lambda m: dimuon.crystalball_pdf(m, 3.09, 0.06, 0.5, 2.9, 3.3), (2.9, 3.3)),
    ("gauss", lambda m: dimuon.gaussian_pdf(m, 3.097, 0.03, 2.9, 3.3), (2.9, 3.3)),
    ("gauss2", lambda m: dimuon.gaussian_pdf(m, 3.05, 0.08, 2.9, 3.3), (2.9, 3.3)),
]


@pytest.mark.parametrize("name,pdf,rng", DENSITY_CASES)
def test_density_integrates_to_one(name, pdf, rng):
    grid = np.linspace(rng[0], rng[1], 20001)
    assert abs(np.trapezoid(pdf(grid), grid) - 1.0) < 1e-4


def test_exponential_flat_limits():
    grid = np.linspace(70, 110, 9)
    for lam in (0.0, 1e-12):
        f = dimuon.exponential_pdf(grid, lam, 70.0, 110.0)
        assert np.allclose(f, 1.0 / 40.0, rtol=1e-9)


def test_voigt_sampling_matches_density():
    rng = np.random.default_rng(10)
    x = dimuon.sample_voigt(rng, 200000, 91.0, 1.5, 0.5 * config.PDG_GAMMA_Z, 70.0, 110.0)
    assert np.all((x >= 70.0) & (x < 110.0))
    edges = np.arange(70.0, 110.0 + 0.25, 0.5)
    counts, _ = np.histogram(x, edges)
    t = np.linspace(0, 1, 21)
    pts = edges[:-1, None] + 0.5 * t[None, :]
    f = dimuon.voigt_pdf(pts.ravel(), 91.0, 1.5, 0.5 * config.PDG_GAMMA_Z, 70.0, 110.0).reshape(-1, 21)
    exp = 200000 * 0.5 * (f.sum(axis=1) - 0.5 * (f[:, 0] + f[:, -1])) / 20
    chi2 = np.sum((counts - exp) ** 2 / exp)
    assert stats.chi2.sf(chi2, len(counts)) > 0.001


def _toy_fit():
    model = dimuon.build_model("z_fixed_width", 70.0, 110.0, 5500)
    truth = {"mass": 91.2, "sigma": 1.8, "lam": -0.03, "Ns": 5000.0, "Nb": 500.0}
    data = dimuon.generate_toy_z(model, truth, 123, fixed_counts=(5000, 500))
    return model, truth, data


def test_fit_recovery_and_yields():
    model, truth, data = _toy_fit()
    res = dimuon.fit_unbinned(data, model)
    assert res["valid"]
    assert abs(res["mass"] - truth["mass"]) < 3 * res["mass_err"]
    assert abs(res["values"]["Ns"] + res["values"]["Nb"] - len(data)) / len(data) < 0.01


def test_binned_and_unbinned_agree():
    model, _, data = _toy_fit()
    ub = dimuon.fit_unbinned(data, model, do_minos=False)
    bn = dimuon.fit_binned(data, model, config.Z_BIN_WIDTH, do_minos=False)
    assert abs(ub["mass"] - bn["mass"]) < 0.05


def test_systematic_combination():
    groups = [0.03, 0.04, 0.0, 0.12, 0.0]
    assert abs(dimuon.total_systematic(groups) - 0.13) < 1e-12
    g = dimuon.group_systematic([0.01, 0.50, -0.02], [True, False, True])
    assert g == pytest.approx(0.02)
    assert dimuon.group_systematic([0.5], [False]) == 0.0


def _bkg_params():
    return [np.log(150.0), -0.8, 0.1, 0.05]


def test_scan_on_background_only_toy():
    counts = dimuon.scan_toy_counts(_bkg_params(), 0)
    masses = np.arange(16.0, 59.5, 2.0)
    out = dimuon.run_scan(counts, 0.0267, masses=masses)
    assert np.all(out["q0"] >= 0)
    assert np.all(np.isfinite(out["z_local"]))


def test_injection_is_found():
    rng = np.random.default_rng(11)
    centres = 0.5 * (dimuon.scan_edges()[1:] + dimuon.scan_edges()[:-1])
    B = dimuon.background_curve(_bkg_params(), centres)
    r = 0.0267
    tpl = dimuon.signal_template(40.0, r)
    S = 6.0 / np.sqrt(np.sum(tpl ** 2 / B))
    counts = rng.poisson(B + S * tpl).astype(float)
    masses = np.arange(16.0, 59.5, 0.5)
    out = dimuon.run_scan(counts, r, masses=masses)
    i = int(np.argmax(out["z_local"]))
    assert abs(out["mass"][i] - 40.0) <= 1.0
    assert out["z_local"][i] > 4.0


def test_global_p_value():
    assert dimuon.global_p_value(3.0, [1.0, 2.0, 3.5, 2.9]) == pytest.approx(0.4)


def test_toy_reproducibility():
    model, truth, _ = _toy_fit()
    a = dimuon.generate_toy_z(model, truth, 5)
    b = dimuon.generate_toy_z(model, truth, 5)
    c = dimuon.generate_toy_z(model, truth, 6)
    assert np.array_equal(a, b)
    assert not (len(a) == len(c) and np.array_equal(a, c))
    s1 = dimuon.scan_toy_counts(_bkg_params(), 3)
    s2 = dimuon.scan_toy_counts(_bkg_params(), 3)
    assert np.array_equal(s1, s2)
