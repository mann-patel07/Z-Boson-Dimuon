"""Z mass fit, systematic variants, toy study and J/psi scale check."""
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

import config
import dimuon

ROOT = Path(__file__).resolve().parent
RES = ROOT / "results"


def log(msg):
    with open(ROOT / "RUN_LOG.md", "a") as f:
        f.write(f"- {time.strftime('%F %T')} {msg}\n")


def fit_variant(df, group, desc, pt_min=config.Z_PT_MIN, mrange=config.Z_RANGE, eta_max=config.ETA_MAX,
                bkg="exp", kind="z_fixed_width", binned=False):
    sel = dimuon.z_selection(df, pt_min=pt_min, mrange=mrange, eta_max=eta_max)
    m = df.loc[sel, "mass"].to_numpy()
    model = dimuon.build_model(kind, mrange[0], mrange[1], len(m), bkg=bkg)
    if binned:
        res = dimuon.fit_binned(m, model, config.Z_BIN_WIDTH, do_minos=False)
    else:
        res = dimuon.fit_unbinned(m, model, do_minos=False)
    # Lowest value of the background shape at the fitted parameters, on 2001 points over x in [-1, 1].
    # This is the polynomial in dimuon.linear_pdf and dimuon.second_order_pdf (their normalization
    # is positive, so it does not change the sign). The exponential background is always positive.
    x = np.linspace(-1.0, 1.0, 2001)
    v = res["values"]
    if bkg == "lin":
        bkg_min = float(np.min(1.0 + v["c"] * x))
    elif bkg == "quad":
        bkg_min = float(np.min(1.0 + v["c1"] * x + v["c2"] * (2.0 * x ** 2 - 1.0)))
    else:
        bkg_min = float("nan")
    return {"group": group, "description": desc, "events": int(len(m)), "mass": res["mass"],
            "error": res["mass_err"], "valid": res["valid"], "at_limit": res["at_limit"],
            "at_limit_params": " ".join(res["at_limit_params"]), "rerun": res["rerun"],
            "bkg_min": bkg_min, "bkg_nonnegative": bool(np.isnan(bkg_min) or bkg_min >= 0.0)}


def toy_job(args):
    i, params, n_mean = args
    model = dimuon.build_model("z_fixed_width", config.Z_RANGE[0], config.Z_RANGE[1], int(n_mean))
    data = dimuon.generate_toy_z(model, params, i)
    model = dimuon.build_model("z_fixed_width", config.Z_RANGE[0], config.Z_RANGE[1], len(data))
    res = dimuon.fit_unbinned(data, model, do_minos=False)
    pull = (res["mass"] - params["mass"]) / res["mass_err"] if res["mass_err"] > 0 else float("nan")
    return {"toy": i, "n": int(len(data)), "mass_fit": res["mass"], "mass_err": res["mass_err"],
            "pull": pull, "valid": res["valid"]}


def main():
    t0 = time.time()
    RES.mkdir(exist_ok=True)
    raw = pd.read_csv(ROOT / config.DATA_FILE)
    df = dimuon.add_recomputed_mass(raw)
    share = float(np.mean(np.abs(df["mass"] - df["M"]) > config.M_DIFF_TOLERANCE))
    log(f"Share of rows where the file's M differs from the recomputed mass by more than "
        f"{config.M_DIFF_TOLERANCE} GeV: {share:.6f} ({int(round(share * len(df)))} of {len(df)}).")
    dimuon.cutflow(df).to_csv(RES / "cutflow.csv", index=False)

    # ---- nominal fit
    sel = dimuon.z_selection(df)
    m = df.loc[sel, "mass"].to_numpy()
    a, b = config.Z_RANGE
    model = dimuon.build_model("z_fixed_width", a, b, len(m))
    nom = dimuon.fit_unbinned(m, model, do_minos=True)
    gof = dimuon.goodness_of_fit(m, model, nom, config.Z_BIN_WIDTH)
    if nom["rerun"]:
        log("Nominal Z fit needed the 10 percent rerun.")
    if not nom["valid"]:
        log("WARNING: nominal Z fit is not valid.")

    # ---- systematic variants
    rows = []
    rows.append(fit_variant(df, "nominal", "nominal fit"))
    rows.append(fit_variant(df, "1 background shape", "linear background", bkg="lin"))
    rows.append(fit_variant(df, "1 background shape", "second-order background", bkg="quad"))
    for rg in config.Z_RANGE_VARIANTS:
        rows.append(fit_variant(df, "2 fit range", f"range {rg[0]:g} to {rg[1]:g} GeV", mrange=rg))
    for pt in config.Z_PT_VARIANTS:
        rows.append(fit_variant(df, "3 muon pt threshold", f"pt > {pt:g} GeV", pt_min=pt))
    rows.append(fit_variant(df, "4 signal model", "width free", kind="z_free_width"))
    rows.append(fit_variant(df, "5 likelihood type", "binned Poisson likelihood", binned=True))
    rows.append(fit_variant(df, "6 detector region (not in total)",
                            f"both muons |eta| < {config.DETECTOR_REGION_ETA:g}",
                            eta_max=config.DETECTOR_REGION_ETA))
    sysdf = pd.DataFrame(rows)
    sysdf["shift_from_nominal"] = sysdf["mass"] - nom["mass"]
    for _, r in sysdf.iterrows():
        if not r["valid"]:
            log(f"Variant '{r['description']}' is not a valid fit and is excluded from the systematic totals.")
        if r["valid"] and not r["bkg_nonnegative"]:
            log(f"Variant '{r['description']}': the background shape is negative at the fitted parameters "
                f"(lowest value {r['bkg_min']:.4f}), so the variant is excluded from the systematic totals.")
        if r["at_limit_params"]:
            log(f"Variant '{r['description']}': parameter(s) at a limit: {r['at_limit_params']}. Fit kept as valid by the Minuit flag.")
        if r["group"] != "nominal" and r["rerun"]:
            log(f"Variant '{r['description']}' needed the 10 percent rerun.")
    sysdf.to_csv(RES / "z_systematics.csv", index=False)

    # ---- width from the free-width fit (same nominal selection as the nominal fit)
    wmodel = dimuon.build_model("z_free_width", a, b, len(m))
    wres = dimuon.fit_unbinned(m, wmodel, do_minos=True)
    # fit_unbinned runs Minos for the mass only, so the width Minos errors and the
    # correlation come from the same fit repeated on a Minuit object that is kept
    wmin = dimuon._make_minuit(dimuon._unbinned_nll(wmodel, m), wmodel, 1.1 if wres["rerun"] else 1.0)
    wmin.migrad()
    wmin.hesse()
    wmin.minos("width")
    assert abs(wmin.values["width"] - wres["values"]["width"]) < 1e-6, "repeated width fit differs from fit_unbinned"
    wme = wmin.merrors["width"]
    width_diff = wres["values"]["width"] - config.PDG_GAMMA_Z
    width_json = {
        "width": wres["values"]["width"], "width_hesse_error": wres["errors"]["width"],
        "width_minos_lower": float(wme.lower), "width_minos_upper": float(wme.upper),
        "width_minos_valid": bool(wme.is_valid),
        "valid": wres["valid"], "at_limit": wres["at_limit"], "at_limit_params": wres["at_limit_params"],
        "reference_width": config.PDG_GAMMA_Z,
        "difference_GeV": width_diff, "difference_in_se": width_diff / wres["errors"]["width"],
        "sigma": wres["values"]["sigma"], "sigma_hesse_error": wres["errors"]["sigma"],
        "correlation_width_sigma": float(wmin.covariance.correlation()["width", "sigma"]) if wres["has_covariance"] else None,
    }
    (RES / "z_width_free.json").write_text(json.dumps(width_json, indent=2))
    groups = {}
    for g in ["1 background shape", "2 fit range", "3 muon pt threshold", "4 signal model", "5 likelihood type"]:
        sub = sysdf[sysdf.group == g]
        counted = (sub["valid"] & sub["bkg_nonnegative"]).tolist()
        groups[g] = dimuon.group_systematic(sub["shift_from_nominal"].tolist(), counted)
    syst = dimuon.total_systematic(list(groups.values()))
    stat = nom["mass_err"]
    total = float(np.sqrt(stat ** 2 + syst ** 2))
    diff = nom["mass"] - config.PDG_MZ
    nominal_json = {
        "n_events": int(len(m)), "range": [a, b], "values": nom["values"], "hesse_errors": nom["errors"],
        "minos": nom["minos"], "valid": nom["valid"], "at_limit": nom["at_limit"],
        "at_limit_params": nom["at_limit_params"], "accurate_covariance": nom["accurate_covariance"],
        "rerun": nom["rerun"], "nll": nom["nll"], "gof": gof,
        "mass": nom["mass"], "stat_error": stat, "group_systematics": groups,
        "indicative_systematic": syst, "total_error": total,
        "reference_mass": config.PDG_MZ, "reference_error": config.PDG_MZ_ERR,
        "difference_GeV": diff, "difference_in_stat_se": diff / stat,
        "difference_in_total_se": diff / total,
    }
    (RES / "z_fit_nominal.json").write_text(json.dumps(nominal_json, indent=2))

    # ---- toys
    params = dict(nom["values"])
    n_mean = int(round(params["Ns"] + params["Nb"]))
    jobs = [(i, params, n_mean) for i in range(config.N_TOYS_Z)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        toys = list(ex.map(toy_job, jobs, chunksize=5))
    tdf = pd.DataFrame(toys)
    tdf.to_csv(RES / "z_toys.csv", index=False)
    good = tdf[tdf.valid]
    toy_summary = {"n_toys": int(len(tdf)), "n_invalid": int((~tdf.valid).sum()),
                   "pull_mean": float(good.pull.mean()), "pull_std": float(good.pull.std(ddof=1)),
                   "pull_mean_err": float(good.pull.std(ddof=1) / np.sqrt(len(good)))}
    (RES / "z_toys_summary.json").write_text(json.dumps(toy_summary, indent=2))

    # ---- J/psi
    sel_j = dimuon.jpsi_selection(df)
    mj = df.loc[sel_j, "mass"].to_numpy()
    ja, jb = config.JPSI_RANGE
    out = {"n_events": int(len(mj)), "range": [ja, jb]}
    scale_rows = []
    for key, kind, label in [("crystal_ball", "jpsi_cb", "J/psi (Crystal Ball)"),
                             ("gaussian", "jpsi_gauss", "J/psi (Gaussian variant)")]:
        jm = dimuon.build_model(kind, ja, jb, len(mj))
        jr = dimuon.fit_unbinned(mj, jm, do_minos=True)
        jg = dimuon.goodness_of_fit(mj, jm, jr, config.JPSI_BIN_WIDTH)
        out[key] = {"values": jr["values"], "hesse_errors": jr["errors"], "minos": jr["minos"],
                    "valid": jr["valid"], "at_limit": jr["at_limit"], "at_limit_params": jr["at_limit_params"],
                    "accurate_covariance": jr["accurate_covariance"], "rerun": jr["rerun"],
                    "gof": jg, "mass": jr["mass"], "mass_err": jr["mass_err"]}
        if not jr["valid"]:
            log(f"J/psi fit '{label}' is not valid.")
        scale_rows.append({"resonance": label, "fitted_mass": jr["mass"], "fit_error": jr["mass_err"],
                           "reference": config.PDG_MJPSI, "relative_scale": (jr["mass"] - config.PDG_MJPSI) / config.PDG_MJPSI,
                           "relative_scale_stat_error": jr["mass_err"] / config.PDG_MJPSI, "valid": jr["valid"]})
    (RES / "jpsi_fit.json").write_text(json.dumps(out, indent=2))
    scale_rows.append({"resonance": "Z (nominal)", "fitted_mass": nom["mass"], "fit_error": stat,
                       "reference": config.PDG_MZ, "relative_scale": diff / config.PDG_MZ,
                       "relative_scale_stat_error": stat / config.PDG_MZ, "valid": nom["valid"]})
    pd.DataFrame(scale_rows).to_csv(RES / "scale_check.csv", index=False)
    log(f"run_z_fit.py finished in {(time.time() - t0) / 60:.1f} min.")
    print(json.dumps(nominal_json, indent=1))


if __name__ == "__main__":
    main()
