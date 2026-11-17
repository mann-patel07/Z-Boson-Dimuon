"""Bump scan with the look-elsewhere effect. Saved toys are skipped on restart."""
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


def toy_job(args):
    i, bkg_p, r = args
    counts = dimuon.scan_toy_counts(bkg_p, i)
    out = dimuon.run_scan(counts, r)
    k = int(np.argmax(out["z_local"]))
    return {"toy": i, "max_z_local": float(out["z_local"][k]), "mass_at_max": float(out["mass"][k]),
            "bkg_fit_valid": bool(out["bkg"]["valid"]),
            "n_invalid_scan_fits": int((~out["valid"]).sum())}


def write_summary(zfit, scan):
    """Headline numbers used by the README (Z fit, systematics, scale check, scan)."""
    jp = json.loads((RES / "jpsi_fit.json").read_text())
    toys = json.loads((RES / "z_toys_summary.json").read_text())
    sc = pd.read_csv(RES / "scale_check.csv").set_index("resonance")
    rows = [("z_mass", zfit["mass"]), ("z_stat_error", zfit["stat_error"]),
            ("z_syst_error", zfit["indicative_systematic"]), ("z_total_error", zfit["total_error"]),
            ("z_reference", zfit["reference_mass"]), ("z_difference_GeV", zfit["difference_GeV"]),
            ("z_difference_in_stat_se", zfit["difference_in_stat_se"]),
            ("z_difference_in_total_se", zfit["difference_in_total_se"]),
            ("z_n_events", zfit["n_events"]), ("z_sigma", zfit["values"]["sigma"]),
            ("z_sigma_err", zfit["hesse_errors"]["sigma"]),
            ("z_gof_chi2", zfit["gof"]["chi2"]), ("z_gof_ndof", zfit["gof"]["ndof"]),
            ("z_gof_p_value", zfit["gof"]["p_value"]),
            ("toy_pull_mean", toys["pull_mean"]), ("toy_pull_std", toys["pull_std"]),
            ("toy_n_invalid", toys["n_invalid"]),
            ("jpsi_mass_cb", jp["crystal_ball"]["mass"]), ("jpsi_mass_cb_err", jp["crystal_ball"]["mass_err"]),
            ("jpsi_gof_p_value_cb", jp["crystal_ball"]["gof"]["p_value"]),
            ("jpsi_mass_gauss", jp["gaussian"]["mass"]),
            ("scale_jpsi_cb", sc.loc["J/psi (Crystal Ball)", "relative_scale"]),
            ("scale_jpsi_gauss", sc.loc["J/psi (Gaussian variant)", "relative_scale"]),
            ("scale_z", sc.loc["Z (nominal)", "relative_scale"])]
    for g, v in zfit["group_systematics"].items():
        rows.append((f"syst_group_{g.split(' ')[0]}", v))
    for k in ["observed_max_z_local", "observed_mass", "p_local", "p_global", "trials_factor",
              "background_gof_p_value", "n_toys", "n_toys_ge_observed", "toy_max_z_median",
              "toy_max_z_95pct", "n_test_masses"]:
        rows.append((f"scan_{k}", scan[k]))
    pd.DataFrame(rows, columns=["key", "value"]).to_csv(RES / "summary.csv", index=False)


def main():
    t0 = time.time()
    RES.mkdir(exist_ok=True)
    df = dimuon.add_recomputed_mass(pd.read_csv(ROOT / config.DATA_FILE))
    sel = dimuon.scan_selection(df)
    counts = dimuon.scan_counts(df.loc[sel, "mass"].to_numpy())
    zfit = json.loads((RES / "z_fit_nominal.json").read_text())
    r = zfit["values"]["sigma"] / zfit["values"]["mass"]
    bkg = dimuon.fit_scan_background(counts)
    (RES / "scan_background_fit.json").write_text(json.dumps(
        {"n_events_in_range": int(counts.sum()), "n_bins": int(len(counts)), "parameters": bkg["p"],
         "parameter_errors": bkg["errors"], "nll": bkg["nll"], "valid": bkg["valid"],
         "chi2": bkg["chi2"], "ndof": bkg["ndof"], "p_value": bkg["p_value"], "r": r}, indent=2))
    ts = time.time()
    obs = dimuon.run_scan(counts, r, bkg=bkg)
    one_scan_s = time.time() - ts
    pd.DataFrame({"test_mass": obs["mass"], "S_fit": obs["S"], "q0": obs["q0"],
                  "z_local": obs["z_local"], "p_local": obs["p_local"]}).to_csv(RES / "scan_local.csv", index=False)
    k = int(np.argmax(obs["z_local"]))
    z_obs = float(obs["z_local"][k])
    p_loc = float(obs["p_local"][k])
    print(f"observed scan took {one_scan_s:.1f} s, max Z_local {z_obs:.2f} at {obs['mass'][k]}")

    toy_path = RES / "scan_toys.csv"
    done = set()
    if toy_path.exists():
        done = set(pd.read_csv(toy_path)["toy"].tolist())
    todo = [i for i in range(config.N_TOYS_SCAN) if i not in done]
    jobs = [(i, bkg["p"], r) for i in todo]
    first = not toy_path.exists()
    with ProcessPoolExecutor(max_workers=4) as ex:
        for k2, row in enumerate(ex.map(toy_job, jobs, chunksize=2), 1):
            pd.DataFrame([row]).to_csv(toy_path, mode="a", header=first, index=False)
            first = False
            if k2 % 25 == 0:
                print(f"{k2}/{len(todo)} toys, {time.time() - t0:.0f} s", flush=True)
    toys = pd.read_csv(toy_path).sort_values("toy")
    assert len(toys) == config.N_TOYS_SCAN and toys["toy"].is_unique
    if not toys["bkg_fit_valid"].all():
        log(f"{int((~toys['bkg_fit_valid']).sum())} toy background fits were not valid.")
    p_glob = dimuon.global_p_value(z_obs, toys["max_z_local"].to_numpy())
    summary = {
        "observed_max_z_local": z_obs, "observed_mass": float(obs["mass"][k]),
        "observed_S": float(obs["S"][k]), "p_local": p_loc, "p_global": p_glob,
        "trials_factor": p_glob / p_loc if p_loc > 0 else float("nan"),
        "n_toys": int(len(toys)), "n_toys_ge_observed": int((toys["max_z_local"] >= z_obs).sum()),
        "toy_max_z_median": float(toys["max_z_local"].median()),
        "toy_max_z_95pct": float(toys["max_z_local"].quantile(0.95)),
        "background_gof_p_value": bkg["p_value"], "background_gof_chi2": bkg["chi2"],
        "background_gof_ndof": bkg["ndof"], "gof_alpha": config.GOF_ALPHA,
        "background_model_poor": bool(bkg["p_value"] < config.GOF_ALPHA),
        "relative_resolution_r": r, "n_test_masses": int(len(obs["mass"])),
        "n_events_scan": int(counts.sum()),
        "n_invalid_scan_fits_observed": int((~obs["valid"]).sum()),
        "n_invalid_scan_fits_toys": int(toys["n_invalid_scan_fits"].sum()),
    }
    (RES / "scan_summary.json").write_text(json.dumps(summary, indent=2))
    write_summary(zfit, summary)
    log(f"run_scan.py finished in {(time.time() - t0) / 60:.1f} min. One observed scan took {one_scan_s:.1f} s.")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
