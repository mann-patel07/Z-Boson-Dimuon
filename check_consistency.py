"""Check that the numbers in README.md equal the values in results/ after the same rounding, and that paths exist."""
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
README = (ROOT / "README.md").read_text()
RES = ROOT / "results"
failures = []
checked = 0

NUM = r"([+-]?\d+(?:\.\d+)?(?:e[+-]?\d+)?)"


def check(label, text_pattern, expected, decimals, scale=1.0):
    """The first number after text_pattern in the README must equal expected, rounded to decimals."""
    global checked
    m = re.search(re.escape(text_pattern) + r"[^0-9+\-]*" + NUM, README)
    checked += 1
    if not m:
        failures.append(f"{label}: pattern '{text_pattern}' not found")
        return
    got = float(m.group(1))
    want = round(expected * scale, decimals)
    if abs(got - want) > 0.5 * 10 ** (-decimals) * 1.0001 and abs(got - want) > 1e-12:
        failures.append(f"{label}: README has {got}, results give {want}")


zf = json.loads((RES / "z_fit_nominal.json").read_text())
ss = json.loads((RES / "scan_summary.json").read_text())
ts = json.loads((RES / "z_toys_summary.json").read_text())
jf = json.loads((RES / "jpsi_fit.json").read_text())
sy = pd.read_csv(RES / "z_systematics.csv")
sc = pd.read_csv(RES / "scale_check.csv").set_index("resonance")
g = zf["gof"]

# headline
check("n events", "I fit the Z boson mass in", zf["n_events"], 0)
check("Z mass", "Z mass from the nominal fit: **", zf["mass"], 4)
check("stat error", "statistical error", zf["stat_error"], 4)
check("syst estimate", "indicative systematic estimate", zf["indicative_systematic"], 4)
check("reference", "Reference value (Particle Data Group):", zf["reference_mass"], 4)
check("difference", "Difference (fit minus reference):", zf["difference_GeV"], 4)
check("diff stat se", "GeV, which is", zf["difference_in_stat_se"], 2)
check("diff total se", "statistical standard errors and", zf["difference_in_total_se"], 2)
check("gof chi2", "chi-square", g["chi2"], 2)
m = re.search(r"chi-square [\d.]+ for (\d+) degrees of freedom, p-value (\d+\.\d+)", README)
checked += 2
if not m:
    failures.append("goodness of fit line not found")
else:
    if int(m.group(1)) != g["ndof"]:
        failures.append(f"ndof: README {m.group(1)}, results {g['ndof']}")
    if abs(float(m.group(2)) - round(g["p_value"], 4)) > 1e-9:
        failures.append(f"gof p: README {m.group(2)}, results {round(g['p_value'], 4)}")

# the agree / disagree sentence must match the numbers
agree = abs(zf["difference_in_total_se"]) <= 2.0
checked += 1
if agree and "do not agree" in README or (not agree and "so the two agree within" in README):
    failures.append("agreement sentence does not match the numbers")

# systematics table rows
for _, r in sy.iterrows():
    row = re.search(r"\| " + re.escape(r["group"]) + r" \| " + re.escape(r["description"]) + r" \| (\d+) \| ([\d.]+) \| ([+-][\d.]+) \|", README)
    checked += 3
    if not row:
        failures.append(f"systematics row missing: {r['description']}")
        continue
    if int(row.group(1)) != int(r["events"]):
        failures.append(f"events for {r['description']}")
    if abs(float(row.group(2)) - round(r["mass"], 4)) > 1e-9:
        failures.append(f"mass for {r['description']}")
    if abs(float(row.group(3)) - round(r["shift_from_nominal"], 4)) > 1e-9:
        failures.append(f"shift for {r['description']}")
for k, v in zf["group_systematics"].items():
    row = re.search(r"\| " + re.escape(k) + r" \| ([\d.]+) \|", README)
    checked += 1
    if not row or abs(float(row.group(1)) - round(v, 4)) > 1e-9:
        failures.append(f"group systematic {k}")
check("total syst", "Indicative systematic estimate:", zf["indicative_systematic"], 4)
check("total error", "Total error (statistical and systematic in quadrature):", zf["total_error"], 4)

# scale check
for name, dec in [("J/psi (Crystal Ball)", 5), ("J/psi (Gaussian variant)", 5), ("Z (nominal)", 5)]:
    row = re.search(r"\| " + re.escape(name) + r" \| ([\d.]+) \| ([\d.]+) \| ([+-][\d.]+) \|", README)
    checked += 3
    if not row:
        failures.append(f"scale row missing: {name}")
        continue
    if abs(float(row.group(1)) - round(sc.loc[name, "fitted_mass"], 5)) > 1e-9:
        failures.append(f"scale fitted mass {name}")
    if abs(float(row.group(2)) - round(sc.loc[name, "reference"], 5)) > 1e-9:
        failures.append(f"scale reference {name}")
    if abs(float(row.group(3)) - round(100 * sc.loc[name, "relative_scale"], 3)) > 1e-9:
        failures.append(f"scale percent {name}")
check("jpsi gof p", "goodness-of-fit p-value", jf["crystal_ball"]["gof"]["p_value"], 4)

# scan
check("scan max z", "Largest local significance:", ss["observed_max_z_local"], 3)
check("scan mass", "at a test mass of", ss["observed_mass"], 1)
check("scan local p", "local p-value", ss["p_local"], 4)
check("scan global p", "Global p-value from", ss["n_toys"], 0)
m = re.search(r"pseudo-experiments: ([\d.]+)\.", README)
checked += 1
if not m or abs(float(m.group(1)) - round(ss["p_global"], 4)) > 1e-9:
    failures.append("global p-value")
check("trials factor", "Trials factor (global p-value divided by local p-value):", ss["trials_factor"], 2)
check("scan gof p", "Background-only goodness-of-fit p-value:", ss["background_gof_p_value"], 4)
check("n test masses", "is moved over", ss["n_test_masses"], 0)
check("pull mean", "(mean", ts["pull_mean"], 3)
check("pull std", "standard deviation", ts["pull_std"], 3)

# paths named in the README
paths = set(re.findall(r"`([\w./-]+\.(?:py|md|csv|json|png|ipynb|txt))`", README))
paths |= set(re.findall(r"\]\(([\w./-]+\.(?:png|md|csv|json))\)", README))
paths |= {"LICENSE"} if "`LICENSE`" in README else set()
paths |= {p for p in re.findall(r"`([\w./-]+/)`", README)}
for p in sorted(paths):
    checked += 1
    if not (ROOT / p).exists():
        failures.append(f"path named in README does not exist: {p}")

if failures:
    print("CONSISTENCY CHECK FAILED")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print(f"consistency check passed: {checked} items checked, README matches results/")
