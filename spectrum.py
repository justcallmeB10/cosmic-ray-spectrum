
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

BASE = Path(r"C:\Users\ssat-\Desktop\GitHub\github.2\cosmic-ray-spectrum")
DATA = BASE / "data"

# --- load data ---
zs = zipfile.ZipFile(DATA / "summary.zip")
za = zipfile.ZipFile(DATA / "auxiliary.zip")
ev = pd.read_csv(zs.open("summary/dataSummarySD1500.csv"))
ex = pd.read_csv(za.open("auxiliary/sd1500exposure.csv"))

# --- selection (Auger notebook: zenith < 60 deg, E >= 2.5 EeV) ---
Emin = 2.5
sel = ev["sd_energy"].notna() & (ev["sd_energy"] >= Emin) & (ev["sd_theta"] < 60)
E = ev.loc[sel, "sd_energy"].values
exposure = ex["sd_exposure"].iloc[-1]          # km^2 sr yr (10% sample)

# --- binning: log10(E/eV) starting at 18.4, width 0.1 ---
edges = 10**np.arange(18.4, 20.3, 0.1) / 1e18  # in EeV
N, _ = np.histogram(E, bins=edges)
dE = np.diff(edges)
Ecen = np.sqrt(edges[:-1] * edges[1:])

# --- flux and Poisson errors ---
J = N / (exposure * dE)
dJ = np.sqrt(N) / (exposure * dE)
mask = N > 0

# --- power-law fits ---
def powerlaw(E, A, gamma):
    return A * E**(-gamma)

def fit_range(lo, hi):
    m = (Ecen >= lo) & (Ecen <= hi) & (N > 0)
    p0 = [J[m][0] * Ecen[m][0]**2.5, 2.5]
    p, cov = curve_fit(powerlaw, Ecen[m], J[m], p0=p0,
                       sigma=dJ[m], absolute_sigma=True)
    return p, np.sqrt(np.diag(cov))

print("events used:", len(E))
print("N first 5:", N[:5])
for lo, hi in [(5, 40), (40, 150)]:
    p, err = fit_range(lo, hi)
    print(f"{lo}-{hi} EeV: gamma = {p[1]:.2f} +/- {err[1]:.2f}")

# --- plot ---
plt.errorbar(Ecen[mask], (Ecen**3 * J)[mask], yerr=(Ecen**3 * dJ)[mask], fmt="o")
plt.xscale("log"); plt.yscale("log")
plt.xlabel("E [EeV]")
plt.ylabel("E³ J [EeV² km⁻² sr⁻¹ yr⁻¹]")
(BASE / "figures").mkdir(exist_ok=True)
plt.savefig(BASE / "figures" / "spectrum_E3J.png", dpi=200, bbox_inches="tight")
plt.show()