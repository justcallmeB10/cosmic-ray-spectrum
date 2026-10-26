
import zipfile
import pandas as pd
from pathlib import Path
BASE = Path(r"C:\Users\ssat-\Desktop\GitHub\github.2\cosmic-ray-spectrum")
DATA = BASE / "data"
zs = zipfile.ZipFile(DATA / "summary.zip")
za = zipfile.ZipFile(DATA / "auxiliary.zip")
ev = pd.read_csv(zs.open("summary/dataSummarySD1500.csv"))
ex = pd.read_csv(za.open("auxiliary/sd1500exposure.csv"))
print(ev.columns.tolist())
print(len(ev))
print(ev.head(3))
print(ex.columns.tolist())
print(ex.head())
import numpy as np
import matplotlib.pyplot as plt
E = ev["sd_energy"].dropna()
print(E.describe())
print(ev.shape, ex.shape)
print(ex["sd_exposure"].sum())
plt.hist(np.log10(E), bins=40)
plt.yscale("log")
plt.xlabel("log10(sd_energy)")
plt.ylabel("events")
plt.show()
print(ex["sd_exposure"].iloc[-1], ex["sd_exposure_all"].iloc[-1])
Emin = 2.5   # EeV, from the Auger notebook
sel = ev["sd_energy"].notna() & (ev["sd_energy"] >= Emin) & (ev["sd_theta"] < 60)
E = ev.loc[sel, "sd_energy"].values
print("events used:", len(E))                   
E = E[E >= Emin]
exposure = ex["sd_exposure"].iloc[-1] 
edges = 10**np.arange(np.log10(Emin), 2.25, 0.1)   
N, _ = np.histogram(E, bins=edges)
dE = np.diff(edges)
Ecen = np.sqrt(edges[:-1] * edges[1:])
J = N / (exposure * dE)
dJ = np.sqrt(N) / (exposure * dE)   
mask = N > 0
plt.errorbar(Ecen[mask], (Ecen**3 * J)[mask], yerr=(Ecen**3 * dJ)[mask], fmt="o")
plt.xscale("log"); plt.yscale("log")
plt.xlabel("E [EeV]"); plt.ylabel("E³ J")
plt.show()
(BASE / "figures").mkdir(exist_ok=True)
plt.xlabel("E [EeV]"); plt.ylabel("E³ J [EeV² km⁻² sr⁻¹ yr⁻¹]")  # check units in the notebook
plt.savefig(BASE / "figures" / "spectrum_E3J.png", dpi=200, bbox_inches="tight")
from scipy.optimize import curve_fit
def powerlaw(E, A, gamma):
    return A * E**(-gamma)
def fit_range(lo, hi):
    m = (Ecen >= lo) & (Ecen <= hi) & (N > 0)
    p0 = [J[m][0] * Ecen[m][0]**2.5, 2.5]
    p, cov = curve_fit(powerlaw, Ecen[m], J[m], p0=p0,
                       sigma=dJ[m], absolute_sigma=True)
    return p, np.sqrt(np.diag(cov))
for lo, hi in [(5, 40), (40, 150)]:
    p, err = fit_range(lo, hi)
    print(f"{lo}-{hi} EeV: gamma = {p[1]:.2f} +/- {err[1]:.2f}")
    plt.errorbar(Ecen[mask], (Ecen**3 * J)[mask], yerr=(Ecen**3 * dJ)[mask], fmt="o")
plt.xscale("log"); plt.yscale("log")
plt.xlabel("E [EeV]")
plt.ylabel("E³ J [EeV² km⁻² sr⁻¹ yr⁻¹]")
(BASE / "figures").mkdir(exist_ok=True)
plt.savefig(BASE / "figures" / "spectrum_E3J.png", dpi=200, bbox_inches="tight")
plt.show()
print(ev["sd_theta"].describe())