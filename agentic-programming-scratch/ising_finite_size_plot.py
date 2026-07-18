import glob
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ising_scan import TC

files = glob.glob("finite_size_result_L*.json")
data = []
for fn in files:
    with open(fn) as f:
        data.append(json.load(f))
data.sort(key=lambda d: d["L"])

sizes = [d["L"] for d in data]
temps = np.array(data[0]["temps"])

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
colors = plt.cm.viridis(np.linspace(0, 1, len(sizes)))
for d, color in zip(data, colors):
    L = d["L"]
    axes[0].errorbar(temps, d["mean_M"], yerr=d["mean_M_err"], fmt="o-", ms=3, capsize=1.5, color=color, label=f"L={L}")
    axes[1].errorbar(temps, d["chi"], yerr=d["chi_err"], fmt="o-", ms=3, capsize=1.5, color=color, label=f"L={L}")
axes[0].axvline(TC, color="gray", ls="--")
axes[1].axvline(TC, color="gray", ls="--")
axes[0].set_xlabel("T"); axes[0].set_ylabel(r"$\langle |M| \rangle$ / spin")
axes[0].set_title("Magnetization"); axes[0].legend(fontsize=8)
axes[1].set_xlabel("T"); axes[1].set_ylabel(r"$\chi$ / spin")
axes[1].set_title("Susceptibility"); axes[1].legend(fontsize=8)
fig.suptitle("Finite-size dependence of the 2D Ising transition")
fig.tight_layout()
fig.savefig("ising_finite_size_overlay.png", dpi=150)

L_arr = np.array(sizes, dtype=float)
chi_max = np.array([max(d["chi"]) for d in data])
chi_max_err = np.array([d["chi_err"][int(np.argmax(d["chi"]))] for d in data])
T_peak = np.array([temps[int(np.argmax(d["chi"]))] for d in data])

log_L = np.log(L_arr)
log_chi = np.log(chi_max)
slope, intercept = np.polyfit(log_L, log_chi, 1)

fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
ax = axes2[0]
ax.errorbar(L_arr, chi_max, yerr=chi_max_err, fmt="o", color="C3")
fit_L = np.linspace(L_arr.min(), L_arr.max(), 50)
ax.plot(fit_L, np.exp(intercept) * fit_L**slope, "--", color="C3", label=f"fit slope = {slope:.2f}")
ax.plot(fit_L, np.exp(intercept) * fit_L**1.75, ":", color="gray", label=r"theory $\gamma/\nu=7/4$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("L"); ax.set_ylabel(r"$\chi_{max}$ / spin")
ax.set_title("Susceptibility peak scaling"); ax.legend()

ax = axes2[1]
ax.plot(1.0 / L_arr, T_peak, "o", color="C0")
fit = np.polyfit(1.0 / L_arr, T_peak, 1)
fit_x = np.linspace(0, (1.0 / L_arr).max(), 50)
ax.plot(fit_x, np.polyval(fit, fit_x), "--", color="C0", label=f"extrapolated Tc = {fit[1]:.4f}")
ax.axhline(TC, color="gray", ls=":", label=f"exact Tc = {TC:.4f}")
ax.set_xlabel("1/L"); ax.set_ylabel(r"$T_{peak}(L)$")
ax.set_title("Pseudo-critical temperature shift"); ax.legend()

fig2.tight_layout()
fig2.savefig("ising_finite_size_scaling.png", dpi=150)

print(f"L values included: {sizes}")
print(f"Fitted chi_max ~ L^{slope:.3f}  (theory: 7/4 = 1.75)")
print(f"Extrapolated Tc(L->inf) = {fit[1]:.4f}  (exact = {TC:.4f})")
print("Saved ising_finite_size_overlay.png and ising_finite_size_scaling.png")
