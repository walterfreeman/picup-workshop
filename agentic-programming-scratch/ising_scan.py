import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ising import random_lattice, metropolis_sweep, total_energy, magnetization

N = 32
N_SPINS = N * N
N_EQUIL = 1000
N_SAMPLE = 4000
SEED = 1
TC = 2.0 / np.log(1.0 + np.sqrt(2.0))  # exact 2D Ising critical temperature
MAX_LAG = 1000
MAX_BLOCKS = 200


def make_checkerboard(n):
    return np.indices((n, n)).sum(axis=0) % 2


def run_chain(T, rng, n_equil, n_sample, n=N, lattice=None):
    if lattice is None:
        lattice = random_lattice(n, rng)
    beta = 1.0 / T
    checkerboard = make_checkerboard(n)
    for _ in range(n_equil):
        metropolis_sweep(lattice, beta, checkerboard, rng)
    energies = np.empty(n_sample)
    mags = np.empty(n_sample)
    for k in range(n_sample):
        metropolis_sweep(lattice, beta, checkerboard, rng)
        energies[k] = total_energy(lattice)
        mags[k] = magnetization(lattice)
    return lattice, energies, mags


def correlation_function(lattice, max_r):
    c = np.empty(max_r + 1)
    for r in range(max_r + 1):
        horiz = np.mean(lattice * np.roll(lattice, -r, axis=1))
        vert = np.mean(lattice * np.roll(lattice, -r, axis=0))
        c[r] = 0.5 * (horiz + vert)
    return c


def autocorr_time(series, c=5.0, max_lag=MAX_LAG):
    """Integrated autocorrelation time via Sokal's self-consistent windowing."""
    n = len(series)
    max_lag = min(max_lag, n // 2)
    x = series - series.mean()
    var = np.mean(x**2)
    if var == 0:
        return 0.5, True
    tau = 0.5
    for t in range(1, max_lag + 1):
        rho_t = np.mean(x[: n - t] * x[t:]) / var
        tau += rho_t
        if t >= c * tau:
            return tau, True
    return tau, False  # window did not converge inside max_lag


def blocked_jackknife(data, block_size, estimator):
    """Leave-one-block-out jackknife error, correct for nonlinear estimators too."""
    n = len(next(iter(data.values())))
    n_blocks = n // block_size
    if n_blocks < 2:
        return estimator(data), np.nan
    trimmed = {k: v[: n_blocks * block_size].reshape(n_blocks, block_size) for k, v in data.items()}
    full_estimate = estimator({k: v.reshape(-1) for k, v in trimmed.items()})
    jack = np.empty(n_blocks)
    mask_template = np.ones(n_blocks, dtype=bool)
    for i in range(n_blocks):
        mask = mask_template.copy()
        mask[i] = False
        reduced = {k: v[mask].reshape(-1) for k, v in trimmed.items()}
        jack[i] = estimator(reduced)
    jack_mean = jack.mean()
    var = (n_blocks - 1) / n_blocks * np.sum((jack - jack_mean) ** 2)
    return full_estimate, np.sqrt(var)


def analyze(T, energies, mags, n_spins=N_SPINS):
    beta = 1.0 / T
    tau_e, conv_e = autocorr_time(energies)
    tau_m, conv_m = autocorr_time(mags)
    tau = max(tau_e, tau_m)
    n = len(energies)
    block_size = max(1, int(round(2 * tau)))
    n_blocks = n // block_size
    if n_blocks > MAX_BLOCKS:
        block_size = n // MAX_BLOCKS
    data = {"E": energies, "M": mags}

    e_val, e_err = blocked_jackknife(data, block_size, lambda d: d["E"].mean() / n_spins)
    m_val, m_err = blocked_jackknife(data, block_size, lambda d: np.abs(d["M"]).mean() / n_spins)
    cv_val, cv_err = blocked_jackknife(data, block_size, lambda d: beta**2 * d["E"].var() / n_spins)
    chi_val, chi_err = blocked_jackknife(
        data, block_size, lambda d: beta * (np.mean(d["M"] ** 2) - np.mean(np.abs(d["M"])) ** 2) / n_spins
    )
    return dict(
        e=e_val, e_err=e_err, m=m_val, m_err=m_err,
        cv=cv_val, cv_err=cv_err, chi=chi_val, chi_err=chi_err,
        tau=tau, conv=(conv_e and conv_m), block_size=block_size, n_blocks=len(energies) // block_size,
    )


def thermalization_check(rng):
    print("Thermalization check at T = Tc (worst case for critical slowing down):")
    hot = random_lattice(N, rng)
    cold = np.ones((N, N), dtype=np.int8)
    _, e_hot, m_hot = run_chain(TC, rng, N_EQUIL, 500, lattice=hot)
    _, e_cold, m_cold = run_chain(TC, rng, N_EQUIL, 500, lattice=cold)
    print(f"  hot start:  <E>/spin={e_hot.mean() / N_SPINS:.4f}  <|M|>/spin={np.abs(m_hot).mean() / N_SPINS:.4f}")
    print(f"  cold start: <E>/spin={e_cold.mean() / N_SPINS:.4f}  <|M|>/spin={np.abs(m_cold).mean() / N_SPINS:.4f}")
    print(f"  -> agreement after {N_EQUIL} equilibration sweeps supports adequate warmup\n")


def main():
    rng = np.random.default_rng(SEED)

    thermalization_check(rng)

    temperatures = np.concatenate([
        np.linspace(1.5, 2.0, 8),
        np.linspace(2.05, 2.5, 20),
        np.linspace(2.6, 3.5, 8),
    ])
    n = len(temperatures)
    mean_E = np.empty(n); mean_E_err = np.empty(n)
    mean_M = np.empty(n); mean_M_err = np.empty(n)
    specific_heat = np.empty(n); specific_heat_err = np.empty(n)
    susceptibility = np.empty(n); susceptibility_err = np.empty(n)
    tau_arr = np.empty(n)

    for i, T in enumerate(temperatures):
        _, energies, mags = run_chain(T, rng, N_EQUIL, N_SAMPLE)
        r = analyze(T, energies, mags)
        mean_E[i], mean_E_err[i] = r["e"], r["e_err"]
        mean_M[i], mean_M_err[i] = r["m"], r["m_err"]
        specific_heat[i], specific_heat_err[i] = r["cv"], r["cv_err"]
        susceptibility[i], susceptibility_err[i] = r["chi"], r["chi_err"]
        tau_arr[i] = r["tau"]
        flag = "" if r["conv"] else "  [tau window unconverged: critical slowing down outran sample length]"
        print(f"T={T:.3f}  <E>={mean_E[i]:.4f}+-{mean_E_err[i]:.4f}  "
              f"<|M|>={mean_M[i]:.4f}+-{mean_M_err[i]:.4f}  "
              f"Cv={specific_heat[i]:.4f}+-{specific_heat_err[i]:.4f}  "
              f"chi={susceptibility[i]:.4f}+-{susceptibility_err[i]:.4f}  "
              f"tau={tau_arr[i]:.1f}  nblocks={r['n_blocks']}{flag}")

    corr_temps = [1.5, TC, 3.5]
    max_r = N // 2
    correlations = {}
    for T in corr_temps:
        lattice, _, _ = run_chain(T, rng, N_EQUIL, N_SAMPLE)
        beta = 1.0 / T
        checkerboard = make_checkerboard(N)
        c_accum = np.zeros(max_r + 1)
        n_snaps = 30
        for _ in range(n_snaps):
            for _ in range(10):
                metropolis_sweep(lattice, beta, checkerboard, rng)
            c_accum += correlation_function(lattice, max_r)
        correlations[T] = c_accum / n_snaps

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    ax = axes[0, 0]
    ax.errorbar(temperatures, mean_M, yerr=mean_M_err, fmt="o-", ms=3, capsize=2)
    ax.axvline(TC, color="gray", ls="--", label=r"$T_c$")
    ax.set_xlabel("T"); ax.set_ylabel(r"$\langle |M| \rangle$ / spin")
    ax.set_title("Magnetization"); ax.legend()

    ax = axes[0, 1]
    ax.errorbar(temperatures, mean_E, yerr=mean_E_err, fmt="o-", ms=3, capsize=2, color="C1")
    ax.axvline(TC, color="gray", ls="--")
    ax.set_xlabel("T"); ax.set_ylabel(r"$\langle E \rangle$ / spin")
    ax.set_title("Energy")

    ax = axes[1, 0]
    ax.errorbar(temperatures, specific_heat, yerr=specific_heat_err, fmt="o-", ms=3, capsize=2, color="C2")
    ax.axvline(TC, color="gray", ls="--")
    ax.set_xlabel("T"); ax.set_ylabel(r"$C_v$ / spin")
    ax.set_title("Specific heat")

    ax = axes[1, 1]
    ax.errorbar(temperatures, susceptibility, yerr=susceptibility_err, fmt="o-", ms=3, capsize=2, color="C3")
    ax.axvline(TC, color="gray", ls="--")
    ax.set_xlabel("T"); ax.set_ylabel(r"$\chi$ / spin")
    ax.set_title("Susceptibility")

    fig.suptitle(f"2D Ising model, {N}x{N} lattice — error bars: blocked jackknife,\n"
                 f"block size set by integrated autocorrelation time")
    fig.tight_layout()
    fig.savefig("ising_observables.png", dpi=150)

    fig3, ax3 = plt.subplots(figsize=(6, 4))
    ax3.plot(temperatures, tau_arr, "o-", ms=3, color="C4")
    ax3.axvline(TC, color="gray", ls="--", label=r"$T_c$")
    ax3.set_xlabel("T"); ax3.set_ylabel(r"integrated autocorrelation time $\tau$ (sweeps)")
    ax3.set_title("Critical slowing down")
    ax3.legend()
    fig3.tight_layout()
    fig3.savefig("ising_autocorrelation_time.png", dpi=150)

    fig2, ax2 = plt.subplots(figsize=(6, 5))
    for T, c in correlations.items():
        label = r"$T_c$" if abs(T - TC) < 1e-9 else f"T = {T}"
        ax2.plot(range(max_r + 1), c, "o-", label=label)
    ax2.set_xlabel("separation r")
    ax2.set_ylabel(r"$\langle s_0 s_r \rangle$")
    ax2.set_title("Spatial spin-spin correlation")
    ax2.legend()
    fig2.tight_layout()
    fig2.savefig("ising_correlation.png", dpi=150)

    print("\nSaved ising_observables.png, ising_autocorrelation_time.png, and ising_correlation.png")


if __name__ == "__main__":
    main()
