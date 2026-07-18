import sys
import json
import time
import numpy as np

from ising_scan import run_chain, analyze, TC  # noqa: F401 (TC re-exported for convenience)

N_EQUIL = 1500
N_SAMPLE = 3000
SEED_BASE = 100

TEMPS = np.concatenate([
    np.linspace(1.6, 2.0, 4),
    np.linspace(2.02, 2.5, 20),
    np.linspace(2.6, 3.2, 4),
]).tolist()


def main():
    L = int(sys.argv[1])
    t0 = time.perf_counter()
    rng = np.random.default_rng(SEED_BASE + L)
    n_spins = L * L
    n = len(TEMPS)
    mean_M = np.empty(n); mean_M_err = np.empty(n)
    cv = np.empty(n); cv_err = np.empty(n)
    chi = np.empty(n); chi_err = np.empty(n)

    for i, T in enumerate(TEMPS):
        # Below Tc, a random start needs ~L^2 sweeps to coarsen into a single
        # ordered domain (curvature-driven growth) -- far more than N_EQUIL
        # for large L. Starting cold (already ordered) sidesteps that entirely.
        lattice = np.ones((L, L), dtype=np.int8) if T < TC else None
        _, energies, mags = run_chain(T, rng, N_EQUIL, N_SAMPLE, n=L, lattice=lattice)
        r = analyze(T, energies, mags, n_spins=n_spins)
        mean_M[i], mean_M_err[i] = r["m"], r["m_err"]
        cv[i], cv_err[i] = r["cv"], r["cv_err"]
        chi[i], chi_err[i] = r["chi"], r["chi_err"]

    record = dict(
        L=L,
        temps=TEMPS,
        mean_M=mean_M.tolist(), mean_M_err=mean_M_err.tolist(),
        cv=cv.tolist(), cv_err=cv_err.tolist(),
        chi=chi.tolist(), chi_err=chi_err.tolist(),
    )
    outfile = f"finite_size_result_L{L}.json"
    with open(outfile, "w") as f:
        json.dump(record, f)

    peak_idx = int(np.argmax(chi))
    elapsed = time.perf_counter() - t0
    print(f"L={L}  chi_max={chi[peak_idx]:.3f}+-{chi_err[peak_idx]:.3f}  "
          f"T_peak={TEMPS[peak_idx]:.3f}  elapsed={elapsed:.1f}s  -> {outfile}")


if __name__ == "__main__":
    main()
