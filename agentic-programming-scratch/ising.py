import numpy as np


def random_lattice(n, rng):
    return rng.choice(np.array([-1, 1], dtype=np.int8), size=(n, n))


def _neighbor_sum(lattice):
    return (
        np.roll(lattice, 1, axis=0) + np.roll(lattice, -1, axis=0)
        + np.roll(lattice, 1, axis=1) + np.roll(lattice, -1, axis=1)
    )


def metropolis_sweep(lattice, beta, checkerboard, rng):
    """One sweep = both checkerboard colors updated once (Metropolis, periodic BCs).

    Same-color sites share no bond, so a color's spins can be flipped
    simultaneously from a single vectorized accept/reject draw.
    """
    for color in (0, 1):
        mask = checkerboard == color
        dE = 2.0 * lattice * _neighbor_sum(lattice)
        accept = (dE <= 0.0) | (rng.random(lattice.shape) < np.exp(-beta * dE))
        flip = mask & accept
        lattice[flip] *= -1
    return lattice


def total_energy(lattice):
    # Only right+down neighbors: each bond <ij> counted exactly once.
    return -np.sum(lattice * (np.roll(lattice, -1, axis=0) + np.roll(lattice, -1, axis=1)))


def magnetization(lattice):
    return lattice.sum()


def run(n=32, temperature=2.269, n_equil_sweeps=1000, n_sweeps=2000, seed=0):
    """Equilibrate then sample. Returns final lattice and per-sweep E, M series."""
    rng = np.random.default_rng(seed)
    lattice = random_lattice(n, rng)
    beta = 1.0 / temperature
    checkerboard = np.indices((n, n)).sum(axis=0) % 2

    for _ in range(n_equil_sweeps):
        metropolis_sweep(lattice, beta, checkerboard, rng)

    energies = np.empty(n_sweeps)
    mags = np.empty(n_sweeps)
    for k in range(n_sweeps):
        metropolis_sweep(lattice, beta, checkerboard, rng)
        energies[k] = total_energy(lattice)
        mags[k] = magnetization(lattice)

    return lattice, energies, mags


if __name__ == "__main__":
    N = 32
    T = 2.269  # ~ critical temperature for the 2D Ising model, J = k_B = 1
    lattice, energies, mags = run(n=N, temperature=T)

    n_spins = N * N
    print(f"N = {N}x{N}, T = {T}")
    print(f"<E>/spin = {energies.mean() / n_spins:.4f}")
    print(f"<|M|>/spin = {np.abs(mags).mean() / n_spins:.4f}")
