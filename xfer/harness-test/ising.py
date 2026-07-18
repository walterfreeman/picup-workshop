"""
2D Ising model simulation via the Metropolis algorithm, with a live
matplotlib visualization of the spin lattice as it evolves.

Usage:
    python3 ising.py                  # interactive animation window
    python3 ising.py --frames 200 --save ising.gif   # headless, save a gif
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def init_lattice(n, rng):
    return rng.choice([-1, 1], size=(n, n))


def delta_energy(lattice, i, j, n):
    s = lattice[i, j]
    neighbors = (
        lattice[(i + 1) % n, j]
        + lattice[(i - 1) % n, j]
        + lattice[i, (j + 1) % n]
        + lattice[i, (j - 1) % n]
    )
    return 2 * s * neighbors


def metropolis_sweep(lattice, n, beta, rng):
    """One sweep = n*n randomly chosen single-spin flip attempts."""
    for _ in range(n * n):
        i = rng.integers(0, n)
        j = rng.integers(0, n)
        dE = delta_energy(lattice, i, j, n)
        if dE <= 0 or rng.random() < np.exp(-beta * dE):
            lattice[i, j] *= -1
    return lattice


def magnetization(lattice):
    return lattice.mean()


def main():
    parser = argparse.ArgumentParser(description="Ising model Metropolis simulation")
    parser.add_argument("--size", type=int, default=100, help="lattice side length")
    parser.add_argument("--temp", type=float, default=2.269, help="temperature (Tc ~= 2.269 for J=1, kB=1)")
    parser.add_argument("--frames", type=int, default=0, help="stop after N frames (0 = run until window closed)")
    parser.add_argument("--save", type=str, default=None, help="save animation to this gif path instead of showing it")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    if args.save:
        import matplotlib
        matplotlib.use("Agg")

    rng = np.random.default_rng(args.seed)
    n = args.size
    beta = 1.0 / args.temp
    lattice = init_lattice(n, rng)

    fig, (ax_lattice, ax_mag) = plt.subplots(1, 2, figsize=(10, 5))
    im = ax_lattice.imshow(lattice, cmap="binary", vmin=-1, vmax=1)
    ax_lattice.set_title(f"Ising lattice, T={args.temp:.3f}")
    ax_lattice.set_xticks([])
    ax_lattice.set_yticks([])

    mag_history = []
    (mag_line,) = ax_mag.plot([], [])
    ax_mag.set_xlim(0, max(args.frames, 100))
    ax_mag.set_ylim(-1, 1)
    ax_mag.set_xlabel("sweep")
    ax_mag.set_ylabel("magnetization")
    ax_mag.axhline(0, color="gray", linewidth=0.5)

    def update(frame):
        metropolis_sweep(lattice, n, beta, rng)
        im.set_data(lattice)
        mag_history.append(magnetization(lattice))
        mag_line.set_data(range(len(mag_history)), mag_history)
        if len(mag_history) > ax_mag.get_xlim()[1]:
            ax_mag.set_xlim(0, len(mag_history) + 50)
        ax_lattice.set_title(f"Ising lattice, T={args.temp:.3f}, sweep {frame + 1}")
        return im, mag_line

    frames = args.frames if args.frames > 0 else None
    anim = FuncAnimation(fig, update, frames=frames, interval=50, blit=False, repeat=False)

    if args.save:
        anim.save(args.save, writer=PillowWriter(fps=20))
        print(f"Saved animation to {args.save}")
    else:
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
