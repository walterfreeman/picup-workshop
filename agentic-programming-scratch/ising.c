#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

static int N;
static int8_t *lattice;

static inline int idx(int i, int j) {
    return i * N + j;
}

static inline double uniform(void) {
    return rand() / (RAND_MAX + 1.0);
}

static void init_lattice(unsigned int seed) {
    srand(seed);
    for (int i = 0; i < N * N; i++) {
        lattice[i] = (int8_t) ((rand() % 2) * 2 - 1);
    }
}

static void metropolis_sweep(double beta) {
    for (int k = 0; k < N * N; k++) {
        int i = rand() % N;
        int j = rand() % N;
        int s = lattice[idx(i, j)];
        int neighbor_sum = lattice[idx((i + N - 1) % N, j)] + lattice[idx((i + 1) % N, j)]
                          + lattice[idx(i, (j + N - 1) % N)] + lattice[idx(i, (j + 1) % N)];
        double dE = 2.0 * s * neighbor_sum;
        if (dE <= 0.0 || uniform() < exp(-beta * dE)) {
            lattice[idx(i, j)] = (int8_t) -s;
        }
    }
}

static double total_energy(void) {
    double e = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int s = lattice[idx(i, j)];
            e -= s * (lattice[idx(i, (j + 1) % N)] + lattice[idx((i + 1) % N, j)]);
        }
    }
    return e;
}

static long magnetization(void) {
    long m = 0;
    for (int i = 0; i < N * N; i++) m += lattice[i];
    return m;
}

int main(void) {
    N = 32;
    double T = 2.269;
    int n_equil = 1000;
    int n_sweeps = 2000;

    lattice = malloc((size_t) N * N * sizeof(int8_t));
    init_lattice(0);

    double beta = 1.0 / T;

    for (int s = 0; s < n_equil; s++) {
        metropolis_sweep(beta);
    }

    double e_sum = 0.0, m_sum = 0.0;
    for (int s = 0; s < n_sweeps; s++) {
        metropolis_sweep(beta);
        e_sum += total_energy();
        m_sum += fabs((double) magnetization());
    }

    int n_spins = N * N;
    printf("N = %dx%d, T = %.3f\n", N, N, T);
    printf("<E>/spin = %.4f\n", e_sum / n_sweeps / n_spins);
    printf("<|M|>/spin = %.4f\n", m_sum / n_sweeps / n_spins);

    free(lattice);
    return 0;
}
