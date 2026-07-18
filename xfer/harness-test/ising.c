/*
 * 2D Ising model, Metropolis algorithm, checkerboard (red-black) update
 * parallelized with OpenMP, rendered live via SDL2.
 *
 * Checkerboard correctness note: partitioning the lattice into two
 * color classes by (i+j)%2 makes every site in one class conditionally
 * independent of the others in that class given the other class fixed
 * (its four neighbors are all the opposite color). So flipping every
 * site of one color in parallel, using only the *pre-sweep* state of
 * the other color, is a valid Gibbs/Metropolis sweep -- not an
 * approximation of the sequential single-spin algorithm.
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp ising.c -o ising -lSDL2 -lm
 *
 * Run:
 *   ./ising [N] [temp] [sweeps_per_frame] [max_frames] [seed]
 *   N                 lattice side length            (default 200)
 *   temp              temperature, Tc ~= 2.269        (default 2.269)
 *   sweeps_per_frame  simulation sweeps per rendered frame (default 1)
 *   max_frames        stop after N frames, 0 = run until window closed (default 0)
 *   seed              RNG seed (default: time-based)
 *
 * Env:
 *   ISING_BENCH=1     skip SDL entirely, print raw sweeps/sec and exit
 */
#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <omp.h>

typedef int8_t spin_t;

static inline uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline uint64_t xorshift64star(uint64_t *state) {
    uint64_t x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    return x * 0x2545F4914F6CDD1DULL;
}

static inline double rng_double(uint64_t *state) {
    return (double)(xorshift64star(state) >> 11) * (1.0 / 9007199254740992.0);
}

static spin_t *lattice;
static int N;
static double exp_table[5]; /* indexed by dE/4 + 2, dE in {-8,-4,0,4,8} */
static uint64_t *thread_state;
static int n_threads;

static void init_lattice(uint64_t seed) {
    uint64_t s = seed;
    for (int i = 0; i < N * N; i++) {
        lattice[i] = (splitmix64(&s) & 1) ? 1 : -1;
    }
}

static void build_exp_table(double temp) {
    double beta = 1.0 / temp;
    for (int k = -2; k <= 2; k++) {
        exp_table[k + 2] = exp(-beta * (4.0 * k));
    }
}

/*
 * One color pass of the checkerboard sweep. color 0 = sites where (i+j)%2==0.
 *
 * NOTE: uses "#pragma omp for", a *worksharing* construct, not
 * "#pragma omp parallel for". This function must be called from inside
 * an already-active "#pragma omp parallel" region -- it does not spawn
 * its own thread team. Spawning/joining a team costs real wall-clock
 * (thread wake + barrier teardown), and paying that cost twice per
 * sweep dominated runtime in early testing at every lattice size tried,
 * making 32 threads *slower* than 1. Opening the team once for the
 * whole run/frame-loop and reusing it via worksharing constructs (this
 * "for", and "single"/"barrier" elsewhere) amortizes that cost instead.
 */
static void sweep_color_worker(int color) {
    #pragma omp for schedule(static)
    for (int i = 0; i < N; i++) {
        uint64_t *state = &thread_state[omp_get_thread_num()];
        int j_start = ((i % 2) == color) ? 0 : 1;
        for (int j = j_start; j < N; j += 2) {
            int idx = i * N + j;
            spin_t s = lattice[idx];
            int up    = lattice[((i - 1 + N) % N) * N + j];
            int down  = lattice[((i + 1) % N) * N + j];
            int left  = lattice[i * N + (j - 1 + N) % N];
            int right = lattice[i * N + (j + 1) % N];
            int neighbor_sum = up + down + left + right;
            int dE = 2 * s * neighbor_sum;
            if (dE <= 0) {
                lattice[idx] = (spin_t)(-s);
            } else {
                double p = exp_table[(dE / 4) + 2];
                if (rng_double(state) < p) lattice[idx] = (spin_t)(-s);
            }
        }
    }
}

/* Must be called from inside an active "#pragma omp parallel" region. */
static void metropolis_sweep_worker(void) {
    sweep_color_worker(0);
    /* "#pragma omp for" has an implicit barrier at loop exit, so every
     * thread has finished writing color 0 before any thread reads it
     * as a neighbor while updating color 1. */
    sweep_color_worker(1);
}

/* Standalone helper for one-off use (outside a persistent parallel region);
 * spawns its own team, so don't call this in a hot per-sweep/per-frame loop. */
static double magnetization(void) {
    long sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N * N; i++) sum += lattice[i];
    return (double)sum / (N * N);
}

static void run_benchmark(int sweeps) {
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    #pragma omp parallel
    {
        for (int k = 0; k < sweeps; k++) metropolis_sweep_worker();
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double secs = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    printf("N=%d threads=%d sweeps=%d total=%.4fs => %.4f ms/sweep, %.1f sweeps/sec\n",
           N, n_threads, sweeps, secs, secs / sweeps * 1000.0, sweeps / secs);
    printf("magnetization=%.4f\n", magnetization());
}

int main(int argc, char **argv) {
    N = argc > 1 ? atoi(argv[1]) : 200;
    double temp = argc > 2 ? atof(argv[2]) : 2.269;
    int sweeps_per_frame = argc > 3 ? atoi(argv[3]) : 1;
    int max_frames = argc > 4 ? atoi(argv[4]) : 0;
    uint64_t seed = argc > 5 ? (uint64_t)strtoull(argv[5], NULL, 10) : (uint64_t)time(NULL);

    /* Empirically, on this machine's shared/virtualized cores, thread
     * counts beyond ~8 lose more to scheduling contention than they gain
     * in parallelism (measured: 8 threads beat both 1 and 32 at every
     * lattice size tried). Cap the default there; OMP_NUM_THREADS still
     * overrides if set, since the right cap is hardware/load dependent. */
    if (!getenv("OMP_NUM_THREADS")) {
        int cap = omp_get_max_threads();
        if (cap > 8) cap = 8;
        omp_set_num_threads(cap);
    }
    n_threads = omp_get_max_threads();
    thread_state = malloc(sizeof(uint64_t) * n_threads);
    for (int t = 0; t < n_threads; t++) {
        uint64_t seedgen = seed + 0x9E3779B9u * (uint64_t)(t + 1);
        thread_state[t] = splitmix64(&seedgen);
    }

    lattice = malloc(sizeof(spin_t) * N * N);
    init_lattice(seed);
    build_exp_table(temp);

    if (getenv("ISING_BENCH")) {
        int bench_sweeps = getenv("ISING_BENCH_SWEEPS") ? atoi(getenv("ISING_BENCH_SWEEPS")) : 200;
        run_benchmark(bench_sweeps);
        free(lattice);
        free(thread_state);
        return 0;
    }

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        return 1;
    }

    int win_size = 800;
    SDL_Window *win = SDL_CreateWindow("Ising model (C + SDL2, OpenMP checkerboard Metropolis)",
                                        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                        win_size, win_size, SDL_WINDOW_SHOWN);
    if (!win) {
        fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        return 1;
    }
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    if (!ren) ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_SOFTWARE);
    SDL_RenderSetLogicalSize(ren, N, N);
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "0"); /* nearest-neighbor */

    SDL_Texture *tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGB24,
                                          SDL_TEXTUREACCESS_STREAMING, N, N);
    uint8_t *pixels = malloc((size_t)N * N * 3);

    int frame = 0;
    int running = 1;
    long mag_sum = 0;
    Uint64 perf_freq = SDL_GetPerformanceFrequency();
    Uint64 fps_t0 = SDL_GetPerformanceCounter();
    int fps_frames = 0;

    /* One thread team for the entire animation loop. Every "#pragma omp
     * single" and "#pragma omp for" below carries an implicit barrier,
     * so all threads stay in lockstep without ever tearing the team down
     * between frames. Only the thread that happens to hit each "single"
     * block touches SDL -- SDL calls are not thread-safe to spray across
     * threads, so all of them are confined to these serial sections. */
    #pragma omp parallel default(shared)
    {
        while (1) {
            #pragma omp single
            {
                SDL_Event e;
                while (SDL_PollEvent(&e)) {
                    if (e.type == SDL_QUIT) running = 0;
                    if (e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_ESCAPE) running = 0;
                }
            } /* implicit barrier: every thread sees the same 'running' */

            if (!running) break;

            for (int k = 0; k < sweeps_per_frame; k++) metropolis_sweep_worker();

            #pragma omp single
            { mag_sum = 0; }

            #pragma omp for reduction(+:mag_sum)
            for (int i = 0; i < N * N; i++) {
                spin_t s = lattice[i];
                uint8_t v = s > 0 ? 255 : 0;
                pixels[i * 3 + 0] = v;
                pixels[i * 3 + 1] = v;
                pixels[i * 3 + 2] = v;
                mag_sum += s;
            }

            #pragma omp single
            {
                SDL_UpdateTexture(tex, NULL, pixels, N * 3);
                SDL_RenderClear(ren);
                SDL_RenderCopy(ren, tex, NULL, NULL);
                SDL_RenderPresent(ren);

                frame++;
                fps_frames++;
                Uint64 now = SDL_GetPerformanceCounter();
                double elapsed = (double)(now - fps_t0) / perf_freq;
                if (elapsed >= 0.5) {
                    double fps = fps_frames / elapsed;
                    double sweeps_per_sec = fps * sweeps_per_frame;
                    char title[256];
                    snprintf(title, sizeof(title),
                             "Ising N=%d T=%.3f | frame %d | %.1f fps | %.0f sweeps/sec | m=%.3f",
                             N, temp, frame, fps, sweeps_per_sec, (double)mag_sum / (N * N));
                    SDL_SetWindowTitle(win, title);
                    fps_t0 = now;
                    fps_frames = 0;
                }

                if (max_frames > 0 && frame >= max_frames) running = 0;
            } /* implicit barrier before the loop re-checks events */
        }
    }

    free(pixels);
    SDL_DestroyTexture(tex);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    free(lattice);
    free(thread_state);
    return 0;
}
