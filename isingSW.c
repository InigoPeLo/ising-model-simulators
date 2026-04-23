/*
 * Modelo de Ising 2D/3D — algoritmo Swendsen-Wang
 *
 * Compilar:
 *   gcc -O3 -march=native -o isingSW isingSW.c -lm
 *
 * Uso:
 *   ./isingSW -L 32 -t 1000 -teq 200 -Tmin 2.0 -Tmax 2.8 -S 20 -d 2
 *   ./isingSW -L 16 -t 500  -teq 100 -Tmin 4.0 -Tmax 5.5 -S 15 -d 3 --ord
 *
 * Salida:
 *   Datos_IsingSW_<d>D_L<L>_<ord>.txt  —  T paso E M
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

//Generador numeros aleatorios xoshiro256 para montecarlo
static uint64_t xs[4];

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}
static inline uint64_t xs_next(void) {
    const uint64_t r = rotl64(xs[1] * 5, 7) * 9;
    const uint64_t t = xs[1] << 17;
    xs[2] ^= xs[0]; xs[3] ^= xs[1];
    xs[1] ^= xs[2]; xs[0] ^= xs[3];
    xs[2] ^= t;     xs[3] = rotl64(xs[3], 45);
    return r;
}
static inline double rand01(void) {
    return (xs_next() >> 11) * (1.0 / (UINT64_C(1) << 53));
}
static void seed_rng(void) {
    FILE *f = fopen("/dev/urandom", "rb");
    if (f && fread(xs, sizeof(xs), 1, f) == 1) {
        fclose(f);
        if (!xs[0] && !xs[1] && !xs[2] && !xs[3]) xs[0] = 1;
        return;
    }
    if (f) fclose(f);
    /* fallback */
    uint64_t s = (uint64_t)time(NULL) ^ (uint64_t)(uintptr_t)&s;
    for (int i = 0; i < 4; i++) {
        s += 0x9e3779b97f4a7c15ULL;
        s = (s ^ (s >> 30)) * 0xbf58476d1ce4e5b9ULL;
        s = (s ^ (s >> 27)) * 0x94d049bb133111ebULL;
        xs[i] = s ^ (s >> 31);
    }
}

//Precalculamos la tabla de vecinos 
static int *neigh = NULL;

static void build_neigh(int L, int dim, int N) {
    neigh = malloc((size_t)N * dim * sizeof(int));
    if (!neigh) { perror("malloc neigh"); exit(EXIT_FAILURE); }
    int L2 = L * L;
    for (int idx = 0; idx < N; idx++) {
        if (dim == 2) {
            int i = idx / L, j = idx % L;
            neigh[idx*2 + 0] = (i+1 < L ? i+1 : 0) * L + j;
            neigh[idx*2 + 1] = i * L + (j+1 < L ? j+1 : 0);
        } else {
            int i = idx / L2, j = (idx / L) % L, k = idx % L;
            neigh[idx*3 + 0] = (i+1 < L ? i+1 : 0) * L2 + j*L + k;
            neigh[idx*3 + 1] = i*L2 + (j+1 < L ? j+1 : 0)*L + k;
            neigh[idx*3 + 2] = i*L2 + j*L + (k+1 < L ? k+1 : 0);
        }
    }
}

//Busqueda de uniones para crecer el cluster
static int *uf_par  = NULL;
static int *uf_rank = NULL;

static inline int uf_find(int i) {
    while (uf_par[i] != i) {
        uf_par[i] = uf_par[uf_par[i]];
        i = uf_par[i];
    }
    return i;
}
static inline void uf_union(int a, int b) {
    a = uf_find(a); b = uf_find(b);
    if (a == b) return;
    if (uf_rank[a] < uf_rank[b]) { int t = a; a = b; b = t; }
    uf_par[b] = a;
    if (uf_rank[a] == uf_rank[b]) uf_rank[a]++;
}

//Generacion del estado inicial
static void init_lattice(int8_t *lat, int N, int ordered) {
    if (ordered) {
        memset(lat, 1, N);
    } else {
        for (int i = 0; i < N; i++)
            lat[i] = (xs_next() & 1) ? 1 : -1;
    }
}

//Calculo de la energia usando la tabla de vecinos para calcular energia total
static double energia(const int8_t *lat, int N, int dim, double J) {
    double E = 0.0;
    for (int idx = 0; idx < N; idx++) {
        int s = lat[idx];
        for (int d = 0; d < dim; d++)
            E -= J * s * lat[neigh[idx*dim + d]];
    }
    return E;
}

//Magnetización por spin
static double magnetizacion(const int8_t *lat, int N) {
    long sum = 0;
    for (int i = 0; i < N; i++) sum += lat[i];
    double m = (double)sum / N;
    return fabs(m);
}

//Barrido de SW excluyendo los pasos de termalización t_eq y midiendo a partir de los ya termalizados t_meas
static void sw_run(int8_t *lat, int dim, int N,
                   int t_eq, int t_meas, double T, double J,
                   double *E_out, double *M_out) {
    const double p_bond = 1.0 - exp(-2.0 * J / T);

    int8_t *cspin = malloc(N);
    if (!cspin) { perror("malloc sw_run"); exit(EXIT_FAILURE); }

    /* Medida del paso 0 (antes de cualquier flip) */
    E_out[0] = energia(lat, N, dim, J);
    M_out[0] = magnetizacion(lat, N);

    int total = t_eq + t_meas;

    for (int t = 0; t < total; t++) {

        /* --- Inicializar UF --- */
        for (int i = 0; i < N; i++) { uf_par[i] = i; uf_rank[i] = 0; }

        /* --- Formar enlaces Swendsen-Wang --- */
        for (int idx = 0; idx < N; idx++) {
            int s = lat[idx];
            for (int d = 0; d < dim; d++) {
                int nb = neigh[idx*dim + d];
                if (s == lat[nb] && rand01() < p_bond)
                    uf_union(idx, nb);
            }
        }

        /* --- Asignar nuevo espín a cada clúster y aplicarlo --- */
        memset(cspin, 0, N);
        for (int idx = 0; idx < N; idx++) {
            int root = uf_find(idx);
            if (!cspin[root])
                cspin[root] = (rand01() < 0.5) ? 1 : -1;
            lat[idx] = cspin[root];
        }

        /* --- Medir solo después de la termalización --- */
        int midx = t - t_eq + 1;
        if (midx > 0) {
            E_out[midx] = energia(lat, N, dim, J);
            M_out[midx] = magnetizacion(lat, N);
        }
    }

    free(cspin);
}

//Argumentos de entrada
typedef struct {
    int    L, t_meas, t_eq, S, dim;
    double Tmin, Tmax, J;
    int    ordered;
} Args;

static void usage(const char *p) {
    fprintf(stderr,
      "\nUso: %s -L <L> -t <medida> [-teq <eq>] -Tmin <T0> -Tmax <T1>"
      " -S <nT> -d <dim> [-J <J>] [--ord]\n\n"
      "  -L     Lado de la red  (N = L^d espines)\n"
      "  -t     Pasos de medida por temperatura\n"
      "  -teq   Pasos de termalización (default: t/4, mín 10)\n"
      "  -Tmin  Temperatura mínima\n"
      "  -Tmax  Temperatura máxima\n"
      "  -S     Número de temperaturas del barrido\n"
      "  -d     Dimensión: 2 o 3\n"
      "  -J     Acoplamiento (default: 1.0)\n"
      "  --ord  Estado inicial ordenado\n\n"
      "Salida:\n"
      "  Datos_IsingSW_<d>D_L<L>_<ord>.txt   T paso E M\n\n",
      p);
    exit(EXIT_FAILURE);
}

static Args parse_args(int argc, char **argv) {
    Args a = {0, 0, -1, 0, 0, 0.0, 0.0, 1.0, 0};
    int gL=0, gt=0, gTn=0, gTx=0, gS=0, gd=0;
    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-L")    && i+1<argc) { a.L      = atoi(argv[++i]); gL=1;  }
        else if (!strcmp(argv[i],"-t")    && i+1<argc) { a.t_meas = atoi(argv[++i]); gt=1;  }
        else if (!strcmp(argv[i],"-teq")  && i+1<argc)   a.t_eq   = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-Tmin") && i+1<argc) { a.Tmin   = atof(argv[++i]); gTn=1; }
        else if (!strcmp(argv[i],"-Tmax") && i+1<argc) { a.Tmax   = atof(argv[++i]); gTx=1; }
        else if (!strcmp(argv[i],"-S")    && i+1<argc) { a.S      = atoi(argv[++i]); gS=1;  }
        else if (!strcmp(argv[i],"-J")    && i+1<argc)   a.J      = atof(argv[++i]);
        else if (!strcmp(argv[i],"-d")    && i+1<argc) { a.dim    = atoi(argv[++i]); gd=1;  }
        else if (!strcmp(argv[i],"--ord"))                a.ordered = 1;
        else { fprintf(stderr,"Argumento desconocido: %s\n", argv[i]); usage(argv[0]); }
    }
    if (!gL||!gt||!gTn||!gTx||!gS||!gd) {
        fprintf(stderr,"Faltan argumentos obligatorios.\n"); usage(argv[0]);
    }
    if (a.dim!=2 && a.dim!=3) { fprintf(stderr,"dim debe ser 2 o 3.\n"); exit(1); }
    if (a.L<=0||a.t_meas<=0||a.S<=0) { fprintf(stderr,"L, t, S deben ser > 0.\n"); exit(1); }
    if (a.Tmin<=0||a.Tmax<a.Tmin)    { fprintf(stderr,"Se requiere 0 < Tmin <= Tmax.\n"); exit(1); }
    if (a.t_eq < 0) a.t_eq = (a.t_meas/4 > 10) ? a.t_meas/4 : 10;
    return a;
}

/* ===========================================================
   Main
   =========================================================== */
int main(int argc, char **argv) {
    Args a = parse_args(argc, argv);
    seed_rng();

    int L   = a.L, dim = a.dim;
    int N   = (dim == 2) ? L*L : L*L*L;
    const char *ord_str = a.ordered ? "True" : "False";

    printf("Swendsen-Wang %dD  L=%d  N=%d  teq=%d  tmeas=%d\n",
           dim, L, N, a.t_eq, a.t_meas);
    printf("Barrido T=[%.3f, %.3f]  S=%d  J=%.2f  inicio=%s\n\n",
           a.Tmin, a.Tmax, a.S, a.J, ord_str);

    build_neigh(L, dim, N);
    uf_par  = malloc(N * sizeof(int));
    uf_rank = malloc(N * sizeof(int));
    int8_t *lat   = malloc(N);
    double *E_out = malloc((a.t_meas + 1) * sizeof(double));
    double *M_out = malloc((a.t_meas + 1) * sizeof(double));
    if (!uf_par||!uf_rank||!lat||!E_out||!M_out) {
        perror("malloc main"); exit(EXIT_FAILURE);
    }

    char fn[256];
    snprintf(fn, sizeof(fn), "Datos_IsingSW_%dD_L%d_%s.txt", dim, L, ord_str);
    FILE *fp = fopen(fn, "w");
    if (!fp) { perror("fopen"); exit(EXIT_FAILURE); }
    fprintf(fp, "# T paso E M\n");

    for (int si = 0; si < a.S; si++) {
        double T = (a.S == 1) ? a.Tmin
                   : a.Tmin + si * (a.Tmax - a.Tmin) / (a.S - 1);

        init_lattice(lat, N, a.ordered);
        sw_run(lat, dim, N, a.t_eq, a.t_meas, T, a.J, E_out, M_out);

        for (int t = 0; t <= a.t_meas; t++)
            fprintf(fp, "%.4f %d %.6g %.6g\n", T, t, E_out[t], M_out[t]);

        double m_avg = 0.0;
        for (int t = 1; t <= a.t_meas; t++) m_avg += M_out[t];
        m_avg /= a.t_meas;

        printf("[OK] T = %.4f  m = %.4f\n", T, m_avg);
    }

    fclose(fp);
    printf("\nCompletado. Datos: %s\n", fn);

    free(neigh); free(uf_par); free(uf_rank);
    free(lat); free(E_out); free(M_out);
    return 0;
}