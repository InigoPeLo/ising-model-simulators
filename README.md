# Ising Model Monte Carlo Simulations (2D & 3D)

This repository contains high-performance C implementations of the 2D and 3D Ising Model using two distinct Monte Carlo algorithms: **Heat Bath** (single-spin flip) and **Swendsen-Wang** (cluster flip). 

These tools are designed for statistical mechanics simulations, allowing you to study phase transitions, energy, and magnetization across a range of temperatures.

## Features

* **Two Algorithms:** Compares local spin updates (Heat Bath) with global cluster updates (Swendsen-Wang) to observe differences in critical slowing down.
* **2D & 3D Support:** Simulates square (2D) and cubic (3D) lattices with periodic boundary conditions.
* **High-Quality RNG:** Uses the robust and extremely fast `xoshiro256**` pseudo-random number generator, seeded via `/dev/urandom`.
* **Optimized:** Implements precomputed probability tables (Heat Bath) and an efficient Union-Find data structure (Swendsen-Wang).

## Included Files

* `isingHB.c` - Heat Bath algorithm implementation.
* `isingSW.c` - Swendsen-Wang algorithm implementation.

### Algorithm Comparison
* **Heat Bath (`isingHB.c`):** A local update algorithm. In each sweep, every spin is updated by sampling directly from its conditional probability distribution based on its neighbors. It is trivial to implement and has low overhead, but suffers from critical slowing down near phase transitions.
* **Swendsen-Wang (`isingSW.c`):** A non-local cluster algorithm. It probabilistically forms bonds between aligned neighboring spins to create clusters, and then flips entire clusters at once. This significantly mitigates critical slowing down compared to single-spin methods.

---

## Compilation

Both scripts require standard C libraries and the math library (`-lm`). For the best performance, compile with high optimization flags:

```bash
# Compile Heat Bath
gcc -O3 -march=native -o isingHB isingHB.c -lm

# Compile Swendsen-Wang
gcc -O3 -march=native -o isingSW isingSW.c -lm

```
## Usage

Both executables share the same command-line interface. 

### Basic Syntax
```bash
./isingHB -L <Size> -t <Measurements> -Tmin <T_start> -Tmax <T_end> -S <Steps> -d <Dimension> [options]
```

### Command-Line Arguments

| Argument | Description | Required |
| :--- | :--- | :---: |
| `-L <int>` | Linear size of the lattice (N = L^d spins). | Yes |
| `-t <int>` | Number of measurement steps per temperature. | Yes |
| `-teq <int>` | Thermalization (equilibration) steps. *(Default: t/4, min 10)* | No |
| `-Tmin <float>` | Minimum temperature for the sweep. | Yes |
| `-Tmax <float>` | Maximum temperature for the sweep. | Yes |
| `-S <int>` | Number of temperature steps in the sweep. | Yes |
| `-d <int>` | Dimensions of the lattice (must be `2` or `3`). | Yes |
| `-J <float>` | Interaction coupling constant. *(Default: 1.0)* | No |
| `--ord` | Start with an ordered initial state (all spins +1). If omitted, starts random. | No |

### Examples

**1. 2D Swendsen-Wang Sweep:**
Simulate a 32x32 lattice in 2D, measuring 1000 times per temperature, with 200 equilibration steps. Sweep from T=2.0 to T=2.8 across 20 points:
```bash
./isingSW -L 32 -t 1000 -teq 200 -Tmin 2.0 -Tmax 2.8 -S 20 -d 2
```

**2. 3D Heat Bath with Ordered Initial State:**
Simulate a 16x16x16 lattice in 3D, measuring 500 times, with 100 equilibration steps. Sweep from T=4.0 to T=5.5 across 15 points, starting from an ordered state:
```bash
./isingHB -L 16 -t 500 -teq 100 -Tmin 4.0 -Tmax 5.5 -S 15 -d 3 --ord
```

---

## Output Format

The programs generate a text file containing the results of the simulation. The filename is automatically generated based on your parameters:
* `Datos_IsingHB_<d>D_L<L>_<ord>.txt`
* `Datos_IsingSW_<d>D_L<L>_<ord>.txt`

**Columns in the output file:**
1. **`T`** - Current Temperature
2. **`paso`** - Measurement Step number
3. **`E`** - Total Energy of the lattice at this step
4. **`M`** - Absolute Magnetization per spin at this step

You can easily plot these files using Python (Matplotlib, Pandas) or Gnuplot to visualize Phase Transitions, Heat Capacity, and Magnetic Susceptibility.
