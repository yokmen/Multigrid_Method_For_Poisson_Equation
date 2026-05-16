# Multigrid Method for the Laplace Equation

Numerical solver in C for a 2D Poisson boundary value problem, comparing three iterative methods built on top of a direct UMFPACK solve:

1. **Two-grid method** with symmetrised Gauss-Seidel smoothing
2. **Multigrid V-cycle** (recursive generalisation of the two-grid scheme)
3. **Conjugate gradient** preconditioned by the two-grid method

> MATH-H401 — Numerical methods — Université Libre de Bruxelles, École Polytechnique de Bruxelles
> Author: Yvahn Sancar — Academic year 2022–2023

## Problem

The physical setup is the electrostatic potential inside a capacitor formed by two rectangular conductors embedded in each other (the inner conductor at −5 V, the outer at +5 V). Assuming invariance along *z*, the problem reduces to a 2D Poisson equation on Ω = [0, 5] × [0, 4] \ [3, 4] × [1, 3]:

```
 -∂²u/∂x² - ∂²u/∂y² = ρ/ε      in Ω
                  u = -5 V     on Γᵢ  (internal Dirichlet)
                  u = +5 V     on Γₑ  (external Dirichlet)
```

with ρ = 0 and ε = 27.7 ε₀.

The PDE is discretised with the second-order central finite-difference stencil, producing a sparse linear system **A u = b** stored in CSR format. UMFPACK (SuiteSparse) provides the reference direct solution and the coarse-grid solve.

## Methods overview

### Two-grid (`twogrid.c`)
Each iteration performs:
- Pre-smoothing — one forward Gauss-Seidel sweep ([gs_prog.c](gs_prog.c))
- Restriction of the residual to the coarse grid via full-weighting, R = ¼ Pᵀ ([restriction.c](restriction.c))
- Coarse-grid solve of A_c u_c = r_c with UMFPACK ([prob_grossier.c](prob_grossier.c))
- Prolongation back to the fine grid via bilinear interpolation ([prolongation.c](prolongation.c))
- Correction u ← u + P u_c
- Post-smoothing — one backward Gauss-Seidel sweep ([gs_reg.c](gs_reg.c))

```
Algorithm 1 — Two-grid method
─────────────────────────────────────────────────────────────────────────────
Initialization : u⁽⁰⁾ = 0
for m = 1, 2, ... do
    Pre-smoothing  (one iteration)        r⁽ᵐ⁾   = b − A u⁽ᵐ⁻¹⁾           (gs_prog.c)
    Restriction on the coarse grid        r_c⁽ᵐ⁾ = R r⁽ᵐ⁾                 (restriction.c)
    Coarse-grid solve                     u_c⁽ᵐ⁾ = A_c⁻¹ r_c⁽ᵐ⁾           (prob_grossier.c)
    Prolongation on the fine grid         P_c⁽ᵐ⁾ = P u_c⁽ᵐ⁾               (prolongation.c)
    Correction                            u_corr⁽ᵐ⁾ = u⁽ᵐ⁾ + P_c⁽ᵐ⁾      (twogrid.c)
    Post-smoothing (one iteration)        A u⁽ᵐ⁾ = u_corr⁽ᵐ⁾ → u⁽ᵐ⁾      (gs_reg.c)
end for
```

```
Algorithm 2 — Forward Gauss-Seidel sweep with one iteration
─────────────────────────────────────────────────────────────────────────────
Initialization : u⁽⁰⁾ = u⁽ᵐ⁻¹⁾ (arbitrary)
    r⁽⁰⁾ = b − A u⁽⁰⁾
    d⁽¹⁾ = tril(A)⁻¹ r⁽⁰⁾
    u⁽¹⁾ = u⁽⁰⁾ + d⁽¹⁾
    r⁽¹⁾ = b − A u⁽¹⁾
```

```
Algorithm 3 — Backward Gauss-Seidel sweep with one iteration
─────────────────────────────────────────────────────────────────────────────
Initialization : u⁽⁰⁾ = u_corr⁽ᵐ⁾
    r⁽⁰⁾ = b − A u⁽⁰⁾
    d⁽¹⁾ = triu(A)⁻¹ r⁽⁰⁾
    u⁽¹⁾ = u⁽⁰⁾ + d⁽¹⁾
```

### Multigrid V-cycle (`multigrid.c`)
Recursive application of the two-grid scheme down to a chosen coarsest level, then prolonged back up. The coarsest problem is solved directly by UMFPACK; every other level is smoothed before/after the coarse-grid correction.

```
Algorithm 4 — Multigrid method (V-cycle)
─────────────────────────────────────────────────────────────────────────────
Initialization :
    u⁽⁰⁾ = 0

for m = 1, 2, ... do
    for k = 0, ..., level do
        if k = 0 then
            Pre-smoothing (one iteration)   r₀⁽ᵐ⁾    = b − A₀ u⁽ᵐ⁻¹⁾
            Restriction                     r_c,0⁽ᵐ⁾ = R₀ r₀⁽ᵐ⁾
        else
            Pre-smoothing (one iteration)   A_k⁽ᵐ⁾ c = r_{k−1}⁽ᵐ⁾  → c_k⁽ᵐ⁾    (*)
            Residual                        r_k⁽ᵐ⁾    = r_{k−1}⁽ᵐ⁾ − A_k c
            Restriction                     r_c,k⁽ᵐ⁾  = R_k r_k⁽ᵐ⁾
            if k = level then
                Coarsest-grid solve         u_coarsest⁽ᵐ⁾ = A_{c,k}⁻¹ r_c,k⁽ᵐ⁾
            end if
        end if
    end for

    for k = level, ..., 0 do
        if k = level then
            Prolongation                    P_c,k⁽ᵐ⁾    = P_k u_coarsest⁽ᵐ⁾
            Correction                      u_corr,k⁽ᵐ⁾ = c_k⁽ᵐ⁾ + P_c,k⁽ᵐ⁾
            Post-smoothing                  A_k⁽ᵐ⁾ u_corr,k⁽ᵐ⁾ = r_k⁽ᵐ⁾ → ũ_k⁽ᵐ⁾
        else
            Prolongation                    P_c,k⁽ᵐ⁾    = P_k ũ_k⁽ᵐ⁾
            Correction                      u_corr,k⁽ᵐ⁾ = c_k⁽ᵐ⁾ + P_c,k⁽ᵐ⁾
            Post-smoothing                  A_k⁽ᵐ⁾ u_corr,k⁽ᵐ⁾ = r_k⁽ᵐ⁾ → ũ_k⁽ᵐ⁾
        end if
    end for
    u⁽ᵐ⁾ = ũ_0⁽ᵐ⁾
end for

(*) Using the zero vector as initial approximation.
```

### Preconditioned conjugate gradient (`conjugate.c`)
Standard CG variant where B⁻¹ rₘ is approximated by one two-grid cycle MG(0, rₘ).

```
Algorithm 5 — Conjugate gradient method (standard variant)
─────────────────────────────────────────────────────────────────────────────
Initialization :
    u₀  = 0 (arbitrary)
    r₀  = b − A u₀
    d₋₁ = 0

for m = 0, 1, ... do
    βₘ = (rₘᵀ B⁻¹ rₘ) / (r_{m−1}ᵀ B⁻¹ r_{m−1})        (β₀ = 0)
    dₘ = B⁻¹ rₘ + βₘ d_{m−1}
    αₘ = (rₘᵀ B⁻¹ rₘ) / (dₘᵀ A dₘ)
    u_{m+1} = uₘ + αₘ dₘ
    r_{m+1} = rₘ − αₘ A dₘ
end for
```

## Results

| Method                          | Behaviour                                                                 |
|---------------------------------|---------------------------------------------------------------------------|
| Direct UMFPACK reference        | ‖r‖ ≈ 2·10⁻¹⁵ (unit roundoff)                                             |
| Two-grid                        | Reaches unit roundoff in ≈ 21 cycles                                      |
| Multigrid V-cycle (levels 2–6)  | Reaches unit roundoff in ≈ 22–25 cycles regardless of level               |
| CG + two-grid preconditioner    | Plateaus at ≈ 10⁻¹¹ in 9–11 iterations (residual stalls — see report §5)  |

The conjugate gradient stagnation is attributed to error introduced when computing vₘ = B⁻¹ rₘ; switching to a flexible CG variant would likely fix it.

## Repository layout

| File                                          | Role                                                          |
|-----------------------------------------------|---------------------------------------------------------------|
| [main.c](main.c)                              | Entry point — pick which method to run                        |
| [prob.c](prob.c) / [prob.h](prob.h)           | Builds the fine-grid CSR matrix A and right-hand side b       |
| [prob_grossier.c](prob_grossier.c)            | Builds the coarse-grid matrix A_c                             |
| [indice.c](indice.c)                          | Index helpers for the 2D ↔ 1D mapping                         |
| [gs_prog.c](gs_prog.c)                        | Forward Gauss-Seidel sweep (pre-smoother)                     |
| [gs_reg.c](gs_reg.c)                          | Backward Gauss-Seidel sweep (post-smoother)                   |
| [restriction.c](restriction.c)                | Full-weighting restriction R = ¼ Pᵀ                           |
| [prolongation.c](prolongation.c)              | Bilinear prolongation P                                       |
| [twogrid.c](twogrid.c)                        | Two-grid iterative scheme                                     |
| [multigrid.c](multigrid.c)                    | Recursive V-cycle multigrid                                   |
| [conjugate.c](conjugate.c)                    | Preconditioned conjugate gradient                             |
| [umfpack.c](umfpack.c)                        | UMFPACK wrapper for direct solves                             |
| [norme.c](norme.c)                            | Residual norm                                                 |
| [time.c](time.c)                              | CPU / wall clock timing helpers                               |
| [plot_potentiel.c](plot_potentiel.c)          | gnuplot output of the potential field                         |
| [Makefile](Makefile)                          | Build configuration                                           |

## Building

The project depends on:
- **SuiteSparse / UMFPACK** (with CHOLMOD, AMD, CAMD, COLAMD, CCOLAMD, METIS-4.0 and `SuiteSparse_config`)
- **BLAS** and **LAPACK** (`-lblas -llapack`)
- **gnuplot** for visualising the potential

> **SuiteSparse is not bundled with this repository — you must install it yourself.** Grab it from <https://people.engr.tamu.edu/davis/suitesparse.html> (or your distro's package manager / `git clone` of the official repo) and build the required modules so that the static libraries land at the paths the Makefile expects:
>
> ```
> SuiteSparse/UMFPACK/Lib/libumfpack.a
> SuiteSparse/CHOLMOD/Lib/libcholmod.a
> SuiteSparse/AMD/Lib/libamd.a
> SuiteSparse/CAMD/Lib/libcamd.a
> SuiteSparse/COLAMD/Lib/libcolamd.a
> SuiteSparse/CCOLAMD/Lib/libccolamd.a
> SuiteSparse/metis-4.0/libmetis.a
> SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a
> ```
>
> Place the `SuiteSparse/` directory next to the project sources, or adjust the `L1…L8` paths at the top of the [Makefile](Makefile) to point at your installation.

Compile with:

```sh
make
```

This produces the `main` executable. To clean:

```sh
make clean
```

## Running

```sh
./main
```

Which method runs is controlled inside [main.c](main.c): by default only the two-grid block is active. Uncomment the relevant section to run the V-cycle multigrid or the preconditioned CG.

Grid parameters (geometry, refinement level, boundary values, number of cycles) are set as constants at the top of `main`:

```c
int cycle = 50;            // max number of iterations
int level = 2;             // refinement level (h = 2^-level / 2000)
double Lx = 0.005, Ly = 0.004;        // domain size
double x_i = 0.003, x_f = 0.004;      // inner conductor bounds (x)
double y_i = 0.001, y_f = 0.003;      // inner conductor bounds (y)
double u_i = -5.0, u_e = 5.0;         // inner / outer Dirichlet values
```

> **Note** — when running the conjugate gradient path, comment out the gnuplot call inside [multigrid.c](multigrid.c) (around line 209) so the preconditioner does not pop up a plot at every iteration.

## Output

Each run prints the residual norm per iteration and produces a gnuplot rendering of the electrostatic potential on the fine grid.
