# Bug: JFNK Jacobian discontinuity from suborbit reclassification

## Summary

When using the implicit theta scheme with JFNK (Jacobian-Free Newton-Krylov)
and `particle_suborbits = true`, Newton's method can diverge because the
finite-difference Jacobian approximation captures an algorithmic discontinuity
rather than the smooth derivative of the nonlinear function.

## Background

The implicit solver uses a Picard fixed-point iteration to advance each
particle. Particles that fail to converge within `max_particle_iterations` are
handed off to a suborbit push (`ImplicitPushXPSubOrbits`), which subdivides the
time step into smaller sub-steps. The transition is tracked by a per-particle
integer attribute `nsuborbits`: particles with `nsuborbits = 1` are pushed
normally; particles with `nsuborbits > 1` use the suborbit path.

The JFNK method approximates the Jacobian-vector product via a finite
difference:

```
J*v ≈ (F(E + ε*v) - F(E)) / ε
```

This requires evaluating the nonlinear function `F` at a slightly perturbed
field `E + ε*v`. For a good approximation, `F` must be smooth with respect to
`E`.

## The bug

In `ImplicitPushXP`, when a particle fails Picard convergence, the code
unconditionally sets `nsuborbits[ip] = 2` and flags the particle for the
suborbit push:

```cpp
// Source/Particles/Pusher/ImplicitPushPX.cpp, line ~631
if (max_iterations > 1 && !convergence) {
    if (nsuborbits) {
        nsuborbits[ip] = 2;  // <-- no guard for linear_stage_of_jfnk
    }
    ...
    amrex::Gpu::Atomic::Add(unconverged_particles_ptr, amrex::Long(1));
}
```

This happens during both the nonlinear evaluation *and* the JFNK linear stage.
The problem: a particle that barely converges in Picard for field `E` may
barely fail for the perturbed field `E + ε*v`. When this happens during the
JFNK evaluation:

1. For `F(E)`: the particle was pushed by `ImplicitPushXP` (normal path,
   `nsuborbits = 1`), depositing current via the standard algorithm.

2. For `F(E + ε*v)`: the same particle fails Picard, gets `nsuborbits = 2`,
   has its weight zeroed, and is re-pushed by `ImplicitPushXPSubOrbits`
   (suborbit path with subdivided time steps), depositing current via a
   different algorithm.

The finite difference `(F(E+ε*v) - F(E))/ε` then captures the difference
between two different deposition algorithms rather than the smooth response
`dF/dE`. This produces a garbage Jacobian approximation, causing Newton to take
wildly wrong steps and diverge.

A secondary effect: the `nsuborbits` change persists after the JFNK evaluation,
contaminating subsequent nonlinear evaluations and making `F(E)` non-repeatable
across Newton iterations.

Note that `ImplicitPushXPSubOrbits` already has the correct guard:

```cpp
// Source/Particles/Pusher/ImplicitPushPX.cpp, line ~1027
if (linear_stage_of_jfnk) { convergence = true; }
```

The same protection was missing from `ImplicitPushXP`.

## The fix

Guard the suborbit reclassification in `ImplicitPushXP` with
`!linear_stage_of_jfnk`:

```cpp
if (max_iterations > 1 && !convergence) {
    if (!linear_stage_of_jfnk) {
        if (nsuborbits) {
            nsuborbits[ip] = 2;
        }
        ...
        amrex::Gpu::Atomic::Add(unconverged_particles_ptr, amrex::Long(1));
    }
}
```

During the JFNK linear stage, a particle that fails Picard simply deposits
current using its approximate (not fully converged) Picard result. This is
acceptable because:

- The perturbation `ε*v` is small, so the Picard iterate is close to what
  the converged result would be.
- What matters for the Jacobian is that `F` is a smooth function of `E`, not
  that each particle is perfectly converged.
- The same algorithmic path is used for both `F(E)` and `F(E + ε*v)`,
  ensuring the finite difference captures the true derivative.

## Why this wasn't caught before

The bug exists on the `development` branch — it is not specific to the
cylindrical mass matrices PR. However, the only test on `development` that
uses `particle_suborbits = true` is a mild 2D Cartesian case
(`inputs_test_2d_theta_implicit_jfnk_vandb_cropping`: 16x16 cells, 10 steps,
electron-positron). In that test, the physics is benign enough that few or no
particles actually fail Picard, so `nsuborbits` never reaches 2 and the
discontinuity is never triggered.

The cylindrical Z-pinch test case exposes the bug because:

- The stagnation physics is extremely stiff (plasma compressed to near-zero
  radius).
- Hundreds of particles fail Picard convergence each step, creating a large
  suborbit population.
- Many particles sit right at the Picard convergence boundary, making them
  sensitive to the small `ε*v` perturbation.

## Test results

With the fix applied to the 1D cylindrical dynamic pinch at step 110861
(restarting from checkpoint 110860, JFNK with `particle_suborbits = true`,
`use_mass_matrices_jacobian = false`):

**Before fix:** Newton diverges immediately at step 110861.

**After fix:** Step 110861 converges with textbook quadratic Newton convergence:

```
Newton: iteration = 0, norm = 7.75e+11 (abs.), 1.00e+00  (rel.)
Newton: iteration = 1, norm = 1.40e+11 (abs.), 1.80e-01  (rel.)
Newton: iteration = 2, norm = 1.16e+10 (abs.), 1.50e-02  (rel.)
Newton: iteration = 3, norm = 3.60e+08 (abs.), 4.65e-04  (rel.)
Newton: iteration = 4, norm = 5.49e+04 (abs.), 7.08e-08  (rel.)
Newton: iteration = 5, norm = 3.08e+00 (abs.), 3.97e-12  (rel.)
Newton: exiting at iteration = 5. Satisfied relative tolerance 1e-08
```

Step 110862 still diverges and requires further investigation (likely related
to increasingly extreme stagnation conditions).
