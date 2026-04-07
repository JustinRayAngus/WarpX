# Bug fixes: JFNK with particle suborbits

Two related bugs were found in the implicit solver's handling of suborbit
particles during the JFNK linear stage. Both are in the non-mass-matrix JFNK
code path (`use_mass_matrices_jacobian = false`, `particle_suborbits = true`).

---

## Bug 1: Suborbit reclassification during JFNK linear stage

### Summary

Newton's method diverges because the finite-difference Jacobian approximation
captures an algorithmic discontinuity rather than the smooth derivative of the
nonlinear function `F(E)`.

### Background

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

### The bug

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
A particle that barely converges in Picard for field `E` may barely fail for
the perturbed field `E + ε*v`. When this happens during the JFNK evaluation:

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

Note that `ImplicitPushXPSubOrbits` already has the correct guard:

```cpp
// Source/Particles/Pusher/ImplicitPushPX.cpp, line ~1027
if (linear_stage_of_jfnk) { convergence = true; }
```

The same protection was missing from `ImplicitPushXP`.

### The fix (commit 6796c9f)

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

---

## Bug 2: Suborbit particles not advanced during JFNK linear stage

### Summary

Particles already flagged as suborbit (`nsuborbits > 1`) from the nonlinear
evaluation are not properly handled during the JFNK linear stage. They are
pushed through `ImplicitPushXP` with a single-step Picard that doesn't converge,
and the Bug 1 fix prevents them from being counted as unconverged, so
`ImplicitPushXPSubOrbits` is never invoked. Their deposited current comes from
a non-converged Picard result, corrupting the Jacobian-vector product.

### How PICNIC handles this

In PICNIC (the predecessor code), suborbit particles live in a separate
container (`m_data_suborbit`). The function `addSubOrbitJ` calls
`advanceSubOrbitParticlesAndSetJ` for these particles at the end of `preRHSOp`,
**regardless** of whether the call is from the linear or nonlinear stage and
whether mass matrices are used. This ensures suborbit particles are always
properly advanced via sub-stepping with their current correctly deposited.

During the PICNIC linear stage, suborbit Picard iteration limits are doubled
for a safety buffer, and if convergence still fails, the suborbit count is
frozen (not increased) -- the particle is simply transferred to a temporary
list and its approximate result is used. This matches the WarpX guard in
`ImplicitPushXPSubOrbits`.

### The bug in WarpX

In `PhysicalParticleContainer::Evolve`, the control flow is:

```cpp
if (evolve_suborbit_particles_only) {       // mass-matrix linear stage
    FindSuborbitParticles(...);              // identifies suborbits, zeros weights
} else {
    ImplicitPushXP(...);                    // pushes ALL particles
}
```

For the non-mass-matrix JFNK linear stage, `evolve_suborbit_particles_only` is
`false`, so all particles go through `ImplicitPushXP` -- including those with
`nsuborbits > 1`. The Bug 1 fix correctly prevents new suborbit classification
during the linear stage, but as a side effect, existing suborbit particles are
also not counted as unconverged. Their weights are not zeroed and
`ImplicitPushXPSubOrbits` is never called for them. They deposit current from
a non-converged single-step Picard result.

### The fix (commit 9bf6ecb)

Call `FindSuborbitParticles` before `ImplicitPushXP` when in the JFNK linear
stage:

```cpp
} else {
    if (implicit_options->linear_stage_of_jfnk) {
        FindSuborbitParticles(pti, offset, np_to_push,
                              num_unconverged_particles,
                              unconverged_indices, saved_weights);
    }
    long const saved_num_unconverged = num_unconverged_particles;
    ImplicitPushXP(pti, ...
                   num_unconverged_particles, unconverged_indices, saved_weights);
    num_unconverged_particles += saved_num_unconverged;
}
```

This ensures:

1. `FindSuborbitParticles` identifies particles with `nsuborbits > 1` and zeros
   their weights before `ImplicitPushXP` runs.
2. `ImplicitPushXP` pushes all particles; suborbit ones have `w = 0` so they do
   not deposit current.
3. `ImplicitPushXPSubOrbits` is triggered (via `num_unconverged_particles > 0`)
   and properly advances the suborbit particles using sub-stepping, depositing
   their current correctly.

This is safe because `ImplicitPushXPSubOrbits` reads the saved particle state
(`x_n`, `ux_n` attributes) directly, overriding any position/velocity changes
from `ImplicitPushXP`.

---

## Why this wasn't caught before

Both bugs exist on the `development` branch -- they are not specific to the
cylindrical mass matrices PR. However, the only test on `development` that uses
`particle_suborbits = true` is a mild 2D Cartesian case
(`inputs_test_2d_theta_implicit_jfnk_vandb_cropping`: 16x16 cells, 10 steps,
electron-positron). In that test, the physics is benign enough that few or no
particles actually fail Picard, so `nsuborbits` never reaches 2 and the bugs
are never triggered.

The cylindrical Z-pinch test case exposes the bugs because:

- The stagnation physics is extremely stiff (plasma compressed to near-zero
  radius).
- Hundreds of particles fail Picard convergence each step, creating a large
  suborbit population.
- Many particles sit right at the Picard convergence boundary, making them
  sensitive to the small `ε*v` perturbation.

## Test results

1D cylindrical dynamic pinch, restarting from checkpoint 110860, JFNK with
`particle_suborbits = true`, `use_mass_matrices_jacobian = false`, 4 MPI ranks:

**Before either fix:** Newton diverges immediately at step 110861.

**After Bug 1 fix only:** Step 110861 converges, but step 110862 diverges
(suborbit particles are not properly advanced during the linear stage).

**After both fixes:** All steps converge:

```
STEP 110861:
Newton: iteration =   0, norm = 7.75323e+11 (abs.), 1.00000e+00 (rel.)
Newton: iteration =   1, norm = 5.30653e+10 (abs.), 6.84428e-02 (rel.)
...
Newton: iteration =  12, norm = 2.67727e+03 (abs.), 3.45310e-09 (rel.)
Newton: exiting at iteration =  12. Satisfied relative tolerance 1e-08

STEP 110862:
Newton: iteration =   0, norm = 3.06359e+11 (abs.), 1.00000e+00 (rel.)
Newton: iteration =   1, norm = 4.19645e+10 (abs.), 1.36978e-01 (rel.)
...
Newton: iteration =  14, norm = 2.70871e+03 (abs.), 8.84162e-09 (rel.)
Newton: exiting at iteration =  14. Satisfied relative tolerance 1e-08

STEP 110863:
Newton: iteration =   0, norm = 9.38006e+10 (abs.), 1.00000e+00 (rel.)
...
Newton: iteration =   7, norm = 4.01454e+01 (abs.), 4.27986e-10 (rel.)
Newton: exiting at iteration =   7. Satisfied relative tolerance 1e-08
```
