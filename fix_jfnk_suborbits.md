# Bug fixes: implicit solver with particle suborbits

Bugs and architectural mismatches were found in the implicit solver's handling
of suborbit particles. Bugs 1 and 2 affect the non-mass-matrix JFNK code path
(`use_mass_matrices_jacobian = false`, `particle_suborbits = true`). Change 3
addresses the mass-matrix Jacobian code path (`use_mass_matrices_jacobian = true`,
`particle_suborbits = true`) by aligning WarpX's behavior with PICNIC.

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

---

## Change 3: Align mass-matrix Jacobian mode with PICNIC (commit 9c48322)

### Summary

When `use_mass_matrices_jacobian = true` and `particle_suborbits = true`, the
solver fails at stagnation (step 110862). Comparison with PICNIC revealed
architectural mismatches in how WarpX handles suborbit particles in the
mass-matrix Jacobian code path.

### Background: how PICNIC handles mass matrices with suborbits

In PICNIC, suborbit (fast-crossing) particles live in a separate container and
are structurally excluded from mass matrix deposition. The parameter
`quasi_freeze_particles_jacobian` (set to `true` when mass matrices are used)
controls the behavior:

- **Linear stage:** No particle push at all for converged particles. The
  Jacobian uses only the mass matrix formula `J = J0 + MM*(E - E0)`.
- **Nonlinear stage:** All particles are pushed (to update positions and
  velocities), but the residual current uses only the current from converged
  particles. Suborbit current is discarded because it is inconsistent with the
  mass matrix linearization.

Note: suborbit particles are already excluded from mass matrix deposition in
WarpX by a guard in `doVillasenorSigmaDeposition` (in
`MassMatricesDeposition.H`):

```cpp
if (nsuborbits && nsuborbits[ip] > 1) { return; }
```

This means no additional filtering is needed for mass matrix deposition itself.

### The problems in WarpX

Three architectural mismatches were identified:

#### 3a. Linear stage pushing suborbit particles

WarpX's linear stage with mass matrices was pushing suborbit particles and
adding their current on top of the mass matrix formula. PICNIC's linear stage
uses only the mass matrix formula with no particle push at all.

**Fix:** Remove the suborbit particle push from the linear stage in
`ImplicitSolver::PreRHSOp`. The Jacobian now relies solely on
`ComputeJfromMassMatrices(true)`:

```cpp
if (m_use_mass_matrices_jacobian && a_from_jacobian) {
    const bool J_from_MM_only = true;
    ComputeJfromMassMatrices( J_from_MM_only );
}
```

#### 3b. Nonlinear residual including suborbit current

WarpX's nonlinear stage used `CumulateJ()` which adds suborbit current to the
total current field. PICNIC uses only converged-particle current for the
nonlinear residual, discarding suborbit contributions.

**Fix:** Add `CopyJ0toJ()` which replaces `current_fp` with
`current_fp_non_suborbit` (converged particles only):

```cpp
else if (m_use_mass_matrices_jacobian && !a_from_jacobian) {
    m_WarpX->PushParticlesandDeposit(...);
    CopyJ0toJ();
}
```

#### 3c. Fast-crossing converged particles not flagged as suborbits

PICNIC's `transferFastParticles` identifies converged particles that cross more
than `max_grid_crossings` cells and moves them to the suborbit container. Their
mass matrix linearization is poor because the particle displacement spans too
many cells. WarpX had no equivalent mechanism.

**Fix:** Add a post-convergence check in `ImplicitPushXP` that flags particles
crossing more than `particle_max_grid_crossings` cells as suborbits:

```cpp
else if (convergence && nsuborbits && !linear_stage_of_jfnk
         && max_crossings > 0) {
    amrex::ParticleReal xp_new = 2.0_prt * xp - xp_n;
    int ncross = static_cast<int>(
        amrex::Math::abs(amrex::Math::floor((xp_new - xyzmin.x) * dinv.x)
                       - amrex::Math::floor((xp_n   - xyzmin.x) * dinv.x)));
    // ... similar for y, z dimensions ...
    if (ncross > max_crossings) {
        nsuborbits[ip] = 2;
        amrex::Gpu::Atomic::Add(unconverged_particles_ptr, amrex::Long(1));
    }
}
```

---

## Why this wasn't caught before

Bugs 1 and 2 exist on the `development` branch -- they are not specific to the
cylindrical mass matrices PR. However, the only test on `development` that uses
`particle_suborbits = true` is a mild 2D Cartesian case
(`inputs_test_2d_theta_implicit_jfnk_vandb_cropping`: 16x16 cells, 10 steps,
electron-positron). In that test, the physics is benign enough that few or no
particles actually fail Picard, so `nsuborbits` never reaches 2 and the bugs
are never triggered.

The Change 3 architectural mismatches also exist on the base branch but are
only visible when both `use_mass_matrices_jacobian = true` and
`particle_suborbits = true`.

The cylindrical Z-pinch test case exposes all issues because:

- The stagnation physics is extremely stiff (plasma compressed to near-zero
  radius).
- Hundreds of particles fail Picard convergence each step, creating a large
  suborbit population.
- Many particles sit right at the Picard convergence boundary, making them
  sensitive to the small `ε*v` perturbation.

---

## Test results

1D cylindrical dynamic pinch, restarting from checkpoint 110860, 4 MPI ranks.

### JFNK (`use_mass_matrices_jacobian = false`, `particle_suborbits = true`)

**Before either fix:** Newton diverges immediately at step 110861.

**After Bug 1 fix only:** Step 110861 converges, but step 110862 diverges
(suborbit particles are not properly advanced during the linear stage).

**After both fixes:** All steps converge through step 111000:

```
STEP 110861:
Newton: iteration =   0, norm = 7.75323e+11 (abs.), 1.00000e+00 (rel.)
Newton: iteration =   1, norm = 5.30653e+10 (abs.), 6.84428e-02 (rel.)
...
Newton: iteration =  12, norm = 2.67727e+03 (abs.), 3.45310e-09 (rel.)
Newton: exiting at iteration =  12. Satisfied relative tolerance 1e-08
```

### Full test matrix (all four combinations)

| `use_mass_matrices_jacobian` | `particle_suborbits` | Base branch | After Bugs 1+2 | After all fixes |
|---|---|---|---|---|
| `false` (JFNK) | `false` | Fails 110861 | Fails 110861 | Fails 110861 |
| `false` (JFNK) | `true` | Fails 110861 | **Converges** to 111000 | **Converges** to 111000 |
| `true` (MM) | `false` | **Converges** to 111000 | **Converges** to 111000 | **Converges** to 111000 |
| `true` (MM) | `true` | Fails 110862 | Fails 110862 | **~110866** then fails |

Notes:
- All runs use `pc_petsc.type = asm/lu` (default). JFNK + suborbits ON
  requires `asm/lu`; it fails with `lu`.
- The JFNK + suborbits OFF case fails because without suborbits, unconverged
  particles corrupt the Jacobian. This is expected behavior, not a bug.
- The MM + suborbits ON case progresses several steps past the previous
  immediate crash (110862 → ~110866) but eventually fails when the mass matrix
  linearization becomes too poor at peak stagnation.
- An epsilon sweep (`newton.jfnk_epsilon` from 1e-3 to 1e-10) and PETSc SNES
  (`implicit_evolve.nonlinear_solver = petsc_snes`) were also tested for the
  MM + suborbits ON case. Neither resolved the failure at step 110866.
