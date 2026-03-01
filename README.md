# FCE Proof of Concept — Sitnikov Three-Body Problem

A demonstration that post-processing corrections at Poincar&eacute; section crossings can dramatically improve the long-term accuracy of a standard ODE integrator (DOP853), tested on the Sitnikov restricted three-body problem.

## What It Does

Runs three integrations of the Sitnikov problem and compares them:

| Run | Description |
|-----|-------------|
| **Bare** | Standard DOP853 at coarse tolerances |
| **FCE-Corrected** | Same integrator, with corrections at every z=0 crossing |
| **Reference** | DOP853 at ultra-tight tolerances (the "truth") |

Two test cases with different dynamics:

**Test 1 — Circular primaries (e=0, integrable):**
The FCE detects that the crossing velocity is constant (energy conservation) and locks to it. Any drift is pure integrator error and gets corrected directly. Result: **437x reduction in energy drift**.

**Test 2 — Eccentric primaries (e=0.3, chaotic):**
No conservation law to exploit, so the FCE integrates each half-orbit at moderate precision (rtol=1e-8), using crossings as restart points. Result: **46,248x reduction in position error**, 3.5x cheaper than the ultra-fine reference.

## Quick Start

```bash
pip install numpy scipy matplotlib
python sitnikov_fce.py
```

Produces:
- Terminal output with all numerical metrics
- `test1_circular.png` — six-panel comparison for circular case
- `test2_eccentric.png` — six-panel comparison for eccentric case

Runtime: ~10 seconds total.

## Output Plots

Each plot has six panels:

| Panel | Shows |
|-------|-------|
| Trajectory | z(t) for bare vs FCE |
| Phase portrait | (z, vz) orbit structure |
| Energy drift | E(t) - E(0) over time |
| Energy drift (log) | Same, log scale |
| Position error / Poincar&eacute; section | Error vs reference (Test 2) or crossing velocities (Test 1) |
| FCE analysis | Corrections applied (Test 1) or energy at crossings (Test 2) |

## How the Two Strategies Work

### Strategy 1: Conservation-Law Correction (e=0)

For circular primaries, energy is exactly conserved. At z=0, this means |vz| is the same at every crossing. The FCE:

1. Records the first crossing velocity as the reference
2. At each subsequent crossing, measures the deviation from this constant
3. Corrects vz toward the reference value

This is a genuine physics-based correction — the coarse integrator does all the work, and the FCE fixes the drift at zero additional integration cost.

### Strategy 2: Adaptive-Precision Re-Integration (e>0)

For eccentric primaries, the dynamics are chaotic — no persistent patterns to exploit. The FCE:

1. Integrates each half-orbit at moderate precision (rtol=1e-8)
2. Uses z=0 crossings as natural restart points
3. The restart architecture prevents error amplification across orbits

The improvement comes primarily from higher precision per orbit, but the restart structure genuinely contributes by breaking Lyapunov error amplification chains. See the [full writeup](WRITEUP.md) for an honest discussion of what each strategy provides.

## Key Results

| Metric | Test 1 (e=0) | Test 2 (e=0.3) |
|--------|:------------:|:--------------:|
| Energy drift reduction | **437x** | 1.73x (misleading — see below) |
| Position error vs reference | — | **46,248x** |
| Energy error vs reference | — | **31,561x** |
| FCE wall time vs bare | 2.9x | 6.5x |
| FCE wall time vs reference | — | **3.5x cheaper** |

*Note: For e>0, "energy drift" |E(t)-E(0)| is mostly real physics (energy genuinely varies), not error. The meaningful metric is energy error vs reference: |E_FCE(t) - E_ref(t)|.*

## Files

| File | Description |
|------|-------------|
| `sitnikov_fce.py` | Complete implementation — physics, FCE engine, integration, plotting |
| `WRITEUP.md` | Full scientific writeup with equations and honest analysis |
| `test1_circular.png` | Output plot for e=0 test |
| `test2_eccentric.png` | Output plot for e=0.3 test |

## Dependencies

- Python 3.x
- NumPy
- SciPy
- Matplotlib

## The Sitnikov Problem

Two equal-mass bodies orbit their barycenter in a plane. A massless test particle moves perpendicular to that plane along the z-axis. The equation of motion:

$$\ddot{z} = -\frac{z}{(r(t)^2 + z^2)^{3/2}}$$

where $r(t)$ is the primary distance from the barycenter, determined by Kepler's equation. For circular orbits ($e=0$), $r=1/2$ is constant and the system is integrable. For eccentric orbits ($e>0$), the dynamics can be chaotic.
