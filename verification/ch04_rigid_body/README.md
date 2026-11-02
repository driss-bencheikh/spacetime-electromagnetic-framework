# Checks of Chapter 4, Classical Rigid-Body Realization (Phase 7)

Each script prints one `PASS <id>: <claim>` or `FAIL <id>: <claim>` line per check and exits 0 iff all pass. Requirements: Python 3, numpy,
scipy, sympy. The run of record is in `reproducibility/phase7_runs/`.

| script | checks | content | method |
|---|---|---|---|
| `rb1_rigid_body.py` | rb1-1 to rb1-8 | conventional rigid-body mechanics in the chapter's conventions: Lie–Poisson bracket on so(3)*, Euler's equations, Casimir and Noether conservation, inertia tensor, stability of steady rotations; the planar body SE(2), its pendulum and its conserved spatial momenta (Section 4.1) | exact (sympy) |
| `rb2_es_map.py` | rb2-1 to rb2-15 | the identification of the closure with a rigid body (Section 4.6, Theorem 4.53): bracket and equations, ranks and Casimirs, the pairing plane as se(2)*, spatial momenta, duality and Lorentz classes, indefinite Lorentz-class inertia, energy bounds, the rotor as the centre of U(3), units, dimension count, solvability; rb2-15 integrates the full closure and records the drift of its invariants | exact; rb2-15 numerical (DOP853, 1e-12) |
| `rb3_orbital.py` | rb3-1 to rb3-13 | the orbital extension (Section 4.7): Hooke's law in a pairing plane, isochrony, Bohlin and Levi-Civita maps, orbit-dependent Kepler coupling, the multi-plane obstruction, the spin map to the circular spin–orbit problem, resonance strengths; rb3-10 to rb3-13 integrate the Kepler benchmark, the Hooke–Kepler image, the circular spin–orbit problem against the closure, and the eccentric case (3:2 capture) | exact; rb3-9 quadrature; rb3-10 to rb3-13 numerical (DOP853, 1e-11 to 1e-12) |

`python3 verification/rigidbody/rb3_orbital.py --write-data` also writes the data of Figure 4.3(b) to `manuscript/figures/ch04/`.
Units of the numerical runs are dimensionless: GM = a = 1 for the orbits (mean motion 1, period 2π) and m = I0 = 1 for the closure. Each numerical
line states its initial data, tolerance and the measured drift of the invariants.
