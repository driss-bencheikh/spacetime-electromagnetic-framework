#!/usr/bin/env python3
"""z5: Phase 3 checks for the classical section e7 (CLASS-01 to CLASS-05, CLASS-07), 1 November 2026.

    python3 verification/ch10_classical_consistency/z5_classical_adjudication.py

Prints one 'PASS z5-<id>: <claim>' or 'FAIL z5-<id>: <claim>' line per claim and exits 0 iff all pass.

Purpose. The checks behind the Phase 3 corrections of sections/e7_classical.tex. It imports nothing from the other verification
scripts. Astronomical inputs are quoted from the references named next to them in the section (Read 2014, Planck 2018,
Rainville et al. 2005, AME2020, Minami and Komatsu 2020, Eskilt and Komatsu 2022); the script does the arithmetic only, and a PASS
line does not validate the observations or the reading they are applied to. Physical constants: scipy.constants (CODATA).
Symbolic statements use sympy in exact arithmetic; the numbers are floating point with the tolerances in the PASS lines.
"""
import math
import sys

import sympy as sp
from scipy import constants as C

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS z5-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL z5-{cid}: {text}")
    sys.stdout.flush()


EV = C.e                      # J per eV
CM3 = 1e-6                    # m^3 per cm^3
GYR = 1e9 * C.Julian_year
MPC = 1e6 * C.parsec
M_SUN = 1.98841e30            # kg, IAU 2015 nominal solar mass (GM_sun / G)

# ------------------------------------------------------------------------- CLASS-01: like with like (z5-1)
u_B = {b: (b * 1e-10) ** 2 / (2 * C.mu_0) / EV * CM3 for b in (3, 6)}          # eV/cm^3 for 3 and 6 microgauss
read_lo, read_hi = 0.2e9, 0.56e9                                                # eV/cm^3, Read 2014 (Table 4, conclusions)
conv = 0.008 * M_SUN / C.parsec ** 3 * C.c ** 2 / EV * CM3                     # Read's conversion 0.008 Msun/pc^3 = 0.3 GeV/cm^3
r_lo, r_hi = u_B[3] / read_hi, u_B[6] / read_lo
T0 = 2.72548                                                                    # K, Fixsen 2009
a_rad = 4 * C.sigma / C.c
u_cmb = a_rad * T0 ** 4 / EV * CM3                                              # eV/cm^3
u0_cmb = 2 / 3 * u_cmb
rho_crit_h2 = 3 * (100e3 / MPC) ** 2 / (8 * math.pi * C.G)                     # kg/m^3 per h^2
u_cdm = 0.120 * rho_crit_h2 * C.c ** 2 / EV * CM3                               # Planck 2018: Omega_c h^2 = 0.120
ok = (math.isclose(u_B[3], 0.2235, rel_tol=2e-3) and math.isclose(u_B[6], 0.894, rel_tol=2e-3)
      and math.isclose(conv, 0.30e9, rel_tol=2e-2) and 3.9e-10 < r_lo < 4.1e-10 and 4.4e-9 < r_hi < 4.5e-9
      and math.isclose(u0_cmb, 0.1737, rel_tol=2e-3) and math.isclose(u_cdm, 1.26e3, rel_tol=5e-3)
      and 1.35e-4 < u0_cmb / u_cdm < 1.45e-4)
report("1", ok, f"CLASS-01: local galactic magnetic energy density {u_B[3]:.3f}-{u_B[6]:.3f} eV/cm^3 (3-6 muG) is "
       f"{r_lo:.2e} to {r_hi:.2e} of the local dark-matter density 0.2-0.56 GeV/cm^3 (Read 2014; his conversion "
       f"0.008 Msun/pc^3 = {conv / 1e9:.3f} GeV/cm^3 reproduced); cosmological means: CMB <u0> = (2/3) a T0^4 = {u0_cmb:.4f} "
       f"eV/cm^3 is {u0_cmb / u_cdm:.2e} of the mean CDM energy density {u_cdm:.0f} eV/cm^3 (Omega_c h^2 = 0.120)")

# ------------------------------------------------------- CLASS-02: smeared fields, commutators, isotropy (z5-2)
x, y, z, a = sp.symbols("x y z a", real=True)
f = sp.exp(-(x ** 2 + y ** 2 + z ** 2) / 2)
g = sp.exp(-((x - a) ** 2 + y ** 2 + z ** 2) / 2)
lim = (-sp.oo, sp.oo)
same = sp.integrate(f * sp.diff(f, x), (x, *lim), (y, *lim), (z, *lim))
shifted = sp.simplify(sp.integrate(g * sp.diff(f, x), (x, *lim), (y, *lim), (z, *lim)))
th, ph = sp.symbols("theta phi", real=True)
khat = [sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)]
avg = lambda expr: sp.simplify(sp.integrate(sp.integrate(expr * sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2 * sp.pi)) / (4 * sp.pi))
iso = all(sp.simplify(avg(sp.KroneckerDelta(i, j) - khat[i] * khat[j]) - sp.Rational(2, 3) * sp.KroneckerDelta(i, j)) == 0
          for i in range(3) for j in range(3))
odd = all(avg(k) == 0 for k in khat)
zz = sp.symbols("z1:7", real=True)
lam = sp.symbols("lambda", positive=True)
zc = [zz[i] + sp.I * zz[i + 3] for i in range(3)]
t_of = lambda w: sp.Abs(sum(c * c for c in w)) ** 2 / (sum(sp.Abs(c) ** 2 for c in w)) ** 2
homog = sp.simplify(t_of([lam * c for c in zc]) - t_of(zc)) == 0
ok = same == 0 and shifted != 0 and sp.simplify(shifted.subs(a, 1)) != 0 and iso and odd and homog
report("2", ok, f"CLASS-02: for one Gaussian test function f, int f d_x f d^3x = {same} (so [E_i(f), B_j(f)] ~ eps_ijk int f d_k f "
       f"vanishes and the six commonly smeared components commute); for f and its translate by a along x the integral is "
       f"{shifted}, nonzero for a != 0; the sphere average of delta_ij - k_i k_j is (2/3) delta_ij and that of k_l is 0 "
       f"(isotropic smeared covariance, no E-B cross term for spherical smearing); t = |z.z|^2/|z|^4 is invariant under z -> lambda z, "
       f"so the law of t does not depend on the variance (exact, sympy)")

# --------------------------------------------------------------- CLASS-04: the Rainville sensitivity budgets (z5-3)
Dsym, Esym = sp.symbols("D E", positive=True)
dm = Esym + Dsym                                          # Delta m c^2 = E + D for a defect D of the initial two-body system
sign_ok = sp.simplify(1 - dm / Esym + Dsym / Esym) == 0   # 1 - Delta m c^2 / E = -D / E
u_MeV = C.physical_constants["atomic mass constant energy equivalent in MeV"][0]
m28, m29 = 27.97692653465, 28.97649466490                 # AME2020 atomic masses (u)
mn_u = C.physical_constants["neutron mass in u"][0]
Sn = (m28 + mn_u - m29) * u_MeV                           # MeV
M0 = (m28 + mn_u) * u_MeV
mu12 = m28 * mn_u / (m28 + mn_u) * u_MeV
central, sigma = 1.4e-7, 4.4e-7                           # D/E = -(1 - Delta m c^2/E) = (1.4 +- 4.4)e-7, Rainville et al. 2005
budgets = {"1 sigma alone": sigma, "central + 1 sigma": central + sigma, "central + 2 sigma": central + 2 * sigma}
res = {k: (v * Sn * 1e6, math.sqrt(2 * v * Sn / mu12)) for k, v in budgets.items()}   # (eV, rad)
ok = (sign_ok and math.isclose(Sn, 8.4736, abs_tol=2e-4)
      and math.isclose(res["1 sigma alone"][0], 3.73, abs_tol=0.01) and 9.0e-5 < res["1 sigma alone"][1] < 9.15e-5
      and math.isclose(res["central + 1 sigma"][0], 4.91, abs_tol=0.01) and 1.03e-4 < res["central + 1 sigma"][1] < 1.05e-4
      and math.isclose(res["central + 2 sigma"][0], 8.64, abs_tol=0.01) and 1.37e-4 < res["central + 2 sigma"][1] < 1.39e-4)
report("3", ok, "CLASS-04: E = Delta m c^2 - D gives 1 - Delta m c^2/E = -D/E (exact), so the central value -1.4e-7 lies on the side "
       f"of a defect; with S_n = {Sn:.4f} MeV and mu_12 = {mu12:.2f} MeV the defect budgets and angles theta = (2D/mu_12)^(1/2) are: "
       + "; ".join(f"{k}: {d:.2f} eV, {t:.2e} rad" for k, (d, t) in res.items())
       + " (sensitivity arithmetic on cited inputs, conditional on the reading and on the reaction assumptions)")

# --------------------------------------------------------------------- CLASS-05: the updated birefringence bound (z5-4)
t_age = 13.787 * GYR                                      # Planck 2018
t_lss = 0.38e6 * C.Julian_year                            # recombination, about 0.38 Myr
H0 = 67.4e3 / MPC
k_mk = math.radians(0.35 + 2 * 0.14) / t_age              # Minami and Komatsu 2020, |central| + 2 sigma
k_ek = math.radians(0.342 + 2 * 0.094) / t_age            # Eskilt and Komatsu 2022, upper error 0.094 deg
k_c = math.radians(0.342) / t_age                         # the rate a constant-kappa model would need for the central value
ok = (2.4e-20 < k_mk < 2.6e-20 and 2.0e-20 < k_ek < 2.2e-20 and 1.3e-20 < k_c < 1.45e-20
      and (t_age - t_lss) / t_age > 0.9999 and 80 < H0 / k_mk < 90 and H0 / k_ek > 95)
report("4", ok, f"CLASS-05: under the coupling model (duality rotation of free light at constant rate in cosmic time), "
       f"|beta| < 0.63 deg (Minami-Komatsu 2020) gives kappa < {k_mk:.2e} 1/s (H0/kappa_max = {H0 / k_mk:.0f}), |beta| < 0.530 deg "
       f"(Eskilt-Komatsu 2022, 0.342 + 2 x 0.094) gives kappa < {k_ek:.2e} 1/s, and the central 0.342 deg corresponds to "
       f"kappa = {k_c:.2e} 1/s; the flight time from last scattering differs from t0 by {t_lss / t_age:.1e} (arithmetic on cited inputs)")

# ------------------------------------------------------- CLASS-07: what the clock and EP tests cannot see (z5-5)
mp, me, c_t, thA, thB, th0 = sp.symbols("m_p m_e c theta_A theta_B theta_0", positive=True)
ratio_common = sp.simplify((c_t * mp) / (c_t * me) - mp / me) == 0
r = lambda thv: 1 / sp.cos(thv)                           # m_g / m_i with m_g = mu, m_i = m_4 = mu cos(theta)
eta = 2 * (r(thA) - r(thB)) / (r(thA) + r(thB))
eta_common = sp.simplify(eta.subs(thB, thA)) == 0
eta_series = sp.simplify(sp.series(eta.subs({thA: th0 * sp.Symbol("u"), thB: th0 * sp.Symbol("v")}), th0, 0, 3).removeO()
                         - th0 ** 2 * (sp.Symbol("u") ** 2 - sp.Symbol("v") ** 2) / 2) == 0
m4A, m4B = sp.symbols("m_4A m_4B", positive=True)        # m_g = m_i = m_4: r = m_4/m_4 = 1 for every body
rA, rB = m4A / m4A, m4B / m4B
eta_m4 = sp.simplify(2 * (rA - rB) / (rA + rB)) == 0
ok = ratio_common and eta_common and eta_series and eta_m4
report("5", ok, "CLASS-07: a common factor c(t) on all masses cancels in m_p/m_e; with inertia m_4 = mu cos(theta) and gravity mu, "
       "the Eotvos parameter of two bodies is (theta_A^2 - theta_B^2)/2 + O(theta^4) and vanishes for a common angle; if inertia "
       "and gravity both couple to m_4 it vanishes identically (exact, sympy)")

print()
print(f"z5: {NPASS[0]} PASS, {len(FAILS)} FAIL" + (f" (failed: {FAILS})" if FAILS else ""))
sys.exit(0 if not FAILS else 1)
