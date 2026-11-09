# Chapter 5 checks — The Clifford Tower and Constituent Duplication

The scripts of this directory are run by `verification/run_all.py`. Each prints one `PASS <id>: <claim>` or `FAIL <id>: <claim>` line per check
and exits with status zero only if every check passes.

Corrections audit (v1.0.1):

```bash
python3 verification/ch05_clifford_tower/a9_blockgamma_certificates.py
```

a9-1 and a9-2, exact integer arithmetic on the 512-dimensional module (about two minutes): the two restricted identities
{Q_mu, Q_nu}|_D = -72 eta_mu nu and (1/4)[G_a, G_b] = Delta_ab|_D in the proof of thm:blockgamma, on the 160-dimensional Dirac block D.
