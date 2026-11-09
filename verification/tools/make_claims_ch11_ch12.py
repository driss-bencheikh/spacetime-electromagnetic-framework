#!/usr/bin/env python3
"""Write the claim-level manifests of Chapters 11 and 12 (Phase 9C, 2 November 2026).

    python3 verification/tools/make_claims_ch11_ch12.py [--runs reproducibility/selected-results/v1.0.1]

For every check of verification/ch11_curvature_gravity and verification/ch12_causality it writes one row of CLAIMS.tsv in that directory: the manuscript
label, the source code, the mathematical assumptions, the method, the expected output, the observed output (the PASS/FAIL line of the most
run record given by --runs, by default reproducibility/selected-results/v1.0.1), the software dependencies, the verification status and the remaining limitations.
Status vocabulary: exactly verified (an identity over the stated general parameters), verified on the stated example, numerically supported,
conditional, obstructed (a no-go that is itself verified), unresolved. A PASS is never read as more than its coded assertion.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
RUNS = ROOT / (sys.argv[sys.argv.index("--runs") + 1] if "--runs" in sys.argv else "reproducibility/selected-results/v1.0.1")
DEPS = "Python 3, sympy (numpy for c1-4, c1-5)"
X = "exactly verified"
META = {
    "ch11_curvature_gravity": [
        ("g1-1", "prop:curvop (conventions)", "g1_curvature_operator.py", "signature (+,-,-,-); frame basis (e^23,e^31,e^12,e^01,e^02,e^03); eps_0123 = +1",
         "exact symbolic", "induced metric = eta; Hodge star = G", X, "fixes the conventions of the whole section"),
        ("g1-2", "prop:curvop", "g1_curvature_operator.py", "R_abcd antisymmetric in each pair; general 36- and 21-entry operators",
         "exact symbolic", "pair symmetry <=> eta-self-adjoint (dim 21); Tr R = R/2; cyclic sum = -(1/2)Tr(GR) eps", X, "algebraic; no metric"),
        ("g1-3", "thm:phase", "g1_curvature_operator.py", "R eta-self-adjoint, theta real", "exact symbolic (trig identities)",
         "Ric(R_theta) formula; Lanczos double dual lemmas", X, "algebraic identity of a given curvature; not a field equation"),
        ("g1-4", "ex:dS", "g1_curvature_operator.py", "constant curvature k", "exact symbolic", "R = k id, R = 12k, Ric(GRG) = -3k g", X,
         "the sign of the induced metric fixes R = +k id"),
        ("g1-5", "cor:phasesel (ii),(iii)", "g1_curvature_operator.py", "R eta-self-adjoint", "exact symbolic",
         "Tr(G R_theta) formula; Q = -(1/8) Gauss-Bonnet density", X, "(i) uses the contracted Bianchi identity (inherited)"),
        ("g1-6", "thm:helmholtz", "g1_curvature_operator.py", "linearisation at g = eta on plane waves; general symmetric h, q",
         "exact symbolic", "antisymmetric part (1/2)cos^2(theta)[(kqk)tr h - (khk)tr q]", "conditional",
         "exact for the linearised symbol; the no-Lagrangian conclusion uses the inherited necessity of the Helmholtz conditions; no action is constructed"),
        ("g1-7", "prop:holst (i)-(iii)", "g1_curvature_operator.py", "orthonormal coframe at a point; F^ab = (1/2)R^ab_cd e^c e^d", "exact symbolic",
         "Palatini = Tr R vol; Holst = -2Tr(GR) vol; phase rotation identity", X, "formal correspondence of Lagrangian densities at a point"),
        ("g1-8", "prop:holst (iii)-(iv)", "g1_curvature_operator.py", "Holst normalisation (1/2)eps eeF - gamma^-1 eeF", "exact symbolic",
         "gamma = -tan(2 theta); Holst 3-forms vanish iff Tr(GR) = 0", X,
         "the variation (connection and coframe equations) is not recomputed; dynamical equivalence with GR is the inherited Holst result"),
        ("g2-1", "prop:kk", "g2_metric_models.py", "flat base; general A_mu(x) and A = Bx dz", "exact symbolic curvature",
         "signature (3,3); R^u_{v mu nu}|Sigma = (dA)_{mu nu}; R^mu_{nu rho sigma}|Sigma = 0", X, "geometric calculation of a given metric"),
        ("g2-2", "thm:ein6", "g2_metric_models.py", "Friedmann base h = diag(1,-t^2,-t^2,-t^2), A = z dt + x^2 dy/2, r > 0", "exact symbolic curvature",
         "all frame components of G^(6) and R^(6) as stated", "verified on the stated example",
         "the general proof is O'Neill's formulas (inherited); geometric calculation, no field equation solved"),
        ("g2-3", "cor:ein6", "g2_metric_models.py", "Einstein-Maxwell Lagrangian a1 R + a2 F^2", "exact symbolic",
         "reduced action pi r_c^2 R - pi r_c^4 F^2/8; delta(sqrt|h| F^2) = 2 sqrt|h| W; -W_00 = (E^2+B^2)/2", X,
         "obstruction: the reduced ratio has the wrong sign; rigidity follows from g2-2"),
        ("g2-4", "ex:nullwave", "g2_metric_models.py", "general H(U,X,Y), a(U)", "exact symbolic curvature",
         "only G_UU = -(1/2)(H_XX + H_YY - l^2 a^2); vacuum solution; negative-energy null wave", X, "a solution of the assumed six-dimensional vacuum equation"),
        ("g2-5", "two-moduli formula (after ex:nullwave)", "g2_metric_models.py", "flat base; alpha, beta, A_z functions of (t,x)", "exact symbolic curvature",
         "G_psipsi - G_chichi formula", "verified on the stated class", "the class restricts the dependence to (t,x)"),
        ("g2-6", "off-slice determinant; ex:landau", "g2_metric_models.py", "general symmetric h; uniform field", "exact symbolic",
         "matrix determinant lemma; R = -2k^2/(1-k^2x^2)^2; circulation -pi b R_c^2", X, "the bound on sum B_i^2 is Cauchy-Schwarz (analytic)"),
        ("g2-7", "thm:ghostvec", "g2_metric_models.py", "flat base; sigma = +-1", "exact symbolic curvature",
         "R^(6) = -(sigma/4) F.F and -sigma r^2 F^2/4", "obstructed",
         "verified no-go for the normal-plane vectors; scalar curvature of given metrics, not a derived reduced theory"),
        ("g2-8", "rem:signconv", "g2_metric_models.py", "Friedmann example", "exact symbolic",
         "R, W odd and G, F^2 even under g -> -g", X, "the conclusion (mismatch convention independent) is a two-line argument"),
        ("g3-1", "prop:parahermitian with prop:twoK", "g3_projector_transport.py", "8x8 realisation", "exact",
         "K = (1/2)[varpi,.]; null eigenspaces; Omega skew", X, "para-Kaehler criterion inherited"),
        ("g3-2", "lem:JN", "g3_projector_transport.py", "8x8 realisation", "exact (sqrt3)", "J_N^2 = -1; commutes with W, e4, e5, e6", X, ""),
        ("g3-3", "holonomy on the bound orbits", "g3_projector_transport.py", "R^{3,3}", "exact symbolic",
         "exp(Phi R_N) formula; hGh^-1 = G iff Phi in 2pi Z; alpha_+(hX) formula", X, ""),
        ("g3-4", "prop:obstruct", "g3_projector_transport.py", "paired-mode projector", "exact",
         "Pi projector of rank 4; rank[J_N,Pi] = 4; [e4e5,Pi]", X, "spin-curvature formula combines with g2-1"),
        ("g3-5", "prop:adapted", "g3_projector_transport.py", "off-diagonal delta; K^2 = 1", "exact on rational test matrices",
         "[[Pi,delta],Pi] = -delta; tangent correction identities", X, "identities linear in delta, tested on explicit matrices"),
        ("g3-6", "adapted generator; expansion of Pi J_N Pi", "g3_projector_transport.py", "model four-plane <W,e4,e5,e6>; axial clock e1", "exact",
         "[J_+ - J_-,Pi] = 0; [(J_+-J_-)/2,e4] = (e6-e5)/sqrt3; blade expansion; six model bivectors fail to commute with Pi", X,
         "correction made: for the axial clock e1, e1e4 commutes with Pi (text clarified)"),
        ("g3-7", "ex:stueckelberg; ex:ghostPi", "g3_projector_transport.py", "paired-mode projector", "exact",
         "T symmetric with Tr T^2 = 1; geodesic condition; trace form signature (3,3)", X,
         "a negative direction is a ghost only once the projector kinetic term is adopted"),
    ],
    "ch12_causality": [
        ("c1-1", "thm:selector", "c1_selector_types.py", "8x8 realisation; six test planes", "exact",
         "(a^b)^2 = -det Gram; types; centraliser dim 7 for nondegenerate planes; orthogonality of coordinate complex structures", X,
         "planes beyond the six tested follow from the identity and Witt; the totally null stratum is c1-12"),
        ("c1-2", "omega_S pairing; thm:cone (i)-(ii)", "c1_selector_types.py", "8x8 realisation", "exact",
         "invariant pairing S+ x S-; m_A span Lambda^2(S+)*; J^2 = -1 on S+", X, ""),
        ("c1-3", "thm:cone (iii)-(v)", "c1_selector_types.py", "general psi in S+; six planes", "exact symbolic",
         "v_C in N-perp, null, image dim 3; definite/indefinite components; v_J = -j", X,
         "nappe-filling is the Klein correspondence (inherited); an algebraic cone is not a causal theory"),
        ("c1-4", "thm:selector3", "c1_selector_types.py", "512x512 signed permutations", "exact (integer); random spinors of S+ numerical",
         "J, Omega, Gamma relations; all 3060 four-planes; 97 = 6 + 91; beta pairings; untwisted square not null on S+",
         "exactly verified (numerical part: 5 samples of S+)", ""),
        ("c1-5", "Dirac current at n = 3 (after thm:selector3)", "c1_selector_types.py", "512x512", "exact matrix identities; 20 random spinors; product states",
         "gamma^0 gamma^k symmetric anticommuting involutions; j future causal; null on product states with a chiral first factor, (2,0,0,0) for |000> + |100>",
         X, "uses the positive inner product, as the text states"),
        ("c1-6", "prop:types", "c1_selector_types.py", "restriction to S+ = sl(4,R)", "exact", "invariants and characteristic polynomials; opposite orientations", X, ""),
        ("c1-7", "thm:degree (i)-(iii)", "c1_selector_types.py", "test elements of the proof", "exact",
         "stated invariants; e = g = 0; B_w invariants and non-simplicity", "exactly verified (inputs of the proof)", "the bounded-below argument is the written proof"),
        ("c1-8", "thm:degree (iv)", "c1_selector_types.py", "F on so(3,3)", "exact symbolic", "F(J) = 0; Hessian rank 1", X, ""),
        ("c1-9", "thm:flat", "c1_selector_types.py", "8x8 realisation", "exact symbolic",
         "n commutes, n^2 = 0; centraliser dim 3; e^{tK} conjugation; invariant forms dim 2, signature (3,3)", "obstructed (verified no-go)",
         "invariant classification, not a dynamical selection mechanism"),
        ("c1-10", "prop:kineticB", "c1_selector_types.py", "trace form = multiple of Killing form", "exact",
         "signature (9,6); Goldstone 2 rotations + 6 boosts; induced (6,2)", "obstructed (verified no-go)",
         "signature of a quadratic kinetic form; ghost freedom of an interacting theory not addressed"),
        ("c1-11", "prop:positivealt", "c1_selector_types.py", "positive form tr(XY^T)/8", "exact symbolic",
         "invariance algebra so3_T + so3_S; Hessian 8, 0^2, 2^12; boosts not symmetries", X, ""),
        ("c1-12", "thm:selector (iii), degenerate stratum", "c1_selector_types.py", "totally null plane <e1+e4, e2+e5>", "exact",
         "C0^2 = 0; centraliser in so(3,3) of dimension 9, against 7 for C = J", X, "one representative of the degenerate strata"),
        ("c2-1", "prop:toy", "c2_dynamics_brackets.py", "free scalar on R^{3,3}", "exact symbolic", "dispersion; energy density; growth beats Sobolev weights", X,
         "Hadamard ill-posedness inherited from the growth"),
        ("c2-2", "prop:orient", "c2_dynamics_brackets.py", "R^{3,3}", "exact", "connecting path with eta = 1", X, ""),
        ("c2-3", "ex:tachyons; 9C.8", "c2_dynamics_brackets.py", "compact surplus times, m R = 3", "exact counting",
         "92, 412, 1652 tachyonic lattice points; stability depends on |k|; 92 tachyonic and 28 non-tachyonic nonzero modes removed", X,
         "free spectrum only; no interacting stability"),
        ("c2-4", "CORE-07 statement (prop:toy discussion)", "c2_dynamics_brackets.py", "non-compact surplus times", "exact",
         "a stable mode with q != 0 is removed by p_u = p_v = 0", X, "constraints sufficient, not minimal"),
        ("c2-5", "prop:constraints; thm:wlgauge (brackets)", "c2_dynamics_brackets.py", "general h^{mu nu}(x), A_mu(x)", "exact symbolic Poisson brackets",
         "first-class algebra with structure functions; {l,H} = 0; H_red", X, "test worldline in a background"),
        ("c2-6", "9C.9 (bracket properties, weak preservation)", "c2_dynamics_brackets.py", "canonical bracket", "exact symbolic",
         "antisymmetry, Leibniz, Jacobi; d p_u/d tau ~ 0", X, "constraint preservation for interacting fields not established"),
        ("c2-7", "thm:wlgauge (gauge invariance)", "c2_dynamics_brackets.py", "general h, A", "exact symbolic", "delta L = 0 identically", X,
         "the gauge principle is chosen; a consistent constrained Hamiltonian theory of interacting fields is not constructed"),
        ("c2-8", "prop:wltau", "c2_dynamics_brackets.py", "general h, A; action with N (H - m^2/2)", "exact symbolic",
         "varying N gives the mass shell 2H = m^2 (derived, not inserted); x' = N h k; d tau_h = |N| m d tau; reparametrisation", X,
         "a universal evolution parameter remains an open hypothesis"),
        ("c2-9", "prop:isometryonly", "c2_dynamics_brackets.py", "polar normal coordinates", "exact symbolic", "l = -p_phi; p_u^2+p_v^2; Bessel modes", X, ""),
    ],
}


def observed(stem_dir):
    lines = {}
    prefix = stem_dir + "__"
    for f in RUNS.glob(prefix + "*.out"):
        for l in f.read_text().splitlines():
            m = re.match(r"^(PASS|FAIL) ([a-z]\d+-\d+):", l)
            if m:
                lines[m.group(2)] = m.group(1)
    return lines


def main():
    for d, rows in META.items():
        obs = observed(d)
        out = ["check\tmanuscript_label\tsource_code\tassumptions\tmethod\texpected_output\tobserved_output\tdependencies\tstatus\tlimitations"]
        for cid, lab, src, ass, meth, exp, stat, lim in rows:
            out.append("\t".join([cid, lab, f"verification/{d}/{src}", ass, meth, exp, obs.get(cid, "not run"), DEPS, stat, lim]))
        (ROOT / "verification" / d / "CLAIMS.tsv").write_text("\n".join(out) + "\n")
        print(d, len(rows), "rows;", sum(1 for r in rows if obs.get(r[0]) == "PASS"), "observed PASS")


if __name__ == "__main__":
    main()
