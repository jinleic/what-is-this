# Convention probe — RESULTS (scratch, non-authoritative). 2026-08-29.

## Setup
m=6, E=F64 (prim 0x43 in log/exp tables), t=3, n=64 = |E| (full-field support),
k_expected = n − mt = 46, D = n − 2t−1 = 57. GF engine validated vs
python-flint fq_default (FQ_ZECH): 500/500 scalar mults, 100/100 poly mults,
8/8 irreducibility verdicts exact.

## Finding 1 (engine bug, fixed before Gate A): Rabin exponent
My first `pis_irreducible` raised Z to the power 2^d (F2-Frobenius) instead of
q^d with q = 2^m (Frobenius of E). Over E = F_{2^m} with m > 1 this is WRONG:
e.g. x²+x+1 has roots 58, 59 in F64 but the 2-Frobenius test called it
irreducible. Fixed to test Z^{q^d} ≡ Z (mod p) and gcd(Z^{q^{d/f}} − Z, p) = 1
for prime f | d; now agrees with flint's is_irreducible on all tested cases
(degree 3 over F64, degree 6 over F64, degree 8 over F4096). Without this fix
Gate A's "irreducible ⇒ squarefree" guard would have been silently violated.

## Finding 2 (construction subtlety): which y_i for interpolating F
The Goppa code is Γ(L,G) = {c ∈ F2^n : Σ c_i a_i^j / G(a_i) = 0, j < t}
(parity rows use 1/G(a_i), NOT 1/λ_i). The GRS embedding concerns the DUAL:
codewords evaluated with multiplier λ′_i = G(a_i)²/Π′(a_i) are F(a_i)·λ′_i
(Apon eq. (1); GIJS Fact 2.3). So for a public row v ∈ F2^n the interpolated
coordinate polynomial f must satisfy f(a_i) = v_i / λ′_i. First probe version
used 1/λ_i with λ_i = G(a_i)/Π′(a_i) (the WRONG GRS embedding, GIJS's other
multiplier λ for GRS_{n−t}) — that interpolates a DIFFERENT ambient code and
the Apon identity (14) need not hold. Gate A uses λ′ (n−2t embedding, D =
n−2t−1) exactly as Apon's Def 4 and Thm 1 require.

## Finding 3 (analytic, decisive for Gate A): F is unique; F′ ≠ 0 forced
The evaluation map {polys deg ≤ D}^k → (E^n)^k, F ↦ (F(a_i))_i, is injective
when D < n, and D = n−2t−1 < n always. So F with y_i = λ′_i F(a_i) is UNIQUE
per public generator matrix — no gauge freedom in choosing F. Then F′ = 0
identically is IMPOSSIBLE for a nonzero row: setting F′(a_i) = 0 in the
identity (15) ΠF′ + Π′F = G²F^{(2)} gives Π′(a_i)F(a_i) = 0, hence F(a_i) = 0
(Pi′(a_i) ≠ 0, since support points are distinct), contradicting y_i ≠ 0 for
rows with a pivot at i. (My first probe claimed "f′ = 0 identically" from a
run with the wrong multiplier per Finding 2; that run is VOID.) The Apon
hypothesis Δ_{p,q} = f_p f_q′ − f_q f_p′ ≠ 0 for some p < q is therefore a
GENUINE genericity condition to machine-check per instance, and Apon's §3.6
nondegeneracy caveat ("If every Δ is zero…") is a real residual degeneracy
risk to be measured (gate B), not a formality.

## Decision carried into src/ (consistent with pre_statement.md)
- Parity rows h_i = a_i^j/G(a_i); binary nullspace (numpy GF(2)) for the code.
- F built by Lagrange interpolation of v_i = y_i/λ′_i at a_i, coordinates
  forced by the pivoted generator rows; guards α–ε machine-checked.
- Support, G, seeds recorded; irreducibility by the corrected Rabin test
  (flint-agreed), and squarefreeness then free.
- If the Δ_{p,q} ≠ 0 guard fails on a seeded instance: log as DEGENERATE
  (finding for gate B), advance seed stream, regenerate — per pre_statement.