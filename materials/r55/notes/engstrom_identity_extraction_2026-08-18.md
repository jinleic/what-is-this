# Engstrom Mixed Identity — Primary-Source Extraction (2026-08-18)

Sources: Engstrom, arXiv:1002.4304v3 (TeX e-print `arxiv3.tex`), "A proof of the
McKay-Radziszowski subgraph counting conjecture"; McKay-Radziszowski 1997 (MR97),
`r55.pdf` (pdftotext -layout). Extraction agent: EngstromLibrarian. All numeric
claims MACHINE-VERIFIED by the extraction agent: the K=4 identity residual is
zero on all 33,867 labeled graphs n <= 6 + 300 random n in 7..10; the
reconstructed MR97 Thm 2.2 display is exact on 36,866 graphs for each m=1..6.

## 1. The normal-form theorem and the K=4 identity

Theorem 2.18 (\label{T2}, arxiv3.tex lines 426-440; Corollary 2.19 lines
481-498): for graphs J1..Jk, J'1..J'l there is a set J of graphs with at most
K = 1 + Σ|V(J_i)| + Σ|V(J'_j)| vertices and constants m_J such that for any G:

    Σ_v Π_i s(J_i, G_v^-) Π_j s(J'_j, G_v^+) = Σ_{J in J} m_J · j(J, G),

with j(J, G) = Σ_{φ inj(V(J),V(G))} (-1)^{|E(φ(J)) ∖ E(G)|} (degree-1 normal
form). Corollary 2.19: J without isolated vertices, m_J(n) polynomials,
|V(J)| + deg m_J(n) ≤ K.

The K=4 identity ("Theorem (Conjecture 1 in McKay and Radziszowski)",
arxiv3.tex lines 79-104), VERBATIM coefficients; n = |V(G)| (GLOBAL order):

    Σ_v [ p1(G_v^+) + p2(G_v^-) + p3(G_v^+, G_v^-) ] = 0,

    p1(X) = n(n−3)s(K1,X) − (n²+2n−6)s(K1,X)² + 3n s(K1,X)³ − 2s(K1,X)⁴
          + 2(n²+n−8)s(K2,X) − 12s(K2,X)² − 12(n−1)s(K1,X)s(K2,X)
          + 12s(K1,X)²s(K2,X) + 72s(C4,X) + 12(n−2)s(K3,X) + 24s(K1,3,X)
          + 24s(P4,X) + 24s(T_{3,1},X) + 12(n+2)s(P3,X)
          − 24s(K1,X)s(P3,X) + 32s(T_{3,2},X);

    p2(Y) = 4s(K2,Y)² − 12s(K1,3,Y) − 8s(C4,Y) − 8s(T_{3,1},Y)
          − 24s(T_{3,2},Y) + 2(n−8)s(P3,Y);

    p3(X,Y) = 4 s(K1,X) s(P3,Y) − 2(n−2) s(K1,X) s(K2,Y)
            + 4 s(K1,X)² s(K2,Y).

Here s(J, H) is the number of induced copies of J in H; T_{3,1} = paw,
T_{3,2} = diamond; every term has total argument-vertex count 4 (the K=4
case). Conjectures 1.1/1.2 (one "difficult" identity per K ≥ 4; uniqueness
modulo easy families) remain OPEN (terminal v3); Engstrom calls this "the
first example of one outside them [the easily described families]".

## 2. MR97 anchors

Thm 2.2 (general T_m identity; display reconstructed from the paper's proof
chain — Lemmas 2.1/2.2 + eq.(2) — and exact-verified): for m ≥ 1,

    Σ_v s(K_m, G_v^-) =
    Σ_v [ (n/m − s(K1, G_v^+) + m − 2) s(K_{m−1}, G_v^+)
          + (m−1) s(K_m, G_v^+)
          + Σ_{j=1}^{m−2} ((1+δ_{j,m−2}) j/(j+1)) s(T_{m−1,j}, G_v^+) ].

T_{m,j}: T_{m,0} = K_m ∪ K1; T_{m,j+1} adds one edge.
Whitney β (Lemma 2.1 VERBATIM): (n−m) s(K_m, G) = Σ_j β_{m,j} s(T_{m,j}, G),
β_{m,j} = m+1 if j = m; 2 if j = m−1; 1 if 0 ≤ j ≤ m−2.

- m ≥ 5 is void at (5,5) [inference from definitions]: G_v^+ is K4-free,
  every G_v^- is K5-free ⇒ both sides collapse to 0.
- Completeness (Lemma 2.3 VERBATIM): "The only identities of degree at most
  6, in which p can be separated as p(G_v^+, G_v^-) = p1(G_v^+) + p2(G_v^-),
  are those of Theorem 2.2 and their linear combinations." The Engstrom row
  is degree-4 and NON-separable (the p3 cross terms), hence outside the
  T-family as a formal identity.

Thm 3.2 (n = 49): both extremal H1, H2 (the two (24,132) census graphs):
s(K2)=132, s(K3)=176, s(K4)=0, s(paw)=1584, s(diamond)=792;
s(K4, complement(H)) = 144 and 138. Per-vertex m=4 row evaluates to 132
exactly on both; 144/138 ≠ 132 gives the contradiction. (Reproduced locally
against the frozen census.)

## 3. Per-vertex mixed row at (5,5), n = 45

F(X, Y) := p1(X) + p2(Y) + p3(X, Y) with n = 45; Σ_v F_v = 0.
X = G_v^+ is a (4,5,d_v)-graph, d_v ∈ {20,…,24}; Y = G_v^- is the complement
of a (4,5,44−d_v)-graph.

Required motif counters (ALL order ≤ 4; NO order-5 induced counts; NO
s(K4,·) at all — unlike the T-family m=4 row):

- X side: K1 (= d_v), K2 (= edges), K3 (= t), P3, C4, K_{1,3}, P4, paw, diamond.
- Y side: K2 (= edges), P3, C4, K_{1,3}, paw, diamond.
- Cross (products, all integer coefficients): d_v·s(P3,Y), d_v·s(K2,Y),
  d_v²·s(K2,Y) — note d_v and s(K2,Y) are STATE-determined in the cone;
  only s(P3,Y) needs a window.

Fast closed forms (all O(ν²) pair aggregation, bitset words):

- induced P3 = Σ_v C(d_v,2) − 3t   (wedges minus triangle triples)
- paw = Σ_{T triangle} Σ_{x∈T} (d_x − 2) − 4·diamond − 12·K4
- 2·C4 + diamond = #{uv ∉ E : c_uv = 2}   (c_uv = common neighbors)
- K_{1,3} = Σ_v C(d_v,3) − paw − 4·K4
- P4 = C(ν,4) − (sum of the other seven 4-vertex induced classes)
- K4(Y) = i4(complement class of Y) — catalog-side i4 windows carry it.

Cost: same class as the landed m4 kernel (O(ν²)/vertex, ν ≤ 24 →
microseconds/vertex; linear catalog streaming over the frozen catalogs).

## 4. Beyond the T-family?

YES, at the level of formal identities: disjoint coordinate support
(C4, K_{1,3}, P3, P4 both sides + 3 mixed product coordinates) vs the
T rows' {d, e, t, paw, diamond, K4} N-side and clique-counts D-side;
MR97 Lemma 2.3 exhausts separable identities of degree ≤ 6, and the mixed
row is non-separable. Engstrom himself: first identity "outside the easily
described families".

CAVEAT [INFERENCE]: after projection onto the frozen (5,5,45) domain
(d ∈ 20..24, K4-free N side, relation system R1–R5 on 4-motifs), coordinate
collapse is not ruled out — whether the mixed row increases the certificate
cone rank is a COMPUTATION, not a formal corollary. At n = 49 the per-vertex
Engstrom values on the four (X, complement-Y) type combos are {144, 432, 0,
288} — all multiples of 144, non-negative, unique zero corner: consistent
with a type-uniform n=49 configuration; it does NOT reproduce the m=4
contradiction there (hence it is a candidate CUT row for n=45, not a
label/replay row).

Also corrected repo fact: with the CORRECT Whitney β (β_{3,2}=2), the
landed universal m4 kernel and MR97's m=4 identity are the same identity
(the ×12 forms coincide after substituting
s(K3∪K1, X) = (d−3)t(X) − paw(X) − 2·diamond(X) − 4·K4(X)).

## 5. Complementology for Y-side windows (Main, verified 2026-08-18)

500 random labeled graphs n in 5..8 PER LINE, naive 4-subset references vs
closed forms (all assertions held; two earlier candidate forms — K13(c)=paw
and C4 self-complementary — were REJECTED by this same test before landing):

    s(P3, c Z)  =   # triples with exactly one edge in Z
                =   Sum_{uv in E(Z)} ( nu - d_u - d_v + c_uv(Z) )
    C4(c Z)     =   2K2(Z)                      (C4 is NOT self-complementary)
    K_{1,3}(c Z)=   s(K3+K1, Z)                 (LANDED m4 kernel counter)
    paw(c Z)    =   s(P3 + K1, Z)
    diamond(c Z)=   s(K2 + 2K1, Z)
    K4(c Z)     =   s(E4, Z) = i4(Z)            (LANDED kernel counter)
    P4(c Z)     =   P4(Z)                       (P4 IS self-complementary)
    e(c Z)      =   C(nu,2) - e(Z)

    2K2(Z)      =   1/2 * Sum_{uv in E} e(G[ common non-neighbors(u,v) ])
    P4(Z)       =   C(nu,4) - sum(the other ten 4-vertex induced classes)

Correction (MixedPlanDrafter, independently re-verified by Main on 600
graphs + a pinned diamond witness): the Z-side claw closed form is the
UNIVERSAL

    s(K_{1,3}, Z) = Sum_v C(d_v, 3) - paw(Z) - 2*diamond(Z) - 4*K4(Z)

(every induced diamond is counted twice from its two degree-3 vertices,
every K4 four times); on the K4-free R(4,5) catalogs the -4*K4 term is 0
by structural assertion, not assumption. An earlier chat-relayed short
form dropping the 2*diamond + 4*K4 terms was WRONG and never loaded into
any artifact.

Hence every Y-side window the p2 row needs (K2/e, P3, C4, K13, paw, diamond)
reduces to Z-side counts computable from the new counters {2K2, K2+2K1,
P3+K1, C4-pair-stats} plus landed kernel counters; no complement graph is
ever materialized for window purposes.
