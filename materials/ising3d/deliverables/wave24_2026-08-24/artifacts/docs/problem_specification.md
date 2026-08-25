# Problem Specification — Exact Solution of the 3D Ising Model

Status: living document. Every convention used anywhere in this repository is fixed here.
Any code or note that deviates must say so explicitly and locally.

## 1. The model

Let `d >= 1` and let `L = (L_1, ..., L_d)` be a tuple of side lengths. The vertex set is

    Lambda = Z_{L_1} x ... x Z_{L_d}      (periodic)      or
    Lambda = {0..L_1-1} x ... x {0..L_d-1} (open / mixed)

with `N = |Lambda| = prod_i L_i` sites. Spins `sigma_r in {-1,+1}`.

Bond set `B(Lambda)`: for each direction `i` and each site `r`, the pair `(r, r+e_i)` is a bond
iff `r_i + 1 < L_i` (open in direction i) or always, with wraparound (periodic in direction i).

Degenerate cases that MUST be handled explicitly (they are a classic source of factor-of-two errors):

* `L_i = 1` with periodic BC in direction `i`: the "bond" `(r, r+e_i)` is a self-loop `(r,r)`.
  Convention used in this repository: **self-loops are dropped**, i.e. a periodic direction of
  length 1 contributes no bonds. This makes `L_i = 1` periodic identical to `L_i = 1` open.
* `L_i = 2` with periodic BC: the two sites are joined by **two distinct parallel bonds**
  (`r -> r+e_i` and `r+e_i -> r+2e_i = r`). Convention: **both bonds are kept**, so a periodic
  direction of length 2 contributes `2 * (N / L_i) * 1 = N` bonds, exactly as for any other
  periodic length. This is the convention under which the transfer matrix `Tr T^{L_i}` is correct.

With these conventions a fully periodic `d`-dimensional lattice with all `L_i >= 2` has exactly
`|B| = d * N` bonds, and every site has coordination number `2d`.

## 2. Hamiltonian, partition function, free energy

    H(sigma) = - sum_{i=1..d} J_i sum_{(r,r') in B_i} sigma_r sigma_{r'} - h sum_r sigma_r

    K_i = beta J_i ,  H_f = beta h

    Z_Lambda(K, H_f) = sum_{sigma in {-1,1}^Lambda}
        exp[ sum_i K_i sum_{(r,r') in B_i} sigma_r sigma_{r'} + H_f sum_r sigma_r ]

    phi(K, H_f) = lim_{|Lambda| -> inf} (1/|Lambda|) log Z_Lambda(K, H_f)      (reduced free energy)

    f = -(1/beta) phi                                                           (physical free energy)

Primary target: `phi(K,K,K,0)` on the simple cubic lattice, `d = 3`.

## 3. Exact combinatorial encodings used in code

For a fixed lattice we always compute the **exact integer joint density of states**

    g(b, m) = # { sigma : B_sat(sigma) = b , M(sigma) = m }

where, writing `n_b = |B|` and `N = |Lambda|`,

    B_sat(sigma) = # { bonds (r,r') : sigma_r sigma_{r'} = +1 }   in {0,...,n_b}
    M(sigma)     = sum_r sigma_r                                   in {-N,...,N}, step 2.

Then, **exactly**,

    sum_{(r,r')} sigma_r sigma_{r'} = 2 B_sat - n_b
    Z = sum_{b,m} g(b,m) exp[K (2b - n_b) + H_f m]                     (isotropic)
      = e^{-K n_b} sum_{b,m} g(b,m) (e^{2K})^b (e^{H_f})^m

For anisotropic couplings we store `g(b_1,...,b_d, m)` with a per-direction satisfied-bond count.

Two standard change of variables:

    v = tanh K,        x = e^{-2K},        u = e^{-4K},       z = e^{-2 H_f}

`Z` is a polynomial with **non-negative integer coefficients** in `(e^{2K}, e^{H_f})` up to the
overall factor `e^{-K n_b}`. All exact-enumeration code returns integer arrays, never floats.

### High-temperature (even subgraph) form

    Z = 2^N (cosh K)^{n_b} sum_{E subset B, E even} v^{|E|}
      = 2^N (cosh K)^{n_b} P(v),          P(v) in Z[v], P(0) = 1

`P(v)` is the even-subgraph generating polynomial of the lattice graph. This is our canonical
exact fingerprint of a finite lattice, and the object from which HT series are extracted.

### Low-temperature form (zero field, d>=2)

    Z = 2 e^{K n_b} sum_{S} x^{|S|},   x = e^{-2K}

sum over "domain wall" configurations = even subgraphs of the dual (for d=2) / closed surfaces
(for d=3), with the leading factor 2 being the two ground states.

## 4. Reference values (benchmarks only, never targets to fit)

| quantity | value | source |
|---|---|---|
| 2D square `K_c` | `ln(1+sqrt 2)/2 = 0.4406867935...` | Onsager 1944, exact |
| 3D sc `K_c` | `0.221654626(5)` | Ferrenberg-Xu-Landau 2018 (to be verified from source) |
| 3D sc `nu` | `0.629971(4)` | conformal bootstrap (to be verified from source) |

## 5. What counts as a solution

See the task statement, Section 2. Operationally, in this repository a claimed exact result must
carry a `results/gates/<name>.json` file recording the pass/fail status of all applicable
validation gates G1-G10, produced by code in `tests/` or `experiments/`.

## 6. Numerical hygiene

* Exact integer arithmetic (Python `int`, numpy `int64` with overflow assertions, or `object`
  dtype) for all combinatorial quantities.
* `mpmath` with explicitly set precision for all transcendental evaluation. Precision is recorded
  in the output artifact.
* Any float comparison in a test carries an explicit tolerance and a justification for it.
