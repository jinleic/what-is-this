# Erdős #389: proved reductions and current construction barriers

All statements use

$$
N=m+k,
\qquad x=m+2k=2N-m,
\qquad
s_p(m,k)=v_p\binom{x}{k}-v_p\binom{N}{m}.
$$

`PROVED` means the argument is given here. `CERTIFIED FINITE` means an exact
script exhausts only its declared finite domain. No novelty claim is made.

## 1. Prime-power floor decomposition — PROVED

For a prime power $q=p^j$, define

$$
\Delta_q(m,k)
=\left\lfloor\frac{x}{q}\right\rfloor
 +\left\lfloor\frac m q\right\rfloor
 -2\left\lfloor\frac Nq\right\rfloor.
$$

Legendre's formula gives

$$
s_p(m,k)=\sum_{j\ge1}\Delta_{p^j}(m,k). \tag{1}
$$

Let $A=\lfloor x/q\rfloor$, $B=\lfloor m/q\rfloor$ and let $r_x,r_m$
be the two residues modulo $q$. If $A+B$ is even, then $\Delta_q=0$. If
$A+B$ is odd, then

$$
\Delta_q=
\begin{cases}
+1,&r_x+r_m<q,\\
-1,&r_x+r_m\ge q.
\end{cases} \tag{2}
$$

Indeed, write $(x+m)/(2q)=(A+B)/2+(r_x+r_m)/(2q)$ and take the floor.
Thus every level contributes exactly $-1$, $0$, or $+1$.

`erdos389.local_floor_contribution`, `slack_floor_contributions`, and 2,000
deterministic random cases in the test suite independently check (1)--(2).

## 2. Exact residue zones — PROVED

### Levels $q\le m$

Put $r_1=m\bmod q$ and $r_2=k\bmod q$. Then

$$
\Delta_q=
\left\lfloor\frac{r_1+2r_2}{q}\right\rfloor
-2\left\lfloor\frac{r_1+r_2}{q}\right\rfloor. \tag{3}
$$

Consequently

$$
\Delta_q=
\begin{cases}
+1,&r_1+r_2<q\le r_1+2r_2,\\
-1,&r_1+r_2\ge q\text{ and }r_1+2r_2<2q,\\
0,&\text{otherwise}.
\end{cases} \tag{4}
$$

A draft derivation incorrectly assigned $+1$ when both
$r_1+r_2\ge q$ and $r_1+2r_2\ge2q$. Equation (3) shows that branch is
$2-2=0$; $(m,k,q)=(5,2,3)$ is an immediate counterexample to the draft.
The incorrect branch is not used anywhere in the code or results.

### Levels $q>m$

Put $t=N\bmod q$. Directly from $x=2N-m$,

$$
\Delta_q=
\begin{cases}
-1,&0\le t<\lceil m/2\rceil,\\
0,&\lceil m/2\rceil\le t<\lceil(q+m)/2\rceil,\\
+1,&\lceil(q+m)/2\rceil\le t<q.
\end{cases} \tag{5}
$$

The exact certificate `data/zone_theorem_m1_8_k2000.json` checks (1),
(3)--(5), and the corollaries below at 330,023 prime-power levels over all
$1\le m\le8$, $1\le k\le2000$.

## 3. Bad-window support and a sharp obstruction — PROVED

Define the top half of the left-binomial numerator interval by

$$
\mathcal B(m,k)=
\{k+\lfloor m/2\rfloor+1,\ldots,k+m\}. \tag{6}
$$

For $q>m$, (5) implies

$$
\Delta_q=-1
\quad\Longleftrightarrow\quad
q\mid w\text{ for the unique possible }w\in\mathcal B(m,k). \tag{7}
$$

Therefore every obstructing prime $p$ either satisfies $p\le m$ or divides an
integer in $\mathcal B(m,k)$. This cuts the large-prime support from all $m$
numerator terms to the top $\lceil m/2\rceil$ terms.

Also, every prime power $q$ with $N<q\le x$ has $\Delta_q=+1$: here
$N\bmod q=N\ge\lceil(q+m)/2\rceil$ is equivalent to $q\le2N-m=x$.

**Single-level obstruction.** If a prime

$$
\max(m,\sqrt x)<p\le N
$$

divides some $w\in\mathcal B(m,k)$, then $s_p(m,k)=-1$ and $(m,k)$ is not a
witness. Equation (7) gives the unique negative level $p$; $p^2>x$ removes
all higher levels. The zone certificate checks 25,536 instances of this
criterion exactly.

## 4. Small-prime controls — PROVED

### Neutral levels

If $q\mid k$, then $\Delta_q=0$. This follows at once from (3), or from (5)
when $q>m$. Hence divisibility by prime powers can delete chosen small-prime
levels without creating a negative contribution.

For every prime $p\le m$, let $a_p$ be minimal with $p^{a_p}>m$, and set

$$
M_m=\prod_{p\le m}p^{a_p}. \tag{8}
$$

If $M_m\mid k$, then adding $m$ to $k$ in base $p$ produces no carry for every
$p\le m$, so

$$
v_p\binom{m+k}{m}=0\qquad(p\le m). \tag{9}
$$

Thus the witness test on this corridor reduces entirely to primes $p>m$.
The bounded analyses have checked 16,826 instances of (9), and the deeper
candidate-family search checked another 112,210 small-prime carry identities.

### Positive forcing levels

If $k\equiv-(m+1)\pmod q$, then $r_1+r_2=q-1$. Equation (3) gives
$\Delta_q=+1$ unless $q\mid m+1$, in which case it gives $0$. For $q>m$ this
is the same good-zone statement, including the neutral boundary $q=m+1$.

Let $J_p=v_p(k+m+1)$ and $L_p=\lfloor\log_p x\rfloor$, with $J_p\le L_p$.
The uncontrolled higher levels contribute at worst $-1$, so

$$
s_p(m,k)\ge2J_p-v_p(m+1)-L_p. \tag{10}
$$

If (10) alone is required to certify every $p\le m$, put
$A=\prod_{p\le m}p^{J_p}\mid N+1$. Since
$p^{L_p}>x/p$,

$$
x^2\ge A^2
\ge\prod_{p\le m}p^{L_p}
>\frac{x^{\pi(m)}}{\prod_{p\le m}p}.
$$

Writing $\vartheta(m)=\sum_{p\le m}\log p$ gives the finite bound

$$
x^{\pi(m)-2}<e^{\vartheta(m)}. \tag{11}
$$

For $m\ge5$, pure all-prime forcing is therefore confined to an effectively
finite $x$-range. A successful general construction must normally use
same-prime compensation across positive and negative levels, not force every
level positive independently.

## 5. Odd-$m$ CRT-safe natural shift — PROVED

For the natural shift $(m,k)\mapsto(m+1,k-1)$,

$$
s'_p=s_p+v_p(m+1)-v_p(x). \tag{12}
$$

Suppose $m$ is odd and $p>m$. Only $p\mid x$ can introduce a negative term in
(12). Put $a=v_p(x)\ge1$. The exceptional possibility $m=1,p=2$ cannot occur
because $x$ is odd, so $p$ is odd. For $q=p^j$, $j\le a$, the congruence
$2N\equiv m\pmod q$ gives

$$
N\bmod q=(q+m)/2,
$$

the lower boundary of the good zone in (5). Hence each of these $a$ levels
contributes $+1$.

For $j>a$, the residue $r=x\bmod q$ is a nonzero multiple of $p^a$, so

$$
r\le q-p^a<q-m.
$$

If this level were bad, then $t=N\bmod q<\lceil m/2\rceil$ and
$r=q+2t-m>q-m$, a contradiction. Thus no higher level is negative and
$s_p\ge a$. Also $p\nmid m+1$, so (12) gives $s'_p\ge0$. Therefore every
obstructing prime for odd $m$ satisfies $p\le m$.

Consequently, if $(m,k)$ is a witness and

$$
\gcd\!\left(m+2k,\prod_{p\le m}p\right)=1, \tag{13}
$$

then $(m+1,k-1)$ is also a witness.

For each odd prime $p\le m$, condition (12) excludes exactly one class

$$
k\equiv-m\,2^{-1}\pmod p.
$$

By CRT, modulo $P_m=\prod_{3\le p\le m}p$ there are exactly
$\prod_{3\le p\le m}(p-1)$ safe classes. This proves that the residue system
itself is consistent; it does **not** prove that a witness lies in one of those
classes.

`data/shift_repair_analysis_m1_50_k50000.json` checks all 5,798 CRT-safe odd
source witnesses in the atlas, and every one shifts successfully. Of 21 failed
odd shifts, only 3 find a later CRT-safe witness within one period $P_m$.
All 21 eventually reach one within the $k\le50{,}000$ rectangle, but the worst
case needs 559 CRT periods. The first safe residue alone is usually not a
witness: CRT avoidance is a valid filter, not a local witness-preserving map.

## 6. Even-$m$ fatal prime-power spikes — PROVED

Let $m$ be even, $(m,k)$ a witness, $x=m+2k$, and let $p>m+1$ have
$a=v_p(x)\ge1$. If

$$
p^{a+1}>x+m, \tag{14}
$$

then the natural shift fails at $p$.

For $j\le a$, $p^j\mid x$ and $x/p^j$ is even, so the level contribution is
zero. For $j>a$, (14) makes all factorial floors zero. Hence $s_p(m,k)=0$,
while the shifted slack is $-a$ because $p\nmid m+1$.

In particular, a prime divisor

$$
p>\max(m+1,\sqrt{x+m})
$$

is fatal. The original bounded certificate identifies 4,578 of the 5,669
failed even-$m$ shifts. Section 10 gives the exact multi-level replacement.

This yields a useful parity dichotomy:

- odd-to-even propagation can eliminate every possible obstruction by making
  $x$ small-prime-rough through finite CRT avoidance;
- even-to-odd propagation is usually destroyed by an isolated large prime or
  prime-power factor of $x$ and needs smoothness/compensation instead.

## 7. A smooth-window sufficient condition — PROVED

If every integer in $\mathcal B(m,k)$ is $m$-smooth, then no prime $p>m$ has
a negative level by (7). Therefore $(m,k)$ is a witness if and only if
$s_p(m,k)\ge0$ for the finite set of primes $p\le m$.

This is a sufficient construction family, not an existence theorem. Long runs
of $m$-smooth integers are highly restrictive, so the condition may generate
only finitely many candidates for fixed $m$; compensated witnesses outside
this family remain possible.

## 8. Local-safe small-prime residue classes — PROVED

### Exact prime-power bands

Fix $p$ and $L\ge1$, and restrict $k$ to

$$
p^L\le m+2k<p^{L+1}.
$$

Every level above $p^L$ is then zero, while each level
$p,p^2,\ldots,p^L$ depends only on $k$ modulo that power by (3)--(5).
Consequently $s_p(m,k)$ depends exactly on $k\bmod p^L$ throughout this band.
If

$$
\mathcal A_{m,p,L}
=\{r\bmod p^L:\sum_{j=1}^L\Delta_{p^j}(m,r)\ge0\},
$$

then $s_p(m,k)\ge0$ if and only if $k\bmod p^L$ lies in
$\mathcal A_{m,p,L}$. On any intersection of fixed prime-power bands, the
complete finite small-prime condition is therefore an exact union of CRT
classes, not merely a sufficient filter. The control artifact verifies 711
values across eight bands; for example
$|\mathcal A_{3,5,3}|=103$ out of $125$ residues.

### Uniform sufficient classes

Fix $p$ and a power $q=p^a>m$. For $0\le r<q$, suppose

$$
r+m<q,
\qquad
\operatorname{carries}_p(r,r+m)\ge
\operatorname{carries}_p(r,m). \tag{15}
$$

Then every $k\equiv r\pmod q$ satisfies $s_p(m,k)\ge0$. Indeed, write
$k=r+qh$. The first inequality in (15) prevents a carry from the low $a$
digits when adding $m$, so the complete left carry count is exactly
$\operatorname{carries}_p(r,m)$. The low $a$ digits of the right addition are
$r+(r+m)$; every carry above them, from adding $h+h$ and the boundary carry,
is additional and nonnegative. The second inequality proves the claim.

Let $\mathcal S_{m,p,a}$ be the residues satisfying (15). It is nonempty:
$r=0$ is always present. For distinct primes the moduli are coprime, so CRT
gives exactly

$$
\prod_{p\le m}|\mathcal S_{m,p,a_p}|
$$

safe classes modulo $\prod_{p\le m}p^{a_p}$, with the exponent $a_p$ chosen
independently for each prime. This strictly extends the single zero-carry
class. At the least exponent $p^{a_p}>m$, the number of carry-free residues is
$\prod_j(p-m_j)$, where $m_j$ are the base-$p$ digits of $m$.

The exact artifact `data/residue_control_analysis_m1_50_k50000.json` classifies
all least-power local-safe residues for $m\le50$ and checks 523,408 extensions
with four different unrestricted high parts. For $m=27$, using zero, one, or
two extra base-$p$ digits at every $p\le27$ gives exact CRT class densities
approximately $0.0046962$, $0.0626106$, and $0.155369$, respectively. These
are densities of sufficient residue classes, not witness densities.

### No fixed CRT class is sufficient — PROVED using Dirichlet's theorem

Fix any modulus $M$, class $c\bmod M$, and one bad-window position
$\lfloor m/2\rfloor< i\le m$. Put

$$
d=\gcd(c+i,M),\qquad
a=(c+i)/d,\qquad
M'=M/d.
$$

Since $\gcd(a,M')=1$, Dirichlet's theorem supplies arbitrarily large primes
$p\equiv a\pmod {M'}$. Choose one with $p>\max(m,2d)$ and set
$k=dp-i$. Then $k\equiv c\pmod M$, while $k+i=dp$ is a bad-window term.
Moreover

$$
x=m+2k=2dp+m-2i<2dp<p^2.
$$

Thus the single-level criterion in Section 3 gives $s_p(m,k)=-1$. Every
residue class of every fixed modulus therefore contains arbitrarily large
non-witnesses. The same holds for every nonempty finite union of fixed CRT
classes. Small-prime CRT can be a pruning or repair tier, but can never by
itself be a sufficient construction.

The control artifact checks four explicit progressions, including
$k=301\equiv1\pmod6$ at $(m,p)=(2,101)$ and
$k=1027\equiv7\pmod {30}$ at $(m,p)=(4,103)$.

## 9. Exact same-prime compensation in the bad window — PROVED

Let $p>m$ divide the unique bad-window term $w=k+i$, where
$\lfloor m/2\rfloor<i\le m$. Put

$$
e=v_p(w),\qquad c=w/p^e,\qquad u_r=c\bmod p^r.
$$

Then

$$
\boxed{\displaystyle
s_p(m,k)=-e+\#\{r\ge1:2u_r>p^r\}.} \tag{16}
$$

For $1\le j\le e$, $p^j\mid w$, so (7) gives one negative contribution.
For $j=e+r$, the residue of $N=m+k$ modulo $p^{e+r}$ is

$$
p^e u_r+(m-i).
$$

There is no wrap because $p^e>m$. Since $p\nmid c$, this residue is never in
the bad zone. Substitution into (5) shows that it is in the strict good zone
exactly when $2u_r>p^r$. Once $p^r>2c$, the inequality is impossible, so the
count is finite. Moreover, $2u_r>p^r$ gives $2c\ge p^r+1$, so
$2w\ge p^{e+r}+p^e>p^{e+r}+m$ and hence $p^{e+r}<2w-m\le x$.
Thus every counted positive level actually occurs before the floor sum ends.
This proves (16), an exact necessary-and-sufficient same-prime compensation
criterion.

There is also an explicit local repair lift. For odd $p>m$,

$$
k+i\equiv-p^e\pmod {p^{2e}} \tag{17}
$$

makes $u_r=p^r-1$ for $1\le r\le e$, producing at least $e$ positive levels
after the $e$ negative levels. For $p=2$, possible here only when $m=1$, the
first prefix is exactly half rather than strictly good; the correct modulus is
$2^{2e+1}$. The naive binary version of (17) fails already at
$(m,k,p)=(1,1,2)$, while $k+1\equiv-2^e\pmod {2^{2e+1}}$ repairs it.

Local-safe conditions (15) and finitely many lifts (17) have coprime
prime-power moduli, so CRT combines them without conflict. This rigorously
repairs every prescribed finite prime set. It does not control new primes
entering the translated bad window; that moving intersection remains the
load-bearing gap.

The finite certificate checks (16) at 50,360 large-prime bad-window events,
440 canonical mirror lifts, and 132 controlled primes in combined CRT systems.

### Compensation-good integers and exact small/large separation — PROVED

For $p^e\Vert w$, define

$$
C_p(w)=
\#\left\{r\ge1:
2\left(\frac{w}{p^e}\bmod p^r\right)>p^r
\right\}.
$$

Call $w$ **$m$-compensation-good** when

$$
C_p(w)\ge v_p(w)
\qquad\text{for every prime }p>m\text{ dividing }w. \tag{18}
$$

Then $(m,k)$ is a witness if and only if both of the following hold:

1. $s_p(m,k)\ge0$ for every prime $p\le m$;
2. every $w\in\mathcal B(m,k)$ is $m$-compensation-good.

Indeed, a prime $p>m$ divides at most one term of $\mathcal B(m,k)$, and
formula (16) says that its slack is exactly $C_p(w)-v_p(w)$. Section 3 says
that a prime $p>m$ not dividing the bad window cannot obstruct. This proves
both directions. Thus the global problem separates exactly into the finite
small-prime tier and a run of $\lceil m/2\rceil$ consecutive
$m$-compensation-good integers. On a fixed prime-power band, Section 8 makes
the first tier an exact finite CRT union; no approximation is involved.

The strict hypothesis $p>m$ is essential: at $(m,k,p)=(5,2,5)$, the prime
divides the bad-window term $5$, but its exact slack is $0$, whereas the local
formula would predict $-1$.

### The full guaranteed first-prefix repair cone — PROVED

The mirror residue in (17) is only one sufficient choice. For odd $p$ and
$e\ge1$, let

$$
H_{p,e}=
\left\{c\bmod p^e:
2(c\bmod p^r)>p^r\text{ for }1\le r\le e
\right\}.
$$

Writing $c=a_0+a_1p+\cdots+a_{e-1}p^{e-1}$, the inequalities are equivalent to

$$
a_0\in\left\{\frac{p+1}{2},\ldots,p-1\right\},
\qquad
a_j\in\left\{\frac{p-1}{2},\ldots,p-1\right\}\quad(1\le j<e).
$$

For $j=0$ this is the first strict upper-half inequality. Inductively, once
the lower prefix is strictly above half, the next digit may equal
$(p-1)/2$; a smaller digit fails even with the largest possible lower prefix,
while every larger digit succeeds. Hence

$$
\boxed{\displaystyle
|H_{p,e}|=\frac{p-1}{2}
\left(\frac{p+1}{2}\right)^{e-1}.} \tag{19}
$$

This is the exact set of cofactor residues whose first $e$ prefixes alone
guarantee repair for every higher lift. Cofactors outside it may still be
repaired by additional higher prefixes.


For every $c\in H_{p,e}$,

$$
k+i\equiv p^e c\pmod {p^{2e}}
$$

has valuation exactly $e$ and supplies at least $e$ positive prefixes, so it
repairs $p$. The mirror lift is the single choice $c\equiv-1\pmod {p^e}$.
For $p=2$, the first strict prefix is impossible for odd $c$, so this cone is
empty; the later-prefix binary repair after (17) remains necessary.

### Retaining all old blockers forces product-scale displacement — PROVED

Suppose distinct primes $p\in S$ are attached to positions $i_p$ of the bad
window at $k_0$, so $p\mid k_0+i_p$, and put $P=\prod_{p\in S}p$. Any $k'$
that retains every attachment satisfies

$$
k'\equiv-i_p\equiv k_0\pmod p
\quad(p\in S),
$$

and therefore $k'\equiv k_0\pmod P$. Consequently every distinct forward
repair retaining all old primes has

$$
\boxed{k'\ge k_0+P.} \tag{20}
$$

This lower bound applies to the entire cone (19), not only to the mirror
residue. It does not restrict a repair that drops an old prime from the
translated bad window.

## 10. Exact natural-shift spike formula — PROVED

Let an odd prime $p>m$ divide $x=m+2k$, put $a=v_p(x)$ and
$y=x/p^a$. Then $p^a>m$ and

$$
\boxed{\displaystyle
s_p(m,k)
=a(m\bmod2)+\sum_{r\ge1}
\left(\left\lfloor\frac{y}{p^r}\right\rfloor\bmod2\right).} \tag{21}
$$

For $j\le a$, division by the odd number $p^j$ preserves the parity of $x$,
which is the parity of $m$; the level contributes one for odd $m$ and zero for
even $m$. For $j=a+r$, the nonzero remainder is a multiple of $p^a$ and is at
most $p^{a+r}-p^a<p^{a+r}-m$, so this level cannot be negative. It contributes
one exactly when $\lfloor y/p^r\rfloor$ is odd. This proves (21).

The natural shift adds $v_p(m+1)-a$ at $p$. Thus odd $m$ gives a
nonnegative shifted slack equal to the sum in (21), recovering Section 5. For
even $m$, the shifted pair fails at $p$ exactly when

$$
\sum_{r\ge1}
\left(\left\lfloor\frac{y}{p^r}\right\rfloor\bmod2\right)
<a-v_p(m+1). \tag{22}
$$

This replaces the sufficient spike condition (14) by a necessary-and-sufficient
finite digit test. For even $m$, since the sum is at most
$\lfloor\log_p y\rfloor$, the shift necessarily fails whenever
$p^{2a}>x$ and $p\nmid m+1$; this strictly strengthens the old sufficient
condition.

### Even-shift spikes are first-term compensation — PROVED

Let the source index $m=2h$ be even and shift to $(m+1,k-1)$. The first term
of the target bad window is

$$
(k-1)+(h+1)=k+h=\frac{x}{2}.
$$

In fact the entire target window is

$$
\mathcal B(m+1,k-1)=\mathcal B(m,k)\cup\{x/2\}.
$$

For an odd source $m$, the corresponding two bad windows are equal. Thus a
natural shift changes the large-prime local conditions only by adjoining
$x/2$ in the even-source case.

For an odd prime $p>m+1$ dividing $x$, put $p^a\Vert x$, $y=x/p^a$, and
$c=y/2$. For every $q=p^r$, write $y=qA+B$ with $0\le B<q$. Because $q$ is
odd and $y$ is even, $A$ and $B$ have the same parity. If $A$ is even, then
$c\bmod q=B/2\le q/2$; if $A$ is odd, then
$c\bmod q=(q+B)/2>q/2$. Therefore

$$
\mathbf 1_{\{2(c\bmod p^r)>p^r\}}
=\left(\left\lfloor\frac{y}{p^r}\right\rfloor\bmod2\right).
$$

Formula (16) for the first bad-window term $x/2$ is consequently identical
to the shifted slack obtained from (21). Large prime-power spike control is
not a separate phenomenon: it is exactly compensation-goodness of the first
target bad-window integer.

The restriction $p>m+1$ is sharp: $p=m+1$ belongs to the target small-prime
tier and the added term $v_p(m+1)$ in (12) must be retained.


In the repaired atlas the exact formula explains 5,620 of 5,669 failed
even-$m$ shifts: the original one-level criterion accounts for 4,578 and the
full digit formula adds 1,042. The remaining 49 even failures have only
small-prime obstructions; all 21 odd failures are likewise small-prime
failures. The artifact checks 24,809 factored shift primes exactly.

## 11. Exact finite evidence and current gaps

`data/zero_carry_corridor_m5_18_t2000.json` exhausts 24,014 candidates
$k=tM_m$. It finds five witnesses: one for $m=5$, three for $m=6$, and one for
$m=7$. The odd witnesses at $m=5,7$ shift successfully; the three even $m=6$
witnesses fail at large primes, exactly as the parity dichotomy predicts.

`data/zero_carry_corridor_m8_12_t10000.json` exhausts 50,000 further candidates
and finds no witness. At the first currently unlisted index $m=27$, the two
partial corridors control respectively $p\le7$ and $p\le11$:
`data/partial_zero_carry_m27_p7_t20000.json` and
`data/partial_zero_carry_m27_p11_t20000.json` each exhaust 20,000 candidates
without a witness. In both corridors **every one of the 20,000 candidates** has
a certified single-level large-prime obstruction from Section 3. Strengthening
small-prime control reduced residual small-prime failures but did not touch the
moving bad-window barrier. These are negative finite results, not obstructions
to other families.

`data/prime_repair_analysis_r64.json` classifies every failed natural shift
whose full radius-64 window stays inside the atlas. All 5,678 original
obstruction sets become nonnegative somewhere in the radius. Nevertheless,
only 492 first joint repairs are witnesses; 3,224 require a later offset
because other primes intervene, and 1,962 have no witness in the radius.
Thus 86.8% of the 3,716 repairable cases exhibit new-prime interference.
Repairing the currently visible primes is not the dominant construction
mechanism; it must be coupled to the moving bad-window sieve.

`data/crt_repair_cover_m1_50_k50000_r64.json` replaces the rough small-prime
filter by the exact band classes above. Every one of the 5,678
boundary-comparable failed shifts reaches the exact $p\le m+1$ tier within the
radius; 5,321 already satisfy it at the failed natural-shift target. Yet the
first tier candidate is a witness only 14 times. The moving window blocks
3,702 before a later witness and blocks another 1,962 with no witness in the
radius. At those first tier candidates the certificate finds 4,759
single-level blockers and 1,847 compensation deficits. Thus the exact
small-prime tier is almost never the binding constraint in this rectangle.

The first unlisted transition is explicit. The published source
$(m,k)=(26,5{,}048{,}891{,}644{,}620)$ is a witness, but its natural-shift
target at $m=27$ fails at
$p=5{,}048{,}891{,}644{,}633$, the first term of its 14-integer bad window.
`data/moving_bad_window_m27_h100000000.json` checks every forward offset
$0\le h<100{,}000{,}000$. Exactly 99,999,997 are rejected by a prime
$p>\sqrt{27+2(k+h)}$ in the translated bad window. Only three offsets survive
that necessary sieve: $24{,}710{,}064$, $25{,}238{,}496$, and
$25{,}238{,}497$. Complete exact witness checks reject them at respectively
14, 6, and 7 primes in the range $27<p\le\sqrt{27+2(k+h)}$; formula (16)
classifies these as uncompensated multi-level obstructions. This is a finite
negative result around one natural starting point, not a lower bound for all
$m=27$ witnesses.

`data/moving_bad_window_m27_h100000000_199999999.json` checks the next
100,000,000 offsets. A single-level prime rejects 99,999,991; complete
Legendre/Kummer checks reject the nine survivors. Seven have only large-prime
compensation deficits, while two also fail at $p=7$. Together the two moving
artifacts certify that the first 200,000,000 forward offsets from the natural
shift target contain no witness.


The three survivors have 27 obstruction primes in total. Repairing any one by
its mirror lift gives 27 candidates within the factorization bound; every one
fails after new primes enter. Simultaneously combining the exact small-prime
tier with all mirror lifts produces CRT solutions of 79, 88, and 157 decimal
digits, far outside that bound. The exact records are in
`data/bad_window_near_miss_m27_h100000000.json`.

`data/compensation_structure_m1_20_k5000_m27_h200000000.json` checks the exact
small/large separation at 100,000 pairs and 587,169 large-prime term factors.
Of 1,294 large-prime-good windows, 190 fail only at the small-prime tier and
1,104 are witnesses. It exhaustively enumerates 19,804 residues for the prefix
cone, checks the even-shift/first-term identity at 19,654 primes, and combines
the moving certificates to show that all 200,000,000 offsets have a certified
large-prime obstruction. Retaining the old blockers in the first three
survivors forces displacements divisible by 69-, 30-, and 35-digit prime
products. The nine later survivors have 113 large blockers and retained-prime
products of 36--79 digits. The larger cone (19) therefore cannot produce a
nearby same-blocker repair; a useful move must drop old primes and control the
new window.


Within the first million offsets, only five windows have as few as one
single-level fatal prime. They have 1, 1, 2, 3, and 7 exact obstruction primes
after full checking. Fourteen one-prime mirror repairs arise from those five
windows; nine remain within the retained $10^{14}$ factorization bound, and
all nine fail because 12--16 new obstruction primes enter the moved window.
The other five mirror candidates and every combined small-prime/large-prime
CRT solution exceed that implementation bound. These exact negative outcomes
are recorded in `data/bad_window_near_miss_m27_h1000000.json`.

The exact local criterion also removes the old residual-sieve restriction
near small $k$. The two contiguous artifacts
`data/compensation_run_m27_k50001_h100000000.json` and
`data/compensation_run_m27_k100050001_h900000000.json` exhaust

$$
50{,}001\le k\le1{,}000{,}050{,}000.
$$

Every one of these $10^9$ candidates has a certified negative local slack at
some prime $p>27$; the large-prime condition rejects every candidate before a
small-prime filter is needed. Across the factored term intervals, the longest
run of consecutive $27$-compensation-good integers has length $9<14$. Combined
with the atlas
rectangle, this proves the finite lower bound

$$
\boxed{k_{27}>1{,}000{,}050{,}000}
$$

for the least possible $m=27$ witness, if one exists.


Open gaps:

1. The exact target is now an intersection between a prime-power-band
   small-prime CRT class and a run of $\lceil m/2\rceil$ consecutive
   $m$-compensation-good integers. No theorem proves that this intersection is
   nonempty.
2. The Dirichlet construction proves that no fixed modulus can be sufficient;
   a successful construction must adapt to the translated bad window.
3. Formula (16) and cone (19) control any prescribed primes, but (20) proves
   that retaining all current blockers is intrinsically product-scale. No
   construction controls the new primes after old blockers are dropped.
4. No construction here reaches a previously unknown $m$; Erdős #389 remains
   open.
