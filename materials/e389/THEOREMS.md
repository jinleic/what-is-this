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

The support statement has a universal one-way form, including \(q\le m\).
Let \(t=N\bmod q\) and \(r=m\bmod q\). Writing \(N=Aq+t\) and \(m=Bq+r\)
in (1) gives

$$
\boxed{\displaystyle
\Delta_q(m,k)=\left\lfloor\frac{2t-r}{q}\right\rfloor,
\qquad
\Delta_q=-1\Longleftrightarrow 2t<r.} \tag{7a}
$$

When this level is negative, \(t<\lceil m/2\rceil\), so

$$
N-t=k+(m-t)\in\mathcal B(m,k),
\qquad q\mid N-t.
$$

Consequently **every** obstructing prime, small or large, divides some integer
in the bad window: if \(s_p<0\), at least one level \(q=p^j\) is negative and
the displayed multiple supplies the term. For \(q\le m\) the converse need
not hold; divisibility of a bad-window term alone does not force that level
negative. The regression suite checks (7a) and its window position over
\(1\le m\le10\), \(1\le k\le100\).

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

### Short-cofactor spike criterion — PROVED

Suppose \(p^e\Vert w\), put \(c=w/p^e\), and assume \(c<p\). In the
large-prime setting \(p>m\) of this section, (16) specializes to

$$
\boxed{\displaystyle
C_p(w)=\mathbf 1_{\{2c>p\}},
\qquad
s_p(m,k)=-e+\mathbf 1_{\{2c>p\}}.} \tag{16a}
$$

Indeed, \(c\bmod p=c\). For every \(r\ge2\), one again has
\(c\bmod p^r=c\), while \(2c<2p\le p^r\), so no prefix after the first can
be positive. This also covers \(p=2\): then \(c=1\), the first prefix is
exactly half, and the strict inequality fails.


In particular, every \(m\)-compensation-good integer \(w\) obeys the necessary
large-prime-power bounds

$$
p^2<2w\quad(e=1),
\qquad
p^{e+1}<w\quad(e\ge2)
$$

for every \(p^e\Vert w\) with \(p>m\). In the second case, goodness forces
\(c>p\); in the first, either \(c>p\) or the strict repair inequality
\(p<2c\) gives the displayed bound.

Consequently a short cofactor can never repair \(e\ge2\). For \(e=1\), it
repairs the prime exactly when \(p<2c\); the boundary \(2c=p\) is bad.

For odd \(p\), there is a sharp threshold without the short-cofactor
assumption. Every counted prefix satisfies \(p^r\le2c-1\), hence

$$
C_p(w)\le\left\lfloor\log_p(2c-1)\right\rfloor.
$$

Therefore compensation of a valuation \(e\) requires

$$
\boxed{\displaystyle
c\ge\frac{p^e+1}{2},
\qquad
w\ge\frac{p^e(p^e+1)}2.} \tag{16c}
$$

The bound is attained. For \(c_0=(p^e+1)/2\) and \(1\le r\le e\),

$$
c_0\bmod p^r=\frac{p^r+1}{2},
$$

so exactly those first \(e\) prefixes are strictly above half and the slack
is zero. Thus \(p^e(p^e+1)/2\) is the least bad-window term whose exact
\(p\)-valuation is \(e\) and can be compensated. This odd-prime threshold
does not apply to \(p=2\), where the strict half-boundary is reachable; the
binary lift after (17) remains the correct control.

There is a sharp natural-shift form. Let \(m\) be even, let an odd
\(p>m+1\) satisfy \(p^a\Vert x=m+2k\), and put \(y=x/p^a\). The first new
bad-window term after shifting is \(w=x/2=p^a(y/2)\). If \(y<2p\), then

$$
\boxed{\displaystyle
s_p(m+1,k-1)=-a+\mathbf 1_{\{y>p\}}.} \tag{16b}
$$

Thus every such spike with \(a\ge2\) is fatal. When \(a=1\), it is fatal
for \(y<p\) and repaired for \(p<y<2p\); equality cannot occur because
\(y\) is even and \(p\) is odd.

More generally, (16c) and the first-term identity give the stronger fatal
range

$$
\boxed{\displaystyle
x<p^{2a}+p^a
\quad\Longrightarrow\quad
s_p(m+1,k-1)<0.} \tag{16d}
$$

Equality is sharp: \(x=p^a(p^a+1)\) makes the first new term the minimal
compensated term in (16c), with shifted slack zero. The regression suite
checks the threshold and the immediately lower cofactor for
\(p\in\{5,7\}\), \(1\le a\le3\).

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

### Half-window translation evicts every large blocker — PROVED

Put $L=\lceil m/2\rceil$ and $I_m=\{\lfloor m/2\rfloor+1,\ldots,m\}$, so
$|I_m|=L$. Suppose $p>m$ divides a term $k+i_p$ of $\mathcal B(m,k)$, with
$i_p\in I_m$. Then

$$
\boxed{p\nmid k+L+i\qquad\text{for every }i\in I_m.} \tag{23}
$$

If $p$ divided $k+L+i$, it would divide the difference $L+i-i_p$. Both
indices lie in $I_m$, so $|i-i_p|\le m-\lfloor m/2\rfloor-1=L-1$ and

$$
1\le L-(L-1)\le L+i-i_p\le L+(L-1)=2\lceil m/2\rceil-1\le m<p,
$$

where the last inequality of the chain uses $2\lceil m/2\rceil-1\le m$, valid
for every $m\ge1$. A nonzero integer smaller than $p$ in absolute value is not
divisible by $p$, which proves (23). Equivalently, $\mathcal B(m,k+L)$ is the
block of $L$ integers immediately following $\mathcal B(m,k)$, and the union of
the two blocks spans at most $m$ consecutive integers, shorter than the gap
between consecutive multiples of any $p>m$.

The offset is also the least uniform one. Fix $1\le\delta\le L-1$, put
$i_p=\lfloor m/2\rfloor+1+\delta$, which lies in $I_m$, take any prime $p>m$,
and set $k=p-i_p\ge1$. Then $p\mid k+i_p$ and

$$
k+\delta+(i_p-\delta)=k+i_p\equiv0\pmod p,
\qquad
i_p-\delta=\lfloor m/2\rfloor+1\in I_m,
$$

so that blocker survives the translation by $\delta$. Hence

$$
\boxed{\min\{\delta\ge1:\ \delta\ \text{evicts every possible blocker}\}
=\lceil m/2\rceil.} \tag{24}
$$

More precisely, a blocker attached at $i_p$ survives exactly the offsets
$\delta\equiv i_p-i\pmod p$ with $i\in I_m$; below $p$ these are the
$\delta\in\{i_p-i:i\in I_m,\ i<i_p\}$.

Eviction therefore costs one half-window step and nothing else. The candidate
$k+L$ has the magnitude of $k$, so it stays inside any factorization bound that
contained the source, and no modulus, product, or CRT system is involved. Two
consequences are immediate.

- The product-scale bound (20) below constrains only repairs that *retain*
  their blockers. Dropping them is free, so (20) is not an obstruction.
- The zero-carry system (20a), its exact class count, and the adaptive
  progression (20b) are not needed in order to evict blockers. Their large
  canonical representatives are a consequence of demanding the *sufficient*
  uniform classes (15) instead of the exact small-prime tier of Section 8,
  which the translated candidate usually already satisfies.

`data/exact_eviction_m1_20_k5000_m27.json` checks (23) at every attached
large prime of $1\le m\le20$, $1\le k\le5{,}000$, and instantiates (24) for
every $m\le200$ and every smaller offset.

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
translated bad window, and by (23) the half-window translation always drops
every old prime. The bound is therefore a property of the retention strategy,
not a barrier to repair.

### Zero-carry CRT eviction of prescribed blockers — PROVED

Put

$$
I_m=\{\lfloor m/2\rfloor+1,\ldots,m\},
\qquad L=|I_m|=\lceil m/2\rceil.
$$

Suppose that each distinct prime \(p\) in a finite set \(S\) satisfies
\(p>m\) and divides one old bad-window term \(k+i_p\), with \(i_p\in I_m\).
Let \(M_m\) be the small-prime zero-carry modulus from (8). The system

$$
\boxed{\displaystyle
K\equiv0\pmod {M_m},
\qquad
K\equiv k+L\pmod p\quad(p\in S)} \tag{20a}
$$

is consistent by CRT. Its first congruence gives the complete small-prime
control (9). For every \(j\in I_m\),

$$
K+j\equiv L+j-i_p\pmod p.
$$

The integer on the right lies between \(1\) and \(m\). Since \(p>m\), it
is nonzero modulo \(p\). Hence no old blocker in \(S\) divides any term of
the new bad window.

The count is also exact. Modulo each \(p\in S\), precisely the \(L\)
classes \(-j\), \(j\in I_m\), put \(p\) into the new window. Therefore,
with the small-prime tier fixed to the zero-carry class, there are exactly

$$
\boxed{\displaystyle
\prod_{p\in S}(p-L)}
$$

eviction classes modulo \(M_m\prod_{p\in S}p\). More generally, replacing
the zero class by all uniform local-safe classes from (15) multiplies this
count by their exact CRT class count.

There is a more useful one-dimensional form. Put \(P=\prod_{p\in S}p\).
Every member of

$$
\boxed{\displaystyle K_t=k+L+tP\qquad(t\ge0)} \tag{20b}
$$

evicts all primes in \(S\). If \(\mathcal C\) is any uniform local-safe
family modulo \(Q\) from (15), then \(\gcd(P,Q)=1\), so \(t\mapsto K_t\bmod Q\)
is a bijection over one period. Exactly \(|\mathcal C|\) values of
\(t\bmod Q\) therefore land in the chosen uniform small-prime-safe family.
This gives an exact adaptive search along the blocker-evicting progression,
not merely an existence statement from a single zero-carry class.

This removes every prescribed old large prime while controlling all
\(p\le m\); unlike (20), it does not retain their attachments. It still
does not control primes newly entering the translated window.

The structure certificate applies the zero-carry construction to all 3,846
large-prime-obstructed source pairs in \(1\le m\le8\), \(1\le k\le500\).
Of the 2,737 canonical candidates with \(m+2K\le10^8\), 103 are witnesses and
2,634 fail only at newly entering large primes; 1,109 candidates exceed that
check bound.

For the 12 single-level-sieve survivors in the two \(m=27\) moving-window
artifacts, the zero-carry representatives have 49--98 digits. The exact
adaptive search `data/adaptive_eviction_m27_h200000000.json` instead scans
(20b) against every least-power local-safe class. All 12 systems reach their
first safe \(t\) between 10 and 706, reducing the candidates to 32--81 digits.
Every candidate still exceeds the retained full-factorization bound. A bounded
exact scan through \(10^6\) performs 12,649,506 trial-prime divisibility checks
and rejects one candidate at \(p=137{,}341\); the other 11 new windows remain
unclassified in that search.

**Superseded.** Both digit ranges are artifacts of the uniform sufficient
classes (15), not of eviction. By (23) the candidate \(k+L\) already evicts
every old blocker, and for all 12 survivor systems it also satisfies the exact
small-prime tier of Section 8 with no search at all. Those candidates have 13
digits, stay inside the retained factorization bound, and are classified
completely in Section 12. The retained artifacts remain valid records of what
the sufficient-class route costs; they are not the best available bound.

An abandoned route is worth recording. Classifying the 32--81-digit candidates
directly needs a general large-integer factorization backend, which this
repository deliberately does not carry. That work was started and then dropped
as unnecessary: the exact-tier translation keeps every candidate at 13 digits,
inside the existing trial factorizer, so no external arithmetic dependency is
required for any statement here.

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

### One integer decides every shift of a witness — PROVED

Goodness is monotone in $m$: being $m$-compensation-good constrains only the
primes $p>m$, and $\{p>m+1\}\subset\{p>m\}$, so every $m$-compensation-good
integer is also $(m+1)$-compensation-good. Combining this with the window
identity above and the separation criterion (18) gives an exact reduction for a
witness source $(m,k)$.

- If $m$ is odd, $\mathcal B(m+1,k-1)=\mathcal B(m,k)$, so the target's large
  tier is inherited and automatically satisfied. Then
  $$
  (m+1,k-1)\ \text{is a witness}
  \iff s_p(m+1,k-1)\ge0\ \text{for every prime }p\le m+1 .
  $$
  This recovers Section 5 without any residue computation.
- If $m$ is even, $\mathcal B(m+1,k-1)=\mathcal B(m,k)\cup\{x/2\}$, so
  $$
  (m+1,k-1)\ \text{is a witness}
  \iff s_p(m+1,k-1)\ge0\ \ (p\le m+1)
  \ \text{and}\ x/2\ \text{is}\ (m+1)\text{-compensation-good}.
  $$

The entire large-prime cost of a natural shift is therefore carried by one
integer, and only in the even-to-odd direction. This is the exact structural
reason for the pair pattern of A375071.

`data/parity_shift_m1_20_k4000_published.json` checks the reduction at all 905
witness sources of $1\le m\le20$, $2\le k\le4{,}000$ and at all 26 published
witnesses, with 2,124 and 463 inherited large-prime monotonicity checks. In the
rectangle, 433 of 435 odd-source shifts are witnesses and the other 2 fail only
at small primes, exactly as the first case allows. Of the 470 even-source
shifts, 34 are witnesses, 378 fail only at the adjoined term, 11 fail only at
small primes, and 47 fail at both.

Every published witness obeys the same split: all 13 odd-source shifts are
witnesses, all 13 even-source shifts fail, and 12 of the published successors
are exactly the natural shift of their predecessor. At $m=26$ the adjoined term
is $x/2=5{,}048{,}891{,}644{,}633$, which is prime, so $c=1$, $e=1$, and
$2c>p$ is impossible: one prime integer accounts for the entire failure of the
published chain at the first unlisted index.

## 11. Powersmooth necessity and an unconditional density ceiling — PROVED

### Every compensated prime power is bounded by the term — PROVED

Let $w$ be $m$-compensation-good and let $p>m$ with $p^e\Vert w$, $e\ge1$.
Threshold (16c) gives $w\ge p^e(p^e+1)/2$, hence

$$
\boxed{p^{v_p(w)}<\sqrt{2w}\qquad\text{for every prime }p>m\text{ dividing }w.}
\tag{25}
$$

In particular a single prime factor $p>m$ with $p^2>2w$ is always fatal: then
$e=1$ and the cofactor $c=w/p<p/2$ satisfies $2c<p$, so (16a) gives
$C_p(w)=0<1=e$. Equivalently, **every bad-window term of a witness is
$\sqrt{2w}$-powersmooth away from the primes $p\le m$.** Section 7 gives the
matching sufficient condition, full $m$-smoothness, so a witness window is
sandwiched between $m$-smooth and $\sqrt{2w}$-powersmooth.

The bound (25) is attained: for odd $p$ the term $w=p^e(p^e+1)/2$ is the least
one with $v_p(w)=e$ that compensates, by the sharpness case of (16c).

### The short-cofactor density is $\log 2$ — PROVED

Call $n$ *short-cofactor* when some prime $p\mid n$ has $p^2>2n$. Such a prime
is unique: two of them would give $pq\mid n$ with $pq>2n$. Writing $n=pc$, the
condition $p^2>2n$ is exactly $p>2c$, and then $c<p$ forces $v_p(n)=1$ and
$p=P(n)$. Hence, counting by cofactor,

$$
S(X)=\#\{n\le X:\ n\ \text{short-cofactor}\}
=\sum_{2c^2<X}\bigl(\pi(X/c)-\pi(2c)\bigr).
$$

The subtracted terms contribute $O(X/\log X)$. For the main sum, the prime
number theorem is uniform on $X/c\ge\sqrt X$, so with $t=\log c$

$$
S(X)=(1+o(1))\,X\!\!\int_0^{\frac12\log X}\!\!\frac{dt}{\log X-t}
=(1+o(1))\,X\log 2 .
$$

Therefore the density of integers passing the necessary condition (25) is
$1-\log 2=0.30685\ldots$, and since a witness requires **every** term of its
bad window to pass, the first term alone gives the unconditional ceiling

$$
\boxed{\#\{k\le X:\ (m,k)\ \text{is a witness}\}\le(1-\log 2+o(1))X}
\tag{26}
$$

for every fixed $m\ge1$. This is the first density statement in this file. It
is not an obstruction: $1-\log 2>0$, so (26) is consistent with infinitely many
witnesses for every $m$. Its content is that a positive proportion of $k$ is
excluded by one term of one exactly computable local condition, and that the
exclusion is quantitative rather than heuristic.

No multiplicative strengthening of (26) is claimed here. Requiring all
$\lceil m/2\rceil$ terms to pass simultaneously is a correlation problem for
large prime factors of consecutive integers, and the available elementary
bounds degrade past two terms.

`data/smooth_density_m1_20_k3000.json` checks (25) at 331,361 large window
prime powers over $1\le m\le20$, $1\le k\le3{,}000$, including 216,747 powers
with $p^{2e}>2w$ that are all fatal, and 1,706 exactly tight minimal terms.
It computes $S(X)$ for $X=10^5,10^6,2\cdot10^6$ by two independent exact
methods, obtaining densities $0.67277$, $0.679867$, and $0.68121$ against the
limit $\log 2=0.693147\ldots$, with the signed deviation shrinking as expected.

### Dominant-prime repair has density zero — PROVED

Let $p^e\Vert w$ with $p>m$ and short cofactor $c=w/p^e<p$. By (16a),
$C_p(w)=\mathbf 1_{\{2c>p\}}$, so $w$ is good at $p$ only when $e=1$ and
$p<2c$, that is

$$
\boxed{w<p^2<2w.} \tag{27}
$$

A dominant prime therefore compensates only inside the narrow window
$p\in(\sqrt w,\sqrt{2w})$. Counting by cofactor, with $n=pc$ and $c<p<2c$,

$$
\#\{n\le X:\ \exists p,\ n=pc,\ c<p<2c\}
=\sum_{c<\sqrt X}
\Bigl(\pi\bigl(\min(2c-1,X/c)\bigr)-\pi(c)\Bigr)
\ll\sum_{c\le\sqrt X}\frac{c}{\log c}
\ll\frac{X}{\log X}.
$$

So apart from a set of density $O(1/\log X)$, an $m$-compensation-good integer
has no prime factor $p>m$ with $p^2>w$ at all. The $\log 2$ proportion of
integers carrying a short-cofactor prime is almost entirely lost rather than
repaired, and (25) fails almost only in its unrepairable direction. This is why
the single-level sieve of Section 3 rejects essentially every candidate in the
moving-window searches before any compensation test is needed.

The certificate counts the narrow window exactly at $X=10^5,10^6,2\cdot10^6$,
obtaining densities $0.06044$, $0.051792$, and $0.0493005$, consistent with the
$1/\log X$ decay, and repeats the smallest count by direct factorization.

## 12. Exact finite evidence and current gaps

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
1,104 are witnesses. It checks the short-cofactor formula at 405,780 factors,
including 363,792 obstructions and 7,173 repeated-prime spikes. It exhaustively
enumerates 19,804 residues for the prefix cone and checks the
even-shift/first-term identity at 19,654 primes.

The same certificate combines the moving-window artifacts to show that all
200,000,000 offsets have a large-prime obstruction. Retaining the old blockers
in the first three survivors forces 69-, 30-, and 35-digit prime products; the
nine later survivors force 36--79 digits. On the smaller certified rectangle,
2,634 of 2,737 exactly checked canonical zero-carry evictions acquire new
large-prime blockers.

`data/exact_eviction_m1_20_k5000_m27.json` replaces both of those repair
directions by the half-window translation (23). It verifies (23) at 587,169
attached large primes over $1\le m\le20$, $1\le k\le5{,}000$ through 3,871,667
divisibility checks, instantiates the sharpness statement (24) for every
$m\le200$ with 9,900 explicit surviving blockers, and checks the residue core
at 676,700 window pairs.

It then classifies the translated candidate exactly. Of the 98,706
large-prime-obstructed sources in that rectangle, all 98,706 translated
candidates evict every old blocker; 1,033 are witnesses, 59,641 fail only at
newly entering large primes, 162 fail only at the small-prime tier, and 37,870
fail at both. On the 9,832 candidates with $k\le500$ the separation-based
classification is confirmed by complete Legendre and Kummer certificates.

For the 12 recorded $m=27$ survivor systems the outcome is decisive. Each
translated candidate $k+14$ has 13 decimal digits, evicts all 140 recorded old
blockers across 1,960 divisibility checks, and satisfies the exact small-prime
tier for every $p\le27$ with no search. All 12 lie inside the retained
factorization bound, so their translated windows are factored completely: 369
large-prime factors, of which 9 to 15 per candidate are exact obstructions.
Every one of the 12 is therefore classified as `new_large_prime_only`. The
earlier 32--81-digit and 49--98-digit candidates and the 11 unclassified cases
are superseded, and only the new-prime side of the problem remains.


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

The same exact criterion now also covers the scale that the shift target lives
at. Four contiguous shards
`data/compensation_run_m27_k5049091644619_h200000000.json`,
`...k5049291644619...`, `...k5049491644619...`, and `...k5049691644619...`
extend the two moving-window artifacts forward by $8\cdot10^8$ candidates. Each
shard rejects all $2\cdot10^8$ of its offsets by an exact negative local slack
at some $p>27$, with no offset ever reaching the small-prime tier, and the
longest compensation-good runs are $8$, $9$, $9$, and $8$. Their union with the
two moving artifacts is contiguous, so

$$
\boxed{\text{no }m=27\text{ witness has }
5{,}048{,}891{,}644{,}619\le k<5{,}049{,}891{,}644{,}619,}
$$

exactly $10^9$ consecutive candidates starting at the published $m=26$
natural-shift target. This is five times the previously certified
neighbourhood, and unlike the earlier sweep every rejection here comes from the
exact compensation criterion rather than the single-level sieve alone. The four
shards used one reduced-priority process each, 717.99--725.21 s of wall time,
for a measured $2.8\cdot10^5$ candidates per second at this scale.

Their combined compensation-good density is $0.123936$, independently
reproducing the $0.123664$ measured in the fresh block below, and the
term-independence model predicts a longest run of $9.03$ against the observed
$8,9,9,8$.

The same artifacts also fix the size of the barrier. Their exact term counts
give a $27$-compensation-good density of $0.133708$ over
$[50{,}015,\,100{,}050{,}027]$ and $0.129223$ over
$[100{,}050{,}015,\,1{,}000{,}050{,}027]$, both far below the proved ceiling
$1-\log 2$ of (26). `data/smooth_density_m1_20_k3000.json` measures the same
density in fresh blocks of 250,000 consecutive integers: $0.142756$ near
$10^6$, $0.124948$ near $10^9$, and $0.123664$ near $5\cdot10^{12}$. The exact
small-prime tier passes for 19,920 of 20,000 consecutive $k$ at the
$5\cdot10^{12}$ scale, so it is not the binding constraint anywhere in this
range.

Run lengths match a term-independence model closely. In the two run artifacts
the model predicts longest runs of $9.01$ and $9.94$ against exact values $8$
and $9$. In the fresh $5\cdot10^{12}$ block the observed spectrum of runs of
length $1$ through $6$ is $19{,}176$, $2{,}369$, $266$, $32$, $6$, $1$ against
model values $19{,}079$, $2{,}374$, $295$, $37$, $4.6$, $0.6$.

One measurement in this session initially contradicted that model: a block
beginning at $5{,}048{,}891{,}644{,}635$ contained a run of $12$, which the
model puts at $2\cdot10^{-6}$ expected occurrences. The cause was sampling
bias, not structure. That block starts inside the bad window of the published
$m=25$ witness $k=5{,}048{,}891{,}644{,}621$, whose 13 window terms are
$25$-compensation-good, hence also $27$-compensation-good because goodness at
$m$ constrains only the primes $p>m$. The certificate now keeps that block as a
labelled control and asserts the full run of 13, next to an unbiased block of
equal length whose longest run is 5. The published data therefore already
supplies 13 of the 14 consecutive good integers that $m=27$ requires, and the
missing one is the term adjoined on the left by the even-to-odd shift of
Section 10.

Extrapolation, **not proved**, places the barrier well past every exhausted
range. Regressing $\log k_m$ on $\lceil m/2\rceil$ over the 26 published
witnesses gives slope $2.1899$ with $R^2=0.989$, an implied per-term density
$0.1119$ consistent with the measured $0.1237$, and a predicted least
$m=27$ witness near $2.8\cdot10^{13}$. The independence model with the measured
densities predicts $5.1\cdot10^{12}$. Both estimates exceed the exhausted
$10^9$ rectangle by three to four orders of magnitude and exceed the
$2\cdot10^8$ offsets swept around the natural-shift target, so the negative
results in this file are exactly what either estimate predicts.

Open gaps:

1. The exact target is a run of $\lceil m/2\rceil$ consecutive
   $m$-compensation-good integers. Section 13 removes the other half of the
   old target: the small-prime tier fails on a set of density zero by (32), so
   existence depends on the run alone, and no band or CRT class can obstruct.
   No theorem proves that such a run exists for a given $m$. Ceiling (26)
   bounds its density from above but cannot establish nonemptiness, and
   Conjecture R of Section 13 is exactly what is missing.
2. The Dirichlet construction proves that no fixed modulus can be sufficient;
   a successful construction must adapt to the translated bad window.
3. Old blockers are no longer part of the difficulty. Translation (23) evicts
   them at zero cost and (24) shows the half-window offset is optimal, while
   formula (16), cone (19), and system (20a) control retained or prescribed
   primes. None of these mechanisms controls the primes newly entering the
   translated window, and the exact classification of 98,706 translated
   candidates shows that this is where every failure now occurs.
4. Necessity (25) shows each window term must be $\sqrt{2w}$-powersmooth away
   from $p\le m$. Proving that $\lceil m/2\rceil$ consecutive integers can meet
   that condition simultaneously is a correlation problem for large prime
   factors of consecutive integers and is not solved here.
5. No construction here reaches a previously unknown $m$; Erdős #389 remains
   open.

## 13. The small-prime tier is asymptotically free — PROVED

Section 9 splits the witness condition exactly into the small-prime tier
$s_p(m,k)\ge0$ for $p\le m$ and a run of $\lceil m/2\rceil$ consecutive
$m$-compensation-good bad-window terms. This section proves that the first
half of the split is not an obstruction: for fixed $m$ the tier holds for all
but a density-zero set of $k$, at an explicit polynomial rate. What remains of
Erdős #389 is exactly the run.

### A digit-sum identity for the slack — PROVED

Write $S_p(n)$ for the base-$p$ digit sum. Kummer's theorem in Legendre form
gives $v_p\binom{a+b}{a}=(S_p(a)+S_p(b)-S_p(a+b))/(p-1)$, so

$$
\boxed{s_p(m,k)=\frac{2S_p(m+k)-S_p(m+2k)-S_p(m)}{p-1}.}\tag{28}
$$

Indeed $s_p=v_p\binom{m+2k}{k}-v_p\binom{m+k}{m}$, the first term equals
$(S_p(k)+S_p(m+k)-S_p(m+2k))/(p-1)$ and the second
$(S_p(m)+S_p(k)-S_p(m+k))/(p-1)$; the two copies of $S_p(k)$ cancel. Once $m$
and $p$ are fixed, (28) makes the slack a statement about the digits of $k$
alone, which is what allows exact counts far beyond any enumeration.

### A carry lower bound — PROVED

Let $\mathrm d_j(n)$ be the $j$-th base-$p$ digit of $n$, let $D$ be the
number of base-$p$ digits of $m$ (so $D=0$ when $m=0$), let

$$
Z=\min\{z\ge0:\ \mathrm d_{D+z}(k)\ne p-1\},\qquad
W=\#\{j>D+Z:\ 2\,\mathrm d_j(k)\ge p\}.
$$

Then

$$
\boxed{s_p(m,k)\ \ge\ W-D-Z.}\tag{29}
$$

Write $c_1$ for the carry count of $m+k$ and $c_2$ for that of $k+(m+k)$, so
that $s_p=c_2-c_1$ by Section 1.

*Upper bound on $c_1$.* A carry leaves position $j$ of $m+k$ only when
$\mathrm d_j(m)+\mathrm d_j(k)+\text{carry in}\ge p$. For $j\ge D$ the first
summand vanishes, so a carry leaves position $j$ only if one enters and
$\mathrm d_j(k)=p-1$. At $j=D+Z$ the digit is not $p-1$, so no carry leaves
it, and by induction none leaves any higher position. Every carry therefore
occurs at one of the $D+Z$ positions below $D+Z$, giving $c_1\le D+Z$.

*Lower bound on $c_2$.* For $j>D+Z$ no carry of $m+k$ enters position $j$ and
$\mathrm d_j(m)=0$, hence $\mathrm d_j(m+k)=\mathrm d_j(k)$. Position $j$ of
$k+(m+k)$ then receives $2\mathrm d_j(k)+\text{carry in}\ge2\mathrm d_j(k)$,
which emits a carry whenever $2\mathrm d_j(k)\ge p$. Hence $c_2\ge W$.

Subtracting gives (29). The inequality is attained: the certificate reports
minimum margin $0$ over its scanned rectangle, so the constants cannot be
removed.

### Exact failure counts and a closed-form bound — PROVED

Put $h=\#\{d:0\le d<p,\ 2d\ge p\}$, so $h=1$ for $p=2$ and $h=(p-1)/2$ for
odd $p$, and

$$
T(n,B)=\sum_{i=0}^{\min(B,n)}\binom{n}{i}h^i(p-h)^{\,n-i}.
$$

By (29) a failure $s_p(m,k)<0$ forces $W\le D+Z-1$. Splitting on the value of
$Z$ and counting digits directly,

$$
\boxed{\#\{k<p^J:\ s_p(m,k)<0\}\ \le\ p^D\Bigl[(p-1)\sum_{z=0}^{J-D-1}
T(J-D-z-1,\ D+z-1)+1\Bigr].}\tag{30}
$$

For $Z=z<J-D$ the digits below $D$ are free ($p^D$ choices), positions
$D,\dots,D+z-1$ are pinned to $p-1$, position $D+z$ has $p-1$ admissible
values, and the $J-D-z-1$ digits above it must carry at most $D+z-1$ high
digits, which is exactly $T(J-D-z-1,D+z-1)$. The trailing $p^D$ covers the
degenerate case $Z=J-D$. Dividing by $p^J$ turns the digit count into a
binomial law:

$$
\frac{\#\{k<p^J:\ s_p<0\}}{p^J}\le
\sum_{z\ge0}\frac{p-1}{p^{z+1}}
\Pr\bigl[\mathrm{Bin}(J-D-z-1,\ h/p)\le D+z-1\bigr]+p^{D-J}.\tag{30a}
$$

### Exponential decay — PROVED

Every prime has $h/p\ge1/3$, with equality at $p=3$ and value $1/2$ at $p=2$.
Let $J\ge16(D+1)$. For $z\le J/16$ the binomial in (30a) has
$n=J-D-z-1\ge\frac78J$ trials, mean $nh/p\ge0.29J$, and threshold
$D+z-1\le J/8$, so Hoeffding's inequality gives
$\Pr\le\exp(-2(nh/p-J/8)^2/n)\le\exp(-0.052J)$. The terms with $z>J/16$ are
dominated by the geometric tail $p^{-J/16}\le2^{-J/16}$, and
$p^{D-J}\le2^{-15J/16}$. Summing the three contributions,

$$
\boxed{\frac{\#\{k<p^J:\ s_p(m,k)<0\}}{p^J}\le4e^{-J/24}
\qquad(J\ge16(D+1)),}\tag{31}
$$

equivalently $\le4X^{-1/(24\log p)}$ at $X=p^J$. The certificate checks (31)
against the exact counts on every reported cell whose hypothesis holds.

### The tier is not the obstruction — PROVED

Let

$$
G_m=\{k\ge1:\ \text{every }w\in\mathcal B(m,k)\text{ is }
m\text{-compensation-good}\}.
$$

Section 9 says $(m,k)$ is a witness exactly when $k\in G_m$ and
$s_p(m,k)\ge0$ for every $p\le m$. By (31) the second condition fails on a set
of density zero,

$$
\#\{k\le X:\ \exists\,p\le m,\ s_p(m,k)<0\}\le
4\pi(m)\,m\,X^{1-1/(24\log m)}=o(X),\tag{32}
$$

so

$$
\#\{k\le X:\ (m,k)\text{ a witness}\}\ \ge\
\#\bigl(G_m\cap[1,X]\bigr)-o(X).
$$

Two consequences. First, if $G_m$ has positive upper density then $(m,k)$ is a
witness for infinitely many $k$. Second, the prime-power band of Section 8 and
the small-prime CRT classes are devices for *constructing* witnesses, not
constraints on their existence: they cannot obstruct, because their failure set
is thin. The triple intersection that Section 12 named as the target collapses
to a single question — does a run of $\lceil m/2\rceil$ consecutive
$m$-compensation-good integers exist?

### The reduction is finitary — PROVED

(32) is an explicit power of $X$, not merely $o(X)$, so the density hypothesis
above can be replaced by an inequality at a *single* scale.

**Corollary.** Fix $m$ and put $c_m=1/(24\log m)$. If for some $X$

$$
\#\bigl(G_m\cap[1,X]\bigr)\;>\;4\pi(m)\,m\,X^{1-c_m},
$$

then $(m,k)$ is a witness for some $k\le X$; if it holds for infinitely many
$X$, for infinitely many $k$.

*Proof.* By (32) the right-hand side bounds the number of $k\le X$ carrying a
tier failure, so some $k\in G_m\cap[1,X]$ carries none, and by Section 9 that
$k$ is a witness. $\square$

This replaces an asymptotic density hypothesis by one finite inequality, and it
is strictly weaker: a run count $\gg X^{1-c}$ with $c<c_m$ suffices, where
positive density demands $\gg X$. The honest cost is where it starts to bite.
For $m=27$, $c_{27}=0.0126$ and $4\pi(m)m=972$, so the criterion is vacuous
below $X=10^{236.3}$ — the bound exceeds $X$ itself — and it beats the measured
run density $(0.123664)^{14}=10^{-12.71}$ only past $X=10^{1241.6}$. It is
therefore a structural sharpening, not a computational one: it says the tier is
free *quantitatively*, and it converts Conjecture R into a target that a single
scale can meet, but no reachable scale meets it by counting alone.

### A conditional resolution

Goodness is monotone in $m$ (Section 10): a $1$-compensation-good integer is
$m$-compensation-good for every $m\ge1$. So define

> **Conjecture R($L$).** For every $L\ge1$ the set of $k$ such that
> $k+1,\dots,k+L$ are all $1$-compensation-good has positive upper density.

Conjecture R($\lceil m/2\rceil$) plus (32) gives infinitely many witnesses for
that $m$; Conjecture R in full gives Erdős #389 for every $m$. The reduction
is unconditional and exact — only R is open.

R is a correlation statement for large prime factors of consecutive integers.
Its first rung $R(1)$ — positive lower density for a single term — **is now
proved in Section 17**; what follows records the two elementary routes tried
here, one of which is genuinely closed and one of which was misread.

*The union bound cannot work asymptotically.* Short-cofactor failures are
pairwise disjoint, since two primes with $p^2,q^2>2w$ would give $pq\mid w$
with $pq>2w$, so their total density is exactly $\log2$ by (26) and is
irreducible. Grouping the level-failure events by their number of usable
levels, a prime with exactly $r$ of them satisfies
$X^{1/(r+2)}<p\le X^{1/(r+1)}$, and Mertens gives $\sum1/p=\log\frac{r+2}{r+1}$
over that range, so the modelled level mass tends to

$$
\sum_{r\ge1}2^{-r}\log\frac{r+2}{r+1}=0.32252\ldots,
$$

for a union total $\log2+0.32252=1.01567>1$. Fixed small primes contribute
nothing in the limit because their level counts grow with $X$. At finite scale
the total dips below one — $0.96981$ at $1.2\cdot10^6$, $0.98660$ at
$1.02\cdot10^7$, $1.00659$ at $10^{10}$ — but there exact counting is stronger
than any bound, so the dip is worthless.

*The obvious repair is false.* One would keep only the level failures that
avoid the short-cofactor class, hoping for a factor $1-\log2$. Exact
classification of every integer in a block refutes this. At $10^6$ the four
disjoint classes have densities $0.144500$ good, $0.629345$ short-cofactor
only, $0.172740$ level-failure only, $0.053415$ both; independence would
predict $0.071745$ for the level-only class, so the observed value is $2.41$
times larger ($2.35$ at $10^7$). The two failure modes are adversely
correlated: a prime factor above $\sqrt{2w}$ leaves a cofactor below
$\sqrt{w/2}$, which has few primes available to fail a level test.

**[RETRACTED — this paragraph previously concluded "the overlap that would
rescue the union bound does not exist". That is false, and Section 17 proves
the opposite. The refutation above is of the *independence* prediction
$0.071745$; the union bound needs only $0.015668$, and the measured overlap is
$0.053415$ — already $3.4$ times enough. The two statements were conflated.]**

Partitioning by the largest prime factor does remove the overlap honestly, but
then requires equidistribution of $p$-smooth cofactors in arithmetic
progressions modulo $p$ — moduli of size about $\sqrt X$ — which is
Bombieri–Vinogradov-strength input and is not attempted here. Section 17 avoids
it by bounding the overlap from below on the *prime* side instead, where a
single-modulus exponential-sum bound suffices.

### Exact evidence

`data/tier_density_m1_30.json` verifies (28) and (29) at 36,400 triples over
$1\le m\le20$, $1\le k\le400$, all $p\le m$, with minimum margin $0$; agrees
with direct enumeration and with a second prefix-decomposition program on 137
cells; and checks (30) on every cell and (31) on the 13 cells whose hypothesis
holds. For $m=27$ the exact tier-failure densities are

| $X$ | union over $p\le27$ | surviving |
|---|---|---|
| $10^6$ | $0.382703$ | $0.617297$ |
| $5.048\,891\,644\,621\cdot10^{12}$ | $0.0132715$ | $0.9867285$ |
| $10^{13}$ | $0.0123698$ | $0.9876302$ |

and the per-prime density at $p=2$ falls $0.0569\to2.86\cdot10^{-4}\to
1.19\cdot10^{-5}\to2.06\cdot10^{-8}\to6.18\cdot10^{-14}\to3.86\cdot10^{-23}$
across $X=2^{20},2^{45},2^{60},2^{90},2^{150},2^{250}$ (ranges of 7, 14,
19, 28, 46 and 76 decimal digits). Measured decay
exponents run from $0.187$ at $p=23$ to $0.306$ at $p=2$.

Two scales are worth separating. At $X=2^{250}$, a 76-digit range, the exact union
of tier-failure densities already lies below the measured
$(0.123664)^{14}=1.96\cdot10^{-13}$ run density, so past that point the tier
cannot be the binding constraint; the same comparison run through the provable
envelope (31) instead of the exact counts needs $X=2^{3343}$, a 1007-digit
range. At the
scale where the least $m=27$ witness is expected, near $5\cdot10^{12}$, the
tier still removes $1.3\%$ of all $k$, so it remains a real filter for search
even though it cannot obstruct existence.


## 14. Congruence-forced powersmooth mass costs modulus — PROVED

Sections 8 and 9 repair prescribed primes by prescribing residues. This
section prices the *mass* such a prescription must carry. Fix a modulus
$M=\prod_p p^{a_p}$ and an integer $k$, and split each bad-window term into
the part that $k\bmod M$ pins and the part it leaves free:

$$
d_i=\prod_{p\mid M}p^{\min(v_p(k+i),\,a_p)},\qquad u_i=\frac{k+i}{d_i}
\qquad(i\in I_m).
$$

Call the pair $(M,k)$ **size-forcing** when the pinned parts alone already meet
the powersmooth threshold of (25) on every term,

$$
u_i<\sqrt{2(k+i)}\qquad\text{for every }i\in I_m,
$$

so that no unforced cofactor can carry a fatal prime whatever it turns out to
be. This is a property of the single pair $(M,k)$; the scope note below
separates it from the stronger reading in which one class must work for
unbounded $k$.

For every size-forcing pair, with $L=\lceil m/2\rceil$,

$$
\boxed{\log M\ \ge\ \frac L4\log\frac k2-\frac L2\bigl(\log L+1\bigr).}
\tag{33}
$$

*Mass needed.* $u_i<\sqrt{2w_i}$ means $d_i>w_i/\sqrt{2w_i}=\sqrt{w_i/2}\ge
\sqrt{k/2}$, so $\sum_i\log d_i>\frac L2\log\frac k2$.

*Mass available.* A window of $L$ consecutive integers contains at most
$1+L/p^j$ multiples of $p^j$, so

$$
\sum_{i\in I_m}\min(v_p(w_i),a_p)=\sum_{j=1}^{a_p}\#\{i\in I_m:p^j\mid w_i\}
\le a_p+\frac{L}{p-1},
$$

giving $\sum_i\log d_i\le\log M+L\sum_{p\mid M}\log p/(p-1)$.

*Correction.* For $p\le L$, Mertens gives
$\sum_{p\le L}\log p/(p-1)\le\log L+1$. For $p>L$ we have $L/(p-1)\le1$, so
those terms contribute at most $\sum_{p\mid M}\log p\le\log M$. Hence
$\frac L2\log\frac k2<2\log M+L(\log L+1)$, which is (33).

### What (33) does and does not say

The hypothesis is satisfiable, so (33) is not vacuous: the certificate builds
size-forcing moduli out of real factorizations for the 18 published witnesses
with $k\le10^9$ and checks (33) on each.

Two readings must be kept apart, and only the second is claimed here.

*A fixed class with $k$ unbounded.* Nothing is left for (33) to do. Section 8
already shows by Dirichlet that for any fixed $M$, any class $c$, and any
single window position, the unforced cofactor takes prime values above
$\sqrt{2w}$ for infinitely many $k\equiv c$, so no fixed class certifies even
one term. Read this way no pair in the class is eventually size-forcing, the
hypothesis of (33) is empty, and the inequality says nothing. That reading is
not used.

*One dyadic range $k\in[N,2N)$.* Here the hypothesis has content, because a
class can be size-forcing at a particular $k$ of the range. Then (33) gives
$\log M\ge\frac L4\log\frac N2-\frac L2(\log L+1)$, so $M\ge N^{L/4-o(1)}$: for
$L\ge5$ the modulus exceeds the range, the class meets $[N,2N)$ at most once,
and it names the witness instead of predicting it. The certificate lists the
exact thresholds at which the right side of (33) passes $\log k$: $k\ge70{,}247$
for $m=27$, $k\ge4{,}708{,}708$ for $m=13$, and none for $L\le4$ (that is
$m\le8$), so cheap forcing survives only at the smallest indices. Applied to
the published witnesses, a size-forcing modulus at
$m=26,\ k=5{,}048{,}891{,}644{,}620$ needs $\log M\ge69.64$ — at least $30$
decimal digits against a $13$-digit target — and every published witness with
$m\ge13$ behaves the same way.

Small moduli are not the business of (33) at all. Section 15 evicts every class
of modulus $M\le\sqrt N/2$ from the block by an elementary construction that
needs no prime localisation, so the surviving question for (33) is only what
happens between $\sqrt N$ and the size-forcing scale $N^{L/4}$. Everything this
file has to say about prime distribution now lives in Section 15; (33) is pure
mass accounting and claims nothing about primes.

What is settled is the accounting: every mechanism of Sections 8 and 9 that
works by prescribing residues — mirror lifts (17), the repair cone (19), the
eviction systems (20a), (20b) — pays $\frac L4\log\frac k2$ in modulus to
force the window by size. The half-window translation (23) escapes because it
prescribes nothing: it moves $k$ by $\lceil m/2\rceil$ and evicts by adjacency.

### Exact evidence

`data/forcing_mass_bound.json` verifies the counting lemma for every residue
class on $1\le m\le30$, $p\in\{2,3,5,7,11\}$, $a\le4$ — 595,020 residues, with
the tightest case leaving slack $0.1$ at $(m,p,a)=(1,11,4)$ — and verifies
$\sum_{p\le P}\log p/(p-1)\le\log P+1$ at every prime $P\le2\cdot10^6$ with
minimum margin $0.856$. For the 18 published witnesses with $k\le10^9$ it also
builds an explicit certifying modulus from the real factorization of every
window term and checks that it dominates (33); at $m=13$, $k=7{,}979{,}077$ the
bracket is $10^7$ from below and $10^{33}$ from above, a factor $4.6$ in the
logarithm, so (33) is not vacuous.


## 15. Every long progression is evicted elementarily — PROVED

Section 14 prices a class that forces its window *by size*. This section
attacks any arithmetic progression directly and shows that no congruence
condition protects its members for long. One theorem does all of the work, and
it uses no prime distribution statement whatsoever — only Bertrand's postulate
and one modular inversion.

### 15.1 The class splitting — exact

Fix $m$, a modulus $M=\prod_p p^{a_p}$ and a class $c\bmod M$. For $i\in I_m$
put

$$
d_i=\prod_{p\mid M}p^{\min(v_p(c+i),\,a_p)},\qquad
q_i=\frac M{d_i},\qquad
a_i\equiv\frac{c+i}{d_i}\pmod{q_i}.
\tag{34}
$$

*The splitting is a property of the class, not of $k$.* If $k\equiv c\pmod M$
and $p^{a_p}\Vert M$ then either $p^{a_p}\mid c+i$, and then $p^{a_p}\mid k+i$,
or $v_p(c+i)<a_p$, and then $v_p(k+i)=v_p(c+i)$. Either way
$\min(v_p(k+i),a_p)=\min(v_p(c+i),a_p)$, so $d_i$ pins the same part of every
term of the class. Writing $u_i=(k+i)/d_i$ for the free cofactor,
$u_i\equiv a_i\pmod{q_i}$ and $\gcd(a_i,q_i)=1$: for $p\mid q_i$ we have
$v_p(c+i)<a_p$, hence $v_p(d_i)=v_p(c+i)$ and $p\nmid u_i$.

Note $d_iq_i=M$ for every $i$. **Choosing a different window position
redistributes the modulus between pinned and free parts; it never lowers the
product.** No route may expect a smaller modulus from the $i$-freedom alone.

### 15.2 Fatality is a factorisation, not a primality test — PROVED

**Lemma.** Let $i\in I_m$ and suppose

$$
k+i=d\,p\,t,\qquad p\ \text{prime},\qquad p>m,\qquad p>2dt .
\tag{35}
$$

Then $k$ is not a witness at $m$.

*Proof.* $p>2dt$ gives $p^2>2dpt=2(k+i)$, so $p^{v_p(k+i)}\ge p>\sqrt{2(k+i)}$,
which violates (25) at the bad-window term $k+i$. $\square$

The content of (35) is that fatality never asks for $u_i$ to *be* prime, only
to *have* a prime factor larger than its cofactor doubled — a condition of
positive density $\log2$ rather than density $1/\log$. That is the entire
reason the theorem below needs no analytic input.

### 15.3 The master theorem — PROVED

**Theorem.** Let $m\ge1$ and let

$$
A=\{\,s+jM\ :\ 0\le j<\Lambda\,\},\qquad
\mathrm{top}(A)=s+(\Lambda-1)M+m,
$$

be any arithmetic progression of positive integers, and let $P_0$ be the least
prime exceeding $\sqrt{2\,\mathrm{top}(A)}$. If

$$
\boxed{\Lambda\ \ge\ P_0,\qquad
P_0=\min\{p\ \text{prime}\ :\ p^2>2\,\mathrm{top}(A)\ \text{and}\ p>m\},}
\tag{36}
$$

then $A$ contains an integer that is **not** a witness at $m$.

*Proof.* Put $p=P_0$ and fix any $i\in I_m$. Apply (34) to the class
$s\bmod M$: $d=d_i$, $q=M/d$, $a\equiv(s+i)/d$, and write $u_0=(s+i)/d$. The
window terms of $A$ at position $i$ are $s+i+jM=d\,(u_0+jq)$.

First, $p>M\ge q$. Indeed $(\Lambda-1)M<\mathrm{top}(A)<p^2/2$ and
$\Lambda\ge p$ give $M<p^2/\bigl(2(p-1)\bigr)\le p$ for $p\ge3$. Hence
$\gcd(q,p)=1$ and we may choose

$$
j\equiv-u_0\,q^{-1}\pmod p,\qquad 0\le j<p\le\Lambda .
$$

Then $p\mid u_0+jq$; set $k=s+jM\in A$, $t=(u_0+jq)/p\ge1$, so that
$k+i=d\,p\,t$. Finally $p^2>2\,\mathrm{top}(A)\ge2(k+i)$, which is precisely
$p>2dt$, and $p>m$ holds by the definition of $P_0$. Lemma 15.2 applies and $k$
is not a witness. $\square$

The condition $p>m$ in (36) is not decoration. It is (25)'s own hypothesis, and
$p^2>2\,\mathrm{top}(A)$ does not imply it at small scales: at $m=5$ with
$A=\{1,\dots,5\}$ the size condition alone offers $p=5$, which (25) says nothing
about. The exhaustive screen of 15.8 rejected an earlier version of this proof
at exactly that pair.

The proof uses exactly two facts about primes: Bertrand's postulate, to know
that $P_0$ exists and is at most $2\sqrt{2\,\mathrm{top}(A)}$, and the
invertibility of $q$ modulo $p$. No prime is asked to lie in a progression, in
a short interval, or anywhere else.

### 15.4 Corollaries

**(a) Effective eviction of a class.** For every $m$, $M$ and class $c\bmod M$
there is a non-witness $k\equiv c\pmod M$ with

$$
\boxed{1\le k<2M\max(2M,m)\qquad\bigl(<4M^2\ \text{once}\ 2M\ge m\bigr).}
\tag{37}
$$

*Proof.* Take $s\in[1,M]$ in the class, $\Lambda=p$ where $p$ is the least
prime above $\max(2M,m)$; then $\mathrm{top}(A)<pM+M+m<p^2/2$ because $p>2M$,
so (36) holds and the evicted $k<s+pM\le2M\max(2M,m)$ by Bertrand. $\square$

**(b) Eviction inside a prescribed block.** If

$$
\boxed{M\cdot P_0(N,m)\le N,\qquad
P_0(N,m)=\min\{p:p^2>4N+2m\ \text{and}\ p>m\},}
\tag{38}
$$

then every class $c\bmod M$ contains a non-witness $k\in[N,2N)$: the class
meets the block in $\Lambda\ge N/M\ge P_0$ terms with
$\mathrm{top}\le2N+m$. Since $P_0<2\max(\sqrt{4N+2m},m)$, the condition holds
for $M\le N/\bigl(2\max(\sqrt{4N+2m},m)\bigr)$, and in practice up to
$M\approx\sqrt N/2$.

The clause $p>m$ is not inherited automatically. When $m^2>4N+2m$ the size
condition alone offers a prime below $m$ — at $m=100$, $N=1000$ it offers $67$ —
and (25) exempts such a prime, so the constructed $k$ would not be evicted while
(38) read without the clause still admits $M\le14$. Every rectangle screened in
15.8 has $N\ge10^6\gg m^2/4$, where the clause is vacuous, so no screen could
have exposed this; the producer has always carried it
(`block_prime` calls `progression_prime`), and a direct test now pins the
$m$-dominated regime.

**(c) The witness set contains no long progression.** Contrapositive of the
master theorem: for every $m$, the set of witnesses contains no arithmetic
progression of length $P_0$ all of whose terms are at most $Y$ — *whatever the
common difference*. Once $2Y\ge m^2$ the size clause governs and that length is
$2\sqrt{2Y}(1+o(1))$; below it the bound is the trivial $m+O(1)$. In particular a class that is entirely
witnesses throughout $[N,2N)$ has $M>N/P_0\sim\sqrt N/2$ and so meets the block
in fewer than $2\sqrt N$ integers, relative density $O(N^{-1/2})$.

**(d) A congruence programme needs $\sqrt N$ classes.** If
$c_1\bmod M_1,\dots,c_r\bmod M_r$ each consist entirely of witnesses at $m$
inside $[N,2N)$, they cover fewer than $2r\sqrt N$ integers of the block, so
covering a proportion $\delta$ of it requires

$$
\boxed{r\ \ge\ \tfrac12\,\delta\sqrt N .}
\tag{39}
$$

(d) is the effective replacement for the qualitative statement of Section 8.
Section 8 says no *single* fixed class works and gives no bound; (39) prices the
entire class-based strategy — mirror lifts (17), the repair cone (19), the
eviction systems (20a), (20b), any CRT repair menu — at $\Omega(\sqrt N)$
classes per block. A construction that wants positive density of witnesses
cannot be a finite congruence system, and one that wants a single witness must
locate it rather than prescribe it.

**(e) A covering criterion — CERTIFIED FINITE per instance.** Eviction by a
single prime uses only $p\mid k+i$, so the pinned parts of (34) cancel. In the
coordinate $x\equiv c\,p^{-1}\pmod M$ the classes evicted by $p$ are exactly
the union of the $L$ intervals $Z_i-i\,p^{-1}$, where
$Z_i=\{z:N\le pz-i<2N\}$ holds $\lfloor N/p\rfloor$ integers. If that union is
all of $\mathbb Z/M\mathbb Z$, then $p$ alone evicts **every** class. Testing it
costs $O(L\log L)$ operations, and the union has measure at most
$L\lfloor N/p\rfloor$, so the criterion reaches

$$
M\ \le\ L\Bigl\lfloor\frac N{P_0}\Bigr\rfloor\ \approx\
\Bigl\lceil\frac m2\Bigr\rceil\cdot\frac{\sqrt N}2 ,
\tag{39b}
$$

a factor $\lceil m/2\rceil$ beyond (38) — and, by the measure bound, no
further. A sufficient deterministic form: if $\delta=p^{-1}\bmod M$ satisfies
$\delta\le\lfloor N/p\rfloor$ and $\lfloor N/p\rfloor+(L-1)\delta\ge M$ then the
$L$ intervals form a single chain that wraps.

Bisected at $N=10^{12}$ with a budget of 2,000 primes, the criterion holds up
to $\theta=0.5145$ for $m=5$, $0.5451$ for $m=13$, $0.5702$ for $m=27$ and $0.5926$ for
$m=51$ — in each case $99.96\%$ of the ceiling (39b) allows, so the criterion is
exhausted to within $0.05\%$ of its own limit. Cost grows only at the edge: at $\theta=0.52$ and $m=27$ the first admissible
prime covers, at $0.55$ the second, at $0.565$ the 22nd, and the certified
edge $0.5702$ the 114th. Forty sampled classes per index were re-verified end to end by
`slack(m,k,p)<0` against the certifying prime. Each success is a theorem for that $(m,N,M)$; none of them is a theorem
for all $N$.

**Remark (self-evicting classes).** If some $d_i$ contains a prime power
exceeding $\sqrt{2(2N+m)}$ then *every* member of the class is fatal at once by
(25), with no construction: the class prescribes its own obstruction.

**Remark (wider windows buy nothing).** Replacing $[N,2N)$ by
$[N,N^{1+\eta})$ lengthens the progression but forces $p>\sqrt{2N^{1+\eta}}$,
so the admissible moduli grow only to $N^{(1+\eta)/2}$: the range must reach
$\approx M^2$, which is (37) again. Unioning adjacent blocks recovers the
initial-segment corollary and no more.

### 15.5 What this replaces

*Against Dirichlet.* Section 8 proves a fixed class contains infinitely many
fatal $k$ by taking the cofactor prime, and names no bound. (37) is the
effective form: the first failure arrives by $4M^2$, with an exhibited witness
of failure.

*Against Linnik.* Routing the same construction through the least prime of
$u\equiv a_i\pmod{q_i}$ gives a fatal $k\ll M^{L_0}$, $L_0\le5.18$ (Xylouris),
and needs the strong form $\pi(x;q,a)\gg x/(\phi(q)\log x)$ to reach a prime
above $2d_i$. Exponent $2$ from Bertrand beats exponent $5.18$ from Linnik.

*Against GRH.* The localisation route asks for a prime in one progression
inside a dyadic interval; GRH gives that for $q\le x^{1/2-\varepsilon}$ and so
evicts $M\le N^{1/2-\varepsilon}$. (38) reaches $\sqrt N/2$ unconditionally,
with a better constant than the conditional route. Nothing a Siegel zero could
damage appears anywhere in Section 15: no class is ever asked to contain a
prime.

### 15.6 The minimal analytic input above $\sqrt N$ — INTERFACE

Above (38) the construction does not break; only its *guarantee* does. Keeping
$i$ and the factorisation $k+i=d_ipt$, eviction of the class from $[N,2N)$ by
the prime $p$ is exactly the statement

$$
a_i\,p^{-1}\bmod q_i\ \in\ J_p:=\Bigl[\tfrac{N+i}{d_ip},\tfrac{2N+i}{d_ip}\Bigr)
\pmod{q_i},\qquad |J_p|=\frac N{d_ip}.
\tag{40}
$$

So the only missing ingredient is:

> **Hypothesis IP$(q,a,P,J)$.** Some prime $p\in(P,2P]$ has
> $a\,p^{-1}\bmod q$ inside the arc $J$.

(38) is the trivial case $|J|\ge q$. Below it, completion turns IP into an
exponential sum: with $H\asymp q$,

$$
\#\{p\sim P:\,ap^{-1}\in J\}
=\frac{|J|}q\,\pi(P)+O\Bigl(\frac{\pi(P)}H+\sum_{0<|h|\le H}\frac1h
\Bigl|\sum_{p\sim P}e_q\bigl(h\,a\,p^{-1}\bigr)\Bigr|\Bigr).
\tag{41}
$$

The main term is independent of $P$: since $|J|=N/(d_iP)$ and $d_iq_i=M$,

$$
\frac{|J|}{q_i}\pi(P)\asymp\frac N{M\log P}.
\tag{42}
$$

**Consequence (the exchange rate).** Suppose the inverse-prime sums admit
$\bigl|\sum_{p\sim P}e_q(hp^{-1})\bigr|\ll\sqrt q\,q^{\varepsilon}$ uniformly
in $h\not\equiv0$ — Weil strength for a complete Kloosterman sum. Then (41) and
(42) evict every class with $N/M\gg M^{1/2+\varepsilon}$, that is

$$
\boxed{M\le N^{2/3-\varepsilon}.}
\tag{43}
$$

More generally a saving $q^{-\delta}$ against the trivial bound $\pi(P)$ evicts
$M\le N^{1/(2(1-\delta))-\varepsilon}$, so $\delta\to\tfrac12$ is what it takes
to reach every class that meets the block at all. The parameter $P$ is free in
$\bigl[\max(\sqrt{4N+2m},m),\,N/d_i\bigr]$ and may be chosen to satisfy whatever
length-versus-modulus hypothesis a given bound requires; it cancels from (42).

This is a strictly weaker demand than the one Section 14's earlier scope note
reached for. It does not ask for a prime in a prescribed progression, nor in a
short interval; it asks for equidistribution of $p^{-1}\bmod q$ over primes of
a dyadic range, in arcs of relative length $N/(PM)$.

### 15.7 The measured barrier

`data/elementary_class_eviction.json` scans primes past (38) at $N=10^{12}$,
$m=27$, 100 uniform classes per exponent, budget 2,000 primes. Writing
$M=N^{\theta}$, the model that treats $p^{-1}\bmod q$ as uniform predicts
$MP_0/(LN)$ primes before a hit. Measured means against that prediction:

| $\theta$ | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|
| predicted | 2.26 | 9.01 | 35.9 | 142.9 | 568.7 |
| measured | 2.11 | 8.81 | 38.2 | 157.5 | 470.9 |

(The last column is censored by the budget, which biases it low.) Every class
tested at $\theta\le0.75$ was evicted with an exact certificate. The residues
$p^{-1}\bmod q$ behave here exactly like independent uniform draws, so the
obstruction above $\sqrt N$ is the missing equidistribution theorem and nothing
arithmetic: there is no conspiracy to find, only a bound to prove.

### 15.8 Exact evidence

`data/elementary_class_eviction.json` checks, without sampling inside each
declared rectangle: (i) on $1\le m\le12$, $1\le k\le300$, every (25) violation
located by real factorisation has negative Legendre slack at that prime, and
none of the 58 witnesses carries one; (ii) the master theorem (36) at its exact
minimal length on all 1,152 progressions with $1\le m\le12$, step $\le24$ and
starts $1,10^4,10^6,10^9$ — and at one shorter length the construction fails in
all 1,152 cases, so its hypothesis is exactly binding; (iii) Corollary 15.4(a)
on all 65,520 classes with $1\le m\le16$, $M\le90$, worst case using $0.672$ of
the bound (37); (iv) Corollary 15.4(b) on all 47,040 classes with
$1\le m\le20$, $M\le48$ at $N=10^6$ and $N=10^{12}$; (v) Corollary 15.4(e)
bisected for $m\in\{5,13,27,51\}$ at $N=10^{12}$, each certified exponent
re-verified on 40 sampled classes end to end; (vi) the sharpness and barrier
measurements above. Every individual eviction anywhere in the artifact is
re-verified by `slack(m,k,p)<0`, independently of the construction that found
it.


## 16. The survivor route, and where the single-position argument is sharp

Section 15 evicts a class by *constructing* a fatal member. There is a second,
independent route: count the members that could survive and show there are
fewer than the class has. It needs no prime localisation either, and the
literature already covers the whole modulus range — what it lacks is one
effective constant.

### 16.1 The survivor reduction — PROVED

For a position $i\in I_m$ write the *survivors* of the class $a\bmod M$ as

$$
S(N,M,a)=\#\bigl\{\,n\in[N,2N):n\equiv a\ (M),\ \text{$n$ satisfies (25)}\,\bigr\},
\tag{44}
$$

where "$n$ satisfies (25)" means $p^{2v_p(n)}\le2n$ for every prime $p>m$
dividing $n$. If

$$
S(N,M,a)\ <\ \#\{n\in[N,2N):n\equiv a\ (M)\}
\tag{45}
$$

then some member violates (25) and the class $a-i\bmod M$ is evicted from
$[N,2N)$. Since the survivors are exactly the $\sqrt{2n}$-powersmooth-above-$m$
integers of the class, (45) is a statement about smooth numbers in one
arithmetic progression, for an individual modulus, with no averaging.

The margin is generous. The Dickman density of integers $n\le x$ with
$P^+(n)\le\sqrt{x}$ is $\rho(2)=1-\log2=0.30685\ldots$, so the expected
survivor count is about $0.307$ of the class, and (45) needs only a bound below
$1$ — a factor $1/\rho(2)=3.26$ of room, or $1.63$ if one insists on the worst
documented deviation mechanism (concentration into an index-two subgroup, which
doubles the count). `LOCALIZATION.md` audits which published bound on
$\Psi(x,y;q,a)$ covers this corner and why none of them closes it: Balog–
Pomerance (1992) and Shiu (1980) both apply to our exact parameters for every
individual modulus, and both carry an unspecified absolute constant.
**This gap is one effective constant, not a missing method.**

### 16.2 The single-position route is sharp at $\sqrt{2N}$ — PROVED

**Proposition.** Let every prime factor of $M$ be at most $m$, and let
$M^2\ge N$. Then every $n\in[N,2N)$ with $M\mid n$ satisfies (25).

*Proof.* Let $p>m$ and $p^v=p^{v_p(n)}$. Then $p\nmid M$, so $p^v$ divides
$j=n/M$, whence $p^v\le j\le(2N-1)/M<2M$, the last step by $M^2\ge N$.
Therefore $p^{2v}\le j^2<j\cdot2M=2n$. $\square$

The threshold is exactly the one Corollary 15.4(b) reaches: eviction is proved
for $M\le N/P_0\approx\sqrt N/2$ and obstructed from $M\ge\sqrt N$, so the
single-position route is sharp to within the factor $2$.

So for $m$-smooth moduli beyond $\sqrt N$ the position $i$ with $M\mid k+i$
is **permanently** non-fatal: no choice of $k$ in the class, no prime and no
prime power can evict through it. Corollary 15.4(b)'s threshold
$M\lesssim\sqrt N/2$ therefore cannot be pushed past $\sqrt N$ by any
single-position argument — the barrier at $\sqrt N$ is a property of the
problem, not an artefact of the construction. The census realises it at
$M=4620=2^2\cdot3\cdot5\cdot7\cdot11$ with $N=10^7$: $M^2=2.13\cdot10^7\ge N$,
and not one of the $2{,}165$ members of the class $0\bmod4620$ violates (25).

The escape is that eviction needs only *one* position, and the proposition pins
only the position divisible by $M$.

**The proposition is the divisibility-side twin of (33).** Read constructively,
16.2 says how to *manufacture* a term satisfying (25): make it divisible by an
$m$-smooth $M$ with $M^2\ge N$. Doing that at every position at once needs
$m$-smooth $M_i\mid k+i$ with $M_i\gtrsim\sqrt k$ for all $i\in I_m$, so CRT
forces $k\gtrsim\prod_iM_i\ge k^{L/2}$, impossible for $L\ge3$. That is the
same verdict (33) reaches by counting mass, obtained here by counting
divisibility: for $m\ge5$ no congruence prescription can supply the whole
window, whichever side one counts from.

### 16.3 The multi-position interface — the weakest input so far

Call an $n$-class *fatal-free* when no member violates (25). A $k$-class $c$
survives every position exactly when $c+i$ is fatal-free for all $i\in I_m$;
since $I_m$ is an interval of $L=\lceil m/2\rceil$ consecutive integers, the
surviving $k$-classes are the starts of cyclic runs of $L$ fatal-free residues.
Hence

$$
\boxed{\text{no }L\text{ consecutive residues}\bmod M\text{ are all fatal-free}
\ \Longrightarrow\ \text{every class}\bmod M\text{ is evicted from }[N,2N).}
\tag{46}
$$

A class containing a prime $n\in(m,2N)$ is never fatal-free — take $p=n$, so
$p^2>2n$ — so (46) follows from: *every interval of $L$ consecutive residues
$\bmod M$ contains one whose class holds a prime of $[N,2N)$.* That is strictly
weaker than everything Section 15.6 asks for:

* it never asks a prescribed class to hold a prime, only one class in each
  window of $L$;
* the union of $L$ consecutive classes is an *interval of residues*, whose
  expected prime count is $L\cdot N/(M\log N)$ — positive up to
  $M\approx LN/\log N$, that is, up to essentially every modulus meeting the
  block more than once;
* it weakens as $L=\lceil m/2\rceil$ grows, and Erdős #389 is open precisely in
  the large-$m$ direction.

The cost is that the residues $c+i$ are consecutive, so the exponential sums
are $\sum_{p\sim N}e(hp/M)$ with $|h|\le M/L$ — Weyl sums over primes with
rational argument, not Kloosterman sums. Vinogradov's $x^{4/5}$ term, and under
GRH the $x^{1/2}$ term, cap that route at $M\lesssim\sqrt N$ again, so (46) is
not yet a theorem. It is the cheapest missing statement in this file, and the
census below shows it holding with room to spare.

### 16.4 Exact evidence

`data/survivor_census.json` builds the exact fatality mask of
$[10^7,2\cdot10^7)$ — $n$ is fatal iff some prime power $p^v\mid n$ with $p>60$
has $p^{2v}>2n$, which is conservative for every $m\le60$ — and censuses every class of 36 moduli in three shapes
(prime, $47$-smooth, primorial) at target exponents
$\theta=\log M/\log N$ from $0.40$ to $0.99$. Every count is exhaustive over the
block.

* Survivor density $0.31223$ against $\rho(2)=0.30685$: ratio $1.0175$.
* Single position ($m=1$): classes evicted up to $\theta=0.85$, with the first
  failures at $\theta=0.523$ — exactly the smooth-modulus classes of 16.2.
* Multi-position: every $k$-class of every tested modulus is evicted up to
  $\theta=0.858$ for $m=3$, $0.900$ for $m=5$, $0.950$ for $m=13$, and
  $0.990$ — the largest tested — for both $m=27$ and $m=51$, with no failure at
  any tested exponent for those two.

At $m=27$, the index whose least witness is unknown, the whole interval between
the proved threshold $\sqrt N/2$ and the trivial ceiling $N$ is therefore
*certified finite* at $N=10^7$: no congruence class of any tested modulus
survives its block. The gap Section 15 leaves is quantitative, and (46) is what
would close it for every $N$.


## 17. The compensation-good integers have positive density — PROVED

Section 13 reduced Erdős #389 to Conjecture R and left R(1) — positive lower
density for a single term — open on two grounds: the first-order union bound
over failure events tends to $\log2+0.322521>1$, and an exact block
classification was read as showing that *"the overlap that would rescue the
union bound does not exist"*. **That second reading was wrong.** This section
replaces it with a proof of R(1).

What the classification actually refuted was the *independence prediction* for
the level-failure-only class. The union bound needs far less than independence:
it overshoots $1$ by $0.015668$, so any lower bound on the overlap above that
number closes R(1). The same census already reported an overlap of $0.053415$
at $10^6$ and $0.063050$ at $10^7$ in the $m=27$ census, and $0.07481$ at
$10^6$ for $m=1$ — three to five times what is required. The
overlap was never too small; it was only never bounded *from below*.

A lower bound looked out of reach because of the direction of counting. Section
13 partitioned by the largest prime factor, which asks for smooth numbers in
progressions to moduli near $\sqrt X$ — Bombieri–Vinogradov strength, correctly
declined there. Counted from the prime side the same overlap is cheap: fix the
level-failing prime $q$ and the small cofactor $t$ and let the large prime $p$
run. The level test becomes a congruence on $p$ modulo $q^r$, a *single* modulus
far below the length of the $p$-range, so Vinogradov's bound for exponential
sums over primes applies to each modulus separately. No average over moduli is
taken, so neither Bombieri–Vinogradov nor GRH enters.

### 17.1 The identity that replaces the union bound

Fix $m\ge1$ and split the failures of (25) at $w$ by Section 13's dichotomy:

* $A$: some prime $p>m$ has $p\mid w$ and $p^2>2w$ (*short cofactor*);
* $B$: some prime $q>m$ has $q^e\|w$, $q^2\le2w$ and $C_q(w/q^e)<e$ (*level*).

$w$ is compensation-good exactly when $w\notin A\cup B$, so on any block

$$
\#\mathrm{good}\,[x,2x)=x-\#A-\#B+\#(A\cap B).
\tag{47}
$$

This is exact. A union bound discards the last term; (47) says that term is
worth exactly as much as the bound overshoots.

### 17.2 Short failures are large prime factors — PROVED

**Lemma 17.2.** $w\in A$ iff some prime $p>m$ divides $w$ with $p^2>2w$; any
such $p$ has $v_p(w)=1$ and $C_p(w)=0$; and $\#A\le(\log2+o(1))x$ on $[x,2x)$.

*Proof.* If $p^2>2w$ and $p^e\|w$ with $e\ge2$ then $p^2\le w<p^2/2$, absurd, so
$e=1$ and the cofactor is $c=w/p<p/2$. Then no level is usable, $C_p(w)=0<1$,
and $p$ fails; conversely a short failure has $p^2>2w$ by definition. For the
count, $p^2>2w\ge2x$ forces $p>\sqrt{2x}$ and $p\mid w<2x$ forces $p<2x$, and
$p$ has at most $x/p+1$ multiples in the block, so
$\#A\le\sum_{\sqrt{2x}<p<2x}(x/p+1)$. Mertens' second theorem gives
$\sum_{\sqrt{2x}<p<2x}1/p=\log\frac{\log2x}{\log\sqrt{2x}}+o(1)=\log2+o(1)$, and
$\pi(2x)=o(x)$. $\square$

Only Mertens is used — no prime number theorem, no primality input.

### 17.3 The level union bound — PROVED

**Lemma 17.3.** $\#B\le\bigl(\Lambda_B+\varepsilon_2(x)+o(1)\bigr)x$ where

$$
\Lambda_B=\sum_{r\ge1}2^{-r}\log\frac{r+2}{r+1}=0.322521\ldots,
\tag{48}
$$

and $\varepsilon_2(x)\to0$ collects the exponents $e\ge2$.

*Proof.* For $e=1$ a prime $q$ with $q\|w$ has
$R_q(w)=\#\{j\ge1:q^j\le2w/q\}=\lfloor\log_q 2w\rfloor-1$ usable levels, and
$C_q(w/q)=0$ demands all $R_q$ of them in the lower half. The admissible
residues form the set $S_R(q)$ of (50) below, of period $q^R$, so counting
$w\in[x,2x)$ with $q\mid w$ and $w/q\in S_R(q)$ gives
$\frac xq2^{-R}\bigl(1+O(q^R q/x)\bigr)$; since $q^{R+1}\le2w$ the relative
error is $O(1/q)$. Grouping by $R=r$, that is
$(2x)^{1/(r+2)}<q\le(2x)^{1/(r+1)}$, Mertens gives $\sum1/q=\log\frac{r+2}{r+1}$
over the band and the sum telescopes to (48); the tail past $r$ terms is at most
$2^{-r}/(r+2)$ because $\log\frac{r+2}{r+1}\le\frac1{r+1}$. For $e\ge2$ the
density of $q^e\|w$ is at most $q^{-e}$ and the failure needs fewer than $e$ of
the $R=\lfloor\log_q2w\rfloor-e$ levels in the upper half, costing at most
$\sum_{i<e}\binom Ri2^{-R}$; every fixed $q$ gains levels as $x$ grows while the
large-$q$ end is killed by $q^{-e}$, so the total is $\varepsilon_2(x)\to0$.
$\square$

### 17.4 The overlap families — PROVED

For $r\ge1$ and $\eta\in(0,\tfrac1{10})$ let $D_r(x)$ be the set of
$w\in[x,2x)$ admitting a factorisation $w=pqt$ with $p,q$ prime, $q\nmid pt$,

$$
p^2>(2w)^{1+2\eta},\qquad q^{r+1}\le2w<q^{r+2},\qquad
q^{r(1+\eta)}\le\frac{2x}{qt},\qquad C_q(pt)=0 .
\tag{49}
$$

**Lemma 17.4.** $D_1(x)$ and $D_2(x)$ are contained in $A\cap B$, are disjoint,
and no element of either has two such factorisations.

*Proof.* $p^2>2w$ puts $w$ in $A$ by Lemma 17.2, and $q^{r+1}\le2w$ with
$r\ge1$ gives $q^2\le2w$, so the failure of $q$ — which is a failure, since
$v_q(w)=1$ and $C_q(w/q)=C_q(pt)=0$ — is a level failure, putting $w$ in $B$.
For uniqueness and disjointness, a band-$r$ prime satisfies
$q>(2w)^{1/(r+2)}\ge(2w)^{1/4}$ for $r\le2$, so two of them together with
$p>(2w)^{1/2}$ would give $w\ge pq_1q_2>(2w)^{1/2+1/4+1/4}=2w$. $\square$

The exponent $1/2+2/(r+2)$ is at least $1$ — so that $(2w)^{1/2+2/(r+2)}\ge2w>w$
and the contradiction closes — exactly for $r\le2$, which is why the
certificate stops at band $2$: bands $r\ge3$ are available but would need a
Bonferroni correction for multiplicity.

### 17.5 The lower-half residue count — PROVED

$$
S_r(q)=\bigl\{\,v\bmod q^r:\ 2\,(v\bmod q^j)\le q^j\ \text{for }j=1,\dots,r\,\bigr\}.
\tag{50}
$$

**Lemma 17.5.** $|S_r(q)|\ge q^r2^{-r}\,(1-1/q)$.

*Proof.* $|S_1(q)|=\lfloor q/2\rfloor+1>q/2$. The conditions with $j<r$ cut out
a $q^{r-1}$-periodic set meeting each period in $|S_{r-1}(q)|$ residues, and the
condition at $j=r$ restricts $v$ to $[0,\lfloor q^r/2\rfloor]$, which contains at
least $\lfloor q/2\rfloor$ whole periods. Hence
$|S_r(q)|\ge\lfloor q/2\rfloor\,|S_{r-1}(q)|\ge\frac{q-1}2\,q^{r-1}2^{-(r-1)}$.
$\square$

### 17.6 Counting the overlap families — PROVED using Vinogradov

**Lemma 17.6.** For $r\in\{1,2\}$ and fixed $\eta>0$,
$\#D_r(x)\ \ge\ \bigl(\lambda_r(\eta)-o(1)\bigr)x$, where

$$
\lambda_r(\eta)=2^{-r}\int_{1/(r+2)}^{1/(r+1)}\frac1\beta\,
\log^{+}\frac{1-\beta}{\max\bigl(\tfrac12+\eta,\ r\beta(1+\eta)\bigr)}\,d\beta .
\tag{51}
$$

*Proof.* Write $w=pqt$ and sum over $q$ and $t$ obeying (49), counting primes
$p\in[x/(qt),2x/(qt))$ with $pt\bmod q^r\in S_r(q)$. Expand the indicator of
$S_r(q)$ — a product of $r$ interval conditions — in additive characters mod
$q^r$: the constant coefficient is $|S_r(q)|/q^r\ge2^{-r}(1-1/q)$ by Lemma 17.5,
and the remaining coefficients have $\ell^1$ norm $O(\log^rq)$. For $h\not\equiv0$
the phase $ht/q^r$ has reduced denominator $q'\in[q,q^r]$ because $q\nmid t$, so
Vinogradov's bound gives
$\sum_{p\le y}e(hpt/q^r)\ll\bigl(y\,q'^{-1/2}+y^{4/5}+\sqrt{yq'}\bigr)\log^4y$
with $y=2x/(qt)$. Summing the three error terms over all admissible $(q,t)$:

* $\sum_{q,t}y\,q^{-1/2}\ll x\log x\sum_{q>x^{1/4}}q^{-3/2}\ll x^{7/8}\log x$;
* $\sum_{q,t}y^{4/5}\ll x^{4/5}\sum_q q^{-4/5}(x^{1/2}/q)^{1/5}\ll x^{9/10}\log\log x$;
* $\sum_{q,t}\sqrt{yq^r}$, the binding term, is for $r=1$ at most
  $\sqrt x\sum_{q\le x^{1/2}}2\sqrt{x^{1/2-\eta}/q}\ll x^{1-\eta/2}/\log x$,
  and for $r=2$, using $t\le x\,q^{-3-2\eta}$ from (49), at most
  $x\sum_{q>x^{1/4}}q^{-1-\eta}\ll x^{1-\eta/4}/\eta$.

All three are $o(x)$ for fixed $\eta>0$; the $\eta$-cuts in (49) exist precisely
to make the third one so. The main term is
$\sum_{q,t}2^{-r}\bigl(\pi(2x/(qt))-\pi(x/(qt))\bigr)$, and with
$q=(2x)^{\beta}$, $t=(2x)^{\tau}$, $p=(2x)^{\alpha}$ the prime number theorem and
Mertens turn it into $2^{-r}x\iint d\alpha\,d\beta/(\alpha\beta)$ over
$\alpha\in[\max(\frac12+\eta,r\beta(1+\eta)),\,1-\beta]$ and $\beta$ in the
band, which is (51). $\square$

### 17.7 The theorem

**Theorem 17.7 (Conjecture R(1)).** For every $m\ge1$,

$$
\liminf_{x\to\infty}\frac{\#\{w\in[x,2x):w\ \text{is }m\text{-compensation-good}\}}{x}
\ \ge\ 1-\log2-\Lambda_B+\lambda_1+\lambda_2\ \ge\ 0.0293 .
\tag{52}
$$

*Proof.* Combine (47) with Lemma 17.2, Lemma 17.3 and
$\#(A\cap B)\ge\#D_1+\#D_2$ from Lemmas 17.4 and 17.6, then let $\eta\to0$ after
$x\to\infty$. $\square$

With $\eta=0.005$ the certified numbers are $\lambda_1\ge0.030515$,
$\lambda_2\ge0.014505$ and $\varepsilon_2\le6\cdot10^{-6}$ at $x=10^{40}$,
against a first-order deficit $\log2+\Lambda_B-1=0.015668$: the overlap supplies
$0.045021$, a factor $2.87$ more than needed. The bound is uniform in
$m\le x^{1/4}$, since every prime it names exceeds $(2x)^{1/4}$; in particular
it holds for $m=1$, which by the monotonicity of Section 10 implies every $m$.

### 17.8 What this settles and what it does not

R(1) is now a theorem, so the first rung of Conjecture R is unconditional and
the union-bound obstruction of Section 13 is dissolved rather than circumvented.

It does not give a good constant: $0.0293$ is far below the measured $0.1237$,
because all three inputs are one-sided. Adding a second Bonferroni term to
$\Lambda_B$ itself — two band primes both failing, worth about $0.015$ — would
raise it, at proof cost and with no new consequence.

**It does not give R(2), and the reason is exact.** For $w$ and $w+1$ the four
events $A_0,B_0,A_1,B_1$ make (47) a four-fold inclusion–exclusion, and the
Bonferroni truncation that would replace it is $S_1-S_2+S_3\ge P(\bigcup)$.
Writing $u=\log2+\Lambda_B$ for the single-term first-order total and
$\pi=d(A\cap B)$, and treating the two coordinates as independent — which the
block census below confirms to within $0.002$ — every term collapses and

$$
S_1-S_2+S_3\;=\;2u-\bigl(u^2+2\pi\bigr)+2\pi u
\;=\;1+(u-1)\bigl(2\pi-(u-1)\bigr).
\tag{53}
$$

So the depth-3 route closes R(2) **iff $\pi<(u-1)/2=0.007834$**. The overlap
certified in 17.6 is $\pi\ge0.045$, and the measured value is $0.087$: too large
by a factor between $5.7$ and $11$. *The same overlap that proves R(1) is what
defeats R(2)*, because at one coordinate it enters as $+\pi$ against a deficit
$u-1$, and at two coordinates as $+2\pi(u-1)$ against a budget of only
$(u-1)^2$. Below the asymptotic regime (53) is negative — the exact block totals
are $0.99234$ at $10^6$ and $0.99435$ at $10^7$, both under $1$ — purely because
$u<1$ there; the finite-scale success is an artefact of $u$ not yet having
converged, and is worthless, exactly as Section 13 said of the same dip.

Even at $\pi<0.007834$ the route would still need the four cross terms of $S_2$
from below: $d(A_0\cap A_1)$ is the density of $w$ with both $w$ and $w+1$
carrying a prime factor above $\sqrt{2w}$, a binary problem of Chen type where
only upper bounds (CRT plus Mertens, giving $\log^22$) are elementary. So R(2)
fails here twice over, once numerically and once structurally, and the pivot is
recorded rather than retried.

What the section does supply is a template: **for these events the overlap is
bounded from below on the prime side, not the smooth side.** Section 13's
partition put the smooth cofactor in the free variable and needed
$\Psi(x,y;q,a)$; putting the large prime there instead needs only a
single-modulus exponential-sum bound. `LOCALIZATION.md` §5 records the same
distinction as the reason Vinogradov's theorem is PARTIAL for input III and
APPLIES here: what matters is $\log q'/\log y$, not the size of $q'$.

### 17.9 Exact evidence

`data/overlap_density_bound.json` certifies each input as a one-sided bound with
an explicit tail — the level sum truncated with tail $2^{-r}/(r+2)$, the
exponent tail with primes past $10^6$ discarded by $2/10^6$, and (51) as a
right-endpoint Riemann sum over a provably decreasing integrand (monotonicity
asserted at every step, not assumed). It then verifies on
$[10^6,10^6+2\cdot10^5)$ and $[10^7,10^7+2\cdot10^5)$, exhaustively, that (47)
holds exactly, that every member of $D_1\cup D_2$ lies in $A\cap B$, that no
member is counted twice, and that the family density ($0.04871$, $0.049755$)
stays below the overlap measured in the same run at $m=1$ ($0.07481$,
$0.08735$). The same artifact records the four-event accounting
behind (53): $S_1,\dots,S_4$ exactly, the verified identity
$S_1-S_2+S_3-S_4=P(\bigcup)$, and the asymptotic depth-3 limit $1.001165$,
asserted to exceed $1$ so that the recorded obstruction cannot go stale.

