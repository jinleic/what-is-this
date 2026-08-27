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
   $m$-compensation-good integers whose start also satisfies the exact
   small-prime tier. No theorem proves that such a run exists for a given $m$.
   Ceiling (26) bounds the density of the target from above but cannot
   establish nonemptiness.
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
