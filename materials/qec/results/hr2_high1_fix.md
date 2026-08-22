# HR2 HIGH-1 — exact replacement text for Theorem G(ii)'s "row space extends" passage

Drop-in text for the false sentence in the proof of Theorem G(ii)
(`reports/paper_pbb_nogo.md:131–134`, `reports/paper_pbb_nogo.tex:173–177`).
Nothing else in the proof changes; the following "Equivalently, via the sector
identity …" sentence and the closing parenthetical (J.1–J.3) stay as they are.

## 1. Text being removed

**md 131–134** (one sentence):

> For the dimension we count generators directly: the row space of $H_Q$ extends that of
> $H_P$ exactly by the $[C\;D]$-parts of first-block combinations taken modulo $S_Z$, so
> $\operatorname{rank}H_Q=\operatorname{rank}H_P+\dim\bar\Delta$, and $k=n-\operatorname{rank}H$
> gives $k_Q=k_P-\dim\bar\Delta$.

**tex 173–177** (same sentence, preceded by a redundant duplicate of the
already-stated fact — drop the duplicate too):

> The pure-$Z$ stabilizer space is $S_Z+\Delta$.  For the dimension we count generators
> directly: the row space of $H_Q$ extends that of $H_P$ exactly by the $[C\;D]$-parts of
> first-block combinations taken modulo $S_Z$, so
> $\operatorname{rank}H_Q=\operatorname{rank}H_P+\dim\bar\Delta$, and
> $k=n-\operatorname{rank}H$ gives $k_Q=k_P-\dim\bar\Delta$.

(The tex line 173 lead-in "The pure-$Z$ stabilizer space is $S_Z+\Delta$." duplicates
tex line 172 / md line 130 verbatim in substance — removing it is part of this patch.)

## 2. Replacement — markdown form

Insert after md line 130 ("…gives pure-$Z$ stabilizer space $S_Z+\Delta$."):

```markdown
For the dimension, project the row space onto the $X$-coordinate: the image is
$\operatorname{rowspace}[A\;B]$ for both codes, while the kernel is $\{0\}\times S_Z$
for $H_P$ and $\{0\}\times(S_Z+\Delta)$ for $H_Q$ — the pure-$Z$ stabilizer space,
by the first part of (ii). Hence
$\operatorname{rank}H_Q=\operatorname{rank}H_P+\dim\bar\Delta$, and
$k=n-\operatorname{rank}H$ gives $k_Q=k_P-\dim\bar\Delta$.
```

## 3. Replacement — LaTeX form

Replaces tex lines 173–177 (keep the trailing "  Equivalently, via the" so line 178 joins unchanged):

```latex
For the dimension, project the row space onto the $X$-coordinate: the image is
$\operatorname{rowspace}[A\;B]$ for both codes, while the kernel is $\{0\}\times S_Z$
for $H_P$ and $\{0\}\times(S_Z+\Delta)$ for $H_Q$ --- the pure-$Z$ stabilizer space,
by the first part of~(ii).  Hence
$\operatorname{rank}H_Q=\operatorname{rank}H_P+\dim\bar\Delta$, and
$k=n-\operatorname{rank}H$ gives $k_Q=k_P-\dim\bar\Delta$.  Equivalently, via the
```

(`\operatorname{rowspace}` is already used in the manuscript: tex lines 127, 211.)

## 4. Why every clause is true

Let $\pi$ be the projection $(x|z)\mapsto x$ onto the $X$-coordinate.

- **Image.** $\operatorname{rowspace}(H_Q)=\{(\lambda[A\;B],\ \lambda[C\;D]+\mu H_Z)\}$
  and $\operatorname{rowspace}(H_P)=\{(\lambda[A\;B],\ \mu H_Z)\}$; in both cases
  $\pi$ maps onto $\operatorname{rowspace}[A\;B]$ (every $u=\lambda[A\;B]$ is the
  image of its own row, and nothing else appears in the first coordinate).
- **Kernel.** $\ker\pi$ meets each row space exactly in its pure-$Z$ elements:
  $\{0\}\times S_Z$ for $H_P$, and $\{0\}\times(S_Z+\Delta)$ for $H_Q$ — which is the
  first part of (ii), proved two sentences earlier; no new fact is invoked.
- **Rank identity.** $\operatorname{rank}=\dim\pi(\text{rowspace})+\dim\ker\pi$ gives
  $\operatorname{rank}H_Q=\operatorname{rank}[A\;B]+\dim(S_Z+\Delta)$ and
  $\operatorname{rank}H_P=\operatorname{rank}[A\;B]+\dim S_Z$; subtracting,
  $\operatorname{rank}H_Q-\operatorname{rank}H_P=\dim(S_Z+\Delta)-\dim S_Z=\dim\bar\Delta$
  with $\bar\Delta=(\Delta+S_Z)/S_Z$. Then $k=n-\operatorname{rank}H$ yields
  $k_Q=k_P-\dim\bar\Delta$.
- No containment of row spaces is claimed or used (the false part of the old text).

**Machine re-verification (this session, GF(2) numpy):** on `phase2_58` (ℓ=m=6),
$\operatorname{rank}H_Q=68=\operatorname{rank}[A\;B]+\dim(S_Z+\Delta)=32+36$,
$\operatorname{rank}H_P=64$, $\dim\bar\Delta=4$, $k_Q=4=k_P-\dim\bar\Delta$; on
`12_6_0193` (ℓ=12, m=6), $\operatorname{rank}H_Q=\operatorname{rank}H_P=132$,
$\dim\bar\Delta=0$, $k_Q=12=k_P$. Meanwhile $\operatorname{rank}[H_Q;H_P]$ is
$92$ and $178$ respectively — confirming the old "extends" clause was false while
the identity above holds.
