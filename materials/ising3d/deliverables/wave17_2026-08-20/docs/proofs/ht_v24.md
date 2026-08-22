# Exact simple-cubic HT coefficient through v^24

## Scope

[COMPUTATION] This is a finite exact coefficient calculation, not a solution of
the three-dimensional Ising model.

## Why the previous route used 242 GB

[LEMMA] The previous broken-bond spin transfer stores a dense array indexed by
all `2^25` spins on a 5x5 section and all 301 broken-bond degrees.  Three
simultaneous int64 arrays therefore require
`3*2^25*301*8 = 242397216768` bytes.  Five
31-bit CRT primes are enough; the modulus count multiplies time, not this peak.
There is no boundary-connectivity state in that implementation.  The real wall
is the cross-section exponent, with polynomial degree as a factor.

## Lower-memory exact extraction

[LEMMA] A graph whose 5x5x5 bounding box first contributes at degree 24 crosses
each of the 12 coordinate cuts positively and evenly.  Its total edge count is
24, so every cut is crossed exactly twice.  The frontier algorithm therefore
tracks only parity on unprocessed edges and ternary cut counters 0,1,2; any
third crossing is discarded.  The frontier has a measured peak of `1727071`
states.  The isolated production run used `352305152` bytes maximum RSS
(`360792736`-byte peak footprint), versus the old 242-GB projection; it read
and wrote zero disk blocks.
The complete coefficient run, including the unchanged old engine on the other
101 canonical boxes, had a process peak RSS of `7305396224` bytes, still below
the workstation's 51.5-GB physical memory.

[LEMMA] Introduce one square-free variable for each cut after restricting its
count to zero or two.  The coefficient containing all 12 variables in the
formal logarithm is exactly the 5x5x5 finite-lattice weight at v^24: any proper
subbox misses a cut, while every full-box degree-24 connected contribution
crosses all cuts twice.

## CRT certificate and coefficient

[THEOREM] The exact 5x5x5 finite-lattice weight is
`8163299968` at v^24.

For each active cut there are at most C(25,2)=300 choices.  Expanding the
square-free logarithm over ordered set partitions gives the strict absolute
bound
`300^12 * sum(k=1..12, S(12,k)*(k-1)!) =
1724666954455386000000000000000000000000`.
The coprime modulus product is `45671921168693645933699105804560590380377589537`, greater than
twice this bound.  Hence centered CRT reconstruction is unique.

[COMPUTATION] Injecting only this minimal-cube weight into the unchanged exact
finite-lattice inversion gives

`[v^24] phi = 2135670379057/8`,

and interaction coefficient `266958797382`.  Every
coefficient through v^22 agrees with the canonical frozen prefix.

## Independent controls

[COMPUTATION] The same cut-frontier/log extractor returns the independently
enumerated minimal-cube weights `16` for 2x2x2,
`9188` for 3x3x3, and
`8655072` for 4x4x4.  The ordinary FLM independently
returns the same three integers.  Disk checkpoint use was zero bytes.
