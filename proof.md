# Position-weighted pressure improvement to 67.3205978423%

## Statement

Let $N(T,2T)$ count nontrivial zeros of the Riemann zeta function with
multiplicity, and let $N_0^s(T,2T)$ count simple zeros on the critical line.
Subject to the same imported analytic interface as
`trmdy/zeta-simple-zeros-673137`, the certificate in this bundle gives

$$
\liminf_{T\to\infty}\frac{N_0^s(T,2T)}{N(T,2T)}
\ge 0.6732059784228011963\ldots
>\frac{6732059784}{10^{10}}.
$$

This is **67.3205978423%**.  The predecessor's tightened value was
67.3137960544%.

## 1. Imported interface

For the already re-optimized seven-term window $v$, write

$$
H(v)\ge H_{\mathrm{cert}}
=\frac{672457041414544284}{10^{18}},
$$

and let $M$ be the Gram matrix of the retained simple-zero atoms.  The
imported stability argument gives

$$
S\ge H_{\mathrm{cert}}N+\mathrm{tr}\,\Psi(M)-o(N),
\qquad S=N_0^s(T,2T),\ N=N(T,2T),
$$

where

$$
\Psi(t)=\begin{cases}(t-1)^2,&0\le t\le2,\\2t-3,&t\ge2.\end{cases}
$$

The window and its kernel are unchanged from the 67.313763% successor.

## 2. New seven-point inequality

For six nonnegative consecutive gaps $g_1,\ldots,g_6$, put
$y_0=0$ and $y_j=g_1+\cdots+g_j$.  Let
$w(x)=k_v(x)^2$.  The new local functional is

$$
F(g)=\sum_{r=1}^6 b_r g_r
 +\sum_{0\le i<j\le6} a_{ij}w(y_j-y_i).
$$

Unlike the predecessor, the linear pressure is **position dependent**.
The exact coefficients are in `candidate.json`.  They satisfy

$$
\sum_{r=1}^6b_r=\frac3{1150}
$$

and, for every pair span $s=1,\ldots,6$,

$$
\sum_{i=0}^{6-s}a_{i,i+s}=2.
$$

The interval verifier proves

$$
F(g_1,\ldots,g_6)\ge
\varepsilon:=\frac{51063}{10^7}=0.0051063
\qquad(g_r\ge0).
$$

The observed floating-point adversarial minimum is
$0.005106397636960345\ldots$, but that observation is not used in the
proof.

### Why position-dependent pressure improves the global bound

Sum the local inequality over all $m-6$ consecutive seven-point windows in
an $m$-point block.  A fixed block gap appears at most once in each of the
six local positions, so its total linear coefficient is at most
$B:=\sum b_r=3/1150$.  A pair spanning $s$ gaps appears at most once at
each admissible position, and its total weight is at most the span capacity
2.  Therefore, with

$$
E_B=2\sum_{i<j}w(y_j-y_i),
$$

one obtains

$$
E_B+B\,\mathrm{span}(B)
\ge A:=\varepsilon(m-6).
$$

Redistributing the six pressure coefficients can raise the local minimum
without increasing $B$, so it is invisible to the global penalty.  This is
the new lever.

## 3. Sharp finite-dimensional Gram profile

For an $m\times m$ positive-semidefinite Gram block $G$ with
$G_{ii}\le1$, define

$$
E=2\sum_{i<j}|G_{ij}|^2,
\qquad \Delta=\mathrm{tr}\,\Psi(G).
$$

Then

$$
\Delta\ge h_m(E),
$$

where

$$
h_m(E)=
\begin{cases}
E,&0\le E\le \dfrac m{m-1},\\
\dfrac Em+2\sqrt{\dfrac{m-1}{m}E}-1,
 &E\ge\dfrac m{m-1}.
\end{cases}
$$

This profile is sharp for every fixed $m$.

### Proof

Let $\lambda_i$ be the eigenvalues, $d_i=G_{ii}$, and
$t=\mathrm{tr}\,G\le m$.  Cauchy--Schwarz gives

$$
E=\sum_i\lambda_i^2-\sum_i d_i^2
\le \sum_i\lambda_i^2-\frac{t^2}{m}=:V.
$$

If all eigenvalues are at most two, then
$\Delta=\sum_i(\lambda_i-1)^2\ge V$.

Otherwise write the $k\ge1$ eigenvalues above two as $2+z_i$, put
$Z=\sum z_i$, and write $s=m-t\ge0$.  Cauchy--Schwarz on the remaining
eigenvalues yields

$$
\Delta\ge k+2Z+\frac{(s+k+Z)^2}{m-k}
\ge1+2Z+\frac{(1+Z)^2}{m-1}.
$$

Also $V\le\Delta+Z^2$.  Eliminating $Z$ gives the displayed formula for
$h_m$.  Equality is attained by a constant-diagonal matrix with spectrum

$$
1+x,\quad 1-\frac{x}{m-1},\ldots,1-\frac{x}{m-1}.
$$

## 4. Block deduction

Set $R=h_m(A)$ and $\eta=R/A$.  Concavity and monotonicity of $h_m$
give

$$
\mathrm{tr}\,\Psi(G_B)
 +\eta B\,\mathrm{span}(B)\ge R-o(1).
$$

The usual shifted-block pinching and averaging then imply

$$
\mathrm{tr}\,\Psi(M)
\ge \frac{R}{m}S
 -\eta B\frac{m-1}{m}N-o(N).
$$

Substitution into the imported interface gives

$$
\frac SN\ge
\frac{mH_{\mathrm{cert}}-\eta B(m-1)}{m-R}.
$$

Scanning integral block lengths selects $m=250$.  Thus

$$
A=\frac{51063}{10^7}(250-6)=1.2459372,
$$

$$
R=h_{250}(A)
=1.2329455538586055771715901086384\ldots,
$$

and

$$
\frac SN\ge
0.6732059784228011963162317558716\ldots.
$$

The enclosed arithmetic is recorded in
`certificates/final-arithmetic.txt`.

## 5. Finite certificate

The bundled verifier uses a deliberately small alternate trust base:

1. `mpmath.iv` interval arithmetic to enclose $w$, $w'$, and $w''$;
2. binary64 lower tables widened with `nextafter`;
3. exhaustive six-dimensional subdivision in C++;
4. a rigorously outward-rounded `long double` interval LDL decomposition to
   prove convexity before applying tangent lower bounds.

The main run records:

- target: $51063/10^7$;
- grid: $1/4000$;
- coarse cells: 56,530;
- initial boxes: 144;
- visited nodes: 2,037,322;
- tangent-pruned boxes: 554,949;
- maximum depth: 73;
- no unresolved terminal cell.

As a hardening check, the same tables prove $F\ge0.00507$ using pure
interval subdivision with tangent pruning disabled: 363,259,330 nodes.

## 6. Research status

The finite inequality has been reproduced in this alternate interval trust
base, and the exact global accounting has been checked independently.  The
result still inherits the analytic inputs of the Anthropic/ainta/successor
research drafts.  Before presenting it as a publication-grade theorem, the
new vector-pressure certificate should be rerun through the repository's
native Arb (`python-flint`) verifier and the finite-dimensional deduction
should receive independent human or Lean review.
