# A 67.3205978423% lower bound for simple zeros of the Riemann zeta function

This repository records and reproducibly verifies

$$
\liminf_{T\to\infty}\frac{N_0^s(T,2T)}{N(T,2T)}
\;\ge\; 0.6732059784228011963\ldots \;>\; \frac{6732059784}{10^{10}},
$$

where $N(T,2T)$ counts nontrivial zeros of $\zeta$ with multiplicity and
$N_0^s(T,2T)$ counts simple zeros on the critical line.

[Proof](proof.md) · [Lineage and prior art](docs/lineage.md) ·
[Exact coefficients](candidate.json) · [Certificates](certificates/)

> [!IMPORTANT]
> This is a same-day research result, not a peer-reviewed theorem, and it is
> not a proof of the Riemann hypothesis. It is a bound on a *proportion*. The
> analytic inputs are imported, not reproved — see
> [Trust boundary](#trust-boundary) for exactly what is and is not
> established here.

## Where this sits

| Result | Bound |
| --- | ---: |
| Claude / Anthropic, Theorem D | 0.672500703679… |
| `ainta/zeta-simple-zeros` | 0.673008527927… |
| `trmdy/zeta-simple-zeros-673137` | 0.673137630699… |
| `tawanerguo-cn/zeta-simple-zeros` (independent branch) | 0.673192911473… |
| `npip99/zeta-zeros` | 0.673195198901… |
| **this repository** | **0.673205978423…** |

Full table with links and routes: [docs/lineage.md](docs/lineage.md).

## The one new idea

The predecessor construction sums a local seven-point inequality over every
window of an $m$-point block. Its linear gap term uses a **uniform** pressure
coefficient. The global accounting charges that term through a single
quantity — the *total* pressure $B=\sum_r b_r$ — because a given block gap
occupies each of the six local positions at most once.

So the six coefficients can be redistributed freely, subject only to holding
their exact total at

$$
\sum_{r=1}^{6} b_r=\frac{3}{1150},
$$

and the global penalty does not notice. The local minimum, however, does.
Redistributing lifts the certified local target from

$$
\varepsilon=\tfrac{1}{200}=0.005
\qquad\longrightarrow\qquad
\varepsilon=\tfrac{51063}{10^7}=0.0051063 ,
$$

which is what produces the improved bound. A sharp finite-$m$ replacement
$h_m$ for the block-defect profile contributes a further $3\times10^{-5}$
percentage points.

This matters because the uniform-pressure branch is already saturated:
`npip99/zeta-zeros` reports the true adversarial minimum of its functional as
$\approx0.0050910$, so no amount of additional search on that construction
reaches $0.0051063$. The redistribution is what opens the door, not more
compute. [docs/lineage.md](docs/lineage.md) isolates the contribution of each
lever arithmetically.

## Independent reproduction

Reproduced on a second platform (Windows 11 / MSYS2 / MinGW-w64, Python
3.12.5, mpmath 1.3.0) against a certificate recorded on Linux 6.18 / Debian
g++ 14.2. Full log: [`certificates/reproduction-windows.txt`](certificates/reproduction-windows.txt).

- All six `mpmath.iv` interval tables were **regenerated from source**, not
  copied, and came out **byte-identical** to the recorded SHA-256 digests.
- The headline search returned `VERIFIED=true` with **every counter matching
  the recorded run exactly** — 2,037,322 nodes, 1,018,733 pruned, 1,018,589
  splits, 1,300,541 convexity proofs, 554,949 tangent prunes, depth 73.
- The tangent-free hardening run ($F\ge0.00507$, no convexity/tangent
  machinery at all) matched exactly too: 363,259,330 nodes, depth 68.
- Re-running at 12 threads instead of 5 gave identical counters, confirming
  the search is thread-count independent.
- `src/check_final_bound.py` output was character-for-character identical to
  [`certificates/final-arithmetic.txt`](certificates/final-arithmetic.txt).

One portability fix was needed to compile outside Debian: `std::fabsl` is not
surfaced by libstdc++'s `<cmath>`, so it is spelled `std::fabs` (which has a
`long double` overload — the same operation). That single line is why this
tree's `verify_positioned.cpp` digest differs from the one in the original
certificate; both digests are recorded.

## Reproduce it yourself

Requires Python 3.10+ with `numpy` and `mpmath`, and a C++20 compiler.

The interval tables are committed for convenience, with their digests recorded
in the certificates. Regenerating them from source is the stronger check, and
it overwrites `tables/` in place — the digests should not move:

```bash
python3 src/build_tables.py --output tables --workers 8
```

```bash
g++ -O3 -std=c++20 -pthread -ffp-contract=off src/verify_positioned.cpp -o verify_positioned
```

```bash
./verify_positioned 51063 10000000 8 tables tangent
```

```bash
./verify_positioned 507 100000 8 tables no-tangent
```

```bash
python3 src/check_final_bound.py
```

Do not compile with `-ffast-math` or without `-ffp-contract=off`; the verifier
depends on IEEE-754 semantics and on FMA contraction being disabled.

Table generation takes a few minutes. The headline search takes about 10
seconds on 5 threads; the tangent-free run about 35 seconds.

## Trust boundary

**Established in this repository.** The finite seven-point inequality
$F\ge51063/10^7$ over all nonnegative gaps, by exhaustive outward-rounded
interval subdivision, in two independent modes (with and without the
convexity/tangent accelerator). The exact rational bookkeeping — pressure
total $3/1150$, all six span capacities equal to 2 — is re-derived and
asserted by the verifier at startup. The final deduction arithmetic, at 80
digits of interval precision.

**Imported, not reproved here.** The analytic interface — the window constant
$H(v)\ge672457041414544284/10^{18}$, the stability argument giving
$S\ge H_{\rm cert}N+\operatorname{tr}\Psi(M)-o(N)$, and the shifted-block
pinching and averaging — comes from Claude/Anthropic's Theorem D and the
`trmdy` re-optimized window. The window and its kernel are unchanged from
that predecessor. If any of those inputs is wrong, this bound falls with them.

**Not done.** No native Arb (`python-flint`) rerun of the new positioned
certificate. No Lean formalization of the new finite-dimensional deduction.
No independent human expert review. The finite-dimensional Gram profile
$h_m$ and the block deduction in [`proof.md`](proof.md) §3–§4 are the parts
most in need of that review; they are ordinary finite linear algebra, but
they are new here and have not been checked by a mathematician.

## Files

| Path | Contents |
| --- | --- |
| [`proof.md`](proof.md) | statement, local-to-global deduction, Gram profile lemma |
| [`docs/lineage.md`](docs/lineage.md) | prior art, and which lever contributes what |
| [`candidate.json`](candidate.json) | every exact rational coefficient |
| `src/build_tables.py` | `mpmath.iv` interval table generator |
| `src/verify_positioned.cpp` | exhaustive branch-and-bound interval verifier |
| `src/check_final_bound.py` | independent high-precision final arithmetic |
| `certificates/` | recorded runs, final arithmetic, second-platform reproduction |

## Attribution

The analytic framework is Claude/Anthropic's. The re-optimized seven-term
window, the weighted seven-point form, and the block-defect approach are from
`ainta/zeta-simple-zeros` and `trmdy/zeta-simple-zeros-673137`. This
repository contributes the position-dependent pressure certificate, the sharp
finite-$m$ profile, and the reproduction above.
