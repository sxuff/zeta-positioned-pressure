# Lineage and prior art

Every entry below bounds the same quantity,

$$
\liminf_{T\to\infty}\frac{N_0^s(T,2T)}{N(T,2T)},
$$

the proportion of simple zeros on the critical line, counted against all
nontrivial zeros with multiplicity. Values are the ones each source
publishes for itself.

| # | Source | Bound | Route |
| --- | --- | ---: | --- |
| 1 | Claude / Anthropic, Theorem D ([`anthropics/zeta-23-lean`](https://github.com/anthropics/zeta-23-lean)) | 0.672500703679… | Montgomery–Taylor window, rank–trace argument |
| 2 | [`ainta/zeta-simple-zeros`](https://github.com/ainta/zeta-simple-zeros) | 0.673008527927… | seven-point stability refinement |
| 3 | [`learademacher/ai-refines-ai-zeta-bound`](https://github.com/learademacher/ai-refines-ai-zeta-bound) | 0.673021361950… | retarget to 191/50000 |
| 4 | [`tekuuu/zeta_sz`](https://github.com/tekuuu/zeta_sz) (GPT-5.6 Pro) | 0.673046504225… | refinement of the ainta argument |
| 5 | [`trmdy/zeta-simple-zeros-673137`](https://github.com/trmdy/zeta-simple-zeros-673137) | 0.673137630699… | re-optimized 7-term window, weighted 7-point, $h(E)=2\sqrt E-1$ |
| 6 | [`tawanerguo-cn/zeta-simple-zeros`](https://github.com/tawanerguo-cn/zeta-simple-zeros) | 0.673192911473… | independent branch: Bellman coboundary, $\cos(1.47s)$, $m=183$ |
| 7 | [`npip99/zeta-zeros`](https://github.com/npip99/zeta-zeros) | 0.673195198901… | raises trmdy's local target to $F\ge509/10^5$, $m=250$ |
| **8** | **this repository** | **0.673205978423…** | **position-dependent pressure, $F\ge51063/10^7$; sharp finite-$m$ profile $h_m$** |

Entries 5, 7 and 8 are the same lineage: 7 and 8 both raise the certified
local minimum $\varepsilon$ of trmdy's weighted seven-point inequality.
Entry 6 is a separate branch and does not share this construction.

A ninth repository, `HRK44/riemann-hypothesis-v1`, advertises 68.220195%, but
it contains only a PDF — no certificate, no verifier — and the figure is
described there as a *finite simple-point value*, not the asymptotic
proportion bounded above. It also exceeds the 0.681828687464… ceiling that
trmdy computes for this certificate class. It is not a comparable claim.

## Why the uniform-pressure branch cannot reach this number

Entry 7 keeps trmdy's *uniform* linear pressure and pushes the certified
local minimum as far as that construction allows. Its own README records the
true adversarial minimum of its functional as $\approx0.0050910$, so
$\varepsilon=0.00509$ is within $10^{-6}$ of saturation: the uniform branch
is finished, at about 0.6731955.

This bundle redistributes the six pressure coefficients while holding their
exact total at $3/1150$. The global span penalty depends only on that total,
so the redistribution is free in the global accounting, and it lifts the
local minimum to $0.005106397636960345\ldots$, certified at
$\varepsilon=51063/10^7$.

Reproducing the arithmetic of the deduction across both profiles isolates
where the gain comes from (percent, best block length $m$ scanned over
$50\le m<4000$):

| $\varepsilon$ | pressure | with $h(E)=2\sqrt E-1$ | with sharp $h_m$ |
| --- | --- | ---: | ---: |
| $1/200$ | uniform | 67.313767% ($m=257$) | — |
| $509/10^5$ | uniform | 67.319524% ($m=250$) | 67.319555% ($m=251$) |
| $51063/10^7$ | positioned | 67.320567% ($m=249$) | **67.320598%** ($m=250$) |

Two things follow. First, the sharp finite-$m$ profile is worth only about
$3\times10^{-5}$ percentage points on its own — the lead is the pressure
redistribution, not the profile. Second, a uniform-pressure argument would
need $\varepsilon\ge0.0051063$ to match this bound, which is above its own
adversarial minimum. It cannot get there without adopting position-dependent
pressure.

## Where this bundle is itself saturated

The certified target $0.0051063$ sits $9.8\times10^{-8}$ below this
construction's own observed adversarial minimum $0.005106397636960345$.
Nearly all of the available slack in the positioned seven-point family has
been taken. Further movement needs a new lever — more points, a different
window, or the independent Bellman-coboundary route of entry 6 — not a
longer search.
