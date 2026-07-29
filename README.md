# gitrialsim

Simulation-guided trial design for gastroenterology and hepatology.

`gitrialsim` helps answer a specific question: *what would the definitive trial of this GI question
actually have to look like?* It holds clinical scenarios as versioned, fully cited parameter packs,
and evaluates candidate designs against them, delegating to established engines wherever they cover
the case and simulating only where they do not.

## What this is not

It is not a general-purpose trial-design package, and it does not try to be. `rpact`, `gsDesign`,
East and FACTS already occupy that space and do it better. This package computes almost no statistics
of its own: sample sizes and group-sequential boundaries come from rpact, and are cross-checked
against gsDesign. What it adds sits in two layers those tools deliberately leave empty, namely the
GI-specific parameterization and the orchestration and reporting around a particular clinical
question.

That constraint is the point. Trial-design software written by someone who is not a biostatistician
is worth exactly as much as its agreement with software that statisticians already trust, so this
package is built to be checkable rather than clever.

## Installation

```r
# install.packages("remotes")
remotes::install_github("tarek-arabi/gi-trial-simulator")
```

R 4.1 or later. Imports `rpact`, `gsDesign`, `yaml` and base `stats`. `ggplot2` and `shiny` are
optional and only needed for figures and the explorer app.

## Thirty seconds

```r
library(gitrialsim)

s <- scenario("ercp_acute_cholangitis")
design_fixed(s)
```

```
<gi_design> fixed
  ercp_acute_cholangitis / mortality_30d
  All-cause mortality at 30 days
    Early ERCP (24 to 48 hours)    0.0658
    Urgent ERCP (within 24 hours)  0.0395
  alpha 0.025 (1-sided)   power 0.900
  maximum n 3,028 total, 1,514 per arm
  engine: rpact::getSampleSizeRates
```

The scenario's event rates are the observed primary-endpoint rates from the largest randomised trial
of ERCP timing in acute cholangitis, which enrolled 304 patients. Powering that same comparison
properly takes roughly 3,000. That gap is the motivating observation for the flagship analysis in
`analysis/ercp_flagship/`.

A group-sequential design trades maximum sample size for expected sample size:

```r
g <- design_group_sequential(s, k = 3)
gs_boundaries(g)
```

```
  analysis information_rate n_cumulative efficacy_z cumulative_alpha_spent
1        1        0.3333333         1022   3.710303           0.0001035057
2        2        0.6666667         2043   2.511427           0.0060483891
3        3        1.0000000         3064   1.993047           0.0249999900
```

Maximum 3,064 against the fixed design's 3,028, but an expected 2,456 under the alternative.

This example has no futility boundary, which is the default. Aim A2 of the flagship analysis adds a
non-binding O'Brien-Fleming futility boundary and so reports a larger maximum, 3,208, and a slightly
larger expected 2,490; buying the option to stop early for futility costs about 4.7 percent on the
maximum. The two sets of numbers describe different designs rather than disagreeing.

## Parameter packs

A pack is a YAML file describing one clinical scenario: its arms, its endpoints with control and
treatment event rates, and the published source every value came from. Packs are data, not code, so
they can be versioned, cited and audited independently of the engine that reads them.

```r
list_packs()
```

Three packs ship with the package, all built only from published aggregate results:

| Pack | Question | Primary source |
|---|---|---|
| `ercp_acute_cholangitis` | Urgent versus early ERCP in acute cholangitis | Jagtap et al., *Gut* 2026, [10.1136/gutjnl-2025-337279](https://doi.org/10.1136/gutjnl-2025-337279) |
| `hrs_terlipressin` | Terlipressin versus midodrine and octreotide in hepatorenal syndrome | Cavallin et al., *Hepatology* 2015;62:567-574, [10.1002/hep.27709](https://doi.org/10.1002/hep.27709) |
| `dph_colonoscopy_sedation` | Diphenhydramine as an adjunct for colonoscopy sedation | Nusrat et al., *Gastrointest Endosc* 2018;88(4):695-702, [10.1016/j.gie.2018.04.2342](https://doi.org/10.1016/j.gie.2018.04.2342) |

Every event rate carries a citation and a `verified` date recording when it was last checked against
its primary source.

Packs are found on a search path, so a pack held outside this repository can extend or override a
shipped one without modifying the package:

```r
options(gitrialsim.pack_paths = "~/my-packs")   # or set GITRIALSIM_PACK_PATH
load_pack("my_scenario")
```

This separation is deliberate. Parameters derived from restricted-use data sources cannot be
published here, and this is how they are used without ever entering the repository.

## What it can do

**Designs.** Fixed and group-sequential designs via rpact, with O'Brien-Fleming and Pocock boundaries
in both classical and alpha-spending forms, and optional futility stopping. A Bayesian adaptive
design with conjugate Beta-binomial posteriors and interim posterior-probability stopping rules,
whose efficacy threshold is calibrated by simulation because it has no closed form. That is the
pattern the FDA's Complex Innovative Trial Design programme contemplates: designs whose operating
characteristics must be established by simulation rather than derived.

**Simulation and reporting.** A vectorised Monte Carlo engine with reproducible L'Ecuyer-CMRG
streams, so results do not depend on the number of parallel workers. Performance measures are
reported following the ADEMP framework of Morris, White and Crowther (*Stat Med* 2019,
[10.1002/sim.8086](https://doi.org/10.1002/sim.8086)), with a Monte Carlo standard error on every
estimate and a helper for justifying the number of replications.

**Virtual cohorts and prognostic adjustment.** Patients with correlated baseline covariates generated
by a Gaussian copula from aggregate marginals, and a PROCOVA-style analysis in which a prognostic
score trained on historical control data enters the trial analysis as a pre-specified covariate. The
package measures the resulting variance reduction and, importantly, checks that the adjustment does
not inflate the type I error rate.

**Emulation and value of information.** A Gaussian-process emulator over simulation output, so the
design space can be explored without re-running Monte Carlo at every point, with posterior standard
deviations so it is visible where the surface is a guess. Expected value of perfect and sample
information, for the question of whether the definitive trial is worth running at all.

**An explorer app.** `run_explorer()` opens a Shiny interface over the same public functions.

## Validation

This is the part that matters most, so it is checked rather than asserted.

Group-sequential boundaries computed through rpact are compared against `gsDesign` across
O'Brien-Fleming and Pocock boundaries, classical and alpha-spending, for two through five analyses.
The largest disagreement across all sixteen configurations is under 1e-6 on the z scale:

```r
benchmark_report()
```

Beyond that cross-engine check, the test suite checks that simulating under the null recovers the
nominal type I error, tests the conjugate posterior probability against exact numerical integration,
and confirms that prognostic-score adjustment preserves the type I error rate.

The simulator is also pinned against rpact's analytic power, and that comparison is worth stating
precisely rather than as a claim of exact agreement. Across five scenarios at 200,000 replications
the simulator sits consistently 0.3 to 0.7 percentage points **above** rpact, always in the same
direction. That is a systematic difference rather than noise, and neither tool is wrong: rpact
evaluates a normal-approximation formula while the simulator simulates the actual discrete binomial
test. The claim the package makes is therefore that the two agree to within about 0.01 on the
probability scale, with the simulator biased slightly high. Details in
[analysis/ercp_flagship/REPORT.md](analysis/ercp_flagship/REPORT.md).

```bash
Rscript -e "pkgload::load_all('.'); testthat::test_local()"
```

Independent verification is welcome and useful. If any quantity here disagrees with rpact, gsDesign,
East or FACTS, please open an issue; disagreement with an established engine is treated as a defect
in this package until shown otherwise.

## Analyses

`analysis/ercp_flagship/` contains a full simulation study of the definitive ERCP-timing trial. Its
ADEMP protocol was written before the study was run, though it entered the repository in the same
commit as the results, so the history does not independently corroborate that ordering; the report
says so in its own header. Protocols for later studies are committed on their own and registered
externally beforehand, which is what makes the claim checkable rather than merely asserted.

## Citation

See `CITATION.cff`, or:

```r
citation("gitrialsim")
```

## License

GPL-3. See `LICENSE.md`.

The maintainer holds sole copyright. Contributions require a contributor licence agreement before
merge; see `CONTRIBUTING.md` for the reasoning. For commercial licensing enquiries, contact
arabi.tarek@gmail.com.
