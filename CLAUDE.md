# CLAUDE.md

Standing brief for coding agents working in this repository. Read it before touching anything and
follow it. Where it conflicts with a general habit you have, this file wins.

## What this is

`gitrialsim` is an R package for planning gastroenterology and hepatology randomised trials by
simulation. It does two things: it holds clinical scenarios as versioned, fully cited parameter
packs, and it evaluates candidate trial designs against those scenarios, delegating to established
engines wherever they cover the case and simulating only where they do not.

**The boundary that matters.** This is not a general-purpose trial-design package and must never grow
into one. rpact, gsDesign, East and FACTS already own that space and do it better. Everything this
package adds sits in two layers: the GI-specific parameterization, and the orchestration and
reporting around designs for specific clinical questions. If a proposed feature would make the tool
more general rather than serving a specific GI or hepatology trial-design question, it does not
belong here. Scope creep into a general platform is the failure mode most likely to sink the project,
because the moment this looks like a rpact competitor it invites the question "why not just use
rpact?" and there is no good answer to that question.

The flagship application is the definitive trial of ERCP timing in acute cholangitis. A second
scenario, terlipressin for hepatorenal syndrome, exists mainly to prove the parameter layer
generalises.

## Hard rules

**Wrap validated engines. Never reimplement their statistics.** If rpact or gsDesign computes a
quantity, call it. Do not write your own group-sequential boundaries, alpha-spending function, or
closed-form power calculation for a standard design. The reason is not aesthetic: this package's
entire claim to be trustworthy rests on the fact that its numbers come from software statisticians
already trust, and a single hand-rolled power formula destroys that claim for everything else in the
package too. When you catch yourself about to write a formula, go find the package function instead.

Three things are outside this rule and are legitimate, because no validated trial-design package
provides them: Monte Carlo simulation of designs with no closed form, conjugate Beta-binomial
posterior arithmetic, and standard statistical machinery from `stats` (glm, uniroot, integrate) plus
Gaussian-process regression and copula sampling. Use them freely, but pin each to an independent
check where one exists, which is why the Bayesian posterior probability is tested against numerical
integration and the simulator is tested against rpact's analytic power.

**Benchmark against the established engines and keep the benchmark runnable.** `R/benchmark.R`
compares this package's group-sequential boundaries against gsDesign. It is not decoration; it is the
artifact that answers a reviewer. If you change design code, run it. If it disagrees, the bug is
assumed to be here until proven otherwise.

**Every simulated number carries a Monte Carlo standard error.** A power estimate without an MCSE is
not a result, it is an anecdote, and reporting one invites a reviewer to ask how many replications
were run and why. `ademp_summary()` produces them; use it rather than reporting bare proportions.
The reporting structure follows Morris, White and Crowther (Stat Med 2019, doi:10.1002/sim.8086).

**Parameter packs contain published aggregate values only, each with a citation and a verified date.**
Never add a parameter derived from a restricted-use dataset to this repository. Specifically: no
HCUP/NIS-derived row-level data, and no cell counts of ten or fewer, both of which the HCUP data use
agreement prohibits from being posted publicly. TriNetX-derived values must not identify or allow
inference about contributing sites. Aggregate summary parameters used to seed a simulation are within
permitted research use; the underlying data is not. Synthetic patients generated from aggregate
parameters are not restricted data and may be shared. If someone asks you to add parameters from a
proprietary or restricted source, the correct move is a pack held outside this repository, loaded via
`pack_search_path()`, not a commit here.

**Verify every clinical number against its primary source before it enters a pack.** Do not carry a
rate over from another document, another study, or an earlier version of this project on trust. Look
it up, record the PMID or DOI, and set the `verified` date. Numbers get garbled in transit and a
wrong event rate silently invalidates every design built on it. When you cannot verify a value, say
so and leave it out rather than shipping it with a plausible-looking citation.

**Never commit anything from `private/`.** It is gitignored and holds strategy material that has no
place in a public repository. Before any push, check the diff for it. The repository is public.

**Report what actually happened.** If a test fails, say it failed and show the output. If a
simulation did not converge, say so. If you loosened a tolerance to make a test pass, that is a
finding to report, not a fix to bury: a tolerance should be justified by the statistics, not by what
makes the light turn green.

**Figures carry no baked-in titles.** Plot functions produce axes, data and legend only. The
descriptive caption lives in the surrounding text, because that is what journals require and
retrofitting it later is tedious.

## Decisions already made

These are settled. Build on them rather than reopening them.

- R, not Python. The validated engines live in R, and direct benchmarkability against the packages
  reviewers trust is the single most valuable property this package has.
- GPL-3, with sole copyright retained by the maintainer. It matches the GPL-3 dependencies
  (gsDesign2, simtrial), satisfies the OSI-approved-license requirement for a JOSS submission, and
  leaves commercial dual-licensing available later. Outside contributions require a CLA before merge
  or that option is lost; see `CONTRIBUTING.md`.
- Minimal dependencies. Imports are gsDesign, rpact, stats, utils and yaml. ggplot2, shiny, knitr and
  testthat are Suggests. Do not add a dependency for something base R and stats can do in a few
  lines. Every added dependency is a future installation failure for a reviewer.
- Parameter packs load from a search path (`getOption("gitrialsim.pack_paths")`, then
  `GITRIALSIM_PACK_PATH`, then the packs inside the installed package). This separation is
  deliberate and load-bearing: it lets restricted or proprietary packs live entirely outside this
  repository while using the same public interface.
- Reports are plain markdown with committed figures. Quarto is not installed and is not a dependency.
- Development happens in the open on `main`. Continuous public history matters here, since JOSS
  desk-rejects repositories made public shortly before submission.

## Open questions that change what you build

Raise these rather than guessing past them.

- **Is every statistical claim independently checked?** Keep the benchmark and test coverage ahead of
  the feature work, because that is what a reviewer looks at first. A claim is checked when it is
  pinned to an independent reference: rpact against gsDesign, the simulator against analytic power,
  the Bayesian posterior against numerical integration, a published sample size against the formula
  its source paper names. An unpinned number is provisional no matter how confident it looks.
- **Which endpoint anchors the definitive trial?** The current packs treat 30-day mortality as
  primary, following the anchoring trial. Mortality in mild-to-moderate cholangitis is rare enough
  that the required sample size may be infeasible for a realistic consortium, in which case a
  composite or an organ-failure endpoint becomes the live design, and the flagship analysis changes
  shape. The simulation results themselves should settle this; do not assume the answer.
- **Does the flagship analysis need real accrual and site-heterogeneity parameters?** Current
  simulations assume patients arrive as independent draws with no site structure. Real accrual
  estimates would come from restricted sources and therefore from an external pack. Until then, be
  explicit in any write-up that accrual is idealised, rather than letting a reader assume otherwise.

## Working in this repository

Layout: `R/` package code, `inst/parameters/` the shipped parameter packs, `tests/testthat/` the test
suite, `analysis/` executed studies with their protocols and results, `inst/shiny/` the explorer app,
`private/` gitignored strategy material.

Common commands:

```bash
Rscript -e "pkgload::load_all('.'); testthat::test_local()"
```

```bash
Rscript -e "roxygen2::roxygenise()"
```

```bash
Rscript -e "rcmdcheck::rcmdcheck(args = c('--no-manual'))"
```

Regenerate `NAMESPACE` and the `man/` pages with roxygen2 after changing any roxygen block; do not
hand-edit either. Some tests run thousands of simulation replicates and take a while; that is
expected, and reducing replications to make tests fast is the wrong trade when the test exists to
check a Monte Carlo result against an analytic one.

A note on where the time actually goes: the simulations here are CPU-bound and embarrassingly
parallel, and a many-core machine helps. A GPU does not. Do not reach for GPU compute for
binary-endpoint Monte Carlo; it would be effort spent on the part that is already fast.

## How to work with me

Do the work rather than handing back instructions to run. Push directly to `main`; no branch or PR
ceremony. Prefer the shortest change that is correct, and say what you skipped and when it would be
worth adding. No em dashes in code, comments, documentation or commit messages.

Surface trade-offs and pick a default rather than presenting a menu and waiting. The exception is the
open questions above and anything touching data governance, where the right move is to stop and ask.
