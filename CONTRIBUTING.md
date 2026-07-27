# Contributing

Thank you for your interest. This project is developed in the open and welcomes issues, bug reports,
validation failures, and parameter-pack corrections.

## Reporting a problem

The most valuable reports are ones that show a number is wrong:

- **A design quantity disagrees with rpact, gsDesign, East or FACTS.** Please open an issue with the
  exact call and both outputs. Disagreement with an established engine is treated as a defect in this
  package until proven otherwise.
- **A parameter pack misstates a published value.** Please cite the primary source, including the
  page or table. Every shipped value carries a `verified` date and is meant to be checkable.
- **A simulation result is not reproducible from the recorded seed.**

## Pull requests

Before a pull request can be merged, its author must agree to a contributor licence agreement
assigning copyright in the contribution to the maintainer, or granting a licence broad enough to
permit relicensing.

This is not bureaucracy for its own sake, and it is worth being straightforward about the reason.
The package is released under the GPL-3. The maintainer holds sole copyright, which keeps open the
option of also offering the software under a separate commercial licence in future. Merging a
contribution whose copyright sits elsewhere would remove that option permanently for everyone,
including for work that has nothing to do with the contribution. If you would rather not sign, please
still open the issue: a well-described defect is more useful than a patch, and the fix can be written
independently.

## Standards for code contributions

- **Do not reimplement statistics that rpact or gsDesign already provide.** Call the package. This is
  the project's central design constraint and the basis of its credibility. Simulation of designs
  with no closed form, conjugate posterior arithmetic, copula sampling, and Gaussian-process
  regression are outside that rule and are welcome.
- Every new numerical function needs a test that would fail if the number were wrong. Tests that only
  check the return type are not sufficient.
- Every simulated result must carry a Monte Carlo standard error.
- Document limitations in the function's own documentation rather than in a commit message.
- No parameter derived from a restricted-use dataset may be added to this repository. See the data
  provenance rules in `CLAUDE.md`.
