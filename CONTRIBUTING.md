# Contributing

Challenges to this work are welcome. The most useful contribution is a report of a suspected error.

## Reporting a suspected error

Open an [issue](https://github.com/driss-bencheikh/spacetime-electromagnetic-framework/issues) and include:

1. **Where.** The version (currently 1.0.1), the theorem, proposition, lemma, corollary, definition or equation number, and the printed
   page number of the PDF in `release/`.
2. **What kind of problem.** One of the categories used in [`docs/errata.md`](docs/errata.md):
   - **Error**: a statement false as written;
   - **Hypothesis**: a missing condition or domain restriction;
   - **Wording**: imprecise or misleading phrasing;
   - **Attribution**: a missing or incomplete reference, including prior work that obtained the same result;
   - **Evidence**: a proof or check weaker than the statement it supports.
3. **Why.** If possible, a counterexample (explicit matrices, vectors or parameter values) or a script that fails. A script should print
   `PASS` or `FAIL` lines like those in `verification/` and run with the packages in `verification/requirements.txt`.

One issue per problem makes each easier to check and to credit.

## What happens next

Each report is checked against the manuscript and, where a computation is involved, by running it. An accepted report is added to
[`docs/errata.md`](docs/errata.md) with its category and credit to the reporter (by name or GitHub handle, as you prefer), and is corrected
in the next version. A report that is not accepted is answered in the issue with the reason.

## Other contributions

Pull requests that fix a script, add a check or correct a typo in the documentation are welcome. Changes to the manuscript text are made
by the author from accepted reports, so please open an issue for those rather than a pull request.
