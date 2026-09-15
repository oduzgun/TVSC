# Temporally Regularized Synthetic Control

R source scripts in development for a chronological synthetic-control workflow. The repository is not an R package and does not contain empirical applications, solver code, prediction code or treatment-effect calculations.

## Current Scope

The active version-2 workflow is:

1. `tvsc.dataprep()` prepares the selected raw panel and now returns integrated raw blocks.
2. `tvsc.split()` constructs expanding-window chronological date assignments from version-2 preparation. Its split object remains schema version 1.
3. `tvsc.vmetrics()` constructs prefix-specific donor-only scaling and predictor metrics.

The older `tvsc.blocks()`, `tvsc.scale()` and `tvsc.metrics()` scripts remain in the repository as version-1 drafts. They are not compatible with version-2 preparation and their old fixtures are legacy evidence, not current integration validation.

## Repository Layout

```text
README.md
LICENSE
R/tvsc.dataprep.R
R/tvsc.split.R
R/tvsc.vmetrics.R
R/tvsc.blocks.R        # legacy version-1 draft
R/tvsc.scale.R         # legacy version-1 draft
R/tvsc.metrics.R       # legacy version-1 draft
tests/test-tvsc.dataprep.R
tests/test-tvsc.split.R
tests/test-tvsc.vmetrics.R
tests/test-tvsc.blocks.R   # legacy version-1 fixture
tests/test-tvsc.scale.R    # legacy version-1 fixture
tests/test-tvsc.metrics.R  # legacy version-1 fixture
docs/tvsc-workflow-paper.tex
docs/tvsc-workflow-strategy.tex
```

No compatibility aliases are provided for earlier `panel.*` names.

## `tvsc.dataprep()`

`tvsc.dataprep()` selects and validates one treated unit, at least two ordered donors and an explicit pre-treatment period. It preserves values as supplied, stores future donor and observed target outcomes separately and constructs full-pre-period raw blocks.

```r
source("R/tvsc.dataprep.R")

panel <- expand.grid(
  state = c("Target", "DonorA", "DonorB"),
  year = 2000:2004,
  stringsAsFactors = FALSE
)
panel$sales <- seq_len(nrow(panel)) * 10
panel$income <- 1000 + seq_len(nrow(panel))

prepared <- tvsc.dataprep(
  data = panel, unit = "state", time = "year", outcome = "sales",
  treated = "Target", donors = c("DonorB", "DonorA"),
  intervention_time = 2003, pre_period = 2000:2002,
  predictors = c("sales", "income")
)

prepared$schema_version
prepared$training
prepared$blocks$X0[["2000"]]
```

The returned `tvsc_prepared` object contains `schema_version = 2L`, `schema`, `training`, `blocks`, `future_donors`, `future_target`, `recipe` and `diagnostics`. The integrated `blocks` object is also version 2 and contains `X1`, `X0`, `Y1`, `Y0`, `periods`, `cutoff`, `dimensions`, `recipe` and diagnostics.

Preparation does not split, scale, impute, aggregate, calibrate metrics, fit donor weights, predict counterfactuals or estimate effects. `recipe$scaling = NULL` means the function applied no package scaling; it does not certify that upstream data were never transformed.

## `tvsc.split()`

`tvsc.split()` accepts version-2 preparation and returns expanding-window date assignments. Window sizes count dates in `schema$pre_period`, not long-panel rows.

```r
source("R/tvsc.dataprep.R")
source("R/tvsc.split.R")

panel <- expand.grid(
  state = c("Target", "DonorA", "DonorB"),
  year = 2000:2009,
  stringsAsFactors = FALSE
)
panel$sales <- seq_len(nrow(panel)) * 10

prepared <- tvsc.dataprep(
  panel, "state", "year", "sales", "Target", c("DonorB", "DonorA"),
  2008, 2000:2007, "sales"
)

splits <- tvsc.split(prepared, initial = 4, horizon = 2, step = 1)
splits$folds[[1]]
splits$diagnostics
```

For eight declared dates, `initial = 4`, `horizon = 2` and `step = 1` produce:

| Fold | Training dates | Validation dates |
| --- | --- | --- |
| 1 | 2000-2003 | 2004-2005 |
| 2 | 2000-2004 | 2005-2006 |
| 3 | 2000-2005 | 2006-2007 |

The output is a `tvsc_splits` object with `schema_version = 1L`, `schema`, `specification`, `folds` and `diagnostics`. It does not modify the prepared object, inspect measurement values, reserve a test set or choose scoring weights.

## `tvsc.vmetrics()`

`tvsc.vmetrics()` combines prefix-specific scaling and predictor-metric construction for one cutoff. It uses the version-2 raw blocks stored by `tvsc.dataprep()`. Do not call the legacy `tvsc.blocks()`, `tvsc.scale()` or `tvsc.metrics()` for this version-2 path.

```r
source("R/tvsc.dataprep.R")
source("R/tvsc.vmetrics.R")

panel <- expand.grid(
  state = c("Target", "DonorA", "DonorB", "DonorC"),
  year = 2000:2005,
  stringsAsFactors = FALSE
)
panel$sales <- seq_len(nrow(panel)) * 10
panel$income <- 100 + seq_len(nrow(panel))

prepared <- tvsc.dataprep(
  panel, "state", "year", "sales", "Target",
  c("DonorC", "DonorB", "DonorA"), 2005, 2000:2004,
  c("sales", "income")
)

metrics <- tvsc.vmetrics(
  prepared, cutoff = 2003,
  predictor_weights = "uniform",
  scaling = "none"
)
metrics$weights
```

The signature is:

```r
tvsc.vmetrics(
  prepared, cutoff, predictor_weights = "uniform",
  scaling = "sd", conversion = "absolute", scales = NULL,
  provenance = NULL, rank_tol = NULL, weight_sum_tol = 1e-8
)
```

`cutoff` must be an actual declared pre-period date with at least three prefix dates. The exact character choices for `predictor_weights` are `"uniform"`, `"constant"` and `"varying"`; named numeric vectors and predictor-by-date matrices are also accepted. Supplied weights are aligned by names, checked for nonnegativity and column sums, and are not renormalized.

Scaling choices are `"sd"`, `"zscore"`, `"none"` and `"custom"`. Learned scaling uses donors only within the selected training prefix. `"sd"` uses one RMS within-date donor sample SD per predictor without centering. `"zscore"` uses each date's donor mean and donor sample SD. `"custom"` requires a finite positive named predictor vector. Reporting outcomes remain in the supplied outcome units.

For empirical `"constant"` and `"varying"` modes, donor-only OLS predicts next-date donor outcomes from the selected transformed features. Constant calibration uses common slopes and date intercepts. Varying calibration fits each date separately and retains the last available metric at the final training date. Rank failure, all-zero slopes and numerical failures raise `tvsc_vmetrics_error` rather than falling back to uniform weights.

The returned `tvsc_vmetrics` object has schema version 1 and contains the prefix scaled blocks, scaling metadata, weights, calibration records, diagnostics and unverified provenance. It does not fit donor weights, choose penalties, score folds, predict outcomes or estimate treatment effects.

## Validation

Run the active version-2 suites from the repository root:

```sh
Rscript --vanilla tests/test-tvsc.dataprep.R
Rscript --vanilla tests/test-tvsc.split.R
Rscript --vanilla tests/test-tvsc.vmetrics.R
```

GitHub Actions runs these three suites and executable examples for `tvsc.dataprep()`, `tvsc.split()` and `tvsc.vmetrics()`, then records `sessionInfo()` as an artifact.

The legacy version-1 fixture files remain present for reference, but they are not current integration checks after preparation moved to schema version 2. Rebuild old preparation objects from the original panel rather than editing version fields by hand.

## Roadmap

| Component | Status |
| --- | --- |
| `tvsc.dataprep()` | Active version-2 source draft with integrated raw blocks. |
| `tvsc.split()` | Active source draft for chronological assignments from version-2 preparation. |
| `tvsc.vmetrics()` | Active source draft for prefix-specific scaling and predictor metrics. |
| `tvsc.blocks()`, `tvsc.scale()`, `tvsc.metrics()` | Legacy version-1 drafts retained for comparison. |
| Donor fitting, prediction, ATT estimation and empirical applications | Planned; not implemented here. |

## License

This repository uses the MIT License. See `LICENSE`.
