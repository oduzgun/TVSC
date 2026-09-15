# Temporally Regularized Synthetic Control

R functions in development for preparing panel data and supporting a chronological synthetic-control workflow.

The project accompanies two manuscripts: a methodological paper on temporally regularized synthetic control and a companion paper on the R workflow. Current source drafts cover **raw panel-data preparation, chronological split construction and cutoff-specific raw predictor blocks**; scaling, estimation and treatment-effect functions are planned.

## Development status

- `panel.dataprep()` is available as an initial source-code draft.
- `panel.split()` is available as an initial source-code draft.
- `panel.blocks()` is available as an initial source-code draft.
- Base-R test fixtures are provided in `tests/test-panel.dataprep.R`.
- Base-R split fixtures are provided in `tests/test-panel.split.R`.
- Base-R block fixtures are provided in `tests/test-panel.blocks.R`.
- This repository is not yet an installable R package or a validated software release.
- No empirical results are supplied by this implementation.

## Available function: `panel.dataprep()`

Prepare a consistently ordered, explicitly selected pre-treatment panel for one target unit and at least two donors. Preserve the supplied measurement values and retain available post-treatment outcomes separately.

The function checks identifiers, dates, panel completeness and required training values. It deliberately does **not** create validation splits, standardize predictors, estimate weights, impute missing values or calculate treatment effects.

### Loading the function

The function uses base R and requires no additional R packages. From the repository root:

```r
source("R/panel.dataprep.R")
```

This is a source-script repository, not an R package. A minimum supported R version has not yet been established through release testing.

### Quick start

The following artificial data illustrate the interface; they are not an empirical application or a simulation study.

```r
source("R/panel.dataprep.R")

panel <- expand.grid(
  state = c("Target", "DonorA", "DonorB"),
  year = 2000:2004,
  stringsAsFactors = FALSE
)
panel$sales <- seq_len(nrow(panel)) * 10
panel$income <- 1000 + seq_len(nrow(panel))

prepared <- panel.dataprep(
  data = panel,
  unit = "state",
  time = "year",
  outcome = "sales",
  treated = "Target",
  donors = c("DonorB", "DonorA"),
  intervention_time = 2003,
  pre_period = 2000:2002,
  predictors = c("sales", "income")
)

prepared$training
prepared$schema$donors
prepared$future_donors
prepared$diagnostics
```

The intended result contains nine training rows: three units observed at three pre-treatment dates. Within each date, the target comes first, followed by `DonorB` and `DonorA` in the supplied donor order. Treat this example as unvalidated until the test command or GitHub Actions workflow succeeds in an R environment.

### Usage

```r
panel.dataprep(
  data, unit, time, outcome, treated, donors,
  intervention_time, pre_period, predictors,
  period_step = 1
)
```

### Arguments

| Argument | Description |
| --- | --- |
| `data` | Nonempty data frame with unique, nonempty column names. |
| `unit` | Name of the unit-identifier column. IDs must be character or plain numeric values; factors are matched using their character labels. |
| `time` | Name of a finite, integer-valued numeric period-index column. Calendar `Date` and date-time objects are not supported in this draft. |
| `outcome` | Name of the outcome column. Required training values must be numeric and finite. |
| `treated` | One target-unit ID, matching the character/numeric type of the unit column. |
| `donors` | At least two distinct donor IDs, excluding the target. The supplied order is preserved in the schema and training panel. |
| `intervention_time` | First treated period, expressed as one integer-valued numeric index. |
| `pre_period` | Explicit sequence of at least three increasing, equally spaced dates, all before treatment. No pre-period or validation split is inferred. |
| `predictors` | Explicit, nonempty vector of unique predictor-column names. May include the outcome; cannot include unit or time identifiers. Required training values must be numeric and finite. |
| `period_step` | Positive integer-valued spacing between successive `pre_period` dates; defaults to `1`. |

The `unit`, `time` and `outcome` arguments must name distinct columns. The outcome is retained and checked even when it is not included among the matching predictors.

### Returned object

A list of class `tvsc_prepared` with the following components:

| Component | Contents |
| --- | --- |
| `schema_version` | Object-schema version, currently `1L`. |
| `schema` | Column roles, predictor names, target and ordered donor IDs, intervention date, pre-period and period spacing. |
| `training` | Selected pre-treatment observations in supplied measurement units, ordered by date and then target/donor order. Contains identifier, outcome and predictor columns. |
| `future_donors` | Available donor observations from the intervention date onward, retaining unit, time and outcome columns. |
| `future_target` | Available observed target outcomes from the intervention date onward, stored separately from training. These are not untreated counterfactuals. |
| `recipe` | Contemporaneous-block declaration, deferred scaling (`NULL`) and an indicator of whether the outcome is a predictor. |
| `diagnostics` | Training-row count, excluded-row count, constant-predictor names and `future_outcomes_validated = FALSE`. |

`scaling = NULL` records that scaling has not been performed; it does not establish a default scaling policy for the future estimator. Predictor matrices and fitted weights are not returned.

### Validation and scope

- Every selected unit must have exactly one row at every declared pre-treatment date. Missing or duplicate training keys cause an error.
- Required training outcome and predictor values must be finite. Categorical predictors require explicit numeric encoding before preparation.
- Constant predictors are retained and flagged, not silently removed.
- Unit IDs and time indices are checked across the entire supplied data frame, including rows outside the selected training panel.
- Future outcome values and future unit-date completeness or uniqueness are **not** validated. Future tables retain input row order and require separate validation before prediction or scoring.
- Future observations may be absent entirely. Future target outcomes are not required to prepare the training panel.
- Rows outside the selected units and retained time windows are excluded. Additional pre-period history outside `pre_period` is not retained for later lag construction.

### Why splitting is separate

Data preparation describes the panel; validation splitting describes an analysis design. The separate `panel.split()` draft creates chronological training and assessment assignments from the prepared object.

The intended workflow permits fitting on the full declared pre-period when tuning choices are supplied. When chronological validation is used, any data-dependent transformations and calibration must instead be learned within the corresponding training prefix—not from the full pre-period before splitting. Those later steps are not implemented here.

## Available function: `panel.split()`

Construct expanding-window chronological splits from a prepared panel. Training and validation assignments are dates, not individual panel rows: all selected units belong to the same window at each assigned date.

The function uses schema metadata from a schema-version-1 `tvsc_prepared` object. It does not modify the prepared object, revalidate the raw panel, transform values, fit weights, choose tuning parameters, score folds or reserve an independent test set.

### Loading the splitter

From the repository root:

```r
source("R/panel.dataprep.R")
source("R/panel.split.R")
```

### Split example

```r
panel <- expand.grid(
  state = c("Target", "DonorA", "DonorB"),
  year = 2000:2009,
  stringsAsFactors = FALSE
)
panel$sales <- seq_len(nrow(panel)) * 10

prepared <- panel.dataprep(
  panel, "state", "year", "sales", "Target", c("DonorB", "DonorA"),
  2008, 2000:2007, "sales"
)

splits <- panel.split(prepared, initial = 4, horizon = 2, step = 1)
splits$folds[[1]]
splits$diagnostics
```

With eight declared pre-treatment dates, `initial = 4`, `horizon = 2` and `step = 1` create three complete folds:

| Fold | Training dates | Validation dates |
| --- | --- | --- |
| 1 | 2000-2003 | 2004-2005 |
| 2 | 2000-2004 | 2005-2006 |
| 3 | 2000-2005 | 2006-2007 |

Validation dates may repeat across folds, and earlier validation dates may become training dates in later folds. This is intentional; folds are chronological assessment assignments, not independent samples.

### Split usage

```r
panel.split(prepared, initial, horizon, step = 1)
```

| Argument | Description |
| --- | --- |
| `prepared` | Unmodified schema-version-1 `tvsc_prepared` object produced by `panel.dataprep()`. |
| `initial` | Required number of dates in the first training window; must be at least three. |
| `horizon` | Required positive number of subsequent validation dates per fold. |
| `step` | Positive number of dates by which the cutoff advances; defaults to `1`. |

`initial`, `horizon` and `step` count dates in `schema$pre_period`, not long-panel rows and not numeric calendar-index units. For example, with biennial period indices, `horizon = 2` means two declared dates, not two calendar years.

### Returned split object

A list of class `tvsc_splits` with the following components:

| Component | Contents |
| --- | --- |
| `schema_version` | Split-object schema version, currently `1L`. |
| `schema` | Copy of the prepared panel schema. This is compatibility metadata, not a checksum of the underlying data. |
| `specification` | Window type, `initial`, `horizon`, `step` and `full_windows_only = TRUE`. |
| `folds` | Ordered list of folds. Each fold contains `id`, `train_index`, `validation_index`, `train_period`, `cutoff`, `validation_period` and `horizons`. |
| `diagnostics` | `fold_count`, date-named `assessment_counts` and `unassessed_period`. |

`train_index` and `validation_index` refer to positions in `schema$pre_period`, not data-frame row numbers. `unassessed_period` includes initial training-only dates and any dates omitted by the stride or complete-window policy.

### Split boundaries

- Every training window starts at the first declared pre-treatment date and expands as the cutoff advances.
- Validation immediately follows the cutoff.
- Only complete validation windows on the requested stride are returned.
- The final validation window is not shortened, and no extra off-stride cutoff is added.
- An error is raised if no complete fold fits.
- Training and validation are disjoint within a fold.
- Repeated assessment dates across folds are retained and reported through `assessment_counts`.
- Downstream fitting must still check schema compatibility and estimate any data-dependent transformations within each training prefix.

## Available function: `panel.blocks()`

Assemble contemporaneous predictor blocks for a chosen training prefix, preserving supplied values. Scaling remains a separate, explicitly specified operation, which may be skipped deliberately.

The function requires `panel.dataprep()` to be loaded. It revalidates required prefix records and numeric values through the existing preparation routine, but it does not modify the prepared object, use later predictor or outcome values, inspect future outcome tables, create lags or summaries, choose predictor weights, fit donor weights, score folds or calculate treatment effects.

### Loading the block builder

From the repository root:

```r
source("R/panel.dataprep.R")
source("R/panel.blocks.R")
```

### Block example

```r
panel <- expand.grid(
  state = c("Target", "DonorA", "DonorB"),
  year = 2000:2007,
  stringsAsFactors = FALSE
)
panel$sales <- seq_len(nrow(panel)) * 10
panel$income <- 1000 + seq_len(nrow(panel))

prepared <- panel.dataprep(
  panel, "state", "year", "sales", "Target", c("DonorB", "DonorA"),
  2006, 2000:2005, c("income", "sales")
)

blocks <- panel.blocks(prepared, cutoff = 2003)
blocks$X1[["2000"]]
blocks$X0[["2000"]]
blocks$Y0
```

### Block usage

```r
panel.blocks(prepared, cutoff)
```

| Argument | Description |
| --- | --- |
| `prepared` | Schema-version-1 `tvsc_prepared` object produced by `panel.dataprep()`, with a contemporaneous recipe and `scaling = NULL`. |
| `cutoff` | Required actual date in `prepared$schema$pre_period`, leaving at least three prefix dates. This is not a row number and not a date-position index. |

Use `tail(prepared$schema$pre_period, 1)` for full-pre-period assembly or a split fold's `cutoff` for fold-specific block construction.

### Returned block object

A list of class `tvsc_blocks` with the following components:

| Component | Contents |
| --- | --- |
| `schema_version` | Block-object schema version, currently `1L`. |
| `schema` | Validated roles, predictor order and donor order, with `pre_period` restricted to the selected prefix. |
| `periods`, `cutoff` | Prefix dates and final included date. |
| `dimensions` | Counts of predictors, donors and periods. |
| `X1` | Date-named list of named length-K target predictor vectors. |
| `X0` | Date-named list of K-by-J donor predictor matrices, with predictor rows and ordered donor columns. A single predictor remains a matrix. |
| `Y1` | Date-named target outcome vector. |
| `Y0` | Date-by-donor outcome matrix in the supplied outcome column's units. |
| `recipe` | Contemporaneous blocks, `scaling = NULL`, and outcome-in-predictors indicator. |
| `diagnostics` | Prefix row count and predictors constant across that prefix. Constant predictors are retained. |

The outcome is returned separately even when it is omitted from matching predictors. It is never silently added to `X1` or `X0`.

### Block boundaries and scaling

- The cutoff must match an actual declared pre-treatment date.
- The selected prefix expands from the first declared pre-treatment date through the cutoff and must contain at least three dates.
- Required unit-date keys and numeric outcome/predictor values are revalidated within the prefix.
- Later pre-treatment values and future tables are excluded from construction and numeric validation.
- Declared donor and predictor ordering are preserved.
- `recipe$scaling = NULL` means these functions applied no package scaling. It does not certify that the researcher supplied original-unit, untransformed values.
- Externally scaled predictor columns are preserved exactly as supplied.
- `Y1` and `Y0` remain in the supplied outcome column's units. To match on a scaled outcome while reporting an original-unit outcome, supply separate columns, for example `outcome = "sales"` and `predictors = c("sales_scaled", "income_scaled")`.
- Externally learned transformations must respect each validation training cutoff. Selecting a prefix here cannot repair leakage from full-period preprocessing.

## Tests

From the repository root, in an environment with R installed:

```sh
Rscript tests/test-panel.dataprep.R
Rscript tests/test-panel.split.R
Rscript tests/test-panel.blocks.R
```

The preparation fixtures cover ordering, incomplete or duplicated training panels, invalid identifiers and dates, nonfinite training values, constant predictors, outcome-only predictors, and separation of training from future outcome values. Split fixtures cover the eight-date example, complete-window and stride rules, date spacing, all-unit membership, unchanged input data, value-independent assignments and invalid arguments. Block fixtures cover exact cells, dimensions and names, row-order invariance, cutoff isolation, constant-predictor diagnostics, a single predictor, omitted outcomes, numeric and factor IDs, non-unit period spacing, malformed inputs, pre-scaled input preservation and transformed-outcome passthrough.

The tests use only base R. The preparation fixtures print `All raw panel preparation checks passed.` when successful. The split fixtures print `All chronological panel split checks passed.` when successful. The block fixtures print `All raw predictor block checks passed.` when successful. The GitHub Actions workflow runs all three test scripts and documented splitter/block examples.

## Roadmap

These names describe the intended interface, not currently callable functions:

| Function or component | Intended responsibility | Status |
| --- | --- | --- |
| `panel.dataprep()` | Validate and retain the raw selected panel. | Documented source draft with base-R fixtures. |
| `panel.split()` | Construct expanding-window training/assessment assignments. | Documented source draft with base-R fixtures. |
| `panel.blocks()` | Assemble contemporaneous, as-supplied predictor and outcome blocks for an explicit cutoff. | Documented source draft with base-R fixtures. |
| Scaling and predictor-weight calibration | Apply declared transformations and construct predetermined predictor metrics. | Planned; design choices remain open. |
| `tvsc.fit()` | Fit simplex donor-weight paths with first- or second-difference regularization. | Planned. |
| A `predict()` method | Continue fitted donor weights and construct counterfactual outcomes. | Planned. |
| `tvsc.att()` | Calculate date-specific gaps and average effects over explicit windows. | Planned. |

The next design checkpoint concerns the scaling interface. Package structure, release automation and distribution details will follow separately.

## Reporting problems

Please include a small reproducible data example, the function call, the complete error message, the expected behavior and the output of `sessionInfo()`. Do not include confidential or identifying data.

## License

This project is licensed under the MIT License. See `LICENSE` for details.
