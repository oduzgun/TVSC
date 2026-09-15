# Temporally Regularized Synthetic Control

R functions in development for preparing panel data and supporting a chronological synthetic-control workflow.

The project accompanies two manuscripts: a methodological paper on temporally regularized synthetic control and a companion paper on the R workflow. The current implementation covers **raw panel-data preparation only**; estimation and treatment-effect functions are planned.

## Development status

- `panel.dataprep()` is available as an initial source-code draft.
- Base-R test fixtures are provided in `tests/test-panel.dataprep.R`.
- This repository is not yet an installable R package or a validated software release.
- No empirical results are supplied by this implementation.

## Available function: `panel.dataprep()`

Prepare a consistently ordered, explicitly selected pre-treatment panel for one target unit and at least two donors. Preserve the original measurement values and retain available post-treatment outcomes separately.

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
| `training` | Selected pre-treatment observations in original measurement units, ordered by date and then target/donor order. Contains identifier, outcome and predictor columns. |
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

Data preparation describes the panel; validation splitting describes an analysis design. A separate, planned `panel.split()` function will create chronological training and assessment assignments from the prepared object.

The intended workflow permits fitting on the full declared pre-period when tuning choices are supplied. When chronological validation is used, any data-dependent transformations and calibration must instead be learned within the corresponding training prefix—not from the full pre-period before splitting. Those later steps are not implemented here.

## Tests

From the repository root, in an environment with R installed:

```sh
Rscript tests/test-panel.dataprep.R
```

The fixtures cover ordering, incomplete or duplicated training panels, invalid identifiers and dates, nonfinite training values, constant predictors, outcome-only predictors, and separation of training from future outcome values.

The tests use only base R and print `All raw panel preparation checks passed.` when successful. No passing-test claim is made until that command or the GitHub Actions workflow succeeds.

## Roadmap

These names describe the intended interface, not currently callable functions:

| Function or component | Intended responsibility | Status |
| --- | --- | --- |
| `panel.dataprep()` | Validate and retain the raw selected panel. | Draft written; R execution pending. |
| `panel.split()` | Construct chronological training/assessment assignments. | Planned. |
| Predictor-block construction and calibration | Build fitting inputs under explicit, cutoff-aware transformation and predictor-weight rules. | Planned; design choices remain open. |
| `tvsc.fit()` | Fit simplex donor-weight paths with first- or second-difference regularization. | Planned. |
| A `predict()` method | Continue fitted donor weights and construct counterfactual outcomes. | Planned. |
| `tvsc.att()` | Calculate date-specific gaps and average effects over explicit windows. | Planned. |

The next implementation checkpoint is to execute and review the raw-preparation tests before extending the interface. Package structure, release automation and distribution details will follow separately.

## Reporting problems

Please include a small reproducible data example, the function call, the complete error message, the expected behavior and the output of `sessionInfo()`. Do not include confidential or identifying data.

## License

This project is licensed under the MIT License. See `LICENSE` for details.
