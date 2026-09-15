#' Construct expanding-window chronological splits for a prepared panel.
#'
#' Assign pre-treatment dates to training and validation windows, keeping all
#' selected units together at each date. Return date assignments and metadata
#' only; do not change the prepared panel, transform values or fit a model.
#'
#' @param prepared An unmodified tvsc_prepared object with schema_version 2L,
#'   returned by tvsc.dataprep(). Checks schema and raw-block structure,
#'   but does not inspect outcome or predictor values or revalidate the panel.
#'   Rebuild older preparation objects from the original panel.
#' @param initial Integer-valued number of dates in the first training window,
#'   at least 3. Required; no automatic window size is chosen.
#' @param horizon Positive integer-valued number of dates in each validation
#'   window immediately following its training cutoff. Required.
#' @param step Positive integer-valued number of dates by which the cutoff
#'   advances between folds. Defaults to 1. Counts dates, not panel rows or
#'   numeric calendar-index units.
#'
#' @return A list of class tvsc_splits containing:
#'   \describe{
#'     \item{schema_version}{Split-object schema version, currently 1L.}
#'     \item{schema}{A copy of the prepared panel's schema, including unit and
#'       predictor roles, donor order, pre-period and period spacing. This is
#'       compatibility metadata, not a checksum of the underlying data.}
#'     \item{specification}{Window type, initial, horizon, step and the
#'       full-window-only policy.}
#'     \item{folds}{An ordered list. Each fold has an id, train_index and
#'       validation_index (positions in schema$pre_period), train_period,
#'       cutoff, validation_period and horizons (1 through horizon). Each
#'       assigned date applies to every selected unit, not individual rows.}
#'     \item{diagnostics}{fold_count; assessment_counts, named by pre-period
#'       date; and unassessed_period, including initial training-only dates
#'       and any dates omitted by the stride or full-window policy.}
#'   }
#'
#' @details Training starts at the first declared pre-period date and expands.
#'   Cutoff positions are initial, initial + step, and so on, provided the
#'   complete subsequent validation window fits within the pre-period. No
#'   shortened final window or extra off-stride cutoff is added. An error is
#'   raised if no full fold can be formed.
#'
#'   Training and validation do not overlap within a fold. Across folds,
#'   validation dates may recur and earlier validation dates may enter later
#'   training prefixes. No independent test set is reserved automatically.
#'   Later scoring must account for the chosen treatment of repeated dates.
#'   All data-dependent transformations and calibration must be fitted within
#'   each training prefix by downstream code, not by this split constructor.
#'
#' @examples
#' panel <- expand.grid(
#'   state = c("Target", "DonorA", "DonorB"),
#'   year = 2000:2009,
#'   stringsAsFactors = FALSE
#' )
#' panel$sales <- seq_len(nrow(panel)) * 10
#' prepared <- tvsc.dataprep(
#'   panel, "state", "year", "sales", "Target", c("DonorB", "DonorA"),
#'   2008, 2000:2007, "sales"
#' )
#' splits <- tvsc.split(prepared, initial = 4, horizon = 2)
#' splits$folds[[1]]$train_period
#' splits$folds[[1]]$validation_period
#' splits$diagnostics
#'
#' @note Revised version-2 preparation compatibility remains unexecuted in R
#'   here; the split output remains version 1. Earlier input-contract evidence:
#'   a user-supplied R console transcript
#'   dated September 15, 2026 shows the split checks passing, including the
#'   eight-date example. The R version and sourced-file revisions were not
#'   supplied; an independent local rerun remains pending.
tvsc.split <- function(prepared, initial, horizon, step = 1) {
  # Check the object contract without reading or transforming panel values.
  fail <- function(message) stop(message, call. = FALSE)
  whole_numeric <- function(value) {
    is.numeric(value) && !is.object(value) && is.null(dim(value)) &&
      length(value) > 0L && all(is.finite(value)) &&
      all(value == floor(value))
  }
  has_fields <- function(value, fields) {
    is.list(value) && !anyDuplicated(names(value)) &&
      all(fields %in% names(value))
  }
  valid_labels <- function(value) {
    is.character(value) && is.null(dim(value)) && length(value) > 0L &&
      !anyNA(value) && all(nzchar(value)) && !anyDuplicated(value)
  }
  valid_ids <- function(value) {
    valid_labels(value) ||
      (is.numeric(value) && !is.object(value) && is.null(dim(value)) &&
       length(value) > 0L && all(is.finite(value)) && !anyDuplicated(value))
  }
  if (!inherits(prepared, "tvsc_prepared") || !is.list(prepared) ||
      !identical(prepared$schema_version, 2L) || !is.list(prepared$schema)) {
    fail("prepared must be a schema-version-2 tvsc_prepared object; rerun tvsc.dataprep() from the original panel.")
  }
  schema <- prepared$schema
  required <- c("unit", "time", "outcome", "predictors", "treated", "donors",
                "intervention_time", "pre_period", "period_step")
  if (!all(required %in% names(schema)) || anyDuplicated(names(schema))) {
    fail("prepared has an incomplete or ambiguous schema; recreate it with tvsc.dataprep().")
  }
  periods <- schema$pre_period
  spacing <- schema$period_step
  intervention <- schema$intervention_time
  if (!whole_numeric(spacing) || length(spacing) != 1L || spacing <= 0 ||
      !whole_numeric(intervention) || length(intervention) != 1L ||
      !whole_numeric(periods) || length(periods) < 3L ||
      any(diff(periods) != spacing) || any(periods >= intervention)) {
    fail("prepared schema must contain regular pre-treatment dates and a valid period_step and intervention_time.")
  }
  roles <- list(schema$unit, schema$time, schema$outcome)
  if (!all(vapply(roles, function(role) {
      valid_labels(role) && length(role) == 1L
    }, logical(1))) || anyDuplicated(unlist(roles, use.names = FALSE)) ||
      !valid_labels(schema$predictors) ||
      any(schema$predictors %in% c(schema$unit, schema$time)) ||
      !valid_ids(schema$treated) || length(schema$treated) != 1L ||
      !valid_ids(schema$donors) || length(schema$donors) < 2L ||
      is.numeric(schema$treated) != is.numeric(schema$donors) ||
      schema$treated %in% schema$donors) {
    fail("prepared has invalid column roles or target/donor metadata; rerun tvsc.dataprep().")
  }
  predictor_names <- unname(schema$predictors)
  donor_names <- as.character(schema$donors)
  date_names <- as.character(periods)
  predictor_count <- length(predictor_names)
  donor_count <- length(donor_names)
  period_count <- length(periods)
  if (!has_fields(prepared, c("schema_version", "schema", "training", "blocks",
      "future_donors", "future_target", "recipe", "diagnostics")) ||
      !is.data.frame(prepared$training) ||
      !identical(names(prepared$training), unname(unique(c(schema$unit, schema$time,
                                                  schema$outcome, schema$predictors)))) ||
      nrow(prepared$training) != period_count * (donor_count + 1L) ||
      !has_fields(prepared$recipe, c("blocks", "scaling", "outcome_in_predictors")) ||
      !identical(prepared$recipe$blocks, "contemporaneous") ||
      !is.null(prepared$recipe$scaling) ||
      !identical(prepared$recipe$outcome_in_predictors,
                 schema$outcome %in% schema$predictors)) {
    fail("prepared has an invalid version-2 training structure or raw recipe; rerun tvsc.dataprep() from the original panel.")
  }
  blocks <- prepared$blocks
  if (!inherits(blocks, "tvsc_blocks") || inherits(blocks, "tvsc_scaled_blocks") ||
      !has_fields(blocks, c("schema_version", "schema", "periods", "cutoff",
                           "dimensions", "X1", "X0", "Y1", "Y0", "recipe", "diagnostics")) ||
      !identical(blocks$schema_version, 2L) || !identical(blocks$schema, schema) ||
      !identical(blocks$periods, periods) ||
      !identical(blocks$cutoff, periods[period_count]) ||
      !identical(blocks$recipe, prepared$recipe) ||
      !identical(blocks$dimensions, list(predictors = predictor_count,
                                         donors = donor_count, periods = period_count)) ||
      !is.list(blocks$X1) || !identical(names(blocks$X1), date_names) ||
      !is.list(blocks$X0) || !identical(names(blocks$X0), date_names) ||
      !is.numeric(blocks$Y1) || !is.null(dim(blocks$Y1)) ||
      !identical(names(blocks$Y1), date_names) ||
      !is.matrix(blocks$Y0) || !is.numeric(blocks$Y0) ||
      !identical(dim(blocks$Y0), c(period_count, donor_count)) ||
      !identical(dimnames(blocks$Y0), list(date_names, donor_names))) {
    fail("prepared has incompatible raw blocks; rerun tvsc.dataprep() from the original panel.")
  }
  for (date_index in seq_along(periods)) {
    target <- blocks$X1[[date_index]]
    donors <- blocks$X0[[date_index]]
    if (!is.numeric(target) || !is.null(dim(target)) ||
        !identical(names(target), predictor_names) ||
        !is.matrix(donors) || !is.numeric(donors) ||
        !identical(dim(donors), c(predictor_count, donor_count)) ||
        !identical(dimnames(donors), list(predictor_names, donor_names))) {
      fail("prepared has incompatible raw block dimensions or ordering; rerun tvsc.dataprep().")
    }
  }

  # Window sizes count dates; require an explicit initial size and horizon.
  if (missing(initial) || missing(horizon)) fail("initial and horizon must be supplied explicitly.")
  arguments <- list(initial = initial, horizon = horizon, step = step)
  for (argument in names(arguments)) {
    value <- arguments[[argument]]
    if (!whole_numeric(value) || length(value) != 1L || value <= 0) {
      fail(paste(argument, "must be one positive integer-valued number of dates."))
    }
  }
  if (initial < 3) fail("initial must contain at least three training dates.")
  period_count <- length(periods)
  if (initial > period_count || horizon > period_count - initial) {
    fail("No full fold fits: initial + horizon must not exceed the number of pre-period dates.")
  }

  # Use only complete windows on the requested stride, never a shifted last fold.
  cutoffs <- seq(from = initial, to = period_count - horizon, by = step)
  folds <- lapply(seq_along(cutoffs), function(fold_id) {
    cutoff_index <- cutoffs[fold_id]
    train_index <- seq_len(cutoff_index)
    validation_index <- as.integer(cutoff_index + seq_len(horizon))
    list(id = fold_id,
         train_index = train_index, validation_index = validation_index,
         train_period = periods[train_index], cutoff = periods[cutoff_index],
         validation_period = periods[validation_index],
         horizons = seq_len(horizon))
  })

  # Expose repeated assessment dates and omissions without choosing scoring weights.
  assessed_indices <- unlist(lapply(folds, function(fold) fold$validation_index),
                            use.names = FALSE)
  counts <- tabulate(assessed_indices, nbins = period_count)
  names(counts) <- as.character(periods)
  structure(list(
    schema_version = 1L,
    schema = schema,
    specification = list(window = "expanding", initial = initial,
                         horizon = horizon, step = step, full_windows_only = TRUE),
    folds = folds,
    diagnostics = list(fold_count = length(folds), assessment_counts = counts,
                       unassessed_period = periods[counts == 0L])
  ), class = "tvsc_splits")
}