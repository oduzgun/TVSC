#' Construct expanding-window chronological splits for a prepared panel.
#'
#' Assign pre-treatment dates to training and validation windows, keeping all
#' selected units together at each date. Return date assignments and metadata
#' only; do not change the prepared panel, transform values or fit a model.
#'
#' @param prepared An unmodified tvsc_prepared object with schema_version 1L,
#'   returned by panel.dataprep(). Only schema metadata are used, not outcome
#'   or predictor values. This function does not revalidate the raw panel.
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
#' prepared <- panel.dataprep(
#'   panel, "state", "year", "sales", "Target", c("DonorB", "DonorA"),
#'   2008, 2000:2007, "sales"
#' )
#' splits <- panel.split(prepared, initial = 4, horizon = 2)
#' splits$folds[[1]]$train_period
#' splits$folds[[1]]$validation_period
#' splits$diagnostics
#'
#' @note Initial base-R implementation.
panel.split <- function(prepared, initial, horizon, step = 1) {
  # Check the object contract without reading or transforming panel values.
  fail <- function(message) stop(message, call. = FALSE)
  whole_numeric <- function(value) {
    is.numeric(value) && !is.object(value) && is.null(dim(value)) &&
      length(value) > 0L && all(is.finite(value)) &&
      all(value == floor(value))
  }
  if (!inherits(prepared, "tvsc_prepared") || !is.list(prepared) ||
      !identical(prepared$schema_version, 1L) || !is.list(prepared$schema)) {
    fail("prepared must be a schema-version-1 tvsc_prepared object from panel.dataprep().")
  }
  schema <- prepared$schema
  required <- c("unit", "time", "outcome", "predictors", "treated", "donors",
                "intervention_time", "pre_period", "period_step")
  if (!all(required %in% names(schema)) || anyDuplicated(names(schema))) {
    fail("prepared has an incomplete or ambiguous schema; recreate it with panel.dataprep().")
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
