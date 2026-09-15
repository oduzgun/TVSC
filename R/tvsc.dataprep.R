#' Prepare raw panel data for temporally regularized synthetic control.
#'
#' Select and validate an explicitly declared pre-treatment panel for one
#' target and at least two donors. Keep measurement values as supplied,
#' order training observations consistently, and retain available post-treatment
#' donor and observed target outcomes in separate tables.
#'
#' @param data A nonempty data frame with unique, nonempty column names.
#'   Predictors may already be transformed; no additional scaling is applied.
#'   Supply the outcome column in the units intended for prediction/reporting.
#' @param unit One character string naming the unit-ID column. IDs must be
#'   nonmissing character or finite plain numeric values. Factor IDs are
#'   matched using their character labels; empty character IDs are rejected.
#' @param time One character string naming a finite, integer-valued numeric
#'   period-index column. Date and date-time objects are not supported.
#' @param outcome One character string naming the outcome column. Required
#'   training values must be plain numeric and finite, even when the outcome
#'   is not included in predictors. Unit, time and outcome columns must differ.
#' @param treated One target ID matching the character/numeric type of unit.
#' @param donors A vector of at least two distinct donor IDs, excluding the
#'   target and matching the character/numeric type of unit. Supplied order
#'   determines donor order in the schema and within each training date.
#' @param intervention_time One integer-valued numeric index identifying the
#'   first treated period, not the last pre-treatment period.
#' @param pre_period An explicit numeric vector of at least three increasing
#'   dates, spaced by period_step and strictly before intervention_time.
#' @param predictors A nonempty character vector of unique predictor-column
#'   names. May include outcome, but not unit or time. Training values must
#'   be plain numeric and finite; encode categorical predictors explicitly.
#' @param period_step One positive integer-valued spacing between successive
#'   pre_period dates. Defaults to 1.
#'
#' @return A list of class tvsc_prepared containing:
#'   \describe{
#'     \item{schema_version}{Object-schema version, currently 1L.}
#'     \item{schema}{Column roles, predictors, target and ordered donor IDs,
#'       intervention_time, pre_period and period_step.}
#'     \item{training}{Selected identifier, outcome and predictor columns,
#'       ordered by date, then target first and donors in supplied order.}
#'     \item{future_donors}{Available post-treatment donor unit, time and
#'       outcome columns, in input row order. May have zero rows.}
#'     \item{future_target}{Available observed post-treatment target unit,
#'       time and outcome columns, in input row order. May have zero rows.}
#'     \item{recipe}{Contemporaneous-block declaration, scaling = NULL and
#'       outcome_in_predictors. Scaling is deferred, not an estimator default.}
#'     \item{diagnostics}{Training and excluded row counts, predictors constant
#'       across all training cells, and future_outcomes_validated = FALSE.}
#'   }
#'
#' @details Each selected unit must occur exactly once at each pre_period date.
#'   Missing or duplicate training keys and nonfinite required training values
#'   cause errors. Constant predictors are retained and flagged, not dropped.
#'   Unit IDs and time indices are checked across all supplied rows; future
#'   outcome values, duplicate future keys and future completeness are not
#'   checked. Future tables require validation before prediction or scoring.
#'
#'   Raw means as supplied, not necessarily never transformed. The function
#'   neither detects nor reverses earlier preprocessing. Keep its provenance
#'   separately; data-dependent preprocessing for chronological validation
#'   must respect each training cutoff even if performed outside this code.
#'
#'   The function does not split, scale, impute, aggregate, construct predictor
#'   matrices or fit donor weights. Only the declared training dates and
#'   selected post-treatment outcome rows are retained; additional pre-period
#'   history is not stored for lag construction. Future target outcomes are
#'   observed treated outcomes, not counterfactuals, and are not required.
#'
#' @examples
#' panel <- expand.grid(
#'   state = c("Target", "DonorA", "DonorB"),
#'   year = 2000:2004,
#'   stringsAsFactors = FALSE
#' )
#' panel$sales <- seq_len(nrow(panel)) * 10
#' panel$income <- 1000 + seq_len(nrow(panel))
#' prepared <- tvsc.dataprep(
#'   data = panel, unit = "state", time = "year", outcome = "sales",
#'   treated = "Target", donors = c("DonorB", "DonorA"),
#'   intervention_time = 2003, pre_period = 2000:2002,
#'   predictors = c("sales", "income")
#' )
#' prepared$training
#' prepared$diagnostics
#'
#' @note Initial implementation using base R only.
tvsc.dataprep <- function(data, unit, time, outcome, treated, donors,
                           intervention_time, pre_period, predictors,
                           period_step = 1) {
  # Local validators reject ambiguous types instead of silently coercing data.
  fail <- function(message) stop(message, call. = FALSE)
  plain_numeric <- function(value) {
    is.numeric(value) && !is.object(value) && is.null(dim(value))
  }
  valid_names <- function(value) {
    is.character(value) && is.null(dim(value)) && length(value) > 0L &&
      !anyNA(value) && all(nzchar(value)) && !anyDuplicated(value)
  }
  normalize_ids <- function(value, argument) {
    if (is.factor(value)) value <- as.character(value)
    if (!(is.character(value) || plain_numeric(value)) ||
        !is.null(dim(value)) || length(value) == 0L || anyNA(value)) {
      fail(paste(argument, "must contain nonmissing character or numeric IDs."))
    }
    if (is.character(value) && any(!nzchar(value))) fail(paste(argument, "contains empty IDs."))
    if (is.numeric(value) && any(!is.finite(value))) fail(paste(argument, "contains nonfinite IDs."))
    unname(value)
  }
  whole_periods <- function(value) {
    plain_numeric(value) && length(value) > 0L &&
      all(is.finite(value)) && all(value == floor(value))
  }

  # Check column roles before selecting any observations.
  if (!is.data.frame(data) || nrow(data) == 0L) fail("data must be a nonempty data frame.")
  if (!valid_names(names(data))) fail("data must have unique, nonempty column names.")
  for (role in list(unit, time, outcome)) {
    if (!valid_names(role) || length(role) != 1L) fail("unit, time and outcome must each name one column.")
  }
  if (anyDuplicated(c(unit, time, outcome))) fail("unit, time and outcome must name distinct columns.")
  if (!valid_names(predictors)) fail("predictors must explicitly name at least one unique predictor column.")
  if (any(predictors %in% c(unit, time))) fail("Unit and time identifiers cannot be matching predictors.")
  required <- unique(c(unit, time, outcome, predictors))
  absent <- setdiff(required, names(data))
  if (length(absent)) fail(paste("Missing columns:", paste(absent, collapse = ", ")))
  data <- as.data.frame(data, optional = TRUE)
  # Match IDs explicitly and preserve the researcher's declared donor order.
  identifiers <- normalize_ids(data[[unit]], "unit column")
  target_id <- normalize_ids(treated, "treated")
  donor_ids <- normalize_ids(donors, "donors")
  if (length(target_id) != 1L) fail("treated must identify exactly one target.")
  if (length(donor_ids) < 2L || anyDuplicated(donor_ids)) fail("donors must identify at least two distinct donors.")
  if (is.numeric(identifiers) != is.numeric(target_id) ||
      is.numeric(identifiers) != is.numeric(donor_ids)) {
    fail("Target and donor IDs must match the unit column's character/numeric type.")
  }
  if (target_id %in% donor_ids) fail("The target cannot be one of its own donors.")
  selected_ids <- c(target_id, donor_ids)
  if (any(!selected_ids %in% identifiers)) fail("Every selected target and donor must occur in data.")
  # Validate the time index and pre-period; do not infer validation splits.
  if (!whole_periods(data[[time]])) fail("time must be a finite integer-valued numeric period index, not calendar dates.")
  if (!whole_periods(intervention_time) || length(intervention_time) != 1L) {
    fail("intervention_time must be one integer-valued first-treated period.")
  }
  if (!whole_periods(period_step) || length(period_step) != 1L || period_step <= 0) {
    fail("period_step must be one positive integer-valued period step.")
  }
  if (!whole_periods(pre_period) || length(pre_period) < 3L ||
      any(diff(pre_period) != period_step)) {
    fail("pre_period must contain at least three increasing dates at period_step spacing.")
  }
  if (any(pre_period >= intervention_time)) fail("Every pre_period date must precede intervention_time.")

  # Restrict training to the declared units and pre-treatment dates.
  selected <- identifiers %in% selected_ids
  training_rows <- which(selected & data[[time]] %in% pre_period)
  training <- data[training_rows, required, drop = FALSE]
  training_ids <- identifiers[training_rows]
  unit_index <- match(training_ids, selected_ids)
  date_index <- match(training[[time]], pre_period)
  # Require one row per unit-date key without averaging or filling missing rows.
  keys <- paste(unit_index, date_index, sep = ":")
  if (anyDuplicated(keys)) {
    position <- which(duplicated(keys))[1L]
    fail(paste("Duplicate training key: unit", training_ids[position],
               "period", training[[time]][position]))
  }
  expected_units <- rep(seq_along(selected_ids), times = length(pre_period))
  expected_dates <- rep(seq_along(pre_period), each = length(selected_ids))
  row_order <- match(paste(expected_units, expected_dates, sep = ":"), keys)
  if (anyNA(row_order)) {
    position <- which(is.na(row_order))[1L]
    fail(paste("Missing training key: unit", selected_ids[expected_units[position]],
               "period", pre_period[expected_dates[position]]))
  }
  # Order by date and then target/donors; retain unscaled measurement values.
  training <- training[row_order, , drop = FALSE]
  rownames(training) <- NULL
  for (column in unique(c(outcome, predictors))) {
    values <- training[[column]]
    if (!plain_numeric(values)) fail(paste("Training column", column, "must be numeric; encode categories explicitly."))
    if (any(!is.finite(values))) {
      position <- which(!is.finite(values))[1L]
      fail(paste("Nonfinite training value in", column, "at unit",
                 training[[unit]][position], "period", training[[time]][position]))
    }
  }

  # Store future outcomes separately; their values and panel keys remain unchecked.
  future_rows <- selected & data[[time]] >= intervention_time
  donor_rows <- future_rows & identifiers %in% donor_ids
  target_rows <- future_rows & identifiers %in% target_id
  future_donors <- data[donor_rows, c(unit, time, outcome), drop = FALSE]
  future_target <- data[target_rows, c(unit, time, outcome), drop = FALSE]
  rownames(future_donors) <- NULL
  rownames(future_target) <- NULL
  # Flag constant predictors without making a scaling or variable-selection decision.
  constant_predictors <- predictors[vapply(predictors, function(column) {
    length(unique(training[[column]])) == 1L
  }, logical(1))]
  # Return the raw panel and metadata, not fitted weights or predictor matrices.
  structure(list(
    schema_version = 1L,
    schema = list(unit = unit, time = time, outcome = outcome,
                  predictors = predictors, treated = target_id, donors = donor_ids,
                  intervention_time = intervention_time, pre_period = pre_period,
                  period_step = period_step),
    training = training,
    future_donors = future_donors,
    future_target = future_target,
    recipe = list(blocks = "contemporaneous", scaling = NULL,
                  outcome_in_predictors = outcome %in% predictors),
    diagnostics = list(training_rows = nrow(training),
                       excluded_rows = nrow(data) - nrow(training) - sum(future_rows),
                       constant_predictors = constant_predictors,
                       future_outcomes_validated = FALSE)
  ), class = "tvsc_prepared")
}
