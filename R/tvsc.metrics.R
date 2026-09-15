#' Construct predetermined predictor metrics for a training prefix.
#'
#' @param blocks A schema-version-1 tvsc_scaled_blocks object. Use
#'   tvsc.scale(raw_blocks, "none") to explicitly choose unscaled predictors.
#' @param predictor_weights Exactly "uniform" (default), a named real numeric
#'   vector covering the predictors, or a real numeric predictor-by-date matrix
#'   with row and column names covering the training prefix exactly.
#' @param tolerance Absolute column-sum tolerance, a finite numeric scalar in
#'   [0, 1). Default 1e-8. Accepted supplied values are never renormalized.
#' @param provenance NULL or a nonempty character scalar describing the source
#'   and allowed training information. This declaration is not verified.
#'
#' @return A schema-version-1 tvsc_metrics object containing a canonical named
#'   predictor-by-date weights matrix, input mode, schema, predictors, periods,
#'   cutoff, dimensions, scaling record, provenance and diagnostics. Diagonal
#'   V matrices can be constructed from columns of weights when fitting.
#'
#' @details Uniform weights equal 1/K. A supplied vector is aligned by predictor
#'   name and repeated across dates; a supplied matrix is aligned by predictor
#'   and date names. Entries must be finite and nonnegative, with positive
#'   column sums within tolerance of one. Zeros are allowed. No clipping,
#'   normalization, calibration, donor fitting or outcome scoring is performed.
#'   Modes describe input shapes, not whether supplied values were estimated.
#'   Even a constant-column matrix has mode supplied_time_varying; diagnostics
#'   separately report exact constancy of its columns.
#'
#'   This constructor checks alignment metadata, not the numerical contents
#'   of X1, X0, Y1 or Y0. It does not recompute scaling or certify that metadata
#'   describe those contents. The eventual fitter must validate its numerical
#'   inputs and compatibility with these metrics. A schema is not a checksum.
#'   No numerical predictor or outcome arrays are copied into the output.
#'
#'   Empirical metrics must be calibrated separately within each training
#'   prefix, before fitting. Matching date labels cannot detect leakage from
#'   a full-pre-period calibration. Extra matrix dates are rejected rather
#'   than silently sliced. Supplied provenance is always marked unverified.
#'   Reuse a prefix's scaled blocks and metrics across orders and penalties.
#'   Effective raw-unit coefficients v[k,t]/scale[k,t]^2 are not normalized
#'   here: weights refer to the chosen scaled predictor coordinates.
#'
#' @examples
#' panel <- expand.grid(state = c("Target", "DonorA", "DonorB"),
#'                      year = 2000:2005, stringsAsFactors = FALSE)
#' panel$sales <- seq_len(nrow(panel)) * 10
#' prepared <- tvsc.dataprep(panel, "state", "year", "sales", "Target",
#'                            c("DonorB", "DonorA"), 2005, 2000:2004, "sales")
#' scaled <- tvsc.scale(tvsc.blocks(prepared, 2003))
#' metrics <- tvsc.metrics(scaled)
#' stopifnot(all(metrics$weights == 1))
#'
#' @note Initial base-R implementation. Empirical calibration is unimplemented.
tvsc.metrics <- function(blocks, predictor_weights = "uniform",
                          tolerance = 1e-8, provenance = NULL) {
  fail <- function(message) stop(message, call. = FALSE)
  has_fields <- function(value, required) {
    is.list(value) && !is.data.frame(value) && !is.null(names(value)) &&
      !anyNA(names(value)) && !anyDuplicated(names(value)) &&
      all(nzchar(names(value))) && all(required %in% names(value))
  }
  plain_numeric <- function(value) {
    typeof(value) %in% c("integer", "double") && !is.object(value) &&
      all(is.finite(value))
  }
  labels_valid <- function(value) {
    is.character(value) && !is.object(value) && is.null(dim(value)) &&
      length(value) > 0L && !anyNA(value) && all(nzchar(value)) &&
      !anyDuplicated(value)
  }
  ids_valid <- function(value) {
    is.null(dim(value)) && !is.object(value) && length(value) > 0L &&
      ((is.character(value) && !anyNA(value) && all(nzchar(value))) ||
         plain_numeric(value)) && !anyDuplicated(value)
  }
  whole_vector <- function(value) {
    plain_numeric(value) && is.null(dim(value)) && length(value) > 0L &&
      all(value == trunc(value))
  }

  # Validate the metadata boundary without reading predictor or outcome values.
  if (!inherits(blocks, "tvsc_scaled_blocks") || !inherits(blocks, "tvsc_blocks") ||
      !has_fields(blocks, c("schema_version", "schema", "periods", "cutoff",
                            "dimensions", "recipe", "X1", "X0", "Y1", "Y0")) ||
      !identical(blocks$schema_version, 1L)) {
    fail("blocks must be a schema-version-1 tvsc_scaled_blocks object.")
  }
  schema <- blocks$schema
  if (!has_fields(schema, c("unit", "time", "outcome", "predictors", "treated",
                           "donors", "intervention_time", "pre_period", "period_step")) ||
      !labels_valid(schema$predictors) || !ids_valid(schema$treated) ||
      length(schema$treated) != 1L || !ids_valid(schema$donors) ||
      length(schema$donors) < 2L ||
      is.character(schema$treated) != is.character(schema$donors) ||
      schema$treated %in% schema$donors) {
    fail("blocks has an incomplete or ambiguous schema.")
  }
  roles <- list(schema$unit, schema$time, schema$outcome)
  if (!all(vapply(roles, function(value) labels_valid(value) && length(value) == 1L,
                  logical(1))) || anyDuplicated(unlist(roles, use.names = FALSE)) ||
      any(schema$predictors %in% c(schema$unit, schema$time))) {
    fail("blocks has invalid column roles.")
  }
  periods <- blocks$periods
  if (!whole_vector(periods) || length(periods) < 3L ||
      !whole_vector(schema$period_step) || length(schema$period_step) != 1L ||
      schema$period_step <= 0 || !all(diff(periods) == schema$period_step) ||
      !whole_vector(schema$intervention_time) || length(schema$intervention_time) != 1L ||
      !all(periods < schema$intervention_time) ||
      !identical(periods, schema$pre_period) || !whole_vector(blocks$cutoff) ||
      length(blocks$cutoff) != 1L || blocks$cutoff != tail(periods, 1L)) {
    fail("blocks must contain aligned regular pre-treatment dates and a matching cutoff.")
  }
  predictor_names <- unname(schema$predictors)
  date_names <- as.character(periods)
  predictor_count <- length(predictor_names)
  period_count <- length(periods)
  dimensions <- blocks$dimensions
  expected_counts <- c(predictors = predictor_count, donors = length(schema$donors),
                       periods = period_count)
  if (!has_fields(dimensions, names(expected_counts)) ||
      !all(vapply(names(expected_counts), function(field) {
        value <- dimensions[[field]]
        whole_vector(value) && length(value) == 1L && value == expected_counts[[field]]
      }, logical(1)))) {
    fail("blocks dimensions must agree with schema and dates.")
  }
  if (!has_fields(blocks$recipe, c("blocks", "scaling", "outcome_in_predictors")) ||
      !identical(blocks$recipe$blocks, "contemporaneous") ||
      !identical(blocks$recipe$outcome_in_predictors,
                 schema$outcome %in% schema$predictors)) {
    fail("blocks must retain a contemporaneous recipe with explicit scaling.")
  }
  scaling <- blocks$recipe$scaling
  formulas <- c(sd = "rms_within_date_sample_sd_no_center",
                zscore = "within_date_mean_and_sample_sd",
                none = "identity_no_arithmetic", custom = "named_fixed_divisors_no_center")
  matrix_aligned <- function(value) {
    plain_numeric(value) && is.matrix(value) &&
      identical(dim(value), c(predictor_count, period_count)) &&
      identical(rownames(value), predictor_names) && identical(colnames(value), date_names)
  }
  if (!has_fields(scaling, c("schema_version", "method", "formula_id", "reference",
                            "periods", "cutoff", "centers", "divisors")) ||
      !identical(scaling$schema_version, 1L) || !labels_valid(scaling$method) ||
      length(scaling$method) != 1L || !scaling$method %in% names(formulas) ||
      !identical(scaling$formula_id, unname(formulas[scaling$method])) ||
      !identical(scaling$periods, periods) || !identical(scaling$cutoff, blocks$cutoff) ||
      !has_fields(scaling$reference, c("mode", "treated", "donors")) ||
      !identical(scaling$reference$treated, schema$treated) ||
      !identical(scaling$reference$donors, schema$donors) ||
      !identical(scaling$reference$mode,
                 if (scaling$method %in% c("sd", "zscore")) "all_units" else "not_estimated") ||
      !matrix_aligned(scaling$centers) || !matrix_aligned(scaling$divisors) ||
      any(scaling$divisors <= 0)) {
    fail("blocks has invalid or misaligned scaling metadata.")
  }
  if (scaling$method != "zscore" &&
      (any(scaling$centers != 0) ||
       !all(scaling$divisors == scaling$divisors[, 1L]) ||
       (scaling$method == "none" && any(scaling$divisors != 1)))) {
    fail("blocks scaling metadata contradicts its method.")
  }
  if (!plain_numeric(tolerance) || !is.null(dim(tolerance)) || length(tolerance) != 1L ||
      tolerance < 0 || tolerance >= 1) {
    fail("tolerance must be a finite numeric scalar in [0, 1).")
  }
  if (!is.null(provenance) && (!labels_valid(provenance) || length(provenance) != 1L)) {
    fail("provenance must be NULL or a nonempty character scalar.")
  }

  # Align supplied weights by exact labels; never select or normalize them.
  if (identical(predictor_weights, "uniform")) {
    mode <- "uniform"
    weights <- matrix(1 / predictor_count, predictor_count, period_count,
                      dimnames = list(predictor_names, date_names))
  } else if (plain_numeric(predictor_weights) && is.null(dim(predictor_weights))) {
    if (!labels_valid(names(predictor_weights)) ||
        !setequal(names(predictor_weights), predictor_names)) {
      fail("Supplied vector names must cover predictors exactly once.")
    }
    mode <- "supplied_constant"
    weights <- matrix(rep(predictor_weights[predictor_names], period_count),
                      predictor_count, period_count, dimnames = list(predictor_names, date_names))
  } else if (plain_numeric(predictor_weights) && is.matrix(predictor_weights)) {
    if (!labels_valid(rownames(predictor_weights)) ||
        !labels_valid(colnames(predictor_weights)) ||
        !setequal(rownames(predictor_weights), predictor_names) ||
        !setequal(colnames(predictor_weights), date_names)) {
      fail("Supplied matrix names must cover predictors and training dates exactly once.")
    }
    mode <- "supplied_time_varying"
    weights <- predictor_weights[predictor_names, date_names, drop = FALSE]
  } else {
    fail("predictor_weights must be exactly uniform, a finite real named vector or a finite real named matrix.")
  }
  column_sums <- colSums(weights)
  if (any(weights < 0) || any(!is.finite(column_sums)) || any(column_sums <= 0) ||
      any(abs(column_sums - 1) > tolerance)) {
    fail("Weights must be nonnegative with positive column sums within tolerance of one.")
  }

  # Retain the declared coordinate system and unverified provenance, not data arrays.
  result <- list(schema_version = 1L, mode = mode, schema = schema,
                 predictors = predictor_names, periods = periods, cutoff = blocks$cutoff,
                 dimensions = list(predictors = predictor_count, periods = period_count),
                 weights = weights, scaling = scaling,
                 provenance = list(description = provenance, verified = FALSE),
                 diagnostics = list(column_sums = column_sums,
                   maximum_unit_sum_error = max(abs(column_sums - 1)),
                   zero_weight_predictors = predictor_names[rowSums(weights) == 0],
                   constant_over_time = all(weights == weights[, 1L]),
                   normalization_applied = FALSE, tolerance = tolerance))
  class(result) <- "tvsc_metrics"
  result
}
