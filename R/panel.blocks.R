#' Assemble as-supplied predictor blocks for a specified training cutoff.
#'
#' Arrange contemporaneous predictors from a prepared panel into target
#' vectors and donor matrices for each date in the selected training prefix.
#' Keep predictor order, donor order and measurement values unchanged.
#' Scaling is a separate operation; no default economic scale is chosen here.
#'
#' @param prepared A schema-version-1 tvsc_prepared object from panel.dataprep(),
#'   with as-supplied training values and its contemporaneous recipe recording
#'   no package scaling. Predictors may have been transformed before preparation.
#' @param cutoff One finite integer-valued numeric date in schema$pre_period,
#'   leaving at least three training dates. Required. This is an actual period
#'   index, not a row number or a position in pre_period. For full-pre-period
#'   assembly, supply the last date in prepared$schema$pre_period.
#'
#' @return A list of class tvsc_blocks containing:
#'   \describe{
#'     \item{schema_version}{Block-object schema version, currently 1L.}
#'     \item{schema}{Validated input roles and donor order, with pre_period
#'       restricted to the selected prefix. No later observations are stored.}
#'     \item{periods, cutoff}{Ordered training dates and the final included date.}
#'     \item{dimensions}{Counts named predictors, donors and periods.}
#'     \item{X1}{Date-named list of named numeric predictor vectors, one per
#'       training date. Each vector has length K in declared predictor order.}
#'     \item{X0}{Date-named list of K-by-J numeric matrices. Rows are predictors;
#'       columns are donors in declared order. K = 1 retains matrix dimensions.}
#'     \item{Y1}{Date-named numeric vector of target training outcomes.}
#'     \item{Y0}{Numeric training-date-by-donor outcome matrix, with date row
#'       names and ordered donor column names. Outcomes remain in supplied units.}
#'     \item{recipe}{Contemporaneous-block declaration, scaling = NULL and
#'       outcome_in_predictors. NULL records no package-applied scaling; it
#'       does not describe or verify preprocessing performed before preparation.}
#'     \item{diagnostics}{training_rows and constant_predictors, recomputed
#'       using only the selected prefix. Constant predictors are retained.}
#'   }
#'
#' @details Requires panel.dataprep() to be loaded. Its validation and ordering
#'   logic is reused on the selected prefix, rather than trusting stored row
#'   order or duplicating the raw-panel validation rules.
#'
#'   Schema metadata and the stored training time column are checked before
#'   selection. Required unit-date keys and numeric outcome/predictor values
#'   are then checked within the prefix only. Later pre-treatment values and
#'   both future tables are excluded from the returned object and from numeric
#'   validation. The outcome is returned even when omitted from predictors;
#'   it is never silently added to the matching blocks.
#'
#'   The function does not modify prepared, create lags or summaries, impute,
#'   center, scale, choose predictor weights, fit donor weights, or calculate
#'   validation scores. The future fitter must explicitly apply or declare its
#'   scaling rule, which may be no additional scaling. Raw blocks are not
#'   automatically complete fitting inputs.
#'
#'   Raw means as supplied, not necessarily in original economic units. This
#'   function neither detects nor reverses earlier transformations. To match
#'   on a scaled outcome while retaining reporting units, supply separate
#'   columns such as outcome = "sales" and predictors = "sales_scaled".
#'   Earlier learned transformations must respect each validation cutoff;
#'   selecting a prefix here cannot undo leakage from full-period scaling.
#'
#' @examples
#' panel <- expand.grid(
#'   state = c("Target", "DonorA", "DonorB"), year = 2000:2007,
#'   stringsAsFactors = FALSE
#' )
#' panel$sales <- seq_len(nrow(panel)) * 10
#' panel$income <- 1000 + seq_len(nrow(panel))
#' prepared <- panel.dataprep(
#'   panel, "state", "year", "sales", "Target", c("DonorB", "DonorA"),
#'   2006, 2000:2005, c("income", "sales")
#' )
#' blocks <- panel.blocks(prepared, cutoff = 2003)
#' blocks$X1[["2000"]]
#' blocks$X0[["2000"]]
#' blocks$Y0
#'
#' @note Initial base-R implementation.
panel.blocks <- function(prepared, cutoff) {
  # Require the raw object contract before using schema fields for selection.
  fail <- function(message) stop(message, call. = FALSE)
  whole_numeric <- function(value) {
    is.numeric(value) && !is.object(value) && is.null(dim(value)) &&
      length(value) > 0L && all(is.finite(value)) &&
      all(value == floor(value))
  }
  if (!inherits(prepared, "tvsc_prepared") || !is.list(prepared) ||
      !identical(prepared$schema_version, 1L) || !is.list(prepared$schema) ||
      !is.data.frame(prepared$training)) {
    fail("prepared must be a schema-version-1 tvsc_prepared object from panel.dataprep().")
  }
  schema <- prepared$schema
  required <- c("unit", "time", "outcome", "predictors", "treated", "donors",
                "intervention_time", "pre_period", "period_step")
  if (!all(required %in% names(schema)) || anyDuplicated(names(schema))) {
    fail("prepared has an incomplete or ambiguous schema; recreate it with panel.dataprep().")
  }
  if (!is.list(prepared$recipe) ||
      !all(c("blocks", "scaling") %in% names(prepared$recipe)) ||
      !identical(prepared$recipe$blocks, "contemporaneous") ||
      !is.null(prepared$recipe$scaling)) {
    fail("panel.blocks requires a contemporaneous raw recipe with scaling = NULL.")
  }

  # The cutoff is a declared date, never an implicitly chosen row or position.
  if (missing(cutoff)) fail("cutoff must be supplied explicitly as a pre-period date.")
  if (!whole_numeric(cutoff) || length(cutoff) != 1L) {
    fail("cutoff must be one finite integer-valued numeric pre-period date.")
  }
  all_periods <- schema$pre_period
  spacing <- schema$period_step
  intervention <- schema$intervention_time
  if (!whole_numeric(spacing) || length(spacing) != 1L || spacing <= 0 ||
      !whole_numeric(intervention) || length(intervention) != 1L ||
      !whole_numeric(all_periods) || length(all_periods) < 3L ||
      any(diff(all_periods) != spacing) || any(all_periods >= intervention)) {
    fail("prepared schema must contain regular pre-treatment dates and valid period metadata.")
  }
  cutoff_index <- match(cutoff, all_periods)
  if (is.na(cutoff_index)) fail("cutoff must match a date in prepared$schema$pre_period.")
  if (cutoff_index < 3L) fail("cutoff must leave at least three training dates.")
  periods <- all_periods[seq_len(cutoff_index)]
  if (!is.character(schema$time) || length(schema$time) != 1L ||
      is.na(schema$time) || !schema$time %in% names(prepared$training) ||
      anyDuplicated(names(prepared$training))) {
    fail("prepared training data must have an unambiguous schema time column.")
  }
  times <- prepared$training[[schema$time]]
  if (!whole_numeric(times)) fail("prepared training time must be a finite integer-valued numeric index.")

  # Reuse the existing validator on the prefix only; do not pass future tables.
  if (!exists("panel.dataprep", mode = "function")) {
    fail("Load panel.dataprep() before calling panel.blocks().")
  }
  prefix <- panel.dataprep(
    data = prepared$training[times %in% periods, , drop = FALSE],
    unit = schema$unit, time = schema$time, outcome = schema$outcome,
    treated = schema$treated, donors = schema$donors,
    intervention_time = intervention, pre_period = periods,
    predictors = schema$predictors, period_step = spacing
  )
  training <- prefix$training
  predictor_names <- unname(prefix$schema$predictors)
  donor_names <- as.character(prefix$schema$donors)
  date_names <- as.character(periods)
  donor_count <- length(donor_names)
  unit_count <- donor_count + 1L
  target_rows <- (seq_along(periods) - 1L) * unit_count + 1L

  # Build date-specific blocks from the validator's target-first ordering.
  X1 <- lapply(target_rows, function(target_row) {
    vapply(predictor_names, function(predictor) training[[predictor]][target_row], numeric(1))
  })
  X0 <- lapply(target_rows, function(target_row) {
    donor_rows <- target_row + seq_len(donor_count)
    block <- t(as.matrix(training[donor_rows, predictor_names, drop = FALSE]))
    storage.mode(block) <- "double"
    dimnames(block) <- list(predictor_names, donor_names)
    block
  })
  names(X1) <- date_names
  names(X0) <- date_names

  # Retain outcomes independently of whether they were chosen as predictors.
  Y1 <- as.numeric(training[[schema$outcome]][target_rows])
  names(Y1) <- date_names
  donor_rows <- unlist(lapply(target_rows, function(target_row) {
    target_row + seq_len(donor_count)
  }), use.names = FALSE)
  Y0 <- matrix(as.numeric(training[[schema$outcome]][donor_rows]),
               nrow = length(periods), ncol = donor_count, byrow = TRUE,
               dimnames = list(date_names, donor_names))

  # Return only prefix data and prefix diagnostics, with scaling still deferred.
  structure(list(
    schema_version = 1L, schema = prefix$schema,
    periods = periods, cutoff = periods[length(periods)],
    dimensions = list(predictors = length(predictor_names),
                      donors = donor_count, periods = length(periods)),
    X1 = X1, X0 = X0, Y1 = Y1, Y0 = Y0,
    recipe = prefix$recipe,
    diagnostics = list(training_rows = nrow(training),
                       constant_predictors = prefix$diagnostics$constant_predictors)
  ), class = "tvsc_blocks")
}
