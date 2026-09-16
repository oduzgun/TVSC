#' Scale training predictor blocks using an explicit panel-specific rule.
#'
#' Transform matching predictors only, preserving reporting outcomes and the
#' input object. Short method names have panel-specific meanings: sd uses one
#' RMS within-date SD per predictor; zscore uses each date's own mean and SD.
#'
#' @param blocks A raw schema-version-1 tvsc_blocks object from panel.blocks().
#'   Its contemporaneous recipe must record scaling = NULL. Externally
#'   transformed values are allowed, but their provenance is not verified.
#' @param scaling One of "sd" (default), "zscore", "none" or "custom".
#'   Names must match exactly. Pooled-level and robust alternatives are deferred.
#' @param scales For custom only, a named plain numeric vector of finite,
#'   strictly positive scales covering every predictor exactly once. Order
#'   need not match the blocks. Must be NULL for all other methods.
#'
#' @return A new object of class c("tvsc_scaled_blocks", "tvsc_blocks").
#'   X1 and X0 contain transformed predictors. Y1, Y0, schema, dates,
#'   dimensions and original raw diagnostics are unchanged. recipe$scaling
#'   records schema_version, method, formula_id, reference units, periods,
#'   cutoff, predictor-by-date centers and divisors, estimated_sds (NULL for
#'   none/custom), zero_variation and zero_scale_fallback masks, zero_policy,
#'   and diagnostics including minimum_positive_applied_scale and fallback_cells.
#'   Centers and divisors are applied parameters, not estimated centers for sd.
#'
#' @details All learned cross-sectional statistics include the target and
#'   every declared donor, using only the supplied training prefix. Sample
#'   SDs use denominator J for J + 1 units. For sd, take the square root of
#'   the mean within-date variances, then divide all dates by that scale,
#'   without centering. RMS evaluation rescales by the largest SD to avoid
#'   unnecessary overflow from squaring. For zscore, subtract the date's
#'   cross-sectional mean and divide by that date's SD. Thus these methods
#'   differ in their date-specific matching coefficients, not just centering.
#'
#'   Exact equality across units at every date permits a scale-1 fallback
#'   for sd. For zscore, exact equality at a date gives its common value as
#'   center and divisor 1, producing zeros. Each fallback is flagged. A zero
#'   or nonfinite computed SD despite unequal values is a numerical error.
#'   Positive scales are not clipped or floored; nonfinite transformed values
#'   are rejected. Finite output is not a certificate of good conditioning.
#'
#'   none records zero centers and unit divisors without arithmetic on X1/X0.
#'   custom divides by the supplied scales without centering. All methods,
#'   including none, mark the output as processed: a second panel.scale()
#'   call is rejected. Start from raw blocks to choose a different method.
#'   This guard does not detect undisclosed upstream transformations.
#'
#'   Rebuild blocks and learned scaling within every validation prefix. Never
#'   scale the full pre-period and then slice it for validation. This function
#'   does not fit weights, calibrate metrics or transform future outcomes.
#'   Predictions must not depend on future target-based z-scores. Supplied
#'   outcome units are never automatically inverted. Scaling alone supplies
#'   no statistical consistency guarantee.
#'
#' @examples
#' panel <- expand.grid(state = c("Target", "DonorA", "DonorB"),
#'                      year = 2000:2005, stringsAsFactors = FALSE)
#' panel$sales <- seq_len(nrow(panel)) * 10
#' prepared <- panel.dataprep(panel, "state", "year", "sales", "Target",
#'                            c("DonorB", "DonorA"), 2005, 2000:2004, "sales")
#' blocks <- panel.blocks(prepared, 2003)
#' scaled <- panel.scale(blocks)
#' scaled$recipe$scaling$divisors
#' unchanged <- panel.scale(blocks, "none")
#' stopifnot(identical(unchanged$X0, blocks$X0), identical(scaled$Y0, blocks$Y0))
#'
#' @note Initial base-R implementation.
panel.scale <- function(blocks, scaling = "sd", scales = NULL) {
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
  whole_vector <- function(value) {
    plain_numeric(value) && is.null(dim(value)) && length(value) > 0L &&
      all(value == trunc(value))
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

  # Validate the raw block boundary before using metadata or performing arithmetic.
  if (!inherits(blocks, "tvsc_blocks") || !has_fields(blocks, c(
      "schema_version", "schema", "periods", "cutoff", "dimensions",
      "X1", "X0", "Y1", "Y0", "recipe", "diagnostics")) ||
      !identical(blocks$schema_version, 1L)) {
    fail("blocks must be a schema-version-1 tvsc_blocks object.")
  }
  if (inherits(blocks, "tvsc_scaled_blocks") ||
      (is.list(blocks$recipe) && !is.null(blocks$recipe$scaling))) {
    fail("Already processed blocks cannot be scaled again; start from raw blocks.")
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
  if (!has_fields(blocks$recipe, c("blocks", "scaling", "outcome_in_predictors")) ||
      !identical(blocks$recipe$blocks, "contemporaneous") ||
      !identical(blocks$recipe$outcome_in_predictors,
                 schema$outcome %in% schema$predictors) ||
      !is.list(blocks$diagnostics)) {
    fail("blocks must retain a contemporaneous raw recipe and diagnostics.")
  }
  predictor_names <- unname(schema$predictors)
  donor_names <- as.character(schema$donors)
  date_names <- as.character(periods)
  predictor_count <- length(predictor_names)
  donor_count <- length(donor_names)
  period_count <- length(periods)
  expected_dimensions <- list(predictors = predictor_count, donors = donor_count,
                              periods = period_count)
  if (!identical(blocks$dimensions, expected_dimensions)) {
    fail("blocks dimensions must agree with the schema.")
  }
  if (!is.list(blocks$X1) || !is.list(blocks$X0) ||
      !identical(names(blocks$X1), date_names) ||
      !identical(names(blocks$X0), date_names)) {
    fail("X1 and X0 must be date-named lists in declared order.")
  }
  for (date_index in seq_len(period_count)) {
    target <- blocks$X1[[date_index]]
    donors <- blocks$X0[[date_index]]
    if (!plain_numeric(target) || !is.null(dim(target)) ||
        !identical(names(target), predictor_names) ||
        !plain_numeric(donors) || !is.matrix(donors) ||
        !identical(dim(donors), c(predictor_count, donor_count)) ||
        !identical(dimnames(donors), list(predictor_names, donor_names))) {
      fail(paste("Invalid finite numeric predictor blocks at date", date_names[date_index]))
    }
  }
  if (!plain_numeric(blocks$Y1) || !is.null(dim(blocks$Y1)) ||
      !identical(names(blocks$Y1), date_names) || !plain_numeric(blocks$Y0) ||
      !is.matrix(blocks$Y0) || !identical(dim(blocks$Y0), c(period_count, donor_count)) ||
      !identical(dimnames(blocks$Y0), list(date_names, donor_names))) {
    fail("Y1 and Y0 must be finite numeric outcomes aligned with the schema.")
  }
  if (!labels_valid(scaling) || length(scaling) != 1L ||
      !scaling %in% c("sd", "zscore", "none", "custom")) {
    fail("scaling must be exactly one of sd, zscore, none or custom.")
  }
  scaling <- unname(scaling)
  if (scaling == "custom") {
    if (!plain_numeric(scales) || !is.null(dim(scales)) ||
        length(scales) != predictor_count || !labels_valid(names(scales)) ||
        !setequal(names(scales), predictor_names) || any(scales <= 0)) {
      fail("custom scales must be finite positive numeric values named once for every predictor.")
    }
    scales <- scales[predictor_names]
  } else if (!is.null(scales)) {
    fail("scales must be NULL unless scaling is custom.")
  }

  # Record applied parameters with a common shape across all four methods.
  centers <- matrix(0, predictor_count, period_count,
                    dimnames = list(predictor_names, date_names))
  divisors <- centers + 1
  fallback <- matrix(FALSE, predictor_count, period_count, dimnames = dimnames(centers))
  estimated_sds <- NULL
  zero_variation <- NULL
  learned <- scaling %in% c("sd", "zscore")
  if (learned) {
    estimated_sds <- centers
    zero_variation <- fallback
    for (date_index in seq_len(period_count)) {
      for (predictor_index in seq_len(predictor_count)) {
        values <- c(blocks$X1[[date_index]][predictor_index],
                    blocks$X0[[date_index]][predictor_index, ])
        identical_values <- all(values == values[1L])
        zero_variation[predictor_index, date_index] <- identical_values
        deviation <- if (identical_values) 0 else stats::sd(values)
        if (!is.finite(deviation) || (!identical_values && deviation <= 0)) {
          fail(paste("Numerical SD failure for predictor", predictor_names[predictor_index],
                     "at date", date_names[date_index]))
        }
        estimated_sds[predictor_index, date_index] <- deviation
        if (scaling == "zscore") {
          center <- if (identical_values) values[1L] else mean(values)
          if (!is.finite(center)) fail("Numerical centering failure.")
          centers[predictor_index, date_index] <- center
          divisors[predictor_index, date_index] <- if (identical_values) 1 else deviation
          fallback[predictor_index, date_index] <- identical_values
        }
      }
    }
    if (scaling == "sd") {
      for (predictor_index in seq_len(predictor_count)) {
        deviations <- estimated_sds[predictor_index, ]
        if (all(zero_variation[predictor_index, ])) {
          fallback[predictor_index, ] <- TRUE
        } else {
          largest <- max(deviations)
          divisor <- largest * sqrt(mean((deviations / largest)^2))
          if (!is.finite(divisor) || divisor <= 0) fail("Numerical RMS SD failure.")
          divisors[predictor_index, ] <- divisor
        }
      }
    }
  } else if (scaling == "custom") {
    divisors[,] <- rep(unname(scales), period_count)
  }

  # none takes a true passthrough path; all other transformations must stay finite.
  result <- blocks
  if (scaling != "none") {
    for (date_index in seq_len(period_count)) {
      target <- (blocks$X1[[date_index]] - centers[, date_index]) / divisors[, date_index]
      donors <- sweep(blocks$X0[[date_index]], 1L, centers[, date_index], "-")
      donors <- sweep(donors, 1L, divisors[, date_index], "/")
      if (any(!is.finite(target)) || any(!is.finite(donors))) {
        fail(paste("Nonfinite transformed predictor values at date", date_names[date_index]))
      }
      result$X1[[date_index]] <- target
      result$X0[[date_index]] <- donors
    }
  }
  formula_ids <- c(sd = "rms_within_date_sample_sd_no_center",
                   zscore = "within_date_mean_and_sample_sd",
                   none = "identity_no_arithmetic", custom = "named_fixed_divisors_no_center")
  result$recipe$scaling <- list(
    schema_version = 1L, method = scaling, formula_id = unname(formula_ids[scaling]),
    reference = list(mode = if (learned) "all_units" else "not_estimated",
                     treated = schema$treated, donors = schema$donors),
    periods = periods, cutoff = blocks$cutoff, centers = centers, divisors = divisors,
    estimated_sds = estimated_sds, zero_variation = zero_variation,
    zero_scale_fallback = fallback,
    zero_policy = if (learned) "unit_divisor_only_for_exact_cross_sectional_equality" else
      if (scaling == "custom") "strictly_positive_supplied_scales" else "identity",
    diagnostics = list(minimum_positive_applied_scale = min(divisors),
                       fallback_cells = sum(fallback), clipping_applied = FALSE)
  )
  class(result) <- c("tvsc_scaled_blocks", "tvsc_blocks")
  result
}
