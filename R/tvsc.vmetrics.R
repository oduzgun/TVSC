#' Construct training-prefix scaling and predictor metrics
#'
#' Source R/tvsc.dataprep.R before this file. No split or legacy scaling
#' operation is required.
#' @param prepared A version-2 tvsc_prepared object with integrated raw blocks.
#' @param cutoff An actual designated pre-period date, with at least three
#'   dates in its training prefix.
#' @param predictor_weights Exactly "uniform", "constant" or "varying", or
#'   supplied weights named by predictor (vector) or predictor/date (matrix).
#' @param scaling Exactly "sd", "zscore", "none" or "custom".
#' @param conversion Exactly "absolute" or "squared"; empirical modes only.
#' @param scales A positive finite named predictor vector for custom only.
#' @param provenance External-input source information, retained unverified.
#' @param rank_tol NULL for max(n, K) times machine epsilon, or a relative
#'   SVD rank tolerance strictly between zero and one.
#' @param weight_sum_tol Absolute supplied-weight sum tolerance in [0, 1).
#' @return A version-1 tvsc_vmetrics object, containing version-2 scaled
#'   blocks and scaling metadata, weights, calibration and diagnostics.
#' @details Learned scales and OLS use donors only. OLS predicts next-date
#'   donor outcomes in original units from the selected transformed features.
#'   Constant calibration uses date intercepts; varying calibration fits each
#'   date and retains the last available metric at the final training date.
#'   Errors inherit tvsc_vmetrics_error and retain stage and diagnostics.
#'   No donor weights, penalties, scores or effects are fitted.
tvsc.vmetrics <- function(prepared, cutoff, predictor_weights = "uniform",
                          scaling = "sd", conversion = "absolute",
                          scales = NULL, provenance = NULL, rank_tol = NULL,
                          weight_sum_tol = 1e-8) {
  if (missing(cutoff)) .tvsc_vm_fail("cutoff must be supplied explicitly.", "prefix")
  blocks <- .tvsc_vm_prefix(prepared, cutoff)
  predictor_names <- unname(blocks$schema$predictors)
  periods <- blocks$periods
  date_names <- as.character(periods)
  predictor_count <- length(predictor_names)
  period_count <- length(periods)
  choice <- function(value, choices) {
    is.character(value) && is.null(dim(value)) && length(value) == 1L &&
      !is.na(value) && value %in% choices
  }
  if (!choice(scaling, c("sd", "zscore", "none", "custom"))) {
    .tvsc_vm_fail("scaling must be exactly sd, zscore, none or custom.", "input")
  }
  if (!choice(conversion, c("absolute", "squared"))) {
    .tvsc_vm_fail("conversion must be exactly absolute or squared.", "input")
  }
  if (!is.null(rank_tol) &&
      (!.tvsc_vm_numeric(rank_tol) || length(rank_tol) != 1L ||
       rank_tol <= 0 || rank_tol >= 1)) {
    .tvsc_vm_fail("rank_tol must be NULL or one finite number strictly between zero and one.", "input")
  }
  if (!.tvsc_vm_numeric(weight_sum_tol) || length(weight_sum_tol) != 1L ||
      weight_sum_tol < 0 || weight_sum_tol >= 1) {
    .tvsc_vm_fail("weight_sum_tol must be one finite number in [0, 1).", "input")
  }
  if (scaling == "custom") {
    if (!.tvsc_vm_numeric(scales) || !.tvsc_vm_labels(names(scales), predictor_names) ||
        any(scales <= 0)) {
      .tvsc_vm_fail("custom scales must be finite positive values named once for every predictor.", "input")
    }
    scales <- scales[predictor_names]
  } else if (!is.null(scales)) {
    .tvsc_vm_fail("scales must be NULL unless scaling is custom.", "input")
  }
  mode <- if (is.character(predictor_weights)) {
    if (!choice(predictor_weights, c("uniform", "constant", "varying"))) {
      .tvsc_vm_fail("predictor_weights must be exactly uniform, constant or varying, or named numeric weights.", "input")
    }
    unname(predictor_weights)
  } else if (is.matrix(predictor_weights)) "supplied_varying" else "supplied_constant"
  empirical <- mode %in% c("constant", "varying")
  outcome_only <- predictor_count == 1L && predictor_names == blocks$schema$outcome
  if (empirical && outcome_only) {
    .tvsc_vm_fail("Outcome-only preparation requires uniform rather than constant or varying calibration.", "input")
  }
  weights <- matrix(1 / predictor_count, predictor_count, period_count,
                    dimnames = list(predictor_names, date_names))
  supplied_check <- NULL
  if (mode == "supplied_constant") {
    if (!.tvsc_vm_numeric(predictor_weights) ||
        !.tvsc_vm_labels(names(predictor_weights), predictor_names)) {
      .tvsc_vm_fail("Supplied constant weights must be finite numeric values named once for every predictor.", "weights")
    }
    weights[,] <- rep(unname(predictor_weights[predictor_names]), period_count)
  } else if (mode == "supplied_varying") {
    if (!typeof(predictor_weights) %in% c("integer", "double") || is.object(predictor_weights) ||
        any(!is.finite(predictor_weights)) ||
        !identical(dim(predictor_weights), c(predictor_count, period_count)) ||
        !.tvsc_vm_labels(rownames(predictor_weights), predictor_names) ||
        !.tvsc_vm_labels(colnames(predictor_weights), date_names)) {
      .tvsc_vm_fail("Supplied weight matrix must have finite values and exactly the predictor and prefix-date names.", "weights")
    }
    weights[,] <- predictor_weights[predictor_names, date_names, drop = FALSE]
  }
  if (mode %in% c("supplied_constant", "supplied_varying")) {
    sums <- colSums(weights)
    supplied_check <- list(tolerance = weight_sum_tol, column_sums = sums,
                           maximum_deviation = max(abs(sums - 1)), normalized = FALSE)
    if (any(weights < 0) || any(!is.finite(sums)) || any(sums <= 0) ||
        any(abs(sums - 1) > weight_sum_tol)) {
      .tvsc_vm_fail("Supplied weights must be nonnegative with positive finite sums within weight_sum_tol of one.",
                    "weights", supplied_check)
    }
  }
  scaled <- .tvsc_vm_scale(blocks, scaling, scales)
  calibration <- list(mode = mode, coefficients = NULL, intercepts = NULL,
                       calibration_dates = periods[FALSE], response_dates = periods[FALSE],
                       conversion = NA_character_, boundary = NA_character_,
                       reason = if (outcome_only && mode == "uniform") "outcome_only_by_design" else
                         if (mode == "uniform") "uniform_requested" else "supplied_weights",
                       rank = NULL)
  if (empirical) {
    fitted <- .tvsc_vm_calibrate(scaled, mode, conversion, rank_tol)
    calibration <- fitted$calibration
    weights <- fitted$weights
  }
  structure(list(
    schema_version = 1L, schema = blocks$schema, periods = periods,
    cutoff = blocks$cutoff, blocks = scaled, weights = weights,
    scaling = scaled$recipe$scaling, calibration = calibration,
    diagnostics = list(supplied_weight_check = supplied_check,
                       weight_column_sums = colSums(weights),
                       scaling = scaled$recipe$scaling$diagnostics,
                       provenance_verified = FALSE),
    provenance = provenance
  ), class = "tvsc_vmetrics")
}

.tvsc_vm_fail <- function(message, stage, diagnostics = list()) {
  stop(structure(list(message = message, call = NULL, stage = stage,
                       diagnostics = diagnostics),
                  class = c("tvsc_vmetrics_error", "error", "condition")))
}

.tvsc_vm_numeric <- function(value) {
  typeof(value) %in% c("integer", "double") && !is.object(value) && is.null(dim(value)) &&
    length(value) > 0L && all(is.finite(value))
}

.tvsc_vm_labels <- function(value, expected) {
  is.character(value) && is.null(dim(value)) && !anyNA(value) &&
    all(nzchar(value)) && !anyDuplicated(value) &&
    length(value) == length(expected) && setequal(value, expected)
}

.tvsc_vm_prefix <- function(prepared, cutoff) {
  fail <- function(message) .tvsc_vm_fail(paste(message, "Rerun tvsc.dataprep() from the original panel."), "prefix")
  fields <- function(value, required) {
    is.list(value) && !anyDuplicated(names(value)) && all(required %in% names(value))
  }
  if (!inherits(prepared, "tvsc_prepared") ||
      !fields(prepared, c("schema_version", "schema", "training", "blocks", "recipe",
                         "future_donors", "future_target", "diagnostics")) ||
      !identical(prepared$schema_version, 2L)) {
    fail("prepared must be a schema-version-2 tvsc_prepared object.")
  }
  schema <- prepared$schema
  if (!fields(schema, c("unit", "time", "outcome", "predictors", "treated", "donors",
                        "pre_period", "intervention_time", "period_step"))) fail("Invalid schema.")
  periods <- schema$pre_period
  spacing <- schema$period_step
  if (!.tvsc_vm_numeric(periods) || length(periods) < 3L || any(periods != floor(periods)) ||
      !.tvsc_vm_numeric(spacing) || length(spacing) != 1L || spacing <= 0 || spacing != floor(spacing) ||
      any(diff(periods) != spacing) || !.tvsc_vm_numeric(schema$intervention_time) ||
      length(schema$intervention_time) != 1L || schema$intervention_time != floor(schema$intervention_time) ||
      any(periods >= schema$intervention_time)) fail("Invalid pre-period metadata.")
  if (!.tvsc_vm_numeric(cutoff) || length(cutoff) != 1L || !cutoff %in% periods) {
    .tvsc_vm_fail("cutoff must be an actual designated pre-period date.", "prefix")
  }
  cutoff_index <- match(cutoff, periods)
  if (cutoff_index < 3L) .tvsc_vm_fail("The training prefix must contain at least three dates.", "prefix")
  prefix_indices <- seq_len(cutoff_index)
  required_columns <- unname(unique(c(schema$unit, schema$time, schema$outcome, schema$predictors)))
  if (!is.data.frame(prepared$training) || !identical(names(prepared$training), required_columns) ||
      !fields(prepared$recipe, c("blocks", "scaling", "outcome_in_predictors")) ||
      !identical(prepared$recipe$blocks, "contemporaneous") || !is.null(prepared$recipe$scaling) ||
      !identical(prepared$recipe$outcome_in_predictors, schema$outcome %in% schema$predictors)) {
    fail("Invalid raw training structure or recipe.")
  }
  if (!is.character(schema$time) || length(schema$time) != 1L || is.na(schema$time) ||
      !.tvsc_vm_numeric(prepared$training[[schema$time]])) fail("Invalid training time index.")
  training <- prepared$training[prepared$training[[schema$time]] %in% periods[prefix_indices], , drop = FALSE]
  reconstructed <- tryCatch(tvsc.dataprep(training, schema$unit, schema$time, schema$outcome,
                                 schema$treated, schema$donors, schema$intervention_time,
                                 periods[prefix_indices], schema$predictors, spacing),
                            error = function(error) fail(conditionMessage(error)))
  predictor_names <- unname(schema$predictors)
  donor_names <- as.character(schema$donors)
  date_names <- as.character(periods)
  predictor_count <- length(predictor_names)
  donor_count <- length(donor_names)
  period_count <- length(periods)
  expected_schema <- reconstructed$schema
  expected_schema$pre_period <- periods
  if (!identical(schema, expected_schema)) fail("Incompatible schema metadata.")
  training_ids <- prepared$training[[schema$unit]]
  if (is.factor(training_ids)) training_ids <- as.character(training_ids)
  if (nrow(prepared$training) != period_count * (donor_count + 1L) ||
      !identical(as.character(training_ids), as.character(rep(c(schema$treated, schema$donors), period_count))) ||
      any(prepared$training[[schema$time]] != rep(periods, each = donor_count + 1L))) fail("Invalid training key ordering.")
  raw <- prepared$blocks
  if (!inherits(raw, "tvsc_blocks") || inherits(raw, "tvsc_scaled_blocks") ||
      !fields(raw, c("schema_version", "schema", "periods", "cutoff", "dimensions",
                     "X1", "X0", "Y1", "Y0", "recipe", "diagnostics")) ||
      !identical(raw$schema_version, 2L) || !identical(raw$schema, schema) ||
      !identical(raw$periods, periods) || !identical(raw$cutoff, periods[period_count]) ||
      !identical(raw$recipe, prepared$recipe) ||
      !identical(raw$dimensions, list(predictors = predictor_count, donors = donor_count, periods = period_count)) ||
      !is.list(raw$X1) || !identical(names(raw$X1), date_names) ||
      !is.list(raw$X0) || !identical(names(raw$X0), date_names) ||
      !is.numeric(raw$Y1) || !is.null(dim(raw$Y1)) || !identical(names(raw$Y1), date_names) ||
      !is.matrix(raw$Y0) || !is.numeric(raw$Y0) ||
      !identical(dim(raw$Y0), c(period_count, donor_count)) ||
      !identical(dimnames(raw$Y0), list(date_names, donor_names))) fail("Incompatible raw blocks.")
  for (date_index in seq_along(periods)) {
    if (!is.numeric(raw$X1[[date_index]]) || !is.null(dim(raw$X1[[date_index]])) ||
        !identical(names(raw$X1[[date_index]]), predictor_names) ||
        !is.matrix(raw$X0[[date_index]]) || !is.numeric(raw$X0[[date_index]]) ||
        !identical(dim(raw$X0[[date_index]]), c(predictor_count, donor_count)) ||
        !identical(dimnames(raw$X0[[date_index]]), list(predictor_names, donor_names))) {
      fail("Incompatible raw block dimensions or ordering.")
    }
  }
  result <- reconstructed$blocks
  if (!identical(raw$X1[prefix_indices], result$X1) ||
      !identical(raw$X0[prefix_indices], result$X0) ||
      !identical(raw$Y1[prefix_indices], result$Y1) ||
      !identical(raw$Y0[prefix_indices, , drop = FALSE], result$Y0)) {
    fail("Raw prefix blocks disagree with the training table.")
  }
  result$schema <- schema
  result
}

.tvsc_vm_scale <- function(blocks, method, scales) {
  predictor_names <- unname(blocks$schema$predictors)
  date_names <- as.character(blocks$periods)
  centers <- matrix(0, length(predictor_names), length(date_names),
                    dimnames = list(predictor_names, date_names))
  divisors <- centers + 1
  fallback <- matrix(FALSE, nrow(centers), ncol(centers), dimnames = dimnames(centers))
  estimated_sds <- NULL
  zero_variation <- NULL
  learned <- method %in% c("sd", "zscore")
  if (learned) {
    estimated_sds <- centers
    zero_variation <- fallback
    for (date_index in seq_along(date_names)) {
      for (predictor_index in seq_along(predictor_names)) {
        values <- blocks$X0[[date_index]][predictor_index, ]
        equal <- all(values == values[1L])
        zero_variation[predictor_index, date_index] <- equal
        deviation <- if (equal) 0 else stats::sd(values)
        if (!is.finite(deviation) || (!equal && deviation <= 0)) {
          .tvsc_vm_fail("Numerical donor SD failure.", "scaling",
                        list(predictor = predictor_names[predictor_index], date = blocks$periods[date_index]))
        }
        estimated_sds[predictor_index, date_index] <- deviation
        if (method == "zscore") {
          center <- if (equal) values[1L] else mean(values)
          if (!is.finite(center)) .tvsc_vm_fail("Numerical donor centering failure.", "scaling")
          centers[predictor_index, date_index] <- center
          divisors[predictor_index, date_index] <- if (equal) 1 else deviation
          fallback[predictor_index, date_index] <- equal
        }
      }
    }
    if (method == "sd") {
      for (predictor_index in seq_along(predictor_names)) {
        deviations <- estimated_sds[predictor_index, ]
        if (all(zero_variation[predictor_index, ])) {
          fallback[predictor_index, ] <- TRUE
        } else {
          largest <- max(deviations)
          divisor <- largest * sqrt(mean((deviations / largest)^2))
          if (!is.finite(divisor) || divisor <= 0) .tvsc_vm_fail("Numerical RMS donor SD failure.", "scaling")
          divisors[predictor_index, ] <- divisor
        }
      }
    }
  } else if (method == "custom") {
    divisors[,] <- rep(unname(scales), ncol(divisors))
  }
  result <- blocks
  if (method != "none") {
    for (date_index in seq_along(date_names)) {
      target <- (blocks$X1[[date_index]] - centers[, date_index]) / divisors[, date_index]
      donors <- sweep(blocks$X0[[date_index]], 1L, centers[, date_index], "-")
      donors <- sweep(donors, 1L, divisors[, date_index], "/")
      if (any(!is.finite(target)) || any(!is.finite(donors))) {
        .tvsc_vm_fail("Nonfinite transformed predictor values.", "scaling", list(date = blocks$periods[date_index]))
      }
      result$X1[[date_index]] <- target
      result$X0[[date_index]] <- donors
    }
  }
  fallback_indices <- which(fallback, arr.ind = TRUE)
  fallback_records <- data.frame(predictor = predictor_names[fallback_indices[, 1L]],
                                 period = blocks$periods[fallback_indices[, 2L]],
                                 stringsAsFactors = FALSE)
  result$recipe$scaling <- list(
    schema_version = 2L, method = method,
    reference = list(mode = if (learned) "donors_only" else "not_estimated", donors = blocks$schema$donors),
    periods = blocks$periods, cutoff = blocks$cutoff, centers = centers, divisors = divisors,
    estimated_sds = estimated_sds, zero_variation = zero_variation,
    zero_scale_fallback = fallback, fallback_records = fallback_records,
    zero_policy = if (learned) "unit_divisor_only_for_exact_donor_equality" else
      if (method == "custom") "strictly_positive_supplied_scales" else "identity",
    diagnostics = list(fallback_cells = sum(fallback), minimum_positive_applied_scale = min(divisors),
                       clipping_applied = FALSE)
  )
  class(result) <- c("tvsc_scaled_blocks", "tvsc_blocks")
  result
}

.tvsc_vm_calibrate <- function(blocks, mode, conversion, rank_tol) {
  predictor_names <- unname(blocks$schema$predictors)
  predictor_count <- length(predictor_names)
  periods <- blocks$periods
  calibration_indices <- seq_len(length(periods) - 1L)
  calibration_dates <- periods[calibration_indices]
  feature_means <- lapply(calibration_indices, function(date_index) rowMeans(blocks$X0[[date_index]]))
  responses <- lapply(calibration_indices, function(date_index) blocks$Y0[date_index + 1L, ])
  response_means <- vapply(responses, mean, numeric(1))
  designs <- lapply(calibration_indices, function(date_index) {
    t(sweep(blocks$X0[[date_index]], 1L, feature_means[[date_index]], "-"))
  })
  centered_responses <- lapply(calibration_indices, function(date_index) responses[[date_index]] - response_means[date_index])
  solve_slopes <- function(design, response, dates) {
    if (any(!is.finite(design)) || any(!is.finite(response))) {
      .tvsc_vm_fail("Nonfinite centered calibration data.", "calibration", list(dates = dates))
    }
    decomposition <- tryCatch(svd(design), error = function(error) {
      .tvsc_vm_fail("Calibration SVD failed.", "calibration", list(dates = dates, cause = conditionMessage(error)))
    })
    singular_values <- decomposition$d
    if (any(!is.finite(singular_values))) {
      .tvsc_vm_fail("Nonfinite calibration singular values.", "calibration", list(dates = dates))
    }
    relative <- if (is.null(rank_tol)) max(dim(design)) * .Machine$double.eps else rank_tol
    threshold <- relative * max(singular_values)
    numerical_rank <- sum(singular_values > threshold)
    diagnostic <- list(requested_tolerance = rank_tol, relative_tolerance = relative,
                       threshold = threshold, singular_values = singular_values,
                       numerical_rank = numerical_rank, required_rank = predictor_count,
                       dates = dates, condition_number = if (numerical_rank < predictor_count) Inf else
                         max(singular_values) / min(singular_values))
    if (numerical_rank != predictor_count) {
      .tvsc_vm_fail("Calibration design is rank deficient under the declared tolerance.", "calibration", diagnostic)
    }
    coefficients <- as.numeric(decomposition$v %*%
                                 (as.numeric(crossprod(decomposition$u, response)) / singular_values))
    names(coefficients) <- predictor_names
    if (any(!is.finite(coefficients))) .tvsc_vm_fail("Nonfinite OLS slopes.", "calibration", diagnostic)
    list(coefficients = coefficients, diagnostic = diagnostic)
  }
  convert <- function(coefficients, dates) {
    largest <- max(abs(coefficients))
    if (largest == 0) .tvsc_vm_fail("All-zero OLS slopes cannot define empirical weights.", "calibration", list(dates = dates))
    relative <- abs(coefficients) / largest
    if (conversion == "squared") relative <- relative^2
    relative / sum(relative)
  }
  intercepts <- stats::setNames(numeric(length(calibration_indices)), as.character(calibration_dates))
  weights <- matrix(0, predictor_count, length(periods), dimnames = list(predictor_names, as.character(periods)))
  if (mode == "constant") {
    fitted <- solve_slopes(do.call(rbind, designs), unlist(centered_responses, use.names = FALSE), calibration_dates)
    coefficients <- fitted$coefficients
    rank <- list(pooled = fitted$diagnostic)
    weights[,] <- rep(convert(coefficients, calibration_dates), length(periods))
    for (date_index in calibration_indices) {
      intercepts[date_index] <- response_means[date_index] - sum(feature_means[[date_index]] * coefficients)
    }
  } else {
    coefficients <- matrix(0, predictor_count, length(calibration_indices),
                            dimnames = list(predictor_names, as.character(calibration_dates)))
    rank <- stats::setNames(vector("list", length(calibration_indices)), as.character(calibration_dates))
    for (date_index in calibration_indices) {
      fitted <- solve_slopes(designs[[date_index]], centered_responses[[date_index]], calibration_dates[date_index])
      coefficients[, date_index] <- fitted$coefficients
      rank[[date_index]] <- fitted$diagnostic
      weights[, date_index] <- convert(fitted$coefficients, calibration_dates[date_index])
      intercepts[date_index] <- response_means[date_index] - sum(feature_means[[date_index]] * fitted$coefficients)
    }
    weights[, length(periods)] <- weights[, length(periods) - 1L]
  }
  if (any(!is.finite(intercepts)) || any(!is.finite(weights))) {
    .tvsc_vm_fail("Nonfinite calibration output.", "calibration", list(dates = calibration_dates))
  }
  list(weights = weights, calibration = list(
    mode = mode, coefficients = coefficients, intercepts = intercepts,
    calibration_dates = calibration_dates, response_dates = periods[-1L],
    conversion = conversion, boundary = if (mode == "varying") "retain_final_available_metric" else "common_metric",
    reason = NULL, rank = rank, reference_donors = blocks$schema$donors
  ))
}
