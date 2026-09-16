source("R/tvsc.dataprep.R")
source("R/tvsc.vmetrics.R")

close_values <- function(actual, expected, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(unname(actual), unname(expected), tolerance = tolerance,
                             check.attributes = FALSE)))
}
expect_error <- function(expression, pattern, stage = NULL) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "tvsc_vmetrics_error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
  if (!is.null(stage)) stopifnot(identical(error$stage, stage))
  invisible(error)
}
donor_ids <- paste0("Donor", 1:5)
periods <- 2000:2005
panel <- expand.grid(unit = c("Target", donor_ids), period = periods,
                      stringsAsFactors = FALSE)
panel$income <- 0
panel$price <- 0
panel$sales <- 0
income_pattern <- c(-2, -1, 0, 1, 2)
price_pattern <- c(2, -1, -2, -1, 2)
for (date_index in seq_along(periods)) {
  donor_rows <- panel$period == periods[date_index] & panel$unit != "Target"
  target_row <- panel$period == periods[date_index] & panel$unit == "Target"
  panel$income[donor_rows] <- 10 * date_index + date_index * income_pattern
  panel$price[donor_rows] <- 30 - date_index + (1 + date_index / 4) * price_pattern
  panel$income[target_row] <- 1000 + date_index
  panel$price[target_row] <- -1000 - date_index
  panel$sales[target_row] <- -date_index
  if (date_index > 1L) {
    previous_rows <- panel$period == periods[date_index - 1L] & panel$unit != "Target"
    panel$sales[donor_rows] <- 100 + date_index - 1 +
      2 * panel$income[previous_rows] - panel$price[previous_rows]
  }
}
prepare <- function(data = panel, predictors = c("income", "price")) {
  tvsc.dataprep(data, "unit", "period", "sales", "Target", rev(donor_ids),
                2006, periods, predictors)
}
prepared <- prepare()
before <- serialize(prepared, NULL)
uniform <- tvsc.vmetrics(prepared, 2005)
stopifnot(identical(before, serialize(prepared, NULL)),
           inherits(uniform, "tvsc_vmetrics"), identical(uniform$schema_version, 1L),
           inherits(uniform$blocks, "tvsc_scaled_blocks"),
           identical(uniform$blocks$schema_version, 2L), identical(uniform$scaling$schema_version, 2L),
           identical(uniform$schema, prepared$schema), identical(uniform$periods, periods),
           identical(uniform$cutoff, 2005L),
           identical(dimnames(uniform$weights), list(c("income", "price"), as.character(periods))),
           all(uniform$weights == 0.5), is.na(uniform$calibration$conversion),
           identical(uniform$blocks$Y1, prepared$blocks$Y1),
           identical(uniform$blocks$Y0, prepared$blocks$Y0),
           identical(uniform$scaling$reference$mode, "donors_only"),
           identical(uniform$scaling$reference$donors, rev(donor_ids)),
           is.null(uniform$lambda), is.null(uniform$donor_weights), is.null(uniform$splits))
expected_sd <- vapply(c("income", "price"), function(predictor) {
  sqrt(mean(vapply(prepared$blocks$X0, function(block) var(block[predictor, ]), numeric(1))))
}, numeric(1))
close_values(uniform$scaling$divisors[, 1], expected_sd)
close_values(uniform$blocks$X1[[1]], prepared$blocks$X1[[1]] / expected_sd)
stopifnot(all(uniform$scaling$centers == 0), !any(uniform$scaling$zero_scale_fallback))
none <- tvsc.vmetrics(prepared, 2005, scaling = "none")
stopifnot(identical(none$blocks$X1, prepared$blocks$X1),
           identical(none$blocks$X0, prepared$blocks$X0),
           is.null(none$scaling$estimated_sds), is.null(none$scaling$zero_variation),
           identical(none$scaling$reference$mode, "not_estimated"))
zscore <- tvsc.vmetrics(prepared, 2005, scaling = "zscore")
for (date_index in seq_along(periods)) {
  donor_block <- prepared$blocks$X0[[date_index]]
  close_values(zscore$scaling$centers[, date_index], rowMeans(donor_block))
  close_values(zscore$scaling$divisors[, date_index], apply(donor_block, 1L, sd))
  close_values(rowMeans(zscore$blocks$X0[[date_index]]), c(0, 0))
  close_values(apply(zscore$blocks$X0[[date_index]], 1L, var), c(1, 1))
}
custom <- tvsc.vmetrics(prepared, 2005, scaling = "custom", scales = c(price = 4, income = 2),
                        provenance = "Declared external divisors")
close_values(custom$blocks$X1[[1]], prepared$blocks$X1[[1]] / c(2, 4))
stopifnot(identical(custom$provenance, "Declared external divisors"),
           isFALSE(custom$diagnostics$provenance_verified), is.null(custom$scaling$estimated_sds))

constant <- tvsc.vmetrics(prepared, 2005, "constant", scaling = "none")
close_values(constant$calibration$coefficients, c(2, -1))
close_values(constant$calibration$intercepts, 101:105)
close_values(constant$weights[, 1], c(2 / 3, 1 / 3))
stopifnot(identical(constant$calibration$calibration_dates, 2000:2004),
           identical(constant$calibration$response_dates, 2001:2005),
           identical(constant$calibration$rank$pooled$required_rank, 2L),
           all(constant$weights == constant$weights[, 1]))
close_values(constant$calibration$rank$pooled$relative_tolerance, 25 * .Machine$double.eps, 0)
varying <- tvsc.vmetrics(prepared, 2005, "varying", scaling = "none")
close_values(varying$calibration$coefficients, matrix(rep(c(2, -1), 5), 2, 5))
close_values(varying$weights, constant$weights)
stopifnot(identical(varying$weights[, 6], varying$weights[, 5]),
           identical(varying$calibration$boundary, "retain_final_available_metric"),
           identical(names(varying$calibration$rank), as.character(2000:2004)))
for (method in c("sd", "zscore", "none", "custom")) {
  for (mode in c("constant", "varying")) {
    result <- tvsc.vmetrics(prepared, 2005, mode, scaling = method,
                            scales = if (method == "custom") c(income = 2, price = 4) else NULL)
    design_blocks <- lapply(result$blocks$X0[1:5], t)
    if (mode == "constant") {
      indicators <- diag(5)[rep(1:5, each = 5), , drop = FALSE]
      design <- cbind(indicators, do.call(rbind, design_blocks))
      response <- as.numeric(t(result$blocks$Y0[2:6, , drop = FALSE]))
      reference <- lm.fit(design, response)$coefficients
      close_values(result$calibration$coefficients, reference[6:7])
      close_values(result$calibration$intercepts, reference[1:5])
    } else {
      for (date_index in 1:5) {
        reference <- lm.fit(cbind(1, design_blocks[[date_index]]),
                             result$blocks$Y0[date_index + 1L, ])$coefficients
        close_values(result$calibration$coefficients[, date_index], reference[-1L])
        close_values(result$calibration$intercepts[date_index], reference[1L])
      }
    }
  }
}
varying_panel <- panel
for (date_index in 1:5) {
  previous <- panel$period == periods[date_index] & panel$unit != "Target"
  response <- panel$period == periods[date_index + 1L] & panel$unit != "Target"
  varying_panel$sales[response] <- 100 + date_index +
    date_index * panel$income[previous] - (date_index + 1) * panel$price[previous]
}
varying_squared <- tvsc.vmetrics(prepare(varying_panel), 2005, "varying",
                                scaling = "none", conversion = "squared")
for (date_index in 1:5) {
  expected <- c(date_index^2, (date_index + 1)^2)
  close_values(varying_squared$weights[, date_index], expected / sum(expected))
}
stopifnot(identical(varying_squared$weights[, 6], varying_squared$weights[, 5]))
scaled_constant <- tvsc.vmetrics(prepared, 2005, "constant")
close_values(scaled_constant$calibration$coefficients, c(2, -1) * expected_sd)

changed_target <- panel
target_rows <- changed_target$unit == "Target"
changed_target$income[target_rows] <- 1e6
changed_target$price[target_rows] <- -1e6
changed_target$sales[target_rows] <- 42
target_result <- tvsc.vmetrics(prepare(changed_target), 2005, "constant")
stopifnot(identical(target_result$scaling, scaled_constant$scaling),
           identical(target_result$calibration, scaled_constant$calibration),
           identical(target_result$weights, scaled_constant$weights))
prefix <- tvsc.vmetrics(prepared, 2003, "varying")
ignored <- prepared
ignored$training$income[ignored$training$period > 2003] <- Inf
ignored$training$sales[ignored$training$period > 2003] <- NA_real_
for (date_index in 5:6) {
  ignored$blocks$X1[[date_index]][] <- NA_real_
  ignored$blocks$X0[[date_index]][] <- Inf
}
ignored$blocks$Y0[5:6, ] <- NA_real_
ignored$blocks$Y1[5:6] <- NA_real_
ignored$future_target <- data.frame(unchecked = NA)
ignored$future_donors <- data.frame(unchecked = Inf)
stopifnot(identical(prefix, tvsc.vmetrics(ignored, 2003, "varying")),
           identical(prefix$periods, 2000:2003), identical(prefix$calibration$response_dates, 2001:2003))
prefix_panel <- panel[panel$period <= 2003, ]
short_prepared <- tvsc.dataprep(prefix_panel, "unit", "period", "sales", "Target", rev(donor_ids),
                               2006, 2000:2003, c("income", "price"))
short_result <- tvsc.vmetrics(short_prepared, 2003, "varying")
stopifnot(identical(prefix$scaling, short_result$scaling),
           identical(prefix$calibration, short_result$calibration),
           identical(prefix$weights, short_result$weights))
stopifnot(identical(none, tvsc.vmetrics(prepare(panel[rev(seq_len(nrow(panel))), ]), 2005, scaling = "none")))

supplied <- tvsc.vmetrics(prepared, 2005, c(price = 0.25, income = 0.75), scaling = "none")
close_values(supplied$weights[, 1], c(0.75, 0.25))
stopifnot(identical(supplied$calibration$reason, "supplied_weights"), is.na(supplied$calibration$conversion))
matrix_weights <- varying_squared$weights[c("price", "income"), as.character(2005:2000)]
supplied_matrix <- tvsc.vmetrics(prepared, 2005, matrix_weights, scaling = "none")
stopifnot(identical(supplied_matrix$weights, varying_squared$weights))
almost <- c(income = 0.75, price = 0.25 + 1e-9)
accepted <- tvsc.vmetrics(prepared, 2005, almost)
stopifnot(identical(accepted$weights[, 1], almost),
           isFALSE(accepted$diagnostics$supplied_weight_check$normalized))
expect_error(tvsc.vmetrics(prepared, 2005, almost, weight_sum_tol = 0), "weight_sum_tol", "weights")
expect_error(tvsc.vmetrics(prepared, 2003, matrix_weights), "prefix-date names", "weights")
for (bad in list(c(income = -1e-15, price = 1), c(income = 0, price = 0),
                 c(income = 1, price = 1))) {
  expect_error(tvsc.vmetrics(prepared, 2005, bad), "Supplied weights", "weights")
}
for (bad in list(c(0.5, 0.5), c(income = 1), c(income = NA_real_, price = 1),
                 c(income = 1 + 1i, price = 0), c(income = 0.5, income = 0.5))) {
  expect_error(tvsc.vmetrics(prepared, 2005, bad), "Supplied constant weights", "weights")
}

flat_panel <- panel
flat_panel$income[flat_panel$unit != "Target"] <- flat_panel$period[flat_panel$unit != "Target"]
flat <- prepare(flat_panel)
for (method in c("sd", "zscore")) {
  result <- tvsc.vmetrics(flat, 2005, scaling = method)
  stopifnot(all(result$scaling$divisors["income", ] == 1),
             all(result$scaling$zero_scale_fallback["income", ]),
             nrow(result$scaling$fallback_records) == 6L)
  close_values(result$blocks$X1[[1]]["income"] - result$blocks$X0[[1]]["income", 1],
                flat$blocks$X1[[1]]["income"] - flat$blocks$X0[[1]]["income", 1])
}
expect_error(tvsc.vmetrics(flat, 2005, "constant"), "rank deficient", "calibration")
one_flat_date <- panel
one_flat_date$income[one_flat_date$unit != "Target" & one_flat_date$period == 2000] <- 10
one_flat <- prepare(one_flat_date)
one_sd <- tvsc.vmetrics(one_flat, 2005)
one_z <- tvsc.vmetrics(one_flat, 2005, scaling = "zscore")
stopifnot(one_sd$scaling$zero_variation["income", 1], !any(one_sd$scaling$zero_scale_fallback),
           one_z$scaling$zero_scale_fallback["income", 1], sum(one_z$scaling$zero_scale_fallback) == 1L)
stopifnot(inherits(tvsc.vmetrics(one_flat, 2005, "constant"), "tvsc_vmetrics"))
rank_error <- expect_error(tvsc.vmetrics(one_flat, 2005, "varying"), "rank deficient", "calibration")
stopifnot(identical(rank_error$diagnostics$dates, 2000L), rank_error$diagnostics$numerical_rank < 2L)
last_flat <- panel
last_flat$income[last_flat$unit != "Target" & last_flat$period == 2004] <- 50
last_error <- expect_error(tvsc.vmetrics(prepare(last_flat), 2005, "varying"), "rank deficient")
stopifnot(identical(last_error$diagnostics$dates, 2004L))
collinear <- panel
collinear$price <- 3 * collinear$income
expect_error(tvsc.vmetrics(prepare(collinear), 2005, "constant"), "rank deficient")
zero_response <- panel
zero_response$sales[zero_response$unit != "Target"] <- 10
expect_error(tvsc.vmetrics(prepare(zero_response), 2005, "constant", scaling = "none"), "All-zero OLS slopes")
expect_error(tvsc.vmetrics(prepare(zero_response), 2005, "varying", scaling = "none"), "All-zero OLS slopes")
override <- tvsc.vmetrics(prepared, 2005, "constant", rank_tol = 1e-8)
stopifnot(identical(override$calibration$rank$pooled$requested_tolerance, 1e-8),
           identical(override$calibration$rank$pooled$relative_tolerance, 1e-8))
near_collinear <- panel
near_collinear$price <- near_collinear$income + 1e-7 * near_collinear$price
near_prepared <- prepare(near_collinear)
near_result <- tvsc.vmetrics(near_prepared, 2005, "constant", scaling = "none")
stopifnot(near_result$calibration$rank$pooled$condition_number > 1e6)
expect_error(tvsc.vmetrics(near_prepared, 2005, "constant", scaling = "none", rank_tol = 1e-4),
              "rank deficient")
huge_slopes <- panel
huge_slopes$sales[huge_slopes$unit != "Target"] <- huge_slopes$sales[huge_slopes$unit != "Target"] * 1e200
huge_result <- tvsc.vmetrics(prepare(huge_slopes), 2005, "constant", scaling = "none", conversion = "squared")
close_values(huge_result$weights[, 1], c(0.8, 0.2))

outcome_only <- prepare(predictors = "sales")
outcome_metrics <- tvsc.vmetrics(outcome_only, 2005)
stopifnot(identical(dim(outcome_metrics$weights), c(1L, 6L)), all(outcome_metrics$weights == 1),
           identical(dim(outcome_metrics$blocks$X0[[1]]), c(1L, 5L)),
           identical(outcome_metrics$calibration$reason, "outcome_only_by_design"))
expect_error(tvsc.vmetrics(outcome_only, 2005, "constant"), "Outcome-only")
expect_error(tvsc.vmetrics(outcome_only, 2005, "varying"), "Outcome-only")
one_feature <- tvsc.vmetrics(prepare(predictors = "income"), 2005, "varying")
stopifnot(all(one_feature$weights == 1), identical(dim(one_feature$calibration$coefficients), c(1L, 5L)))
numeric_panel <- panel
numeric_panel$unit <- match(numeric_panel$unit, c("Target", donor_ids))
numeric_prepared <- tvsc.dataprep(numeric_panel, "unit", "period", "sales", 1,
                                 c(6, 5, 4, 3, 2), 2006, as.numeric(periods), c("income", "price"))
close_values(tvsc.vmetrics(numeric_prepared, 2005, "constant")$weights, scaled_constant$weights)
factor_panel <- panel
factor_panel$unit <- factor(factor_panel$unit)
stopifnot(identical(tvsc.vmetrics(prepare(factor_panel), 2005), uniform))
spaced_panel <- panel
spaced_panel$period <- 2000 + 2 * (spaced_panel$period - 2000)
spaced <- tvsc.dataprep(spaced_panel, "unit", "period", "sales", "Target", rev(donor_ids),
                        2012, seq(2000, 2010, 2), c("income", "price"), period_step = 2)
spaced_metrics <- tvsc.vmetrics(spaced, 2006, "varying")
close_values(spaced_metrics$weights, prefix$weights)
stopifnot(identical(spaced_metrics$calibration$response_dates, c(2002, 2004, 2006)))

expect_error(tvsc.vmetrics(prepared), "supplied explicitly")
for (bad in list(NULL, NA_real_, Inf, "2005", TRUE, 2, 2006, c(2003, 2005), matrix(2005))) {
  expect_error(tvsc.vmetrics(prepared, bad), "actual designated")
}
expect_error(tvsc.vmetrics(prepared, 2001), "at least three")
for (bad in c("empirical_constant", "empirical_varying", "vary", "Uniform")) {
  expect_error(tvsc.vmetrics(prepared, 2005, bad), "exactly uniform")
}
expect_error(tvsc.vmetrics(prepared, 2005, scaling = "SD"), "scaling must be exactly")
expect_error(tvsc.vmetrics(prepared, 2005, conversion = "abs"), "conversion must be exactly")
expect_error(tvsc.vmetrics(prepared, 2005, scales = c(income = 2, price = 4)), "scales must be NULL")
for (bad in list(NULL, c(income = 0, price = 1), c(income = Inf, price = 1), c(1, 2))) {
  expect_error(tvsc.vmetrics(prepared, 2005, scaling = "custom", scales = bad), "custom scales")
}
for (bad in list(NA_real_, Inf, -1, 0, 1, "1e-8", c(1e-8, 1e-9), 1i)) {
  expect_error(tvsc.vmetrics(prepared, 2005, rank_tol = bad), "rank_tol must be")
}
for (bad in list(NULL, NA_real_, Inf, -1, 1, "1e-8", c(1e-8, 1e-9), 1i)) {
  expect_error(tvsc.vmetrics(prepared, 2005, weight_sum_tol = bad), "weight_sum_tol must be")
}
invalid <- prepared
invalid$schema_version <- 1L
expect_error(tvsc.vmetrics(invalid, 2005), "schema-version-2")
invalid <- prepared
invalid$blocks$schema_version <- 1L
expect_error(tvsc.vmetrics(invalid, 2005), "Incompatible raw blocks")
invalid <- prepared
invalid$blocks$X0[[1]] <- invalid$blocks$X0[[1]][, 5:1, drop = FALSE]
expect_error(tvsc.vmetrics(invalid, 2005), "dimensions or ordering")
invalid <- prepared
invalid$blocks$X1[[1]][1] <- invalid$blocks$X1[[1]][1] + 1
expect_error(tvsc.vmetrics(invalid, 2005), "disagree with the training table")
invalid <- prepared
invalid$recipe$scaling <- list(method = "none")
expect_error(tvsc.vmetrics(invalid, 2005), "raw training structure or recipe")
invalid <- prepared
invalid$training$income[1] <- NA_real_
expect_error(tvsc.vmetrics(invalid, 2005), "Nonfinite training value")
tiny_scale <- c(income = .Machine$double.xmin, price = 1)
expect_error(tvsc.vmetrics(prepared, 2005, scaling = "custom", scales = tiny_scale), "Nonfinite transformed")
overflow_panel <- panel
overflow_panel$income[overflow_panel$unit != "Target"] <- rep(c(-1e308, 1e308, 0, 1e308, -1e308), 6)
expect_error(tvsc.vmetrics(prepare(overflow_panel), 2005), "Numerical donor SD failure")

isolated <- new.env(parent = baseenv())
sys.source("R/tvsc.dataprep.R", envir = isolated)
sys.source("R/tvsc.vmetrics.R", envir = isolated)
stopifnot(identical(isolated$tvsc.vmetrics(prepared, 2005, "constant"), scaled_constant),
           !exists("tvsc.split", envir = isolated, inherits = FALSE),
           !exists("tvsc.scale", envir = isolated, inherits = FALSE))
cat("All prefix scaling and predictor-metric checks passed.\n")