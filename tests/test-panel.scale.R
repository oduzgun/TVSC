source("R/panel.dataprep.R")
source("R/panel.split.R")
source("R/panel.blocks.R")
source("R/panel.scale.R")

expect_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}
expect_equal <- function(actual, expected) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-12)))
}

# Each date has three units with hand-calculated cross-sectional SDs.
panel <- expand.grid(state = c("Target", "DonorA", "DonorB"),
                     year = 2000:2005, stringsAsFactors = FALSE)
date_position <- panel$year - 1999
unit_offset <- match(panel$state, c("Target", "DonorA", "DonorB")) - 2
panel$sales <- seq_len(nrow(panel)) * 10
panel$signal <- 100 + 10 * (date_position - 1) + unit_offset * date_position
panel$intermittent <- 50 + unit_offset * c(0, 2, 0, 4, 8, 16)[date_position]
panel$shared <- 1000 + date_position
panel$fixed <- 7
predictors <- c("signal", "intermittent", "shared", "fixed")
prepare <- function(data = panel, matching = predictors) {
  panel.dataprep(data, "state", "year", "sales", "Target", c("DonorB", "DonorA"),
                 2005, 2000:2004, matching)
}
prepared <- prepare()
blocks <- panel.blocks(prepared, 2003)
before <- serialize(blocks, NULL)
scaled <- panel.scale(blocks)
stopifnot(identical(before, serialize(blocks, NULL)),
          identical(class(scaled), c("tvsc_scaled_blocks", "tvsc_blocks")),
          identical(scaled, panel.scale(blocks, "sd")))
for (field in setdiff(names(blocks), c("X1", "X0", "recipe"))) {
  stopifnot(identical(scaled[[field]], blocks[[field]]))
}
stopifnot(identical(scaled$recipe$blocks, blocks$recipe$blocks),
          identical(scaled$recipe$outcome_in_predictors, blocks$recipe$outcome_in_predictors))
metadata <- scaled$recipe$scaling
expected_sds <- rbind(signal = 1:4, intermittent = c(0, 2, 0, 4),
                      shared = rep(0, 4), fixed = rep(0, 4))
colnames(expected_sds) <- as.character(2000:2003)
expect_equal(metadata$estimated_sds, expected_sds)
expected_scales <- c(signal = sqrt(7.5), intermittent = sqrt(5), shared = 1, fixed = 1)
expect_equal(metadata$divisors, matrix(rep(expected_scales, 4), 4, 4,
                                      dimnames = dimnames(expected_sds)))
stopifnot(all(metadata$centers == 0), all(metadata$zero_scale_fallback[3:4, ]),
          !any(metadata$zero_scale_fallback[1:2, ]),
          identical(metadata$periods, blocks$periods), metadata$cutoff == 2003,
          identical(metadata$reference$donors, c("DonorB", "DonorA")),
          identical(metadata$reference$treated, "Target"),
          metadata$reference$mode == "all_units", metadata$diagnostics$fallback_cells == 8,
          metadata$diagnostics$minimum_positive_applied_scale == 1,
          isFALSE(metadata$diagnostics$clipping_applied))
for (date_index in seq_along(blocks$periods)) {
  expect_equal(scaled$X1[[date_index]], blocks$X1[[date_index]] / expected_scales)
  expect_equal(scaled$X0[[date_index]], sweep(blocks$X0[[date_index]], 1, expected_scales, "/"))
}

# zscore uses its own divisor at each date; exactly equal cells become zero.
zscored <- panel.scale(blocks, "zscore")
zmeta <- zscored$recipe$scaling
expect_equal(zmeta$estimated_sds, expected_sds)
expect_equal(unname(zmeta$centers["signal", ]), c(100, 110, 120, 130))
expect_equal(zscored$X1[[1]], c(signal = -1, intermittent = 0, shared = 0, fixed = 0))
expect_equal(zscored$X0[[1]], matrix(c(1, 0, 0, 0, 0, 0, 0, 0), 4, 2,
                                    dimnames = list(predictors, c("DonorB", "DonorA"))))
stopifnot(identical(zscored$Y1, blocks$Y1), identical(zscored$Y0, blocks$Y0),
          identical(zmeta$zero_scale_fallback, zmeta$zero_variation),
          all(zmeta$divisors[zmeta$zero_variation] == 1),
          !identical(zmeta$divisors, metadata$divisors))
for (date_index in seq_along(blocks$periods)) {
  for (predictor in predictors) {
    values <- c(zscored$X1[[date_index]][predictor], zscored$X0[[date_index]][predictor, ])
    if (zmeta$zero_variation[predictor, date_index]) {
      stopifnot(all(values == 0))
    } else {
      expect_equal(mean(values), 0)
      expect_equal(stats::sd(values), 1)
    }
  }
}

# none is exact passthrough; custom scales are aligned by names, not position.
unchanged <- panel.scale(blocks, "none")
stopifnot(identical(unchanged$X1, blocks$X1), identical(unchanged$X0, blocks$X0),
          is.null(unchanged$recipe$scaling$estimated_sds),
          all(unchanged$recipe$scaling$centers == 0),
          all(unchanged$recipe$scaling$divisors == 1))
custom_scales <- c(fixed = 2, shared = 100, intermittent = 10, signal = 20)
custom <- panel.scale(blocks, "custom", custom_scales)
for (date_index in seq_along(blocks$periods)) {
  expect_equal(custom$X1[[date_index]], blocks$X1[[date_index]] / custom_scales[predictors])
}
stopifnot(identical(custom$Y0, blocks$Y0), identical(custom$Y1, blocks$Y1),
          is.null(custom$recipe$scaling$estimated_sds))
for (processed in list(scaled, zscored, unchanged, custom)) {
  expect_error(panel.scale(processed), "Already processed")
}
forged <- scaled
class(forged) <- "tvsc_blocks"
expect_error(panel.scale(forged), "Already processed")

# Common shifts change levels but neither within-date SDs nor matching gaps.
shifted <- panel
shifted$signal <- shifted$signal + 10000 * (shifted$year - 1999)^2
shift_blocks <- panel.blocks(prepare(shifted), 2003)
shift_sd <- panel.scale(shift_blocks)
shift_zscore <- panel.scale(shift_blocks, "zscore")
expect_equal(shift_sd$recipe$scaling$divisors, metadata$divisors)
expect_equal(shift_zscore$X1, zscored$X1)
expect_equal(shift_zscore$X0, zscored$X0)
weights <- c(0.25, 0.75)
for (date_index in seq_along(blocks$periods)) {
  raw_gap <- blocks$X1[[date_index]] - drop(blocks$X0[[date_index]] %*% weights)
  sd_gap <- scaled$X1[[date_index]] - drop(scaled$X0[[date_index]] %*% weights)
  shifted_gap <- shift_sd$X1[[date_index]] - drop(shift_sd$X0[[date_index]] %*% weights)
  z_gap <- zscored$X1[[date_index]] - drop(zscored$X0[[date_index]] %*% weights)
  expect_equal(sd_gap, raw_gap / metadata$divisors[, date_index])
  stopifnot(isTRUE(all.equal(shifted_gap, sd_gap, tolerance = 1e-9)))
  expect_equal(z_gap, raw_gap / zmeta$divisors[, date_index])
}

# Prefix isolation: later corruption/removal cannot enter learned parameters.
later <- prepared
later$training$signal[later$training$year > 2003] <- Inf
later$future_target <- "not a scaling input"
later$future_donors <- NULL
stopifnot(identical(scaled, panel.scale(panel.blocks(later, 2003))))
later$training <- later$training[later$training$year <= 2003, ]
stopifnot(identical(scaled, panel.scale(panel.blocks(later, 2003))))
splits <- panel.split(prepared, 3, 1)
for (fold in splits$folds) {
  prefix <- panel.scale(panel.blocks(prepared, fold$cutoff))
  expect_equal(unname(prefix$recipe$scaling$divisors["signal", 1]),
               sqrt(mean(seq_along(fold$train_period)^2)))
}
full <- panel.scale(panel.blocks(prepared, 2004))
stopifnot(full$recipe$scaling$divisors["signal", 1] != metadata$divisors["signal", 1])

# Single predictors, supplied transformations and numeric IDs retain their shape.
single <- panel.blocks(prepare(matching = "signal"), 2003)
for (method in c("sd", "zscore", "none", "custom")) {
  result <- panel.scale(single, method, if (method == "custom") c(signal = 20) else NULL)
  stopifnot(identical(dim(result$X0[[1]]), c(1L, 2L)),
            identical(dim(result$recipe$scaling$divisors), c(1L, 4L)),
            identical(names(result$X1[[1]]), "signal"))
}
pre_scaled <- panel
pre_scaled$signal <- (pre_scaled$signal - 10) / 20
pre_blocks <- panel.blocks(prepare(pre_scaled), 2003)
stopifnot(identical(panel.scale(pre_blocks, "none")$X1, pre_blocks$X1),
          identical(panel.scale(pre_blocks)$Y0, blocks$Y0))
numeric_panel <- panel
numeric_panel$state <- match(numeric_panel$state, c("Target", "DonorA", "DonorB"))
numeric_prepared <- panel.dataprep(numeric_panel, "state", "year", "sales", 1, c(3, 2),
                                  2005, 2000:2004, "signal")
numeric_scaled <- panel.scale(panel.blocks(numeric_prepared, 2003))
stopifnot(identical(numeric_scaled$recipe$scaling$reference$donors, c(3, 2)),
          identical(colnames(numeric_scaled$X0[[1]]), c("3", "2")))
named_blocks <- panel.blocks(prepare(matching = c(alias = "signal")), 2003)
stopifnot(identical(rownames(panel.scale(named_blocks)$recipe$scaling$divisors), "signal"))

# No clipping of small positive scales; floating-point failures are explicit.
tiny <- single
for (date_index in seq_along(tiny$periods)) {
  tiny$X1[[date_index]][] <- -1e-100
  tiny$X0[[date_index]][, ] <- c(1e-100, 0)
}
tiny_scaled <- panel.scale(tiny)
stopifnot(all(tiny_scaled$recipe$scaling$divisors < 1e-99),
          !any(tiny_scaled$recipe$scaling$zero_scale_fallback))
expect_equal(unname(tiny_scaled$X1[[1]]), -1)
underflow <- tiny
underflow$X1[[1]][] <- -1e-200
underflow$X0[[1]][, ] <- c(1e-200, 0)
expect_error(panel.scale(underflow), "Numerical SD failure")
overflow <- tiny
overflow$X1[[1]][] <- -1e200
overflow$X0[[1]][, ] <- c(1e200, 0)
expect_error(panel.scale(overflow), "Numerical SD failure")
stopifnot(identical(panel.scale(overflow, "none")$X1, overflow$X1))
expect_error(panel.scale(single, "custom", c(signal = 1e-320)), "Nonfinite transformed")

# Unsupported methods and malformed scale vectors never trigger partial matching.
for (bad in list(NULL, NA_character_, "s", "z-score", "mad", "iqr", "minmax", "maxabs",
                 "pooled_sd", c("sd", "none"), TRUE, 1, matrix("sd"))) {
  expect_error(panel.scale(blocks, bad), "scaling must be exactly")
}
for (bad in list(NULL, rep(1, 4), c(signal = 1),
                 c(signal = 1, intermittent = 1, shared = 1, extra = 1),
                 c(signal = 1, signal = 1, shared = 1, fixed = 1),
                 c(signal = 0, intermittent = 1, shared = 1, fixed = 1),
                 c(signal = -1, intermittent = 1, shared = 1, fixed = 1),
                 c(signal = Inf, intermittent = 1, shared = 1, fixed = 1),
                 c(signal = NA_real_, intermittent = 1, shared = 1, fixed = 1),
                 setNames(rep("1", 4), predictors),
                 setNames(rep(1 + 1i, 4), predictors), matrix(1, 4, 1))) {
  expect_error(panel.scale(blocks, "custom", bad), "custom scales must be")
}
for (method in c("sd", "zscore", "none")) {
  expect_error(panel.scale(blocks, method, custom_scales), "scales must be NULL")
}

# Validate schema, date alignment, recipes, shapes and finite values at entry.
expect_error(panel.scale(list()), "schema-version-1")
invalid <- blocks
invalid$schema_version <- 2L
expect_error(panel.scale(invalid), "schema-version-1")
invalid <- blocks
invalid$schema$donors <- NULL
expect_error(panel.scale(invalid), "incomplete or ambiguous")
invalid <- blocks
invalid$periods <- rev(invalid$periods)
expect_error(panel.scale(invalid), "regular pre-treatment")
invalid <- blocks
invalid$cutoff <- 2002
expect_error(panel.scale(invalid), "regular pre-treatment")
invalid <- blocks
invalid$recipe$scaling <- NULL
invalid$recipe <- invalid$recipe[names(invalid$recipe) != "scaling"]
expect_error(panel.scale(invalid), "contemporaneous raw recipe")
invalid <- blocks
invalid$dimensions$donors <- 3L
expect_error(panel.scale(invalid), "dimensions must agree")
invalid <- blocks
names(invalid$X1) <- rev(names(invalid$X1))
expect_error(panel.scale(invalid), "date-named lists")
invalid <- blocks
invalid$X0[[1]] <- invalid$X0[[1]][, 2:1, drop = FALSE]
expect_error(panel.scale(invalid), "Invalid finite numeric predictor")
invalid <- blocks
invalid$X1[[1]][1] <- NA_real_
expect_error(panel.scale(invalid, "none"), "Invalid finite numeric predictor")
invalid <- blocks
invalid$Y0[1, 1] <- Inf
expect_error(panel.scale(invalid), "finite numeric outcomes")
invalid <- blocks
invalid$X1[[1]] <- factor(invalid$X1[[1]])
expect_error(panel.scale(invalid), "Invalid finite numeric predictor")
invalid <- blocks
invalid$X1[[1]] <- invalid$X1[[1]] + 1i
expect_error(panel.scale(invalid), "Invalid finite numeric predictor")

cat("All panel scaling checks passed.\n")
