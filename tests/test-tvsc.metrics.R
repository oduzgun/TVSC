source("R/tvsc.dataprep.R")
source("R/tvsc.split.R")
source("R/tvsc.blocks.R")
source("R/tvsc.scale.R")
source("R/tvsc.metrics.R")

expect_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}
expect_equal <- function(actual, expected) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-12)))
}

# Build an explicit prefix; uniform weights do not depend on its data values.
panel <- expand.grid(state = c("Target", "DonorA", "DonorB"),
                     year = 2000:2005, stringsAsFactors = FALSE)
panel$sales <- seq_len(nrow(panel)) * 10
panel$income <- 100 + seq_len(nrow(panel))^2
panel$price <- rep(c(2, 5, 8), 6)
predictors <- c("sales", "income", "price")
prepared <- tvsc.dataprep(panel, "state", "year", "sales", "Target",
                           c("DonorB", "DonorA"), 2005, 2000:2004, predictors)
raw <- tvsc.blocks(prepared, 2003)
scaled <- tvsc.scale(raw)
before <- serialize(scaled, NULL)
uniform <- tvsc.metrics(scaled)
expected <- matrix(1 / 3, 3, 4, dimnames = list(predictors, as.character(2000:2003)))
expect_equal(uniform$weights, expected)
stopifnot(identical(before, serialize(scaled, NULL)),
          identical(class(uniform), "tvsc_metrics"), uniform$schema_version == 1L,
          uniform$mode == "uniform", identical(uniform$schema, scaled$schema),
          identical(uniform$scaling, scaled$recipe$scaling),
          identical(uniform$periods, scaled$periods), uniform$cutoff == 2003,
          uniform$diagnostics$constant_over_time,
          !uniform$diagnostics$normalization_applied,
          !uniform$provenance$verified, is.null(uniform$provenance$description),
          !any(c("X1", "X0", "Y1", "Y0", "future_target") %in% names(uniform)))

# Named vectors and matrices are aligned without changing supplied values.
supplied <- c(price = 0.3, sales = 0.5, income = 0.2)
constant <- tvsc.metrics(scaled, supplied, provenance = "Externally fixed example weights")
expected_constant <- matrix(rep(supplied[predictors], 4), 3, 4, dimnames = dimnames(expected))
stopifnot(identical(constant$weights, expected_constant),
          identical(constant$weights, tvsc.metrics(scaled, supplied[predictors])$weights),
          constant$mode == "supplied_constant", constant$diagnostics$constant_over_time,
          constant$provenance$description == "Externally fixed example weights",
          !constant$provenance$verified)
varying <- matrix(c(0.6, 0.3, 0.1, 0.4, 0.4, 0.2, 0.2, 0.5, 0.3, 0.5, 0.2, 0.3),
                  3, 4, dimnames = dimnames(expected))
dynamic <- tvsc.metrics(scaled, varying[3:1, 4:1])
stopifnot(identical(dynamic$weights, varying), dynamic$mode == "supplied_time_varying",
          !dynamic$diagnostics$constant_over_time)
constant_matrix <- tvsc.metrics(scaled, expected_constant)
stopifnot(constant_matrix$mode == "supplied_time_varying",
          constant_matrix$diagnostics$constant_over_time)
zero <- tvsc.metrics(scaled, c(sales = 1, income = 0, price = 0))
stopifnot(identical(zero$diagnostics$zero_weight_predictors, c("income", "price")))

# Diagonal metric arithmetic matches the predictor-weighted loss.
donor_weights <- c(0.25, 0.75)
for (date_index in seq_along(scaled$periods)) {
  residual <- scaled$X1[[date_index]] - drop(scaled$X0[[date_index]] %*% donor_weights)
  diagonal <- diag(dynamic$weights[, date_index], nrow = length(predictors))
  expect_equal(drop(crossprod(residual, diagonal %*% residual)),
               sum(dynamic$weights[, date_index] * residual^2))
}
for (method in c("sd", "zscore", "none", "custom")) {
  processed <- tvsc.scale(raw, method,
                          if (method == "custom") c(sales = 10, income = 20, price = 2) else NULL)
  metrics <- tvsc.metrics(processed, supplied)
  stopifnot(identical(metrics$weights, constant$weights),
            identical(metrics$scaling, processed$recipe$scaling))
}
metadata_only <- scaled
metadata_only$X1 <- "not inspected"
metadata_only$X0 <- NULL
metadata_only["X0"] <- list(NULL)
metadata_only$Y1 <- Inf
metadata_only$Y0 <- "not inspected"
stopifnot(identical(tvsc.metrics(metadata_only), uniform))

# Every fold has its own exact date domain; extra dates are not silently sliced.
splits <- tvsc.split(prepared, 3, 1)
for (fold in splits$folds) {
  prefix <- tvsc.scale(tvsc.blocks(prepared, fold$cutoff))
  metrics <- tvsc.metrics(prefix)
  stopifnot(identical(colnames(metrics$weights), as.character(fold$train_period)))
  prefix_weights <- matrix(1 / 3, 3, length(fold$train_period),
                           dimnames = list(predictors, as.character(fold$train_period)))
  expect_equal(tvsc.metrics(prefix, prefix_weights)$weights, metrics$weights)
}
full <- tvsc.scale(tvsc.blocks(prepared, 2004))
expect_error(tvsc.metrics(scaled, tvsc.metrics(full)$weights), "training dates exactly once")

# Single-predictor dimensions and numeric donor identities are retained.
numeric_panel <- panel
numeric_panel$state <- match(numeric_panel$state, c("Target", "DonorA", "DonorB"))
single_prepared <- tvsc.dataprep(numeric_panel, "state", "year", "sales", 1,
                                  c(3, 2), 2005, 2000:2004, c(alias = "sales"))
single <- tvsc.scale(tvsc.blocks(single_prepared, 2003))
for (input in list("uniform", c(sales = 1),
                   matrix(1, 1, 4, dimnames = list("sales", as.character(2000:2003))))) {
  metrics <- tvsc.metrics(single, input)
  stopifnot(identical(dim(metrics$weights), c(1L, 4L)), all(metrics$weights == 1),
            identical(metrics$schema$donors, c(3, 2)))
}

# Numerical tolerance is explicit and never changes admitted weights.
near <- c(sales = 1 + 5e-9, income = 0, price = 0)
near_metrics <- tvsc.metrics(scaled, near)
stopifnot(identical(unname(near_metrics$weights[, 1]), unname(near[predictors])),
          near_metrics$diagnostics$maximum_unit_sum_error > 0)
expect_error(tvsc.metrics(scaled, near, tolerance = 0), "within tolerance")
expect_error(tvsc.metrics(scaled, c(sales = 0, income = 0, price = 0), tolerance = 0.99),
             "positive column sums")
for (bad in list(NULL, NA_real_, Inf, -1, 1, TRUE, "0.01", c(0, 0.01), matrix(0), 1i)) {
  expect_error(tvsc.metrics(scaled, tolerance = bad), "tolerance must be")
}
for (bad in list("", NA_character_, c("source", "other"), 1, list(source = "external"))) {
  expect_error(tvsc.metrics(scaled, provenance = bad), "provenance must be")
}
for (bad in list(NULL, "u", "constant", TRUE, list(varying),
                 c(sales = Inf, income = 0, price = 0),
                 c(sales = NA_real_, income = 0, price = 0),
                 supplied + 1i, array(1, c(3, 4, 1)))) {
  expect_error(tvsc.metrics(scaled, bad), "predictor_weights must be")
}
for (bad in list(unname(supplied), c(sales = 1),
                 c(sales = 1, sales = 0, price = 0), c(sales = 1, income = 0, extra = 0))) {
  expect_error(tvsc.metrics(scaled, bad), "vector names")
}
for (bad in list(c(sales = -1e-12, income = 1, price = 0), supplied * 2,
                 c(sales = 1e308, income = 1e308, price = 1e308))) {
  expect_error(tvsc.metrics(scaled, bad), "within tolerance")
}
for (bad in list(unname(varying), varying[, 1:3], varying[c(1, 1, 3), ],
                 varying[, c(1, 1, 3, 4)], t(varying))) {
  expect_error(tvsc.metrics(scaled, bad), "matrix names")
}
invalid_weights <- varying
invalid_weights[1, 2] <- -0.01
expect_error(tvsc.metrics(scaled, invalid_weights), "nonnegative")

# Invalid metadata is rejected without pretending to authenticate source data.
expect_error(tvsc.metrics(raw), "tvsc_scaled_blocks")
expect_error(tvsc.metrics(list()), "tvsc_scaled_blocks")
invalid <- scaled
invalid$schema_version <- 2L
expect_error(tvsc.metrics(invalid), "schema-version-1")
invalid <- scaled
invalid$schema$donors <- NULL
expect_error(tvsc.metrics(invalid), "incomplete or ambiguous")
invalid <- scaled
invalid$periods <- rev(invalid$periods)
expect_error(tvsc.metrics(invalid), "regular pre-treatment")
invalid <- scaled
invalid$cutoff <- 2002
expect_error(tvsc.metrics(invalid), "matching cutoff")
invalid <- scaled
invalid$dimensions$predictors <- 2L
expect_error(tvsc.metrics(invalid), "dimensions must agree")
invalid <- scaled
invalid$recipe$scaling <- NULL
stopifnot(!"scaling" %in% names(invalid$recipe))
expect_error(tvsc.metrics(invalid), "contemporaneous recipe with explicit scaling")
invalid <- scaled
invalid$recipe["scaling"] <- list(NULL)
stopifnot("scaling" %in% names(invalid$recipe), is.null(invalid$recipe$scaling))
expect_error(tvsc.metrics(invalid), "scaling metadata")
invalid <- scaled
invalid$recipe$scaling$reference$donors <- rev(invalid$schema$donors)
expect_error(tvsc.metrics(invalid), "scaling metadata")
invalid <- scaled
invalid$recipe$scaling$divisors[1, 1] <- 0
expect_error(tvsc.metrics(invalid), "scaling metadata")
invalid <- scaled
invalid$recipe$scaling$divisors[1, 1] <- invalid$recipe$scaling$divisors[1, 1] * 2
expect_error(tvsc.metrics(invalid), "contradicts its method")
invalid <- scaled
invalid$recipe$scaling$centers[1, 1] <- 1
expect_error(tvsc.metrics(invalid), "contradicts its method")
invalid <- scaled
colnames(invalid$recipe$scaling$divisors) <- rev(colnames(invalid$recipe$scaling$divisors))
expect_error(tvsc.metrics(invalid), "scaling metadata")

cat("All panel metrics checks passed.\n")
