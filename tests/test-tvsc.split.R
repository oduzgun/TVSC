source("R/tvsc.dataprep.R")
source("R/tvsc.split.R")

panel <- expand.grid(state = c("Target", "DonorA", "DonorB"),
                     year = 2000:2009, stringsAsFactors = FALSE)
panel$sales <- seq_len(nrow(panel)) * 10
prepare <- function(data = panel) {
  tvsc.dataprep(data, "state", "year", "sales", "Target",
                 c("DonorB", "DonorA"), 2008, 2000:2007, "sales")
}
expect_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}

# Reproduce the agreed eight-date example and preserve the original object.
prepared <- prepare()
before <- serialize(prepared, NULL)
splits <- tvsc.split(prepared, initial = 4, horizon = 2)
stopifnot(inherits(splits, "tvsc_splits"), identical(splits$schema, prepared$schema),
          identical(before, serialize(prepared, NULL)), length(splits$folds) == 3L,
          identical(splits$folds[[1]]$train_period, 2000:2003),
          identical(splits$folds[[1]]$validation_period, 2004:2005),
          identical(splits$folds[[2]]$train_period, 2000:2004),
          identical(splits$folds[[2]]$validation_period, 2005:2006),
          identical(splits$folds[[3]]$validation_period, 2006:2007),
          identical(unname(splits$diagnostics$assessment_counts), c(0L, 0L, 0L, 0L, 1L, 2L, 2L, 1L)),
          identical(splits$diagnostics$unassessed_period, 2000:2003))

# Dates apply to all units, and horizons are measured from each training cutoff.
for (fold in splits$folds) {
  training <- prepared$training[prepared$training$year %in% fold$train_period, ]
  validation <- prepared$training[prepared$training$year %in% fold$validation_period, ]
  stopifnot(all(table(training$year) == 3L), all(table(validation$year) == 3L),
            max(fold$train_period) == fold$cutoff,
            all(fold$validation_period > fold$cutoff),
            length(intersect(fold$train_index, fold$validation_index)) == 0L,
            all(fold$validation_period < prepared$schema$intervention_time),
            identical(fold$horizons, 1:2),
            identical(fold$train_period, prepared$schema$pre_period[fold$train_index]))
}

# Full-window and stride rules also handle a single fold and unassessed tails.
single <- tvsc.split(prepared, initial = 6, horizon = 2)
stopifnot(length(single$folds) == 1L, identical(single$folds[[1]]$validation_period, 2006:2007))
stride <- tvsc.split(prepared, initial = 3, horizon = 2, step = 2)
stopifnot(length(stride$folds) == 2L,
          identical(stride$folds[[2]]$validation_period, 2005:2006),
          identical(stride$diagnostics$unassessed_period, c(2000:2002, 2007L)))
large_step <- tvsc.split(prepared, initial = 3, horizon = 1, step = 100)
stopifnot(length(large_step$folds) == 1L)
gaps <- tvsc.split(prepared, initial = 3, horizon = 1, step = 2)
stopifnot(identical(gaps$diagnostics$unassessed_period, c(2000:2002, 2004L, 2006L)))

# Numeric index spacing is separate from the number of dates per window.
spaced_panel <- panel
spaced_panel$year <- 2000 + 2 * (spaced_panel$year - 2000)
spaced <- tvsc.dataprep(spaced_panel, "state", "year", "sales", "Target",
                        c("DonorB", "DonorA"), 2016, seq(2000, 2014, 2), "sales",
                        period_step = 2)
spaced_splits <- tvsc.split(spaced, 4, 2)
stopifnot(identical(spaced_splits$folds[[1]]$validation_period, c(2008, 2010)),
          identical(spaced_splits$folds[[1]]$horizons, 1:2),
          all((spaced_splits$folds[[1]]$validation_period - spaced_splits$folds[[1]]$cutoff) /
              spaced$schema$period_step == 1:2))

# Splits depend on declared metadata, not row order or outcome values.
stopifnot(identical(splits, tvsc.split(prepare(panel[rev(seq_len(nrow(panel))), ]), 4, 2)))
changed <- panel
changed$sales <- changed$sales + 100
changed$sales[changed$year >= 2008] <- NA_real_
stopifnot(identical(splits, tvsc.split(prepare(changed), 4, 2)),
          identical(splits, tvsc.split(prepare(panel[panel$year < 2008, ]), 4, 2)))

# Reject missing, nonscalar, nonnumeric, nonfinite and infeasible window inputs.
expect_error(tvsc.split(prepared, horizon = 2), "supplied explicitly")
expect_error(tvsc.split(prepared, initial = 4), "supplied explicitly")
for (bad in list(NULL, NA_real_, NaN, Inf, -1, 0, 1.5, "4", TRUE, c(3, 4), matrix(4))) {
  expect_error(tvsc.split(prepared, bad, 2), "initial must be one positive")
  expect_error(tvsc.split(prepared, 4, bad), "horizon must be one positive")
  expect_error(tvsc.split(prepared, 4, 2, bad), "step must be one positive")
}
expect_error(tvsc.split(prepared, 2, 2), "at least three")
expect_error(tvsc.split(prepared, 8, 1), "No full fold fits")
expect_error(tvsc.split(prepared, 4, 5), "No full fold fits")
short <- tvsc.dataprep(panel, "state", "year", "sales", "Target",
                       c("DonorB", "DonorA"), 2008, 2000:2002, "sales")
expect_error(tvsc.split(short, 3, 1), "No full fold fits")

# Reject unsupported object versions and invalid date metadata.
expect_error(tvsc.split(list(), 4, 2), "tvsc_prepared")
invalid <- prepared
invalid$schema_version <- 1L
expect_error(tvsc.split(invalid, 4, 2), "schema-version-2")
expect_error(tvsc.split(invalid, 4, 2), "original panel")
invalid$schema_version <- 3L
expect_error(tvsc.split(invalid, 4, 2), "schema-version-2")
invalid <- prepared
invalid$schema$donors <- NULL
expect_error(tvsc.split(invalid, 4, 2), "incomplete or ambiguous schema")
invalid <- prepared
invalid$schema$pre_period <- rev(invalid$schema$pre_period)
expect_error(tvsc.split(invalid, 4, 2), "regular pre-treatment dates")
invalid <- prepared
invalid$schema$pre_period[1] <- NA_real_
expect_error(tvsc.split(invalid, 4, 2), "regular pre-treatment dates")
invalid <- prepared
invalid$schema$intervention_time <- 2007
expect_error(tvsc.split(invalid, 4, 2), "regular pre-treatment dates")
invalid <- prepared
invalid$schema$period_step <- 0
expect_error(tvsc.split(invalid, 4, 2), "regular pre-treatment dates")

stopifnot(identical(prepared$schema_version, 2L),
          identical(splits$schema_version, 1L))
invalid <- prepared
invalid$blocks <- NULL
expect_error(tvsc.split(invalid, 4, 2), "version-2 training structure")
invalid <- prepared
invalid$blocks$schema_version <- 1L
expect_error(tvsc.split(invalid, 4, 2), "incompatible raw blocks")
invalid <- prepared
invalid$blocks$periods <- rev(invalid$blocks$periods)
expect_error(tvsc.split(invalid, 4, 2), "incompatible raw blocks")
invalid <- prepared
invalid$blocks$schema$donors <- rev(invalid$blocks$schema$donors)
expect_error(tvsc.split(invalid, 4, 2), "incompatible raw blocks")
invalid <- prepared
invalid$blocks$X0[[1]] <- invalid$blocks$X0[[1]][, 2:1, drop = FALSE]
expect_error(tvsc.split(invalid, 4, 2), "dimensions or ordering")
invalid <- prepared
invalid$blocks$X1[[1]] <- unname(invalid$blocks$X1[[1]])
expect_error(tvsc.split(invalid, 4, 2), "dimensions or ordering")
invalid <- prepared
invalid$blocks$X0[[1]] <- drop(invalid$blocks$X0[[1]])
expect_error(tvsc.split(invalid, 4, 2), "dimensions or ordering")
invalid <- prepared
invalid$blocks$Y0 <- invalid$blocks$Y0[, 2:1, drop = FALSE]
expect_error(tvsc.split(invalid, 4, 2), "incompatible raw blocks")
invalid <- prepared
invalid$blocks$dimensions$predictors <- 2L
expect_error(tvsc.split(invalid, 4, 2), "incompatible raw blocks")
invalid <- prepared
invalid$recipe$scaling <- list(method = "none")
expect_error(tvsc.split(invalid, 4, 2), "raw recipe")
invalid <- prepared
invalid$schema$donors <- c("DonorA", "DonorA")
expect_error(tvsc.split(invalid, 4, 2), "target/donor metadata")
invalid <- prepared
invalid$training <- invalid$training[-1, ]
expect_error(tvsc.split(invalid, 4, 2), "training structure")
values_ignored <- prepared
values_ignored$training$sales[] <- NA_real_
values_ignored$blocks$X1 <- lapply(values_ignored$blocks$X1, function(block) {
  block[] <- NA_real_
  block
})
values_ignored$blocks$X0 <- lapply(values_ignored$blocks$X0, function(block) {
  block[] <- Inf
  block
})
values_ignored$blocks$Y1[] <- NaN
values_ignored$blocks$Y0[] <- NA_real_
values_ignored$future_target$sales[] <- NA_real_
stopifnot(identical(splits, tvsc.split(values_ignored, 4, 2)))
cat("All chronological split checks for version-2 preparation passed.\n")