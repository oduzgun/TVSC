source("R/panel.dataprep.R")
source("R/panel.split.R")
source("R/panel.blocks.R")

panel <- expand.grid(state = c("Target", "DonorA", "DonorB"),
                     year = 2000:2007, stringsAsFactors = FALSE)
panel$sales <- seq_len(nrow(panel)) * 10
panel$income <- 1000 + seq_len(nrow(panel))
panel$constant <- ifelse(panel$year <= 2003, 1, panel$year)
prepare <- function(data = panel, predictors = c("income", "sales", "constant")) {
  panel.dataprep(data, "state", "year", "sales", "Target",
                 c("DonorB", "DonorA"), 2006, 2000:2005, predictors)
}
expect_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}

# Check exact hand-calculated cells, orientation, names and unchanged raw input.
prepared <- prepare()
before <- serialize(prepared, NULL)
blocks <- panel.blocks(prepared, 2003)
stopifnot(inherits(blocks, "tvsc_blocks"), identical(before, serialize(prepared, NULL)),
          identical(blocks$periods, 2000:2003),
          identical(blocks$schema$pre_period, 2000:2003),
          identical(blocks$dimensions, list(predictors = 3L, donors = 2L, periods = 4L)),
          identical(names(blocks$X1), as.character(2000:2003)),
          identical(names(blocks$X0), as.character(2000:2003)),
          identical(blocks$X1[[1]], c(income = 1001, sales = 10, constant = 1)),
          identical(blocks$X0[[1]], matrix(c(1003, 30, 1, 1002, 20, 1), 3, 2,
                    dimnames = list(c("income", "sales", "constant"), c("DonorB", "DonorA")))),
          identical(blocks$Y1, setNames(c(10, 40, 70, 100), as.character(2000:2003))),
          identical(blocks$Y0, matrix(c(30, 20, 60, 50, 90, 80, 120, 110), 4, 2,
                    byrow = TRUE, dimnames = list(as.character(2000:2003), c("DonorB", "DonorA")))),
          is.null(blocks$recipe$scaling), blocks$diagnostics$training_rows == 12L,
          identical(blocks$diagnostics$constant_predictors, "constant"),
          is.null(blocks$future_donors), is.null(blocks$future_target))

# Recover every original cell across dates, not only the first block.
for (date in blocks$periods) {
  target_row <- panel$state == "Target" & panel$year == date
  date_name <- as.character(date)
  for (predictor in prepared$schema$predictors) {
    stopifnot(blocks$X1[[date_name]][[predictor]] == panel[[predictor]][target_row])
    for (donor in prepared$schema$donors) {
      donor_row <- panel$state == donor & panel$year == date
      stopifnot(blocks$X0[[date_name]][predictor, donor] == panel[[predictor]][donor_row])
    }
  }
}

# Predictor labels come from column names, not optional names on the argument vector.
named_predictors <- panel.blocks(prepare(predictors = c(label = "income", other = "sales")), 2003)
stopifnot(identical(names(named_predictors$X1[[1]]), c("income", "sales")),
          identical(rownames(named_predictors$X0[[1]]), c("income", "sales")))

# Reorder stored rows and columns without changing any returned block.
reordered <- prepared
reordered$training <- reordered$training[rev(seq_len(nrow(reordered$training))),
                                         rev(names(reordered$training)), drop = FALSE]
stopifnot(identical(blocks, panel.blocks(reordered, 2003)))

# Ignore later values and future tables, including their absence or corruption.
changed <- prepared
later <- changed$training$year > 2003
changed$training$income[later] <- Inf
changed$training$sales[later] <- NA_real_
changed$training$constant[later] <- -999
changed$future_donors <- NULL
changed$future_target <- "not a fitting input"
stopifnot(identical(blocks, panel.blocks(changed, 2003)))
changed$training <- changed$training[!later, ]
stopifnot(identical(blocks, panel.blocks(changed, 2003)))

# Support both the complete pre-period and a cutoff returned by panel.split().
full <- panel.blocks(prepared, tail(prepared$schema$pre_period, 1))
stopifnot(length(full$X0) == 6L, identical(full$periods, 2000:2005),
          length(full$diagnostics$constant_predictors) == 0L)
splits <- panel.split(prepared, initial = 4, horizon = 2)
stopifnot(identical(blocks, panel.blocks(prepared, splits$folds[[1]]$cutoff)))

# Do not add the outcome to predictor blocks when it is explicitly omitted.
without_outcome <- panel.blocks(prepare(predictors = "income"), 2003)
stopifnot(length(without_outcome$X1[[1]]) == 1L,
          identical(dim(without_outcome$X0[[1]]), c(1L, 2L)),
          identical(rownames(without_outcome$X0[[1]]), "income"),
          identical(without_outcome$Y1, blocks$Y1),
          identical(without_outcome$Y0, blocks$Y0),
          isFALSE(without_outcome$recipe$outcome_in_predictors))
outcome_only <- panel.blocks(prepare(predictors = "sales"), 2003)
stopifnot(identical(dim(outcome_only$X0[[1]]), c(1L, 2L)),
          identical(unname(outcome_only$X0[[1]][1, ]), unname(blocks$Y0[1, ])),
          isTRUE(outcome_only$recipe$outcome_in_predictors))

# Preserve pre-scaled matching columns and keep the reporting outcome separate.
# Fixed illustrative constants avoid learning scales from later observations.
scaled_panel <- panel
scaled_panel$sales_scaled <- (scaled_panel$sales - 50) / 20
scaled_panel$income_scaled <- (scaled_panel$income - 1000) / 100
scaled_prepared <- panel.dataprep(
  scaled_panel, "state", "year", "sales", "Target", c("DonorB", "DonorA"),
  2006, 2000:2005, c("sales_scaled", "income_scaled")
)
scaled_blocks <- panel.blocks(scaled_prepared, 2003)
stopifnot(identical(scaled_blocks$Y1, blocks$Y1),
          identical(scaled_blocks$Y0, blocks$Y0),
          is.null(scaled_blocks$recipe$scaling),
          isFALSE(scaled_blocks$recipe$outcome_in_predictors))
for (date in scaled_blocks$periods) {
  date_name <- as.character(date)
  target_row <- scaled_panel$state == "Target" & scaled_panel$year == date
  for (predictor in scaled_prepared$schema$predictors) {
    stopifnot(scaled_blocks$X1[[date_name]][[predictor]] == scaled_panel[[predictor]][target_row])
    for (donor in scaled_prepared$schema$donors) {
      donor_row <- scaled_panel$state == donor & scaled_panel$year == date
      stopifnot(scaled_blocks$X0[[date_name]][predictor, donor] == scaled_panel[[predictor]][donor_row])
    }
  }
}

# A transformed outcome is also passed through, never automatically inverted.
transformed_outcome <- panel.dataprep(
  scaled_panel, "state", "year", "sales_scaled", "Target", c("DonorB", "DonorA"),
  2006, 2000:2005, "income_scaled"
)
transformed_blocks <- panel.blocks(transformed_outcome, 2003)
stopifnot(identical(transformed_blocks$Y1, (blocks$Y1 - 50) / 20),
          identical(transformed_blocks$Y0, (blocks$Y0 - 50) / 20))

# Numeric and factor IDs preserve declared donor order; spacing is not a row index.
numeric_panel <- panel
numeric_panel$state <- match(numeric_panel$state, c("Target", "DonorA", "DonorB"))
numeric_prepared <- panel.dataprep(numeric_panel, "state", "year", "sales", 1,
                                  c(3, 2), 2006, 2000:2005, "sales")
numeric_blocks <- panel.blocks(numeric_prepared, 2002)
stopifnot(identical(colnames(numeric_blocks$X0[[1]]), c("3", "2")),
          identical(numeric_blocks$schema$donors, c(3, 2)))
factor_panel <- panel
factor_panel$state <- factor(factor_panel$state)
stopifnot(identical(blocks, panel.blocks(prepare(factor_panel), 2003)))
spaced_panel <- panel
spaced_panel$year <- 2000 + 2 * (spaced_panel$year - 2000)
spaced <- panel.dataprep(spaced_panel, "state", "year", "sales", "Target",
                        c("DonorB", "DonorA"), 2012, seq(2000, 2010, 2), "sales", 2)
stopifnot(identical(panel.blocks(spaced, 2004)$periods, c(2000, 2002, 2004)))
expect_error(panel.blocks(spaced, 2005), "must match a date")

# Fail explicitly for invalid cutoffs and malformed object metadata.
expect_error(panel.blocks(prepared), "supplied explicitly")
for (bad in list(NULL, NA_real_, NaN, Inf, "2003", TRUE, c(2002, 2003), 2003.5, matrix(2003))) {
  expect_error(panel.blocks(prepared, bad), "finite integer-valued numeric pre-period date")
}
expect_error(panel.blocks(prepared, 2001), "at least three")
expect_error(panel.blocks(prepared, 2006), "must match a date")
expect_error(panel.blocks(prepared, 4), "must match a date")
expect_error(panel.blocks(list(), 2003), "tvsc_prepared")
invalid <- prepared
invalid$schema_version <- 2L
expect_error(panel.blocks(invalid, 2003), "schema-version-1")
invalid <- prepared
invalid$schema$donors <- NULL
expect_error(panel.blocks(invalid, 2003), "incomplete or ambiguous schema")
invalid <- prepared
invalid$schema$pre_period <- rev(invalid$schema$pre_period)
expect_error(panel.blocks(invalid, 2003), "regular pre-treatment dates")
invalid <- prepared
invalid$recipe$scaling <- "standardize"
expect_error(panel.blocks(invalid, 2003), "contemporaneous raw recipe")

# Revalidate required prefix keys and values, including a non-predictor outcome.
invalid <- prepared
invalid$training <- invalid$training[-1, ]
expect_error(panel.blocks(invalid, 2003), "Missing training key")
invalid <- prepared
invalid$training <- rbind(invalid$training, invalid$training[1, ])
expect_error(panel.blocks(invalid, 2003), "Duplicate training key")
invalid <- prepared
invalid$training$income[1] <- NA_real_
expect_error(panel.blocks(invalid, 2003), "Nonfinite training value")
invalid <- prepare(predictors = "income")
invalid$training$sales[1] <- Inf
expect_error(panel.blocks(invalid, 2003), "Nonfinite training value")
invalid <- prepared
invalid$training$income <- factor(invalid$training$income)
expect_error(panel.blocks(invalid, 2003), "must be numeric")

cat("All raw predictor block checks passed.\n")
