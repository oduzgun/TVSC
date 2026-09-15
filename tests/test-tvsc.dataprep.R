source("R/tvsc.dataprep.R")
panel <- expand.grid(state = c("Target", "DonorB", "DonorA"),
                     year = 2000:2004, stringsAsFactors = FALSE)
panel$sales <- seq_len(nrow(panel)) * 10
panel$income <- 1000 + seq_len(nrow(panel))
panel$constant <- 1
prepare <- function(data = panel) {
  tvsc.dataprep(data, "state", "year", "sales", "Target",
                 c("DonorA", "DonorB"), 2003, 2000:2002,
                 c("sales", "income", "constant"))
}
expect_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}
prepared <- prepare()
stopifnot(identical(prepared$schema_version, 2L),
          inherits(prepared$blocks, "tvsc_blocks"),
          identical(prepared$blocks$schema_version, 2L),
          identical(prepared$blocks$schema, prepared$schema),
          identical(prepared$blocks$recipe, prepared$recipe),
          identical(prepared$blocks$periods, 2000:2002),
          identical(prepared$blocks$cutoff, 2002L),
          identical(prepared$blocks$dimensions,
                    list(predictors = 3L, donors = 2L, periods = 3L)),
          identical(prepared$blocks$diagnostics$constant_predictors, "constant"),
          identical(names(prepared$blocks$X1), as.character(2000:2002)),
          identical(names(prepared$blocks$X0), as.character(2000:2002)))
for (period in prepared$schema$pre_period) {
  date_name <- as.character(period)
  target_row <- panel$state == "Target" & panel$year == period
  expected_target <- vapply(prepared$schema$predictors, function(predictor) {
    panel[[predictor]][target_row]
  }, numeric(1))
  stopifnot(identical(prepared$blocks$X1[[date_name]], expected_target),
            identical(unname(prepared$blocks$Y1[date_name]), panel$sales[target_row]))
  for (donor in prepared$schema$donors) {
    donor_row <- panel$state == donor & panel$year == period
    expected_donor <- vapply(prepared$schema$predictors, function(predictor) {
      panel[[predictor]][donor_row]
    }, numeric(1))
    stopifnot(identical(prepared$blocks$X0[[date_name]][, donor], expected_donor),
              identical(unname(prepared$blocks$Y0[date_name, donor]), panel$sales[donor_row]))
  }
}
stopifnot(identical(unname(prepared$blocks$X0[[1]]["sales", ]), c(30, 20)),
          identical(unname(prepared$blocks$Y1), c(10, 40, 70)))
stopifnot(inherits(prepared, "tvsc_prepared"), nrow(prepared$training) == 9L,
          identical(prepared$schema$donors, c("DonorA", "DonorB")),
          identical(prepared$training$state, rep(c("Target", "DonorA", "DonorB"), 3)),
          identical(prepared$diagnostics$constant_predictors, "constant"),
          is.null(prepared$recipe$scaling), is.null(prepared$splits),
          isFALSE(prepared$diagnostics$future_outcomes_validated))
reordered <- prepare(panel[rev(seq_len(nrow(panel))), ])
stopifnot(identical(prepared$training, reordered$training),
          identical(prepared$blocks, reordered$blocks))
changed <- panel
changed$sales[changed$year >= 2003] <- NA_real_
changed$income[changed$year >= 2003] <- Inf
stopifnot(identical(prepared$training, prepare(changed)$training),
          identical(prepared$blocks, prepare(changed)$blocks))
without_future <- prepare(panel[panel$year < 2003, ])
stopifnot(identical(prepared$training, without_future$training),
          identical(prepared$blocks, without_future$blocks),
          nrow(without_future$future_donors) == 0L, nrow(without_future$future_target) == 0L)
without_target <- prepare(panel[!(panel$year >= 2003 & panel$state == "Target"), ])
stopifnot(nrow(without_target$future_target) == 0L, nrow(without_target$future_donors) == 4L,
          identical(prepared$blocks, without_target$blocks))
expect_error(prepare(rbind(panel, panel[1, ])), "Duplicate training key")
expect_error(prepare(panel[-1, ]), "Missing training key")
invalid <- panel
invalid$income[1] <- NA_real_
expect_error(prepare(invalid), "Nonfinite training value")
invalid <- panel
invalid$income <- factor(invalid$income)
expect_error(prepare(invalid), "must be numeric")
invalid <- panel
invalid$year[1] <- 2000.5
expect_error(prepare(invalid), "integer-valued numeric period index")
future_duplicates <- prepare(rbind(panel, panel[panel$year == 2004, ]))
stopifnot(identical(prepared$training, future_duplicates$training),
          identical(prepared$blocks, future_duplicates$blocks),
          isFALSE(future_duplicates$diagnostics$future_outcomes_validated))
single <- tvsc.dataprep(panel, "state", "year", "sales", "Target",
                        c("DonorA", "DonorB"), 2003, 2000:2002, "sales")
stopifnot(identical(single$schema$predictors, "sales"), isTRUE(single$recipe$outcome_in_predictors),
          identical(dim(single$blocks$X0[[1]]), c(1L, 2L)),
          identical(names(single$blocks$X1[[1]]), "sales"),
          identical(single$blocks$Y1, prepared$blocks$Y1),
          identical(single$blocks$Y0, prepared$blocks$Y0))
expect_error(tvsc.dataprep(panel, "state", "year", "sales", "Target",
                           c("Target", "DonorA"), 2003, 2000:2002, "sales"), "own donors")
expect_error(tvsc.dataprep(panel, "state", "year", "sales", "Target",
                           c("DonorA", "DonorB"), 2003, c(2000, 2002, 2004), "sales"), "period_step spacing")
expect_error(tvsc.dataprep(panel, "state", "year", "sales", "Target",
                           c("DonorA", "DonorB"), 2002, 2000:2002, "sales"), "precede intervention_time")
expect_error(tvsc.dataprep(panel, "state", "year", "sales", "Target",
                           c("DonorA", "DonorB"), 2003, 2000:2002, c("sales", "sales")), "unique predictor column")
numeric_panel <- panel
numeric_panel$state <- match(numeric_panel$state, c("Target", "DonorA", "DonorB"))
numeric_prepared <- tvsc.dataprep(numeric_panel, "state", "year", "sales", 1,
                                 c(3, 2), 2003, 2000:2002, "sales")
stopifnot(identical(numeric_prepared$schema$donors, c(3, 2)),
          identical(colnames(numeric_prepared$blocks$X0[[1]]), c("3", "2")),
          identical(unname(numeric_prepared$blocks$Y0),
                    unname(prepared$blocks$Y0[, c("DonorB", "DonorA")])))
expect_error(tvsc.dataprep(numeric_panel, "state", "year", "sales", "1",
                           c("2", "3"), 2003, 2000:2002, "sales"), "character/numeric type")
covariates <- tvsc.dataprep(panel, "state", "year", "sales", "Target",
                           c("DonorA", "DonorB"), 2003, 2000:2002,
                           c("constant", "income"))
stopifnot(identical(names(covariates$blocks$X1[[1]]), c("constant", "income")),
          identical(rownames(covariates$blocks$X0[[1]]), c("constant", "income")),
          isFALSE(covariates$recipe$outcome_in_predictors),
          identical(covariates$blocks$Y1, prepared$blocks$Y1))
transformed <- panel
transformed$income <- transformed$income / 100
transformed$sales <- transformed$sales / 10
transformed_prepared <- prepare(transformed)
stopifnot(identical(unname(transformed_prepared$blocks$X1[[1]]["income"]), 10.01),
          identical(unname(transformed_prepared$blocks$Y1), c(1, 4, 7)),
          is.null(transformed_prepared$blocks$recipe$scaling))
factor_panel <- panel
factor_panel$state <- factor(factor_panel$state)
stopifnot(identical(prepare(factor_panel)$blocks, prepared$blocks))
spaced_panel <- panel
spaced_panel$year <- 2000 + 2 * (spaced_panel$year - 2000)
spaced <- tvsc.dataprep(spaced_panel, "state", "year", "sales", "Target",
                       c("DonorA", "DonorB"), 2006, c(2000, 2002, 2004),
                       "sales", period_step = 2)
stopifnot(identical(names(spaced$blocks$X0), c("2000", "2002", "2004")),
          identical(unname(spaced$blocks$Y0), unname(single$blocks$Y0)))
full <- tvsc.dataprep(panel, "state", "year", "sales", "Target",
                     c("DonorA", "DonorB"), 2005, 2000:2004,
                     c("sales", "income", "constant"))
prefix <- tvsc.dataprep(full$training, "state", "year", "sales", "Target",
                       c("DonorA", "DonorB"), 2005, 2000:2002,
                       c("sales", "income", "constant"))
changed_training <- full$training
changed_training$sales[changed_training$year > 2002] <- NA_real_
changed_training$income[changed_training$year > 2002] <- Inf
changed_prefix <- tvsc.dataprep(changed_training, "state", "year", "sales", "Target",
                               c("DonorA", "DonorB"), 2005, 2000:2002,
                               c("sales", "income", "constant"))
stopifnot(identical(prefix$blocks, changed_prefix$blocks),
          identical(prefix$blocks$X0, full$blocks$X0[1:3]))
isolated <- new.env(parent = baseenv())
sys.source("R/tvsc.dataprep.R", envir = isolated)
isolated_prepared <- isolated$tvsc.dataprep(panel, "state", "year", "sales", "Target",
                                          c("DonorA", "DonorB"), 2003, 2000:2002,
                                          c("sales", "income", "constant"))
stopifnot(identical(isolated_prepared, prepared),
          !exists("tvsc.blocks", envir = isolated, inherits = FALSE))
cat("All version-2 raw panel preparation checks passed.\n")