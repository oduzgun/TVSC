source("R/panel.dataprep.R")
panel <- expand.grid(state = c("Target", "DonorB", "DonorA"),
                     year = 2000:2004, stringsAsFactors = FALSE)
panel$sales <- seq_len(nrow(panel)) * 10
panel$income <- 1000 + seq_len(nrow(panel))
panel$constant <- 1
prepare <- function(data = panel) {
  panel.dataprep(data, "state", "year", "sales", "Target",
                 c("DonorA", "DonorB"), 2003, 2000:2002,
                 c("sales", "income", "constant"))
}
expect_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}
prepared <- prepare()
stopifnot(inherits(prepared, "tvsc_prepared"), nrow(prepared$training) == 9L,
          identical(prepared$schema$donors, c("DonorA", "DonorB")),
          identical(prepared$training$state, rep(c("Target", "DonorA", "DonorB"), 3)),
          identical(prepared$diagnostics$constant_predictors, "constant"),
          is.null(prepared$recipe$scaling), is.null(prepared$splits),
          isFALSE(prepared$diagnostics$future_outcomes_validated))
reordered <- prepare(panel[rev(seq_len(nrow(panel))), ])
stopifnot(identical(prepared$training, reordered$training))
changed <- panel
changed$sales[changed$year >= 2003] <- NA_real_
changed$income[changed$year >= 2003] <- Inf
stopifnot(identical(prepared$training, prepare(changed)$training))
without_future <- prepare(panel[panel$year < 2003, ])
stopifnot(identical(prepared$training, without_future$training),
          nrow(without_future$future_donors) == 0L, nrow(without_future$future_target) == 0L)
without_target <- prepare(panel[!(panel$year >= 2003 & panel$state == "Target"), ])
stopifnot(nrow(without_target$future_target) == 0L, nrow(without_target$future_donors) == 4L)
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
          isFALSE(future_duplicates$diagnostics$future_outcomes_validated))
single <- panel.dataprep(panel, "state", "year", "sales", "Target",
                        c("DonorA", "DonorB"), 2003, 2000:2002, "sales")
stopifnot(identical(single$schema$predictors, "sales"), isTRUE(single$recipe$outcome_in_predictors))
expect_error(panel.dataprep(panel, "state", "year", "sales", "Target",
                           c("Target", "DonorA"), 2003, 2000:2002, "sales"), "own donors")
expect_error(panel.dataprep(panel, "state", "year", "sales", "Target",
                           c("DonorA", "DonorB"), 2003, c(2000, 2002, 2004), "sales"), "period_step spacing")
expect_error(panel.dataprep(panel, "state", "year", "sales", "Target",
                           c("DonorA", "DonorB"), 2002, 2000:2002, "sales"), "precede intervention_time")
expect_error(panel.dataprep(panel, "state", "year", "sales", "Target",
                           c("DonorA", "DonorB"), 2003, 2000:2002, c("sales", "sales")), "unique predictor column")
numeric_panel <- panel
numeric_panel$state <- match(numeric_panel$state, c("Target", "DonorA", "DonorB"))
numeric_prepared <- panel.dataprep(numeric_panel, "state", "year", "sales", 1,
                                 c(3, 2), 2003, 2000:2002, "sales")
stopifnot(identical(numeric_prepared$schema$donors, c(3, 2)))
expect_error(panel.dataprep(numeric_panel, "state", "year", "sales", "1",
                           c("2", "3"), 2003, 2000:2002, "sales"), "character/numeric type")
cat("All raw panel preparation checks passed.\n")
