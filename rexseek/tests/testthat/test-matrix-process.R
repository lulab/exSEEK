testthat::context('Testing matrix process')
if (basename(getwd()) == 'testthat') setwd('../..')

testthat::test_that('read_mat()', {
	mat_raw <- read_mat('inst/extdata/scirep_sequential_qc.txt', n_max = 1000)

	testthat::expect_true(is.matrix(mat_raw))
	testthat::expect_identical(typeof(mat_raw), 'integer')
	testthat::expect_identical(ncol(mat_raw), 191L)
});

testthat::test_that('filter_low()', {
	mat <- read_mat('inst/extdata/scirep_sequential_qc.txt', n_max = 1000) %>% filter_low()

	testthat::expect_identical(dim(mat), c(100L, 191L))
});


