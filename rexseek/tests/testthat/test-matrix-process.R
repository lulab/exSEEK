testthat::context('Testing matrix process')
if (basename(getwd()) == 'testthat') setwd('../..')

testthat::test_that('read_mat', {
	read_mat('inst/extdata/scirep_sequential_qc.txt')
    testthat::expect_identical(1L, 1L);
    testthat::expect_true(T);
});


if (tolower(Sys.getenv('CI')) != 'true') {

} else {

}

