testthat::context('Testing foobar')
if (basename(getwd()) == 'testthat') setwd('../..')

testthat::test_that('something', {
    testthat::expect_identical(1L, 1L);
    testthat::expect_true(T);
});


if (tolower(Sys.getenv('CI')) != 'true') {

} else {
    
}

