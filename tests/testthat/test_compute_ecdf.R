test_that('we require at least minSamples', {

    x <- c(1,2,3,4,5,6,7,8,9)
    f <- compute_ecdf(x, 0)
    expect_true(is.na(f(4)))

    x <- c(1,2,3,4,5,6,7,8,9,10)
    f <- compute_ecdf(x, 0)
    expect_equal(f(4), 0.6666666, tolerance = 1e-7)

    x <- c(1,2,3,4,5,6,7,8,9,10)
    f <- compute_ecdf(x, 0, minSamples = 11)
    expect_true(is.na(f(4)))

})

test_that('we require at least minSamples with minV', {

    x <- c(1,2,3,4,5,6,7,8,9,10)
    f <- compute_ecdf(x, 2)
    expect_true(is.na(f(4)))

})
