expression <- read.table(text = "
                         transcript s1 s2 s3
                         A1 1 1 5
                         A2 1 1 5
                         A3 1 1 5
                         B1 3 NA 5
                         B2 2 NA 5
                         B3 1 NA 5
                         C1 NA 1 5
                         C2 NA 2 5
                         C3 NA 3 5
                         D1 3 1 5
                         D2 2 2 5
                         D3 1 3 5
                         ", header = TRUE, stringsAsFactors = FALSE)

test_that('output is as expected', {

    long <- tpm2long(expression, measured)

    # size is expected
    expect_true(is.data.frame(long))
    expect_equal(ncol(long), 3)
    expect_equal(nrow(long), nrow(expression) * (ncol(expression) - 1))

    # everyone gets what they should
    for (s in c('s1','s2','s3')) {
        expect_identical(filter(long, sample == s)$measured,
                         expression[,s])
    }

})
