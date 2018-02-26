expression <- read.table(text = "
                       transcript sample tpmNormal tpmTumor
                         A1 s1 1 1
                         A2 s1 1 1
                         A3 s1 1 1
                         A1 s2 3 NA
                         A2 s2 2 NA
                         A3 s2 1 NA
                         A1 s3 NA 1
                         A2 s3 NA 2
                         A3 s3 NA 3
                         A1 s4 3 1
                         A2 s4 2 2
                         A3 s4 1 3
                         ", header = TRUE, stringsAsFactors = FALSE)
tx2gene <- read.table(text = "
                      transcript gene
                      A1 A
                      A2 A
                      A3 A
                      B1 B
                      ", header = TRUE, stringsAsFactors = FALSE)
dPSI <- calculate_dPSI(expression, tx2gene)
pdPSI <- score_dPSI(dPSI, minSamples = 2)

test_that('output is as expected', {

    expect_true(is.data.frame(pdPSI))
    expect_equal(ncol(pdPSI), 9)

})

test_that('dPSI p-value calculation is correct', {

    # transcript A1
    for (tx in c('A1','A2','A3')) {
        f <- compute_ecdf(subset(dPSI, transcript == tx)$psiNormal, 0, 2)
        expect_equal(subset(pdPSI, transcript == tx)$pdPSI,
                     f(subset(pdPSI, transcript == tx)$dPSI))
    }

})
