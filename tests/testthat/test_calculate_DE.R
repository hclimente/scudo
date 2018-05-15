geneExpression <- read.table(text = "
                                    gene sample tpmCtrl tpmCase
                                    A s1 2 3
                                    A s2 3 1
                                    A s3 1 10
                                    A s4 NA 2
                                    A s5 3 NA
                                    A s6 1 4
                                    ", header = TRUE, stringsAsFactors = FALSE)

DE <- calculate_DE(geneExpression)

test_that('output is as expected', {

    # we have the all the information for all every gene in every sample
    expect_equal(ncol(DE), 5)
    expect_equal(nrow(DE), nrow(geneExpression))

    # all samples are present
    expect_true(all(geneExpression$sample %in% DE$sample))

    # samples have the right value-types
    expect_false(any(is.na(DE$tpmCtrl)))
    expect_equal(DE$sample[is.na(DE$tpmCase)], 's5')

    # delta is the right difference
    expect_equal(DE$DE, DE$tpmCtrl - DE$tpmCase)
    expect_equal(is.na(DE$tpmCase), is.na(DE$DE))

})

test_that('output is compatible with score_delta', {

    pDE <- score_delta(DE, gene, tpmCtrl, DE, minSamples = 2)
    expect_true(any(!is.na(pDE[,'p'])))

})
