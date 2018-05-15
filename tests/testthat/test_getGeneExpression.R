xpr <- readTxExpression('ctrlExpression.tsv', 'caseExpression.tsv') %>%
    getGeneExpression('tx2gene.csv')

test_that('output is as expected', {

    expect_true(is.data.frame(xpr))
    expect_equal(ncol(xpr), 4)
    expect_equal(colnames(xpr), c('gene', 'sample', 'tpmCtrl', 'tpmCase'))

})

test_that('patients are matched', {

    matched <- xpr$sample %in% c('S1', 'S2')
    expect_false(xpr$tpmCtrl[matched] %>% is.na %>% any)
    expect_false(xpr$tpmCase[matched] %>% is.na %>% any)

    onlyCase <- xpr$sample %in% c('S3', 'S4')
    expect_true(xpr$tpmCtrl[onlyCase] %>% is.na %>% all)
    expect_false(xpr$tpmCase[onlyCase] %>% is.na %>% any)

    onlyCtrl <- xpr$sample %in% c('S5', 'S6')
    expect_false(xpr$tpmCtrl[onlyCtrl] %>% is.na %>% any)
    expect_true(xpr$tpmCase[onlyCtrl] %>% is.na %>% all)

})

test_that('sum is correct', {

    expect_equal(xpr$tpmCtrl[xpr$sample == 'S1'], c(15.999 * 3, 15.999 * 2, 15.999))
    expect_equal(xpr$tpmCtrl[xpr$sample == 'S2'], c(7.999 * 3, 7.999 * 2, 7.999))
    expect_equal(xpr$tpmCtrl[xpr$sample == 'S5'], c(31.999 * 3, 31.999 * 2, 31.999))
    expect_equal(xpr$tpmCtrl[xpr$sample == 'S6'], c(3.999 * 3, 3.999 * 2, 3.999))

    expect_equal(xpr$tpmCase[xpr$sample == 'S1'], c(1.999 * 3, 1.999 * 2, 1.999))
    expect_equal(xpr$tpmCase[xpr$sample == 'S2'], c(3.999 * 3, 3.999 * 2, 3.999))
    expect_equal(xpr$tpmCase[xpr$sample == 'S3'], c(31.999 * 3, 31.999 * 2, 31.999))
    expect_equal(xpr$tpmCase[xpr$sample == 'S4'], c(7.999 * 3, 7.999 * 2, 7.999))

})
