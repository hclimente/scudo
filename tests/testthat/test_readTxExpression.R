xpr <- readTxExpression('ctrlExpression.tsv', 'caseExpression.tsv')

test_that('output is as expected', {

    expect_true(is.data.frame(xpr))
    expect_equal(ncol(xpr), 4)
    expect_equal(colnames(xpr), c('transcript', 'sample', 'tpmCtrl', 'tpmCase'))

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

test_that('log2 is reverted', {

    expect_true(all(xpr$tpmCtrl[xpr$sample == 'S1'] == 15.999))
    expect_true(all(xpr$tpmCtrl[xpr$sample == 'S2'] == 7.999))
    expect_true(all(xpr$tpmCtrl[xpr$sample == 'S5'] == 31.999))
    expect_true(all(xpr$tpmCtrl[xpr$sample == 'S6'] == 3.999))

    expect_true(all(xpr$tpmCase[xpr$sample == 'S1'] == 1.999))
    expect_true(all(xpr$tpmCase[xpr$sample == 'S2'] == 3.999))
    expect_true(all(xpr$tpmCase[xpr$sample == 'S3'] == 31.999))
    expect_true(all(xpr$tpmCase[xpr$sample == 'S4'] == 7.999))

})
