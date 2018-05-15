compute_isoform_switches('ctrlExpression.tsv', 'caseExpression.tsv', 'tx2gene.csv', 'out')

test_that('output is as expected', {

    expect_true(file.exists('out'))

})
