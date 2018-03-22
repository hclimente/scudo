normalGeneExpression <- read.table(text = "
                                   gene s1 s2 s3 s4
                                   A 1 2 3 1
                                   B 2 3 1 1
                                   C 3 1 2 1
                                   D 1 0 3 1
                                   E 2 0 1 1
                                   ", header = TRUE, stringsAsFactors = FALSE)

tumorGeneExpression <- read.table(text = "
                                  gene s1 s2 s3 s5
                                  A 7 8 9 1
                                  B 8 9 7 1
                                  C 9 7 8 1
                                  D 7 0 9 1
                                  E 8 0 7 1
                                  ", header = TRUE, stringsAsFactors = FALSE)

normalizedCounts <- normalize_counts(normalGeneExpression, tumorGeneExpression)

nSamples <- colnames(normalGeneExpression[,-1])
tSamples <- colnames(tumorGeneExpression[,-1])
shared <- intersect(nSamples, tSamples)
nExclusive <- setdiff(nSamples, tSamples)
tExclusive <- setdiff(tSamples, nSamples)

test_that('output is as expected', {

    # we have the all the information for all every gene in every sample
    nrows <- nrow(tumorGeneExpression)*(length(shared)+length(nExclusive)+length(tExclusive))
    expect_equal(ncol(normalizedCounts), 4)
    expect_equal(nrow(normalizedCounts), nrows)

    # all samples are present
    expect_true(all(shared %in% normalizedCounts$sample))
    expect_true(nExclusive %in% normalizedCounts$sample)
    expect_true(tExclusive %in% normalizedCounts$sample)

    # samples have the right value-types
    expect_true(filter(normalizedCounts, sample %in% nExclusive)$tumorExpression %>% is.na %>% all)

})

test_that('normalization is correct', {


})
