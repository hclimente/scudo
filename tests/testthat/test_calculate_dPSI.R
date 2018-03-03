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

test_that('output is as expected', {

    expect_true(is.data.frame(dPSI))
    expect_equal(ncol(dPSI), 8)

})

test_that('PSI calculation is correct', {

    # normal
    expect_equal(subset(dPSI, sample == 's1' & transcript == 'A1')$psiNormal, 1/3)
    expect_equal(subset(dPSI, sample == 's1' & transcript == 'A2')$psiNormal, 1/3)
    expect_equal(subset(dPSI, sample == 's1' & transcript == 'A3')$psiNormal, 1/3)
    expect_equal(subset(dPSI, sample == 's2' & transcript == 'A1')$psiNormal, 1/2)
    expect_equal(subset(dPSI, sample == 's2' & transcript == 'A2')$psiNormal, 1/3)
    expect_equal(subset(dPSI, sample == 's2' & transcript == 'A3')$psiNormal, 1/6)
    expect_equal(subset(dPSI, sample == 's3' & transcript == 'A1')$psiNormal, 1/2)
    expect_equal(subset(dPSI, sample == 's3' & transcript == 'A2')$psiNormal, 1/3)
    expect_equal(subset(dPSI, sample == 's3' & transcript == 'A3')$psiNormal, 1/6)
    expect_equal(subset(dPSI, sample == 's4' & transcript == 'A1')$psiNormal, 1/2)
    expect_equal(subset(dPSI, sample == 's4' & transcript == 'A2')$psiNormal, 1/3)
    expect_equal(subset(dPSI, sample == 's4' & transcript == 'A3')$psiNormal, 1/6)

    # tumor
    expect_equal(subset(dPSI, sample == 's1' & transcript == 'A1')$psiTumor, 1/3)
    expect_equal(subset(dPSI, sample == 's1' & transcript == 'A2')$psiTumor, 1/3)
    expect_equal(subset(dPSI, sample == 's1' & transcript == 'A3')$psiTumor, 1/3)
    expect_true(is.na(subset(dPSI, sample == 's2' & transcript == 'A1')$psiTumor), NA)
    expect_true(is.na(subset(dPSI, sample == 's2' & transcript == 'A2')$psiTumor), NA)
    expect_true(is.na(subset(dPSI, sample == 's2' & transcript == 'A3')$psiTumor), NA)
    expect_equal(subset(dPSI, sample == 's3' & transcript == 'A1')$psiTumor, 1/6)
    expect_equal(subset(dPSI, sample == 's3' & transcript == 'A2')$psiTumor, 1/3)
    expect_equal(subset(dPSI, sample == 's3' & transcript == 'A3')$psiTumor, 1/2)
    expect_equal(subset(dPSI, sample == 's4' & transcript == 'A1')$psiTumor, 1/6)
    expect_equal(subset(dPSI, sample == 's4' & transcript == 'A2')$psiTumor, 1/3)
    expect_equal(subset(dPSI, sample == 's4' & transcript == 'A3')$psiTumor, 1/2)

})

test_that('dPSI calculation is correct', {

    expect_equal(subset(dPSI, sample == 's1' & transcript == 'A1')$dPSI, 0)
    expect_equal(subset(dPSI, sample == 's1' & transcript == 'A2')$dPSI, 0)
    expect_equal(subset(dPSI, sample == 's1' & transcript == 'A3')$dPSI, 0)
    expect_true(is.na(subset(dPSI, sample == 's2' & transcript == 'A1')$dPSI), NA)
    expect_true(is.na(subset(dPSI, sample == 's2' & transcript == 'A2')$dPSI), NA)
    expect_true(is.na(subset(dPSI, sample == 's2' & transcript == 'A3')$dPSI), NA)
    expect_equal(subset(dPSI, sample == 's3' & transcript == 'A1')$dPSI, 1/2 - 1/6)
    expect_equal(subset(dPSI, sample == 's3' & transcript == 'A2')$dPSI, 0)
    expect_equal(subset(dPSI, sample == 's3' & transcript == 'A3')$dPSI, 1/6 - 1/2)
    expect_equal(subset(dPSI, sample == 's4' & transcript == 'A1')$dPSI, 1/2 - 1/6)
    expect_equal(subset(dPSI, sample == 's4' & transcript == 'A2')$dPSI, 0)
    expect_equal(subset(dPSI, sample == 's4' & transcript == 'A3')$dPSI, 1/6 - 1/2)

})
