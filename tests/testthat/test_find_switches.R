tx_pdPSI <- read.table(text = "
                    transcript gene sample dPSI p tpmCtrl tpmCase
                    A1 A s1 0.1 0.01 1 1
                    A2 A s1 0.2 0.01 1 1
                    A3 A s1 -0.1 0.01 1 1
                    A1 A s2 -0.1 0.01 1 1
                    A2 A s2 0.2 0.01 1 1
                    A3 A s2 0.1 0.01 1 1
                    ", header = TRUE, stringsAsFactors = FALSE)

gene_pDE <- read.table(text = "
                  gene sample p
                  A s1 0.2
                  A s2 0.2
                  ", header = TRUE, stringsAsFactors = FALSE)

test_that('dPSI conditions the switches', {

    tmp <- mutate(tx_pdPSI, dPSI = c(0.04, 0.1, -0.1, -0.1, 0.1, 0.04))

    expect_equal(filter(tmp, dPSI != 0.1) %>%
                 find_switches(gene_pDE, deltaPSI = 0.04),
                 read.table(text = "
                            gene ctrl case samples
                            A A3 A1 s1
                            A A1 A3 s2
                            ", header = TRUE, stringsAsFactors = FALSE))

    expect_equal(nrow(find_switches(filter(tmp, dPSI != 0.1), gene_pDE)), 0)

})

test_that('pdPSI conditions the switches', {

    tmp <- mutate(tx_pdPSI, p = c(0.01, 0.1, 0.01, 0.01, 0.1, 0.01))

    expect_equal(find_switches(tmp, gene_pDE),
                 read.table(text = "
                            gene ctrl case samples
                            A A1 A3 s2
                            A A3 A1 s1
                            ", header = TRUE, stringsAsFactors = FALSE))

    expect_equal(find_switches(tmp, gene_pDE, pdPSI = 0.1),
                 read.table(text = "
                            gene ctrl case samples
                            A A1 A2 s2
                            A A3 A2 s1
                            ", header = TRUE, stringsAsFactors = FALSE))

})

test_that('minExpression conditions the switches', {

    tmp <- mutate(tx_pdPSI, tpmCtrl = c(1, 1, 1, 0.05, 1, 1),
                            tpmCase = c(1, 0.05, 1, 1, 1, 1))

    expect_equal(find_switches(tmp, gene_pDE),
                 read.table(text = "
                            gene ctrl case samples
                            A A3 A1 s1
                            ", header = TRUE, stringsAsFactors = FALSE))

    expect_equal(find_switches(tmp, gene_pDE, minExpression = 0.05),
                 read.table(text = "
                            gene ctrl case samples
                            A A1 A2 s2
                            A A3 A2 s1
                            ", header = TRUE, stringsAsFactors = FALSE))

})

test_that('pDE conditions the switches', {

    tmp <- mutate(gene_pDE, p = c(0.01, 0.001))

    expect_equal(find_switches(tx_pdPSI, tmp, pDE = 0.001),
                 read.table(text = "
                            gene ctrl case samples
                            A A1 A2 s2
                            A A3 A2 s1
                            ", header = TRUE, stringsAsFactors = FALSE))

    expect_equal(find_switches(tx_pdPSI, tmp),
                 read.table(text = "
                            gene ctrl case samples
                            A A3 A2 s1
                            ", header = TRUE, stringsAsFactors = FALSE))

})

test_that('we aggreggate patients with matching switches', {

    tmp <- mutate(tx_pdPSI, dPSI = c(0.04, 0.1, -0.1, 0.04, 0.1, -0.1))

    expect_equal(find_switches(tmp, gene_pDE, deltaPSI = 0.04),
                 read.table(text = "
                            gene ctrl case samples
                            A A3 A2 s1,s2
                            ", header = TRUE, stringsAsFactors = FALSE))

})
