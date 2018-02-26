#!/usr/bin/env Rscript

# cd /Users/hclimente/projects/spada/spada/tests/data/calculate_switches
# /Users/hclimente/projects/spada/scripts/calculate_switches.R tx_normal tx_tumor gn_normal gn_tumor tx2gene switches.tsv

#' @importFrom readr cols read_csv read_tsv write_tsv
#' @importFrom dplyr left_join rename select %>%
compute_isoform_switches <- function(normalTxExpressionFile, tumorTxExpressionFile,
                                     normalGeneExpressionFile, tumorGeneExpressionFile, tx2geneFile, outfile) {

    print("Differential transcript-expression")
    tx2gene <- read_csv(tx2geneFile, col_types = 'cc')

    txDPSI <-left_join(read_tsv(normalTxExpressionFile, col_types = cols(transcript = "c", .default = "d")) %>% tpm2long(tpmNormal),
                       read_tsv(tumorTxExpressionFile, col_types = cols(transcript = "c", .default = "d")) %>% tpm2long(tpmTumor),
                       by = c("transcript", "sample")) %>%
        select(gene, transcript, sample, tpmTumor, tpmNormal) %>%
        calculate_dPSI(tx2gene) %>%
        score_delta(transcript, psiNormal, dPSI) %>%
        rename(pdPSI = p)

    print("Calculate differential gene-expression")
    normalGeneExpression <- read_tsv(normalGeneExpressionFile,
                                     col_types = cols(gene = 'c', .default = 'i'))
    tumorGeneExpression <- read_tsv(tumorGeneExpressionFile,
                                    col_types = cols(gene = 'c', .default = 'i'))
    gnDiffExpr <- calculate_diff_expression(normalGeneExpression, tumorGeneExpression)

    print("Compute switches")
    left_join(txDPSI, gnDiffExpr, by = "gene") %>%
        find_switches %>%
        write_tsv(outfile)

}
