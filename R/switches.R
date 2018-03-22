#' Find isoform switches
#'
#' @description Find the isoform switches that match all of the following
#' conditions: 1) genes are not differentially expressed; 2) abs(deltaPSI) is
#' > 0.05 and more extreme than 1% of the observed variance between controls; 3)
#' involved transcripts are sufficiently expressed (TPM > 0.1).
#' @param tx_pdPSI data.frame with transcript-level annotation per sample, as
#' outputted by \code{score_delta}. Must contain all of the following columns:
#' \code{dPSI} indicating the deltaPSI and \code{p} scoring how extreme it is;
#' \code{tpmNormal} and \code{tpmTumor} indicating the expression in the two
#' conditions; and \code{transcript}, \code{gene} and \code{sample}.
#' @param gene_pDE data.frame with gene-level annotation per sample, as
#' outputted by \code{score_delta}. Must contain all of the following columns:
#' \code{p} indicating if the gene is differentially expressed; and \code{gene}
#' and \code{sample}.
#' @param pDE Minimum p-value to consider a gene differentially expressed.
#' @param pdPSI Minimum p-value to consider a scplicing change significant.
#' @param deltaPSI Minimum PSI change to consider a scplicing change significant.
#' @param minExpression Minimum expression of a transcript to consider it
#' expressed.
#' @return A data.frame with the measured isoform switches.
#' @importFrom dplyr filter group_by summarize ungroup %>%
#' @export
find_switches <- function(tx_pdPSI, gene_pDE, pDE = 0.01, pdPSI = 0.01,
                          deltaPSI = 0.05, minExpression = 0.1) {

    tx_pdPSI <- rename(tx_pdPSI, p_dPSI = p)
    gene_pDE <- rename(gene_pDE, p_DiffEExpression = p)

    left_join(tx_pdPSI, gene_pDE, by = c("gene", "sample")) %>%
        # remove differentially expressed genes
        filter(p_DiffEExpression >= pDE) %>%
        # remove low deltaPSI
        filter(p_dPSI <= pdPSI & abs(dPSI) >= deltaPSI) %>%
        # remove lowly expressed transcript
        filter((dPSI > 0 & tpmTumor >= minExpression) |
               (dPSI < 0 & tpmNormal >= minExpression)) %>%
        group_by(gene, sample) %>%
        # remove unelligible genes
        filter(n() > 1 & any(dPSI > 0) & any(dPSI < 0)) %>%
        summarize(normal = transcript[which.min(dPSI)],
                  tumor = transcript[which.max(dPSI)]) %>%
        ungroup %>%
        group_by(gene, normal, tumor) %>%
        summarize(samples = paste(sample, collapse = ','))

}

# cd /Users/hclimente/projects/spada/spada/tests/data/calculate_switches
# /Users/hclimente/projects/spada/scripts/calculate_switches.R tx_normal tx_tumor gn_normal gn_tumor tx2gene switches.tsv
#' @importFrom readr cols read_csv read_tsv write_tsv
#' @importFrom dplyr left_join rename select %>%
#' @export
compute_isoform_switches <- function(normalTxExpressionFile,
                                     tumorTxExpressionFile,
                                     normalGeneExpressionFile,
                                     tumorGeneExpressionFile,
                                     tx2geneFile, outfile) {

    print("Differential transcript-expression")
    tx2gene <- read_csv(tx2geneFile, col_types = 'cc')

    colt <- cols(transcript = "c", .default = "d")
    normalTxExpression <- read_tsv(normalTxExpressionFile, col_types = colt) %>%
        tpm2long(tpmNormal)
    tumorTxExpression <- read_tsv(tumorTxExpressionFile, col_types = colt) %>%
        tpm2long(tpmTumor)

    tx_pdPSI <-left_join(normalTxExpression, tumorTxExpression,
                         by = c("transcript", "sample")) %>%
        select(gene, transcript, sample, tpmTumor, tpmNormal) %>%
        calculate_dPSI(tx2gene) %>%
        score_delta(transcript, psiNormal, dPSI)

    rm(normalTxExpression, tumorTxExpression)

    print("Calculate differential gene-expression")
    colt <- cols(gene = "c", .default = "i")
    normalGeneExpression <- read_tsv(normalGeneExpressionFile, col_types = colt)
    tumorGeneExpression <- read_tsv(tumorGeneExpressionFile, col_types = colt)

    gene_pDE <- calculate_dExpr(normalGeneExpression, tumorGeneExpression) %>%
        score_delta(transcript, normalExpression, DE)

    rm(normalGeneExpression, tumorGeneExpression)

    print("Compute switches")
    find_switches(tx_pdPSI, gene_pDE) %>%
        write_tsv(outfile)

}
