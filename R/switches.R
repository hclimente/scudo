#' Find isoform switches
#'
#' @description Find the isoform switches that match all of the following
#' conditions: 1) genes are not differentially expressed; 2) abs(deltaPSI) is
#' > 0.05 and more extreme than 1% of the observed variance between controls; 3)
#' involved transcripts are sufficiently expressed (TPM > 0.1).
#' @param tx_pdPSI data.frame with transcript-level annotation per sample, as
#' outputted by \code{score_delta}. Must contain all of the following columns:
#' \code{dPSI} indicating the deltaPSI and \code{p} scoring how extreme it is;
#' \code{tpmCtrl} and \code{tpmCase} indicating the expression in the two
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
        filter((dPSI > 0 & tpmCase >= minExpression) |
               (dPSI < 0 & tpmCtrl >= minExpression)) %>%
        group_by(gene, sample) %>%
        # remove unelligible genes
        filter(n() > 1 & any(dPSI > 0) & any(dPSI < 0)) %>%
        summarize(ctrl = transcript[which.min(dPSI)],
                  case = transcript[which.max(dPSI)]) %>%
        ungroup %>%
        group_by(gene, ctrl, case) %>%
        summarize(samples = paste(sample, collapse = ','))

}

#' Compute isoform switches
#'
#' @description Top-level function to find isoform switches from transcript- and
#' gene-level expression files.
#'
#' @param ctrlExpressionFile Tab-separated table with normalized
#' transcript-level expression for the control samples, in log2(TPM + 0.001).
#' @param caseExpressionFile Tab-separated table with normalized
#' transcript-level expression for the case samples, in log2(TPM + 0.001).
#' @param tx2geneFile  Tab-separated table with the assignation of transcripts
#' to genes (one per line).
#' @param outfile Output file.
#' @param ... arguments passed to the \code{\link{find_switches}} function.
#' @importFrom readr cols read_csv read_tsv write_tsv
#' @importFrom dplyr left_join rename select %>%
#' @export
compute_isoform_switches <- function(ctrlExpressionFile, caseExpressionFile,
                                     tx2geneFile, outfile, ...) {

    print("Differential transcript-expression")
    txExpression <- readTxExpression(ctrlTxExpression, caseTxExpression)

    tx_pdPSI <- txExpression %>%
        calculate_dPSI(tx2gene) %>%
        score_delta(transcript, psiCtrl, dPSI)

    print("Calculate differential gene-expression")
    gene_pDE <- getGeneExpression(txExpression, tx2geneFile) %>%
        calculate_DE %>%
        score_delta(gene, ctrlExpression, DE)
    rm(txExpression)

    print("Compute switches")
    find_switches(tx_pdPSI, gene_pDE, ...) %>%
        write_tsv(outfile)

}
