#' Calculate deltaPSI
#'
#' @description Calculates the deltaPSI between two conditions. When a subject
#' has matching samples for the two conditions, that value is taken. Else, the
#' deltaPSI is compared to the median of the controls.
#' @param expression data.frame with the expression information in long format,
#' with column names 'gene','sample','tpmCtrl' and 'tpmCase'.
#' @param tx2gene data.frame with two columns, named "transcript" and "gene",
#' which contains the transcript-gene mapping.
#' @return A data.frame with the PSIs per sample and condition, and the
#' corresponding deltaPSI between conditions.
#' @importFrom dplyr group_by mutate select summarize ungroup %>%
#' @importFrom stats median
#' @export
calculate_dPSI <- function(expression, tx2gene) {

    dPSI <- expression %>%
        left_join(tx2gene, by = "transcript") %>%
        group_by(gene, sample) %>%
        mutate(psiCtrl = .data$tpmCtrl/sum(.data$tpmCtrl),
               psiCase = .data$tpmCase/sum(.data$tpmCase)) %>%
        ungroup

    ctrlMedianPsi <- dPSI %>%
        group_by(transcript) %>%
        summarize(psiCtrlMedian = median(.data$psiCtrl, na.rm = T))

    left_join(dPSI, ctrlMedianPsi, by = "transcript") %>%
        mutate(psiCtrl = ifelse(is.na(.data$psiCtrl),
                                  .data$psiCtrlMedian, .data$psiCtrl),
               dPSI = .data$psiCtrl - .data$psiCase) %>%
        select(-psiCtrlMedian)

}

#' Measure p of delta
#'
#' @description Assigns a probability of an observed difference between
#' conditions as function of the differences observed between the controls.
#' @param df data.frame containing the differences.
#' @param probe Unquoted name of the biomarker column (e.g. transcript).
#' @param ctrl Unquoted name of the expression in the control condition column
#' (e.g. tpmCtrl).
#' @param delta Unquoted name of the actual difference between conditions column
#' (e.g. dPSI).
#' @param minValue Minimum value to consider a biomarker expressed.
#' @param minSamples Minimum number of samples to consider the differences
#' informative. If the number of samples with a valid value (i.e. >= minValue)
#' is < minSamples, NAs will be returned.
#' @return A copy of df with an extra column \code{p} indicating the probability
#' of the difference.
#' @importFrom dplyr enquo group_by mutate ungroup %>%
#' @export
score_delta <- function(df, probe, ctrl, delta, minValue = 0, minSamples = 10) {

    probe <- enquo(probe)
    ctrl <- enquo(ctrl)
    delta <- enquo(delta)

    df %>%
        group_by(!!probe) %>%
        mutate(p = compute_ecdf(!!ctrl, minValue, minSamples)(!!delta)) %>%
        ungroup

}

#' Calculate sample-wise differential expression
#'
#' @description Computes the gene expression difference between ctrl and case
#' conditions of the same sample. If no ctrl sample is available, computes the
#' difference with the median expression.
#' @param geneExpression data.frame containing the gene expression in long
#' format.
#' @return A data.frame with the difference in expression \code{DE} between
#' the conditions per gene and sample.
#' @importFrom dplyr everything left_join full_join select %>%
#' @importFrom stats median
#' @export
calculate_DE <- function(geneExpression) {

    ctrlMedian <- geneExpression %>%
        group_by(gene) %>%
        summarize(ctrlMedian = median(.data$ctrlExpression, na.rm = T))

    left_join(geneExpression, ctrlMedian, by = 'gene') %>%
        mutate(ctrlExpression = ifelse(is.na(.data$ctrlExpression),
                                  .data$ctrlMedian, .data$ctrlExpression),
               DE = .data$ctrlExpression - .data$caseExpression) %>%
        select(-ctrlMedian)

}
