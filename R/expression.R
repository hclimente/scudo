#' Convert TPM to long format
#'
#' @description Converts TPM, wide-format matrix to a long format tibble.
#' @param wide Wide format data.frame with TPM values. First column contains the
#' transcript id, and succesive columns are the different samples.
#' @param colname Unquoted name of the desired expression column name.
#' @return A data.frame with three columns: transcript, sample and expression.
#' @importFrom dplyr enquo
#' @importFrom tidyr gather
#' @export
tpm2long <- function(wide, colname) {

    colname <- enquo(colname)
    gather(wide, key = sample, value = !!colname, -1)

}

#' Calculate deltaPSI
#'
#' @description Calculates the deltaPSI between two conditions. When a subject
#' has matching samples for the two conditions, that value is taken. Else, the
#' deltaPSI is compared to the median of the controls.
#' @param expression data.frame with the expression information in long format,
#' as outputed by \code{tpm2long}.
#' @param tx2gene data.frame with two columns, named "transcript" and "gene",
#' which contains the transcript-gene mapping.
#' @return A data.frame with the PSIs per sample and condition, and the
#' corresponding deltaPSI between conditions.
#' @importFrom dplyr group_by mutate select summarize ungroup %>%
#' @export
calculate_dPSI <- function(expression, tx2gene) {

    dPSI <- expression %>%
        left_join(tx2gene, by = "transcript") %>%
        group_by(gene, sample) %>%
        mutate(psiNormal = .data$tpmNormal/sum(.data$tpmNormal),
               psiTumor = .data$tpmTumor/sum(.data$tpmTumor)) %>%
        ungroup

    normalMedianPsi <- dPSI %>%
        group_by(transcript) %>%
        summarize(psiNormalMedian = median(.data$psiNormal, na.rm = T))

    left_join(dPSI, normalMedianPsi, by = "transcript") %>%
        mutate(psiNormal = ifelse(is.na(.data$psiNormal),
                                  .data$psiNormalMedian, .data$psiNormal),
               dPSI = .data$psiNormal - .data$psiTumor) %>%
        select(-psiNormalMedian)

}

#' Measure p of delta
#'
#' @description Assigns a probability of an observed difference between
#' conditions as function of the differences observed between the controls.
#' @param df data.frame containing the differences.
#' @param probe Unquoted name of the biomarker column (e.g. transcript).
#' @param ctrl Unquoted name of the expression in the control condition column
#' (e.g. tpmNormal).
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

#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom dplyr left_join select %>%
#' @importFrom tibble tibble
calculate_diff_expression <- function(normalGeneExpression, tumorGeneExpression) {

    # get sample names
    normalSamples <- colnames(normalGeneExpression)
    tumorSamples <- colnames(tumorGeneExpression)

    # generate expression matrix
    geneExpression <- left_join(tumorGeneExpression, normalGeneExpression,
                                by = 'gene')
    genes <- geneExpression$gene
    geneExpression <- select(geneExpression, -gene) %>% as.matrix
    rownames(geneExpression) <- genes

    # normalize expression and calculate log2
    y <- DGEList(counts = geneExpression)
    y <- calcNormFactors(y)

    logNormGeneExpression <- cpm(y, normalized.lib.sizes=TRUE) %>%
        log2
    logNormGeneExpression[logNormGeneExpression == -Inf] <- NA

    # calculate differential expression
    tibble(gene = rownames(logNormGeneExpression),
           pDE = apply(logNormGeneExpression, 1, function(x) {
               wilcox.test(x[normalSamples], x[tumorSamples])$p.value
           }))

}

#' Find isoform switches
#'
#' @description Find the isoform switches that match all of the following
#' conditions: 1) genes are not differentially expressed; 2) abs(deltaPSI) is > 0.05
#' and more extreme than 1% of the observed variance between controls; 3)
#' involved transcripts are sufficiently expressed (TPM > 0.1).
#' @param txInfo data.frame with transcript-level annotation per subject. Must
#' contain all of the following columns: \code{pDE} indicating if the gene is
#' differentially expressed; \code{dPSI} indicating the deltaPSI and
#' \code{pdPSI} scoring how extreme it is; \code{tpmNormal} and \code{tpmTumor}
#' indicating the expression in the two conditions.
#' @return A data.frame with the measured isoform switches.
#' @importFrom dplyr filter group_by summarize ungroup %>%
#' @export
find_switches <- function(txInfo) {

    txInfo %>%
        # remove differentially expressed genes
        filter(pDE > 0.01) %>%
        # remove low deltaPSI
        filter(pdPSI < 0.01 & abs(dPSI) > 0.05) %>%
        # remove lowly expressed transcript
        filter((dPSI > 0 & tpmTumor > 0.1) | (dPSI < 0 & tpmNormal > 0.1)) %>%
        group_by(gene, sample) %>%
        # remove unelligible genes
        filter(n() > 1 & any(dPSI > 0) & any(dPSI < 0)) %>%
        summarize(normal = transcript[which.min(dPSI)],
                  tumor = transcript[which.max(dPSI)]) %>%
        ungroup %>%
        group_by(gene, normal, tumor) %>%
        summarize(samples = paste(samples, collapse = ','))

}
