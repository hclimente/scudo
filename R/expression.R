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
#' @importFrom stats median
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

#' Calculate sample-wise differential expression
#'
#' @description Computes the gene expression difference between normal and tumor
#' conditions of the same sample. If no normal sample is available, computes the
#' difference with the median expression. Prior to that, performs a
#' normalization of the data.
#' @param normalGeneExpression data.frame containing the gene expression in wide
#' format. A column \code{gene} indicates the gene, and the rest of the column
#' are the CPM for each sample.
#' @param tumorGeneExpression Equivalent to normalGeneExpression with tumor
#' data.
#' @return A data.frame with the difference in expression \code{DE} between
#' the conditions per gene and sample.
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom dplyr everything left_join full_join select %>%
#' @importFrom stats median
#' @export
calculate_DE <- function(geneExpression) {

    normalMedian <- geneExpression %>%
        group_by(gene) %>%
        summarize(normalMedian = median(.data$normalExpression, na.rm = T))

    left_join(geneExpression, normalMedian, by = 'gene') %>%
        mutate(normalExpression = ifelse(is.na(.data$normalExpression),
                                  .data$normalMedian, .data$normalExpression),
               DE = .data$normalExpression - .data$tumorExpression) %>%
        select(-normalMedian)

}

#' @export
normalize_counts <- function(normalGeneExpression, tumorGeneExpression) {

    # generate expression matrix
    geneExpression <- left_join(normalGeneExpression, tumorGeneExpression,
                                by = 'gene')
    genes <- geneExpression[,'gene']
    geneExpression <- select(geneExpression, -gene) %>% as.matrix
    rownames(geneExpression) <- genes

    # get sample names
    border <- ncol(normalGeneExpression)
    normalSamples <- colnames(geneExpression)[1:border - 1]
    tumorSamples <- colnames(geneExpression)[border:(border - 2 + ncol(tumorGeneExpression))]

    # normalize expression and calculate log2
    y <- DGEList(counts = geneExpression)
    y <- calcNormFactors(y)

    logNormGeneExpression <- cpm(y, normalized.lib.sizes = TRUE) %>% log2
    logNormGeneExpression[logNormGeneExpression == -Inf] <- NA

    normalGeneExpression <- data.frame(logNormGeneExpression[,normalSamples]) %>%
        mutate(gene = genes) %>%
        select(gene, everything() ) %>%
        tpm2long(normalExpression) %>%
        mutate(sample = gsub('\\.x$', '', sample))
    tumorGeneExpression <- data.frame(logNormGeneExpression[,tumorSamples]) %>%
        mutate(gene = genes) %>%
        select(gene, everything() ) %>%
        tpm2long(tumorExpression) %>%
        mutate(sample = gsub('\\.y$', '', sample))

    full_join(normalGeneExpression, tumorGeneExpression, by = c('gene','sample'))
}
