#' @importFrom tidyr gather
#' @export
tpm2long <- function(wide, colname) {

    gather(wide, key = sample, value = !!colname, -1)

}

#' @importFrom dplyr group_by mutate select summarize ungroup %>%
#' @export
calculate_dPSI <- function(expression, tx2gene) {

    dPSI <- expression %>%
        left_join(tx2gene, by = "transcript") %>%
        group_by(gene, sample) %>%
        mutate(psiNormal = tpmNormal/sum(tpmNormal),
               psiTumor = tpmTumor/sum(tpmTumor)) %>%
        ungroup

    normalMedianPsi <- dPSI %>%
        group_by(transcript) %>%
        summarize(psiNormalMedian = median(psiNormal, na.rm = T))

    left_join(dPSI, normalMedianPsi, by = "transcript") %>%
        mutate(psiNormal = ifelse(is.na(psiNormal), psiNormalMedian, psiNormal),
               dPSI = psiNormal - psiTumor) %>%
        select(-psiNormalMedian)

}

#' @importFrom dplyr group_by mutate ungroup %>%
#' @export
score_dPSI <- function(dPSI, minValue = 0, minSamples = 10) {

    pdPSI <- by(dPSI$psiNormal, dPSI$transcript, compute_ecdf,
                minValue, minSamples)

    dPSI %>%
        group_by(transcript) %>%
        mutate(pdPSI = pdPSI[[unique(transcript)]](dPSI)) %>%
        ungroup

}

#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom dplyr left_join select %>%
#' @importFrom tibble tibble
#' @export
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

#' @importFrom dplyr filter group_by summarize ungroup %>%
#' @export
find_switches <- function(txInfo) {

    txInfo %>%
        # remove differentially expressed genes
        filter(pDE > 0.01) %>%
        # remove low deltaPSI
        filter(pdPSI < 0.01 & dPSI > 0.05) %>%
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
