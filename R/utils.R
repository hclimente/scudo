#' Get empirical distribution
#'
#' @description  Get the empirical distribution of the differences between a
#' vector
#' @param x numerical vector
#' @param minValue minimum value to consider that measure
#' @param minSamples Minimum number of valid samples.
#' @return \itemize{
#' \item {function returning NA when there are less than 10 valid values (ie non
#' NA and higher than threshold)}
#' \item {function returning 1 to any value higher than 0 when all the
#' differences are 0.}
#' \item {return ecdf otherwise, counting that the minimum p-value will be 1/n+1
#' }}
#' @importFrom stats ecdf
#' @export
compute_ecdf <- function(x, minValue, minSamples){

    v <- as.numeric(x)

    # discard cases with less than 10 valid cases
    if ( sum(!is.na(v)) < minSamples | sum(v[!is.na(v)] >= minValue) < minSamples){
        f <- uninformative_output
    } else {
        diffMatrix <- abs(outer(v,v,"-"))
        subtraction <- diffMatrix[lower.tri(diffMatrix, diag = FALSE)]

        if (all(subtraction[!is.na(subtraction)]==0)){
            f <- uninformative_input
        } else {
            f <- ecdf(subtraction)
        }
    }

    return(f)
}

#' Pseudo-ecdf when all values are 0
#'
#' @description ecdf function where anything higher than 0 gets a p-value of 1.
#' @param x numerical vector
#' @return Numerical vector with 0 when the input is 0, and 1 when it is > 0.
uninformative_input <- function(x) {

    p <- rep(NA, length(x))
    p[x == 0] <- 0
    p[x > 0] <- 1
    return(p)

}

#' Pseudo-ecdf when input is insufficient
#'
#' @description ecdf function where anything returns NA.
#' @param x numerical vector.
#' @return Numerical vector with NAs.
uninformative_output <- function(x){

    return(rep(NA, length(x)))

}

#' Read expression
#'
#' @param ctrlExpressionFile File containing the transcript expression for
#' the control samples, measured in log2(TPMs + 0.001).
#' @param caseExpressionFile File containing the transcript expression for
#' the case samples, measured in log2(TPMs + 0.001).
#' @return Data frame with transcript expression in TPM in long format.
#' @importFrom dplyr full_join mutate select
#' @importFrom readr cols read_tsv
#' @importFrom tidyr gather
read_tx_expression <- function(ctrlExpressionFile, caseExpressionFile) {

    colt <- cols(transcript = "c", .default = "d")
    ctrlTxExpression <- read_tsv(ctrlExpressionFile, col_types = colt) %>%
        gather(key = sample, value = expression, -1) %>%
        mutate(tpmCtrl = -0.001 + 2^(expression))
    caseTxExpression <- read_tsv(caseExpressionFile, col_types = colt) %>%
        gather(key = sample, value = expression, -1) %>%
        mutate(tpmCase = -0.001 + 2^(expression))

    full_join(ctrlTxExpression, caseTxExpression,
              by = c('transcript', 'sample')) %>%
        select(transcript, sample, tpmCtrl, tpmCase)

}

#' Read expression
#'
#' @param txExpression Data frame with expression in TPM in long format.
#' @param tx2geneFile File containing genes-transcript correspondence, one
#' per line.
#' @return Data frame with gene expression in TPM in long format.
get_gene_expression <- function(txExpression, tx2geneFile) {

    tx2gene <- read_csv(tx2geneFile, col_types = 'cc')

    inner_join(txExpression, tx2gene, by = "transcript") %>%
        select(-transcript) %>%
        group_by(gene, sample) %>%
        summarise_all(funs(sum))

}
