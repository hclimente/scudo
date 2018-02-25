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
#' @export
compute_ecdf <- function(x, minValue, minSamples = 10){

    v <- as.numeric(x)

    # discard cases with less than 10 valid cases
    if ( sum(!is.na(v)) < minSamples | sum(v[!is.na(v)] >= minValue) < minSamples){
        return(function(x){return(rep(NA, length(x)))})
    } else {
        diffMatrix <- abs(outer(v,v,"-"))
        subtraction <- diffMatrix[lower.tri(diffMatrix, diag = FALSE)]

        if (all(subtraction[!is.na(subtraction)]==0)){
            return(function(x){
                p <- rep(NA,length(x))
                p[x==0] <- 0
                p[x>0] <- 1
                return(p)})
        } else {
            return(ecdf(subtraction))
        }
    }
}
