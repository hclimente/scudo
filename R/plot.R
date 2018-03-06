#' @importFrom ggplot2 aes geom_histogram geom_point ggplot ggplotGrob labs
#' scale_fill_manual
#' @importFrom grid textGrob
#' @importFrom gridExtra grid.arrange
plot_transcripts <- function(tx_dPSI, switch, what) {

    gene <- switch$gene
    nTx <- switch$normal
    tTx <- switch$tumor

    palette <- c("Tumor no switch"="#7fc97f",
                 "Normal"="#beaed4",
                 "Tumor switch"="#fdc086")

    samples <- strsplit(switch$samples, ',', fixed = T)
    tx_dPSI <- tx_dPSI %>%
        filter(transcript %in% c(nTx, tTx)) %>%
        mutate(switched = ifelse(sample %in% samples, 'Yes', 'No'))

    if (what == 'expression') {
        scatter <- ggplot(tx_dPSI,
                          aes(x = tpmNormal, y = tpmTumor, color = switched)) +
            geom_point()
    } else if (what == 'psi') {
        scatter <- ggplot(tx_dPSI,
                          aes(x = psiNormal, y = psiTumor, color = switched)) +
            geom_point()
    }

    nHist <- ggplot(tx_dPSI,
                    aes(x = tpmNormal, fill = switched)) +
                    geom_histogram(alpha = 0.75) +
                    labs(x = 'Normal isoform TPM', y = '') +
                    scale_fill_manual(values = palette)

    tHist <- ggplot(tx_dPSI,
                    aes(x = tpmTumor, fill = switched)) +
        geom_histogram(alpha = 0.75) +
        labs(x = 'Tumor isoform TPM', y = '') +
        scale_fill_manual(values = palette)

    x <- ggplotGrob(scatter +
                        theme(legend.position = 'bottom') +
                        labs(color = '') +
                        guides(color =
                                   guide_legend(nrow = 3, byrow = TRUE,
                                        override.aes = list(shape = 15, size = 8)),
                               shape =
                                   guide_legend(nrow = 2, byrow = TRUE,
                                        override.aes = list(size = 8))))$grobs

    legend <- x[[which(sapply(x, function(y) y$name) == "guide-box")]]

    grid.arrange(tHist, scatter, legend, nHist,
                 ncol = 2, nrow = 2, widths = c(1.5,5), heights = c(5,1.5),
                 top = textGrob(paste(gene, nTx, tTx),
                              gp = gpar(fontsize=20,font=3)))

}

#' @importFrom ggplot2 aes geom_point ggplot labs
plot_genes <- function(gene_pDE, switch, what) {

    gene <- switch$gene
    nTx <- switch$normal
    tTx <- switch$tumor

    samples <- strsplit(switch$samples, ',', fixed = T)
    gene_pDE <- gene_pDE %>%
        filter(gene == gene) %>%
        mutate(switched = ifelse(sample %in% samples, 'Yes', 'No'))

    ggplot(gene_pDE, aes(x = log2(normalExpression), y = log2(tumorExpression),
                              color = switched)) +
        geom_point() +
        labs(x = "log2(normal gene cpm)", y = "log2(tumor gene cpm)",
             color = 'Switched')

}
