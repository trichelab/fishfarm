#' extract clonal prevalence from a Canopy run for plotting
#'
#' In order to plot Canopy configs nicely, we use timescape.
#' This function massages the Canopy output for timescape. 
#' 
#' @param output_tree     Canopy output
#' @param tp_order        timepoint ordering
#'
#' @return                a data.frame 
#' 
#' @import Canopy
#' @import timescape
#' 
#' @export
#'
getClonalPrev <- function(output_tree, tp_order) {
  
  n <- ncol(output_tree$VAF)
  k <- max(output_tree$edge)
  kP <- nrow(output_tree$P)
  k0 <- k - kP # internal-only nodes
  prevalence <- data.frame(output_tree$P)
  rownames(prevalence) <- seq_len(nrow(prevalence))
  sampleNames <- colnames(output_tree$P)

  clonal_prev <- data.frame( 
    timepoint = do.call(c, lapply(sampleNames, rep, k)),
    clone_id = rep(seq_len(k), n), 
    clonal_prev = do.call(c, lapply(prevalence, function(x) c(x, rep(0, k0)))),
    stringsAsFactors = FALSE
  )
  rownames(clonal_prev) <- with(clonal_prev, 
                                paste0(timepoint, "_clone", clone_id))
  clonal_prev$timepoint <- factor(clonal_prev$timepoint, levels = tp_order)
  return(clonal_prev[order(clonal_prev$timepoint, clonal_prev$clone_id), ])
  
}
