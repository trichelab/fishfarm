#' extract clonal mutations from a Canopy run and its VAFs 
#'
#' In order to plot Canopy configs nicely, we use timescape.
#' This function massages Canopy output into mutation data.
#' Note: if genelist has impacts, coords, breakpoints, etc, they are added!
#' 
#' @param output_tree     canopy output
#' @param VAF             variant allele frequencies
#' @param genelist        coordinates for the genes (default is to autopopulate)
#'
#' @return a data.frame
#' 
#' @import Canopy
#' @import timescape
#' @import Homo.sapiens
#'
#' @export
#'
assembleMutations <- function(output_tree, VAF, genelist=NULL) {

  # genelist could be preannotated, hint hint
  if (is.null(genelist)) { 
    genelist <- sort(genes(Homo.sapiens))
    names(genelist) <- mapIds(Homo.sapiens,names(genelist),"SYMBOL","ENTREZID")
    genelist <- subset(genelist, 
                       seqnames(genelist) %in% paste0("chr", c(1:22,"X","Y")) &
                         !is.na(names(genelist)))
  }
  clonalmut <- output_tree$clonalmut
  mutgenes <- genelist[intersect(names(genelist), unlist(clonalmut))]
  mutgenes$GENEID <- NULL # CharacterList breaks data.frame building
  if (length(mutgenes) == 0) {
    stop("None of your clonal mutations are in your gene list. Aborting.")
  }

  # template for fishplot annotation
  mutations <- data.frame(symbol=c(), 
                          chrom=c(), 
                          coord=c(), 
                          clone_id=c(), 
                          timepoint=c(), 
                          type=c(),
                          fraction=c(),
                          VAF=c()) 
  for (nm in names(mcols(mutgenes))) mutations[,nm] <- c() 

  # horrible kludge below; swap out mutgenes for actual coords/impact ASAP!
  addmut <- function(sym, cid, timept, VAF, type="NS", mutations, mutgenes) {
    entry <- data.frame(
      symbol=sym, 
      chrom=as.character(seqnames(mutgenes[sym])),
      coord=start(mutgenes[sym]),
      clone_id=cid,
      timepoint=timept,
      type=type,
      fraction=paste0(round(VAF*100, 1), "%"),
      VAF=VAF
    )
    for (nm in names(mcols(mutgenes))) entry[,nm] <- mcols(mutgenes[sym])[[nm]]
    rbind(entry, mutations)
  }

  # this is also disgusting 
  for (cid in seq_along(clonalmut)) {
    for (timept in colnames(VAF)) {
      subVAF <- subset(VAF, VAF[, timept] > 0) 
      clonal <- intersect(clonalmut[[cid]], rownames(subVAF))
      clonal <- intersect(clonal, names(mutgenes))
      if (length(clonal) > 0) {
        for (sym in clonal) {
          if (subVAF[sym, timept] > 0) { 
            # message(sym, " has VAF ", subVAF[sym, timept], " in clone ", cid)
            mutations <- addmut(sym, cid, timept, subVAF[sym, timept], 
                                mutations=mutations, mutgenes=mutgenes)
          }
        }
      }
    }
  }
  return(mutations)

}
