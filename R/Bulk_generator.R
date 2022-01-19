#' Pseudo Bulk generator from single cell data
#'
#' @param sce numeric matrix or data.frame containing the single cell data. columns : cells, rows : genes
#' @param phenoData data.frame containing cellID and cellType as character
#' @param Num.mixtures integer number of pseudo bulk to generate (default = 1000)
#' @param nb_CT_random whether the cells types present in pseudo bulk should be selected randomly or not. If not, all the cell types will be kept (default = TRUE)
#' @param pool.size integer number of single cell per mixture (default = 100)
#' @param min.percentage minimum percentage of cells to put in a mixture (if seleted) (default = 1)
#' @param max.percentage maximum percentage of cells to put in a mixture (if seleted) (default = 99)
#' @param seed random seed (default = 24)
#' @param forbidden_CT character vector indicating cellType not to put in the mixture.
#'
#' @return list with T a data.frame of the Num.mixture pseudo bulk and P a data.frame of the cellType proportion used to build the mixtures.
#' @import dplyr
#' @importFrom stats runif
#' @importFrom utils ?
#' @export
#' @examples Bulk_generator(sc_counts, phenoData)

Bulk_generator <- function(sce, phenoData, Num.mixtures = 1000, nb_CT_random = TRUE, pool.size = 100, min.percentage = 1, max.percentage = 99, seed = 24, forbidden_CT = NULL){

  CT = unique(phenoData$cellType)
  ?stopifnot(length(CT) >= 2)

  set.seed(seed)

  cell.distribution = data.frame(table(phenoData$cellType),stringsAsFactors = FALSE)
  colnames(cell.distribution) = c("CT","max.n")

  Tissues = list()
  Proportions = list()

  for(y in 1:Num.mixtures){
    #Only allow feasible mixtures based on cell distribution
    while(!exists("P")){
      if (!is.null(forbidden_CT)){
        CT = CT[!(CT %in% forbidden_CT)]
      }
      if (nb_CT_random){
        num.CT.mixture = sample(x = 2:length(CT),1)
      }
      else {
        num.CT.mixture = length(CT)
      }

      selected.CT = sample(CT, num.CT.mixture, replace = FALSE)

      P = runif(num.CT.mixture, min.percentage, max.percentage)
      P = round(P/sum(P), digits = log10(pool.size))  #sum to 1
      P = data.frame(CT = selected.CT, expected = P, stringsAsFactors = FALSE)

      missing.CT = CT[!CT %in% selected.CT]
      missing.CT = data.frame(CT = missing.CT, expected = rep(0, length(missing.CT)), stringsAsFactors = FALSE)

      P = rbind.data.frame(P, missing.CT)
      potential.mix = merge(P, cell.distribution)
      potential.mix$size = potential.mix$expected * pool.size

      if( !all(potential.mix$max.n >= potential.mix$size) | sum(P$expected) != 1){
        rm(list="P")
      }

    }
    # Using info in P to build T simultaneously
    chosen_cells <- sapply(which(P$expected != 0), function(x){

      n.cells = P$expected[x] * pool.size
      chosen = sample(phenoData$cellID[phenoData$cellType == P$CT[x]],
                      n.cells)

      chosen
    }) %>% unlist()

    T <- apply(sce[,colnames(sce) %in% chosen_cells], 1, function(x) sum(x)) %>% as.data.frame()
    colnames(T) = paste0("mix",y)

    row.names(P) = P$CT
    P = P[CT,]
    P = data.frame(x = P$expected)
    colnames(P) = paste0("mix",y)
    row.names(P) = CT

    Tissues[[y]] <- T
    Proportions[[y]] <- P

    rm(list=c("T","P","chosen_cells","missing.CT"))
  }

  P = do.call(cbind.data.frame, Proportions)
  T = do.call(cbind.data.frame, Tissues)

  return(list(T = T, P = P))
}
