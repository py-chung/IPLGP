#' EQF plot
#'
#' Generate an EQF plot based on the result of the permutation process
#' used to detect the QTL hotspot.
#'
#' @param result list. The data list of the output from LOD.QTLdetect(),
#' EQF.permu(), or Qhot.EQF().
#' @param plot.all logical. When set to TRUE, it directs the function to
#' output one figure of the EQF values over the bins.
#' @param plot.chr logical. When set to TRUE, it instructs the function to
#' output the figures of the EQF values over the bins for each chromosome.
#' @param plot.main logical of character. When set to TRUE, it will use the
#' default title on the plot. When set to FALSE, it will be no title on the
#' plot. Users can also use this argument to set their own title.
#'
#' @return
#'
#' One or several EQF plots.
#'
#' @export
#'
#' @references
#'
#' Wu, P.-Y., M.-.H. Yang, and C.-H. KAO 2021 A Statistical Framework
#' for QTL Hotspot Detection. G3: Genes, Genomes, Genetics: jkab056. <doi: 10.1093/g3journal/jkab056>
#'
#' @seealso
#' \code{\link[QTLEMM]{LOD.QTLdetect}}
#' \code{\link[QTLEMM]{EQF.permu}}
#' \code{\link[QTLEMM]{Qhot.EQF}}
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "LODexample.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' EQF.plot(LOD.QTLdetect.result)
#' EQF.plot(EQF.permu.result)
EQF.plot <- function(result, plot.all = TRUE, plot.chr = TRUE, plot.main = TRUE){

  dat <- result
  name0 <- names(dat)
  if(length(name0) == 6){
    datatest <- name0 != c("detect.QTL.number", "QTL.matrix", "EQF.matrix",
                           "linkage.QTL.number", "LOD.threshold", "bin")
  } else if(length(name0) == 9){
    datatest <- name0 != c("EQF.matrix", "bin", "LOD.threshold", "cluster.number", "cluster.id",
                           "cluster.matrix", "permu.matrix.cluster", "permu.matrix.Q", "EQF.threshold")
  } else if(length(name0) == 10){
    datatest <- name0 != c("EQF.matrix", "bin", "bin.size", "EQF.trait", "EQF.detect", "EQF.nondetect",
                           "cluster.matrix", "permu.matrix.cluster", "permu.matrix.Q", "EQF.threshold")
  } else {
    stop("Input data error, please input the original output data of LOD.QTLdetect, EQF.permu ,or Qhot.EQF.", call. = FALSE)
  }

  if(TRUE %in% (datatest)){
    stop("Input data error, please input the original output data of LOD.QTLdetect, EQF.permu ,or Qhot.EQF.", call. = FALSE)
  }

  if(!plot.all[1] %in% c(0,1) | length(plot.all) > 1){plot.all <- TRUE}
  if(!plot.chr[1] %in% c(0,1) | length(plot.chr) > 1){plot.chr <- TRUE}
  if(length(plot.main) > 1){plot.main <- TRUE}

  EQF <- dat$EQF.matrix
  clumatrix <- dat$cluster.matrix
  thre <- dat$LOD.threshold
  eqfthre <- dat$EQF.threshold
  bin <- dat$bin
  bin.size <- dat$bin.size
  nc <- nrow(bin)
  lcr <- bin[, 2]
  ncr <- c()
  for(i in 1:nc){
    ncr[i] <- sum(bin[1:i, 2])
  }
  cr0 <- c()
  for(i in 1:nc){
    cr0 <- c(cr0, rep(i, bin[i, 2]))
  }
  eqf.all <- apply(EQF, 2, sum)

  gate <- ceiling(length(cr0)/nc/2)
  x0 <- 1:length(cr0)+(cr0-1)*gate
  xn <- ncr+(0:(nc-1))*gate-lcr/2

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  yli <- max(eqf.all)*1.2
  xli <- max(x0)/40
  yax <- 2
  if(yli > 20){
    yax <- 5
  }
  if(plot.all & nc>1){
    graphics::par(mfrow = c(1, 1))
    graphics::par(mai = c(1, 1, 1, 1))

    if(plot.main==TRUE){
      if(is.null(bin.size)){
        ma1 <- paste("LOD threshold =", thre)
        ma <- "EQF plot from LOD data"
      } else {
        ma1 <- paste("bin size =", bin.size)
        ma <- "EQF plot from flanking marker data"
      }
      if(length(clumatrix) > 0){
        ma2 <- paste("  # of group =", nrow(clumatrix))
        graphics::par(mai = c(1, 1, 1, 1.5))
      } else {
        ma2 <- NULL
      }
      masb <- paste(ma1, ma2)
    } else if(plot.main==FALSE) {
      ma <- NULL
      masb <- NULL
    } else {
      ma <- plot.main
      masb <- NULL
    }

    plot(x0, eqf.all, type = "h", ylab = "EQF", xlab = "chromosome", main = ma,
         xaxt = "n", ylim = c(0, yli), yaxt = "n",  cex.main = 1.5, cex.lab = 1.2, axes = FALSE)
    graphics::mtext(masb, side = 3, line = 0.5, cex = 1.2)
    graphics::axis(side = 1, pos = -yli/35, at = xn, labels = 1:nc, cex.axis = 1.2, tick = FALSE)
    lse <- 4000/max(x0)
    graphics::segments(x0, rep(-yli/35, length(x0)), x0, rep(-yli/70, length(x0)), lwd = lse)
    graphics::segments(-xli, -yli/35, max(x0)+xli, -yli/35)
    graphics::segments(-xli, -yli/35, -xli, yli)
    graphics::segments(-xli, yli, max(x0)+xli, yli)
    graphics::segments(max(x0)+xli, -yli/35, max(x0)+xli, yli)
    if(length(eqfthre) > 0){
      for(j in 1:nrow(eqfthre)){
        graphics::axis(2, eqfthre[j, 1], las = 2)
        graphics::axis(4, paste(rownames(eqfthre)[j], " (", eqfthre[j, 2], ")", sep = ""), at = eqfthre[j], las = 2)
      }
      graphics::abline(h = eqfthre[, 1], col = "red")
    } else {graphics::axis(2, seq(0, yli, 5), las = 2)}
  }

  if(plot.chr | (nc == 1 & plot.all)){
    if(nc >= 4){
      graphics::par(mai = c(0.8, 0.8, 0.8, 0.8))
      if(length(clumatrix) > 0){
        graphics::par(mai = c(0.8, 0.8, 0.8, 1.2))
      }
    } else {
      graphics::par(mai = c(1, 1, 1, 1.5))
      if(length(clumatrix) > 0){
        graphics::par(mai = c(1, 1, 1, 1))
      }
    }

    if(plot.main==TRUE){
      if(is.null(bin.size)){
        ma1 <- paste(", LOD threshold =", thre)
        ma <- "EQF plot from LOD data"
      } else {
        ma1 <- paste(", bin size =", bin.size)
        ma <- "EQF plot from flanking marker data"
      }
    } else if(plot.main==FALSE) {
      ma <- NULL
      ma1 <- NULL
    } else {
      ma <- plot.main
      ma1 <- NULL
    }

    for(i in 1:nc){
      eqf <- eqf.all[cr0 == i]
      plot(eqf, type = "h", ylab = "EQF", xlab = "position (bin)", main = ma,
           ylim = c(0, max(eqf.all)*1.2), yaxt = "n", cex.main = 1.5, cex.lab = 1.2)
      masb <- paste("chr ", i, ma1, sep = "")
      graphics::mtext(masb, side = 3, line = 0.3, cex = 1.2)
      if(length(eqfthre)>0){
        for(k in 1:nrow(eqfthre)){
          graphics::axis(2, eqfthre[k, 1], las = 2.5)
          graphics::axis(4, paste(rownames(eqfthre)[k], " (", eqfthre[k, 2], ")", sep = ""), at = eqfthre[k], las = 2.5)
        }
        graphics::abline(h = eqfthre[, 1], col = "red")
      } else {graphics::axis(2, seq(0, yli, 5), las = 2.5)}
    }
  }
}


