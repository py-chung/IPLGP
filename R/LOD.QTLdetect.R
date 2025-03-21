#' QTL Detect by LOD
#'
#' Detect QTL by the likelihood of odds (LOD) matrix.
#'
#' @param LOD matrix. The LOD matrix, which is a t*p matrix, where t is
#' the number of traits and p is the number of bins on the chromosomes.
#' Missing values should be denoted as NA in the matrix.
#' @param bin matrix. An n*2 matrix that represents the number of bins on
#' each chromosome, where n is the number of chromosomes. The first column
#' denotes the chromosome number, and the second column denotes the number
#' of bins on that chromosome. It's important to ensure that chromosomes
#' are divided in order.
#' @param thre numeric. The LOD threshold. Any LOD score under this
#' threshold will be calculated as 0.
#' @param QTLdist numeric. The minimum distance (in bins) among different
#' linked significant QTL.
#' @param console logical. Determines whether the process of the algorithm
#' will be displayed in the R console or not.
#'
#' @return
#' \item{detect.QTL.number}{The number of detected QTL in each trait.}
#' \item{QTL.matrix}{The QTL position matrix. Where the elements 1
#' donates the position of QTL; elements 0 donate the bins whose LOD
#' score is under the LOD threshold; other positions are shown as NA.}
#' \item{EQF.matrix}{The matrix denotes the EQF value of each bin.}
#' \item{linkage.QTL.number}{The linkage QTL number of all detected
#' QTL. In other words, it is the table that denote how many QTL are
#' on one chromosome.}
#' \item{LOD.threshold}{The LOD threshold used in this analysis.}
#' \item{bin}{The bin information matrix used in this analysis.}
#'
#' @export
#'
#' @references
#'
#' Wu, P.-Y., M.-.H. Yang, and C.-H. KAO 2021 A Statistical Framework
#' for QTL Hotspot Detection. G3: Genes, Genomes, Genetics: jkab056. <doi: 10.1093/g3journal/jkab056>
#'
#' @seealso
#' \code{\link[QTLEMM]{EQF.permu}}
#' \code{\link[QTLEMM]{EQF.plot}}
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "LODexample.RDATA", package = "QTLEMM"))
#' dim(LODexample) # 100 traits, 633 bins on chromosome
#'
#' # run and result
#' result <- LOD.QTLdetect(LODexample, bin, thre = 3, QTLdist = 10)
#' result$detect.QTL.number
#'
LOD.QTLdetect <- function(LOD, bin, thre = 3, QTLdist = 20, console = TRUE){

  if(is.null(LOD) | is.null(bin)){
    stop("Input data is missing, please cheak and fix.", call. = FALSE)
  }

  datatry <- try(LOD*LOD, silent=TRUE)
  if(class(datatry)[1] == "try-error" | class(LOD)[1] != "matrix"){
    stop("LOD data error, please cheak your LOD data.", call. = FALSE)
  }

  bintest <- c(ncol(bin) != 2, NA %in% bin, bin[,1] != sort(bin[,1]), sum(bin[, 2]) != ncol(LOD))
  datatry <- try(bin*bin, silent=TRUE)
  if(class(datatry)[1] == "try-error" | TRUE %in% bintest){
    stop("bin data error, or the number of bin does not match the LOD data.", call. = FALSE)
  }

  datatry <- try(bin*bin, silent=TRUE)
  if(class(datatry)[1] == "try-error" | ncol(bin) != 2 | NA %in% bin | sum(bin[,2]) != ncol(LOD)){
    stop("bin data error or the number of bins is not match the LOD data, please cheak and fix.", call. = FALSE)
  }

  if(!is.numeric(thre) | length(thre) > 1 | min(thre) < 0){
    stop("Parameter thre error, please input a positive number.", call. = FALSE)
  }

  if(!is.numeric(QTLdist) | length(QTLdist) > 1 | min(QTLdist) < 0 | max(QTLdist) > min(bin[, 2])){
    stop("Parameter QTLdist error, please input a positive integer. Or the number is too big, please input a smaller number.", call. = FALSE)
  }

  if(!console[1] %in% c(0,1) | length(console) > 1){console <- TRUE}

  nt <- nrow(LOD)
  ns <- ncol(LOD)
  nc <- nrow(bin)

  cr0 <- c()
  for(i in 1:nc){
    cr0 <- c(cr0, rep(i, bin[i, 2]))
  }

  lcr <- bin[, 2]
  ncr <- c()
  for(i in 1:nc){
    ncr[i] <- sum(bin[1:i, 2])
  }
  colnames(bin) <- c("chromosome", "# of bins")

  LOD[is.na(LOD)] <- 0
  LOD[LOD < thre] <- 0

  det <- matrix(0, nt, ns)
  cat("step\tprocess \n"[console])
  t0 <- Sys.time()
  for(j in 1:nt){
    cat1 <- paste("detect", paste(j, nt, sep = "/"), "\n", sep = "\t")
    if(Sys.time()-t0 > 1){
      cat(cat1[console])
      t0 <- Sys.time()
    }

    k1 <- c()
    for(i in 1:nc){
      x <- LOD[j, cr0 == i]
      l <- lcr[i]
      k0 <- rep(0, l)
      for(m in 1:l){
        x0 <- x[m]
        if(x0 == 0){det[m] <- 0
        } else {
          if(m %in% c(1:QTLdist)){
            maxm <- m+QTLdist
            if(maxm > l){maxm <- l}
            if(x0 < max(x[1:(maxm)])){
              k0[m] <- NA
            } else {
              k0[m] <- 1
            }
          } else if(m %in% c((l-QTLdist+1):l)){
            if(x0 < max(x[(m-QTLdist):l])){
              k0[m] <- NA
            } else {
              k0[m] <- 1
            }
          } else {
            if(x0<max(x[(m-QTLdist):(m+QTLdist)])){
              k0[m] <- NA
            } else {
              k0[m] <- 1
            }
          }
        }
      }
      k1 <- c(k1, k0)
    }
    det[j,] <- k1
  }
  cat(cat1[console])

  summ <- function(x){sum(x, na.rm = TRUE)}
  qtldetect <- apply(det, 1, summ)

  link <- matrix(0, nt, nc)
  cat("step\tprocess \n"[console])
  t0 <- Sys.time()
  for(j in 1:nt){
    cat1 <- paste("linkage", paste(j, nt, sep = "/"), "\n", sep = "\t")
    if(Sys.time()-t0 > 1){
      cat(cat1[console])
      t0 <- Sys.time()
    }
    for(i in 1:nc){
      x <- det[j, cr0 == i]
      a <- table(x)
      if(length(a) == 2){link[j, i] <- a[2]}
    }
  }
  cat(cat1[console])
  link <- table(link)

  EQF <- matrix(0, nt, ns)
  t0 <- Sys.time()
  for(i in 1:nt){
    cat1 <- paste("EQF caculating", paste(i, nt, sep = "/"), "\n", sep = "\t")
    if(Sys.time()-t0 > 1){
      cat(cat1[console])
      t0 <- Sys.time()
    }
    QTL <- det[i,]
    if(!1 %in% QTL){next}
    lod <- LOD[i,]
    ka <- c()
    for(j in 1:nc){
      QTLcr <- QTL[cr0 == j]
      n <- length(QTLcr)
      k0 <- rep(0, n)
      if(!1 %in% QTLcr){
        ka <- c(ka, k0)
        next
      }
      QTLlo <- which(QTLcr == 1 & !is.na(QTLcr))
      QTLci <- lod[cr0 == j]
      for(k in 1:length(QTLlo)){
        pii <- QTLlo[k]
        QTLci1 <- rep(-1, length(QTLci))
        QTLci1[QTLci>(QTLci[pii]-1) & QTLci>thre] <- 1
        QTLci2 <- QTLci1*c(-1, (QTLci1[-n]))
        QTLci3 <- which(QTLci2 == -1)
        QTLci4 <- sort(c(QTLci3, pii))
        QTLci5 <- QTLci4-c(QTLci3, (n+1))
        ci0 <- min(which(QTLci5<0))
        ci1 <- QTLci3[ci0-1]
        ci2 <- n
        if(ci0 <= length(QTLci3)){ci2 <- QTLci3[ci0]-1}
        sdi <- ((ci2-ci1)/(2*1.96))
        k1 <- stats::pnorm((1:n)+0.5, pii, sdi)-stats::pnorm((1:n)-0.5, pii, sdi)
        if(length(k1) == 0){
          k1 <- 0
        }
        k0 <- k0+k1
      }
      ka <- c(ka, k0)
    }
    EQF[i,] <- ka
  }
  cat(cat1[console])

  EQF.all <- apply(EQF, 2, sum)

  return(list(detect.QTL.number = qtldetect, QTL.matrix = det, EQF.matrix = EQF,
              linkage.QTL.number = link, LOD.threshold = thre, bin = bin))
}
