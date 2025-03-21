#' EQF Permutation
#'
#' The EQF matrix cluster permutation process for QTL hotspot detection.
#'
#' @param LOD.QTLdetect.result list. The data list of the output from
#' LOD.QTLdetect().
#' @param ptime integer. The permutation times.
#' @param alpha numeric. The type 1 error rate of detecting the hotspot.
#' @param Q logical. When set to TRUE, the function will additionally carry
#' out the permutation of the Q method as the control group, which will be
#' indicated as 'B' in the output.
#' @param console logical. Determines whether the process of the algorithm
#' will be displayed in the R console or not.
#'
#' @return
#' \item{EQF.matrix}{The matrix denotes the EQF value of each bin.}
#' \item{bin}{The bin information matrix used in this analysis.}
#' \item{LOD.threshold}{The LOD threshold used in this analysis.}
#' \item{cluster.number}{The number of QTLs in each cluster group.}
#' \item{cluster.id}{The serial number of traits in each cluster group.}
#' \item{cluster.matrix}{The new EQF matrix after the clustering process.}
#' \item{permu.matrix.cluster}{The permutation result of the clustering
#' method, which has been sorted by order.}
#' \item{permu.matrix.Q}{The permutation result of the Q method, which has
#' been sorted by order.}
#' \item{EQF.threshold}{The EQF threshold is calculated from the
#' permutation process.}
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
#' \code{\link[QTLEMM]{EQF.plot}}
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "LODexample.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' result <- EQF.permu(LOD.QTLdetect.result, ptime = 50)
#' result$cluster.number
EQF.permu <- function(LOD.QTLdetect.result, ptime = 1000, alpha = 0.05, Q = TRUE, console = TRUE){

  datatest <- names(LOD.QTLdetect.result) != c("detect.QTL.number", "QTL.matrix", "EQF.matrix",
                                               "linkage.QTL.number", "LOD.threshold", "bin")
  if(TRUE %in% (datatest) | length(datatest) != 6){
    stop("Input data error, please input the original output data of LOD.QTLdetect.", call. = FALSE)
  }

  if(!is.numeric(ptime) | length(ptime) > 1 | min(ptime) < 0){
    stop("Parameter ptime error, please input a positive integer.", call. = FALSE)
  }

  if(!is.numeric(alpha) | length(alpha) > 1 | min(alpha) < 0 | max(alpha) > 1){
    stop("Parameter alpha error, please input a positive number between 0 and 1.", call. = FALSE)
  }

  if(!Q[1] %in% c(0,1) | length(Q) > 1){Q <- TRUE}

  if(!console[1] %in% c(0,1) | length(console) > 1){console <- TRUE}

  dat <- LOD.QTLdetect.result
  detect <- dat$QTL.matrix
  EQF <- dat$EQF.matrix
  nt <- nrow(detect)
  ns <- ncol(detect)
  thre <- dat$LOD.threshold
  bin <- dat$bin
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

  clu.det <- detect
  clu.det[is.na(clu.det)] <- 0
  rownames(clu.det) <- 1:nt
  det.all <- apply(clu.det, 1, sum)
  clu.det <- clu.det[det.all > 0,]

  clu.eqf <- EQF
  clu.eqf[is.na(clu.det)] <- 0
  rownames(clu.eqf) <- 1:nt
  clu.eqf <- clu.eqf[det.all > 0,]

  qtlind <- as.numeric(rownames(clu.det))
  cluster <- list()

  maxeqf <- 1
  add.group <- c()
  cat("cluster group \n"[console])
  while(class(clu.eqf)[[1]] == "matrix"){
    peqf <- c()
    clu.eqf1 <- apply(clu.eqf, 2, sum)
    peqf <- c(peqf, which.max(clu.eqf1))
    if(max(clu.eqf1) == 0){
      break
    }
    group <- which(clu.det[, peqf] == 1)
    if(length(group) == 0){
      clu.eqf[, peqf] <- 0
      next
    }
    add <- 1
    while(add > 0){
      group1 <- group
      group.det <- clu.det[group1,]
      if(length(group1) > 1){
        clu.det1 <- apply(group.det, 2, sum)
      } else {clu.det1 <- group.det}
      group.pos <- which(clu.det1 > 0)
      group.pos <- group.pos[!group.pos %in% peqf]
      peqf <- c(peqf, group.pos)
      if(length(group.pos) > 0){
        for(i in 1:length(group.pos)){
          group2 <- which(clu.det[, group.pos[i]] == 1)
          group1 <- c(group1, group2[!group2 %in% group1])
        }
      }
      add <- length(group1)-length(group)
      group <- group1
    }
    ind <- sort(as.numeric(rownames(clu.det)[group]))
    cluster[[maxeqf]] <- ind
    add.group <- c(add.group, ind)
    clu.det <- clu.det[-group,]
    clu.eqf <- clu.eqf[-group,]
    if(maxeqf%%10 == 0){
      cat(paste(maxeqf, "\n")[console])
    }
    maxeqf <- maxeqf+1
  }
  cat(paste(maxeqf, "\n")[console])

  cluster[[maxeqf]] <- qtlind[!qtlind %in% add.group]
  Lagend <- c()
  for(i in 1:maxeqf){
    Lagend[i+2] <- length(cluster[[i]])
  }
  Lagend[1] <- maxeqf
  Lagend[2] <- sum(Lagend[3:(maxeqf+2)])
  names(Lagend) <- c("ngroup", "nQTL", (paste("# g", 1:maxeqf, sep = "")))

  clumatrix <- matrix(0, maxeqf, ns)
  for(i in 1:maxeqf){
    submatrix <- EQF[cluster[[i]],]
    if(class(submatrix)[[1]] == "matrix"){
      clumatrix[i,] <- apply(submatrix, 2, sum)
    } else {clumatrix[i,] <- submatrix}
  }

  eqf.premutation <- function(eqfmatrix, ptime = ptime, console = console){
    t0 <- Sys.time()
    cat("permutation time \n"[console])
    n <- ncol(eqfmatrix)
    g <- nrow(eqfmatrix)
    out <- matrix(0, ptime, n)

    for(i in 1:ptime){
      p.matrix <- matrix(0, g, n)
      for(j in 1:g){
        eqfpos <- eqfmatrix[j, eqfmatrix[j,] > 0]
        k <- sample(1:n, length(eqfpos))
        p.matrix[j, k] <- eqfpos
      }
      out[i,] <- sort(apply(p.matrix, 2, sum), decreasing = TRUE)
      cat1 <- paste(paste(i, ptime, sep = "/"), "\n", sep = "\t")
      if(Sys.time()-t0 > 5){
        cat(cat1[console])
        t0 <- Sys.time()
      }
    }
    cat(paste(paste(ptime, ptime, sep = "/"), "\n", sep = "\t")[console])
    return(out)
  }

  cat("cluster group permutation \n"[console])
  permu.clu <- eqf.premutation(clumatrix, ptime, console)

  sortEQF <- function(x, a = ceiling(ptime*alpha)){
    return(sort(x, decreasing = TRUE)[a])
  }

  thre.clu <- apply(permu.clu, 2, sortEQF)

  a1 <- thre.clu[1]
  a2 <- thre.clu[2]
  a5 <- thre.clu[5]
  a20 <- thre.clu[20]
  eqfthre <- c(a20, a5, a2, a1)
  names(eqfthre) <- c(20, 5, 2, 1)

  permu.q <- NULL
  if(Q){
    cat("Q-method permutation \n"[console])
    q0 <- EQF[apply(EQF, 1, sum) > 0,]
    permu.q <- eqf.premutation(q0, ptime, console)
    thre.q <- apply(permu.q, 2, sortEQF)
    b1 <- thre.q[1]
    beta <- max(which(b1<thre.clu))
    eqfthre <- c(b1, a20, a5, a2, a1)
    names(eqfthre) <- c(paste("B:", beta), 20, 5, 2, 1)
  }

  eqfthre <- round(eqfthre, 1)

  eqf.all <- apply(EQF, 2, sum)
  real.spot.num <- c()
  for(k in 1:length(eqfthre)){
    real.spot.num[k] <- length(which(eqf.all > eqfthre[k]))
  }

  eqfthre <- cbind(eqfthre, real.spot.num)

  return(list(EQF.matrix = EQF, bin = bin, LOD.threshold = thre, cluster.number = Lagend, cluster.id = cluster,
              cluster.matrix = clumatrix, permu.matrix.cluster = permu.clu, permu.matrix.Q = permu.q,
              EQF.threshold = eqfthre))
}
