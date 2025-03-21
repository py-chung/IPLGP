#' EQF Matrix Conversion
#'
#' Convert the QTL flanking marker data to EQF matrix. And the EQF matrix 
#' cluster permutation process can be further carried out to detect QTL 
#' hotspots.
#' 
#' @param DataQTL data.frame. A data frame with 5 columns for QTL information.
#' The columns represent the serial number of QTLs, the trait names, the
#' chromosome numbers, the left flanking marker positions(in cM) of QTLs, and
#' the right flanking marker positions(in cM) of QTLs.
#' @param Datachr vector. The length of each chromosome(in cM).
#' @param bin.size numeric. The bin size(in cM) for QTL detection. If the 
#' distance of flanking marker of a QTL is lass than the bin size, it will be
#' mark in the EQF matrix and will participate in the cluster grouping process.
#' @param permu logical. When set to TRUE, the function will carry out the 
#' cluster grouping process and cluster group permutation. 
#' @param ptime integer. The permutation times.
#' @param alpha numeric. The type 1 error rate of detecting the hotspot.
#' @param Q logical. When set to TRUE, the function will additionally carry
#' out the permutation of the Q method as the control group, which will be
#' indicated as 'B' in the output.
#' @param console logical. Determines whether the process of the algorithm
#' will be displayed in the R console or not.
#' 
#' @return
#' \item{EQF.matrix}{The matrix denotes the EQF value of each bin for every 
#' QTL in this database.}
#' \item{bin}{The bin information matrix whose first column denotes the 
#' chromosome number and the second column denotes the number of bins on 
#' that chromosome.}
#' \item{bin.size}{The bin size set in this analysis.}
#' \item{EQF.trait}{The matrix denotes the EQF value of each bin for every 
#' trait of this database.}
#' \item{EQF.detect}{The matrix denotes the EQF value of each bin for the 
#' trait that have the QTL detected in the set bin size.}
#' \item{EQF.nondetect}{The matrix denotes the EQF value of each bin for the 
#' trait that have no QTL detected in the set bin size.}
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
#' \code{\link[QTLEMM]{Qhot}}
#' \code{\link[QTLEMM]{EQF.plot}}
#' 
#' @examples
#' # load the example data
#' load(system.file("extdata", "QHOTEQFexample.RDATA", package = "QTLEMM"))
#' 
#' #' # run and result
#' result <- Qhot.EQF(QTL.example, chr.example, bin.size = 2, permu = TRUE, 
#' ptime = 100, alpha = 0.05, Q = FALSE)
Qhot.EQF <- function(DataQTL, Datachr, bin.size = 0.5, permu = TRUE, 
                     ptime = 1000, alpha = 0.05, Q = TRUE, console = TRUE){
  
  DataQTL <- data.frame(DataQTL)
  if(ncol(DataQTL) != 5 | TRUE %in% is.na(DataQTL)){
    stop("Input data DataQTL error, please check and fix.", call. = FALSE)
  }
  colnames(DataQTL) <- c("X", "Trait", "chr", "L", "R")
  
  if(!is.numeric(DataQTL[, 4]) | !is.numeric(DataQTL[, 5])){
    stop("Input data DataQTL error, please check and fix.", call. = FALSE)
  }
  
  if(length(table(DataQTL[, 1])) != nrow(DataQTL)){
    DataQTL <- DataQTL[order(DataQTL[, 3], DataQTL[, 4]),]
    DataQTL[, 1] <- 1:nrow(DataQTL)
  }
  
  datatry <- try(Datachr%*%(1:max(DataQTL[,3])), silent=TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% Datachr){
    stop("Input data Datachr error, please check and fix.", call. = FALSE)
  }
  
  if(!is.numeric(bin.size) | length(bin.size) > 1 | min(bin.size) < 0){
    stop("Parameter bin.size error, please input a positive number.", call. = FALSE)
  }
  
  if(!permu[1] %in% c(0,1) | length(permu) > 1){permu <- TRUE}
  
  if(!is.numeric(ptime) | length(ptime) > 1 | min(ptime) < 0){
    stop("Parameter ptime error, please input a positive integer.", call. = FALSE)
  }
  
  if(!is.numeric(alpha) | length(alpha) > 1 | min(alpha) < 0 | max(alpha) > 1){
    stop("Parameter alpha error, please input a positive number between 0 and 1.", call. = FALSE)
  }
  
  if(!Q[1] %in% c(0,1) | length(Q) > 1){Q <- TRUE}
  
  if(!console[1] %in% c(0,1) | length(console) > 1){console <- TRUE}
  
  lcr <- ceiling(Datachr/bin.size)
  nc <- length(Datachr)
  cr0 <- c()
  for(i in 1:nc){
    cr0 <- c(cr0, rep(i, lcr[i]))
  }
  ncr <- c()
  for(i in 1:nc){
    ncr[i] <- sum(lcr[1:i])
  }
  n0 <- nrow(DataQTL)
  nbin <- sum(lcr)
  
  scr <- cumsum(lcr)
  EQF <- matrix(0, n0, nbin)
  tq <- c()
  cat("detection step \n"[console])
  t0 <- Sys.time()
  for(i in 1:n0){
    crq <- DataQTL[i, 3]
    gqj <- (floor(DataQTL[i, 4]/bin.size)+1):ceiling(DataQTL[i, 5]/bin.size)
    if((DataQTL[i,5]-DataQTL[i, 4])<bin.size){
      gqj <- mean(gqj)
      EQF[i, scr[crq]-lcr[crq]+gqj] <- 1
      tq <- rbind(tq, c(i, DataQTL[i, 2]))
    } else {
      pii <- mean(gqj)
      sdi <- (length(gqj)/(2*1.96))
      k1 <- stats::pnorm((1:lcr[crq])+0.5, pii, sdi)-stats::pnorm((1:lcr[crq])-0.5, pii, sdi)
      EQF[i, (scr[crq]-lcr[crq]+1):scr[crq]] <- k1/sum(k1)
    }
    ti <- as.numeric(floor((Sys.time()-t0)))
    cat(paste(paste(i, n0, sep = "/"), "\n")[ti*console])
    t0 <- t0+(Sys.time()-t0)*ti
  }
  cat(paste(paste(n0, n0, sep = "/"), "\n")[console])
  rownames(EQF) <- DataQTL[, 2]
  EQFq <- EQF[as.numeric(tq[, 1]),]
  eqf.all <- apply(EQF, 2, sum)
  
  tt <- table(rownames(EQF))
  eqf.t <- matrix(0, length(tt), nbin)
  row.names(eqf.t) <- names(tt)
  for(i in 1:length(tt)){
    et <- EQF[rownames(EQF)==names(tt)[i],]
    if(is.null(dim(et))){
      ett <- et
    } else {
      ett <- apply(et, 2, sum) 
    }
    eqf.t[row.names(eqf.t)==names(tt)[i],] <- ett
  }
  
  eqf.c <- NULL
  permu.clu <- NULL
  permu.q <- NULL
  eqfthre <- NULL
  
  if(permu){
    eqf.tnq <- eqf.t[!row.names(eqf.t)%in%tq[, 2],]
    eqf.c0 <- eqf.tq <- eqf.t[row.names(eqf.t)%in%tq[, 2],]
    
    detect <- c()
    for(i in 1:nbin){
      detect[i] <- length(which(EQF[, i]==1))
    }
    od <- order(detect, decreasing = TRUE)
    od <- od[1:length(detect[detect>0])]
    
    tp0 <- matrix(0, length(od), nrow(eqf.c0))
    wn <- which(row.names(eqf.c0) %in% (row.names(EQF)[EQF[, od[1]]==1]))
    tp0[1,wn] <- 1
    cat("cluster grouping \n"[console])
    tgn <- 1
    for(i in 2:length(od)){
      wn <- which(row.names(eqf.c0) %in% (row.names(EQF)[EQF[, od[i]]==1]))
      tp0[i,wn] <- 1
      tgn <- tgn+1
      for(j in 1:i){
        if(sum(tp0[i,]*tp0[j,])>0 & j<i){
          tp0[i,] <- tp0[i,]+tp0[j,]
          tp0[i,tp0[i,]>0] <- 1
          tp0[j,] <- 0
          tgn <- tgn-1
        }
      }
      tgni <- floor(tgn/10)
      cat(paste(tgn, "groups", "\n")[tgni*console])
    }
    cat(paste(tgn, "groups", "\n")[console])
    tp0 <- tp0[apply(tp0,1,sum)>0,]
    
    if(is.null(dim(tp0))){
      eqf.c0 <- apply(eqf.tq[tp0,], 2, sum)
    } else {
      eqf.c0 <- matrix(0, nrow(tp0), nbin)
      for(i in 1:nrow(eqf.c0)){
        tp1 <- which(tp0[i,]==1)
        if(length(tp1)>1){
          eqf.c0[i,] <- apply(eqf.tq[tp1,], 2, sum)
        } else {
          eqf.c0[i,] <- eqf.tq[tp1,]
        }
      }
    }
    eqf.c <- rbind(eqf.c0, eqf.tnq)
    
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
    permu.clu <- eqf.premutation(eqf.c, ptime = ptime, console = console)
    
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
    
    if(Q){
      q0 <- EQF
      cat("Q-method permutation \n"[console])
      permu.q <- eqf.premutation(q0, ptime = ptime, console = console)

      
      thre.q <- apply(permu.q, 2, sortEQF)
      b1 <- thre.q[1]
      beta <- max(which(b1<thre.clu))
      eqfthre <- c(b1, a20, a5, a2, a1)
      names(eqfthre) <- c(paste("B:", beta), 20, 5, 2, 1)
      
      eqfthre <- round(eqfthre, 2)
    }
    
    real.spot.num <- c()
    for(k in 1:length(eqfthre)){
      real.spot.num[k] <- length(which(eqf.all > eqfthre[k]))
    }
    
    eqfthre <- cbind(eqfthre, real.spot.num)
  }
  
  bin <- cbind(1:nc, lcr)
  
  return(list(EQF.matrix = EQF, bin = bin, bin.size = bin.size, EQF.trait = eqf.t, EQF.detect = eqf.tq, 
              EQF.nondetect = eqf.tnq, cluster.matrix = eqf.c, permu.matrix.cluster = permu.clu, 
              permu.matrix.Q = permu.q, EQF.threshold = eqfthre))

}