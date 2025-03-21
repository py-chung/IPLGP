#' Generate D Matrix
#'
#' Generate the genetic design matrix of specified QTL number and effects.
#'
#' @param nQTL integer. The number of QTLs.
#' @param type character. The population type of the dataset. Includes
#' backcross (type="BC"), advanced intercross population (type="AI"), and
#' recombinant inbred population (type="RI"). The default value is "RI".
#' @param a integer or vector. A integer or vector to determines which
#' additive effects of QTLs will be considered in this design matrix. If
#' a=TRUE, the additive effect of all QTLs will be considered. If
#' a=FALSE, no additive effect will be considered.
#' @param d integer or vector. A integer or vector to determines which
#' dominant effects of QTLs will be considered in this design matrix. If
#' d=TRUE, the dominant effect of all QTLs will be considered.If
#' d=FALSE, no dominant effect will be considered.
#' @param aa vector or matrix. The additive-by-additive interaction. Two
#' format can be used in this parameter. One format is a vector, where
#' every two elements indicate a combination of additive-by-additive
#' interaction. The other format is a 2*i matrix, where i is the number
#' of combinations of interaction, and each column indicates the two
#' interacting QTLs. Additionally, if aa=TRUE, all combinations of
#' additive-by-additive interaction will be considered. If aa=FALSE, no
#' additive-by-additive interaction will be considered.
#' @param dd vector or matrix. The dominant-by-dominant interaction. The
#' format is the same as that in aa.
#' @param ad vector or matrix. The additive-by-dominant interaction. The
#' format is the same as that in aa. Note that in each pair of QTLs, the
#' first element indicates the additive effect, and the second element
#' indicates the dominant effect.
#'
#' @return
#' The genetic design matrix, where the elements represent the coded
#' variables of the QTL effects. It is a g*p matrix, where g is the number
#' of possible QTL genotypes, and p is the number of effects in the MIM
#' model.
#'
#' @note
#' For the 'type' parameter, if type="BC", the design matrix exclusively
#' contains additive effects and additive-by-additive interactions. However,
#' if type="AI" or type="RI", it encompasses additive and dominance effects
#' along with all interactions.
#'
#' For instance, when aa=c(1,3,2,4,5,6), it denotes that the interaction
#' between QTL1 and QTL3, the interaction between QTL2 and QTL4, and that
#' between QTL5 and QTL6 will be considered in the design matrix. Furthermore,
#' the matrix format can be expressed as aa=matrix(c(1,3,2,4,5,6),2,3).
#' Similarly, parameters DD and AD are also expressed in the same format.
#'
#' @export
#'
#' @references
#'
#' KAO, C.-H. and Z.-B. ZENG 1997 General formulas for obtaining the maximum
#' likelihood estimates and the asymptotic variance-covariance matrix in QTL
#' mapping when using the EM algorithm. Biometrics 53, 653-665. <doi: 10.2307/2533965.>
#'
#' KAO, C.-H., Z.-B. ZENG and R. D. TEASDALE 1999 Multiple interval mapping
#' for Quantitative Trait Loci. Genetics 152: 1203-1216. <doi: 10.1093/genetics/152.3.1203>
#'
#' @examples
#' D.make(4, d = c(1,3,4), aa = c(1,2,2,3), dd = c(1,3,1,4), ad = c(1,2,2,1,2,3,3,4))
#'
#' aa <- matrix(c(1,2,3,4,4,5), 2, 3)
#' aa
#' D.make(5, type = "BC", a = c(1,3,4,5), aa = aa)
D.make <- function(nQTL, type = "RI", a = TRUE, d = TRUE, aa = FALSE, dd = FALSE, ad = FALSE){

  if(!is.numeric(nQTL) | length(nQTL) != 1 | min(nQTL) < 0){
    stop("Parameter nQTL error, please input a positive integer.", call. = FALSE)
  }

  if(!type[1] %in% c("AI","RI","BC") | length(type) > 1){
    stop("Parameter type error, please input AI, RI, or BC.", call. = FALSE)
  }

  if(type == "BC"){
    d <- 0
    dd <- 0
    ad <- 0
  }

  if(!FALSE %in% (c(a, d, aa, dd, ad) == 0)){
    stop("Please set as least one effect in the design matrix.", call. = FALSE)
  }


  if(!0 %in% a){
    if("TRUE" %in% as.character(a)){a <- 1:nQTL}
    a <- c(a)
    datatry <- try(a%*%a, silent = TRUE)
    if(class(datatry)[1] == "try-error" | min(a) < 1 | max(a) > nQTL | length(a) > nQTL){
      stop("Effects setting error, please check the parameters of effects.", call. = FALSE)
    }
  } else {a <- NULL}

  if(!0 %in% d){
    if("TRUE" %in% as.character(d)){d <- 1:nQTL}
    d <- c(d)
    datatry <- try(d%*%d, silent = TRUE)
    if(class(datatry)[1] == "try-error" | min(d) < 1 | max(d) > nQTL | length(d) > nQTL){
      stop("Effects setting error, please check the parameters of effects.", call. = FALSE)
    }
  } else {d <- NULL}

  if(type == "BC"){
    D.matrix <- matrix(c(0.5, -0.5), 2, 1)
    row.names(D.matrix) <- c(2, 1)
    colnames(D.matrix) <- "a1"
    if(nQTL > 1){
      for(i in 2:nQTL){
        D.matrix <- rbind(cbind(0.5, D.matrix), cbind(-0.5, D.matrix))
      }
      colnames(D.matrix) <- paste("a", 1:nQTL, sep = "")
      D2 <- matrix(0, 2^nQTL, nQTL)
      for(i in 1:nQTL){
        D2[, i] <- rep(rep(c(2, 1), each = 2^(nQTL-i)), 2^(i-1))
      }
      row.names(D.matrix) <- apply(D2, 1, paste, collapse = "")
    }
  } else {
    D.matrix <- matrix(c(1, 0, -1, -0.5, 0.5, -0.5), 3, 2)
    row.names(D.matrix) <- c(2, 1, 0)
    colnames(D.matrix) <- c("a1", "d1")
    if(nQTL > 1){
      for(i in 2:nQTL){
        D.matrix <- rbind(cbind(1, -0.5, D.matrix), cbind(0, 0.5, D.matrix), cbind(-1, -0.5, D.matrix))
      }
      colnames(D.matrix) <- paste(rep(c("a", "d"), nQTL), rep(1:nQTL, each=2), sep = "")
      D2 <- matrix(0, 3^nQTL, nQTL)
      for(i in 1:nQTL){
        D2[, i] <- rep(rep(c(2, 1, 0), each = 3^(nQTL-i)), 3^(i-1))
      }
      row.names(D.matrix) <- apply(D2, 1, paste, collapse = "")
    }
  }

  if(nQTL > 1){
    if(!0 %in% aa){
      if(length(dim(aa)) > 0){
        if(dim(aa)[1] == 2){
          aa <- c(aa)
        } else {
          stop("Epistasis effects setting error, please check the parameters of epistasis effects.", call. = F)
        }
      }
      if("TRUE" %in% as.character(aa) | "T" %in% as.character(aa)){
        aa <- utils::combn(nQTL, 2)
      } else {
        datatry <- try(aa%*%aa, silent = TRUE)
        if(class(datatry)[1] == "try-error" | min(aa) < 1 | max(aa) > nQTL | length(aa) > choose(nQTL, 2)*2 | length(aa)%%2 == 1){
          stop("Epistasis effects setting error, please check the parameters of epistasis effects.", call. = F)
        }
        aa <- matrix(aa, 2, length(aa)/2)
      }
      for(j in 1:ncol(aa)){
        D1 <- matrix(D.matrix[, colnames(D.matrix) == paste("a", aa[1, j], sep = "")]*
                       D.matrix[, colnames(D.matrix) == paste("a", aa[2, j], sep = "")], nrow(D.matrix), 1)
        colnames(D1) <- paste("a", aa[1, j], ":", "a", aa[2, j], sep = "")
        D.matrix <- cbind(D.matrix, D1)
      }
    } else {aa <- NULL}
    if(type != "BC"){
      if(!0 %in% dd){
        if(length(dim(dd)) > 0){
          if(dim(dd)[1] == 2){
            dd <- c(dd)
          } else {
            stop("Epistasis effects setting error, please check the parameters of epistasis effects.", call. = F)
          }
        }
        if("TRUE" %in% as.character(dd) | "T" %in% as.character(dd)){
          dd <- utils::combn(nQTL, 2)
        } else {
          datatry <- try(dd%*%dd, silent = TRUE)
          if(class(datatry)[1] == "try-error" | min(dd) < 1 | max(dd) > nQTL | length(dd) > choose(nQTL, 2)*2 | length(dd)%%2 == 1){
            stop("Epistasis effects setting error, please check the parameters of epistasis effects.", call. = FALSE)
          }
          dd <- matrix(dd, 2, length(dd)/2)
        }
        for(j in 1:ncol(dd)){
          D1 <- matrix(D.matrix[, colnames(D.matrix) == paste("d", dd[1, j], sep = "")]*
                         D.matrix[, colnames(D.matrix) == paste("d", dd[2, j], sep = "")], nrow(D.matrix), 1)
          colnames(D1) <- paste("d", dd[1, j], ":", "d", dd[2, j], sep = "")
          D.matrix <- cbind(D.matrix, D1)
        }
      } else {dd <- NULL}
      if(!0 %in% ad){
        if(length(dim(ad)) > 0){
          if(dim(ad)[1] == 2){
            ad <- c(ad)
          } else {
            stop("Epistasis effects setting error, please check the parameters of epistasis effects.", call. = F)
          }
        }
        if("TRUE" %in% as.character(ad) | "T" %in% as.character(ad)){
          ad <- t(utils::combn(nQTL, 2))
          ad <- rbind(ad, cbind(ad[, 2], ad[, 1]))
          ad <- t(ad[order(ad[, 1], ad[, 2]), ])
        } else {
          datatry <- try(ad%*%ad, silent = TRUE)
          if(class(datatry)[1] == "try-error" | min(ad) < 1 | max(ad) > nQTL | length(ad) > choose(nQTL, 2)*4 | length(ad)%%2 == 1){
            stop("Epistasis effects setting error, please check the parameters of epistasis effects.", call. = FALSE)
          }
          ad <- matrix(ad, 2, length(ad)/2)
        }
        for(j in 1:ncol(ad)){
          D1 <- matrix(D.matrix[, colnames(D.matrix) == paste("a", ad[1, j], sep = "")]*
                         D.matrix[, colnames(D.matrix) == paste("d", ad[2, j], sep = "")], nrow(D.matrix), 1)
          colnames(D1) <- paste("a", ad[1, j], ":", "d", ad[2, j], sep = "")
          D.matrix <- cbind(D.matrix, D1)
        }
      } else {ad <- NULL}
    }
    a <- unique(a)
    a <- sort(a)
    a.cut <- (1:nQTL)[!(1:nQTL) %in% a]

    d <- unique(d)
    d <- sort(d)
    d.cut <- (1:nQTL)[!(1:nQTL) %in% d]
    if(type == "BC"){d.cut <- NULL}

    DD <- D.matrix

    if(length(a.cut) > 0 & length(d.cut) > 0){
      DD <- D.matrix[, -c((a.cut*2-1), (d.cut*2))]
    } else {
      if(length(a.cut) > 0){
        d.cut <- NULL
        if(type == "BC"){
          DD <- D.matrix[, -(a.cut)]
        } else {DD <- D.matrix[, -(a.cut*2-1)]}
      }
      if(length(d.cut) > 0){
        a.cut <- NULL
        DD <- D.matrix[, -(d.cut*2)]
      }
    }

    if(!is.matrix(DD)){
      cut0 <- c((a.cut*2-1), (d.cut*2))
      cut0 <- cut0[!is.na(cut0)]
      DD <- matrix(DD, nrow(D.matrix), 1)
      colnames(DD) <- colnames(D.matrix)[setdiff((1:ncol(D.matrix)), cut0)]
      row.names(DD) <- row.names(D.matrix)
    } else if (ncol(DD) == 0){
      DD <- 1
    }
    D.matrix <- DD

  } else {
    if(sum(d, na.rm = T) == 0){
      DD <- matrix(D.matrix[, 1], nrow(D.matrix), 1)
      colnames(DD) <- "a1"
      row.names(DD) <- row.names(D.matrix)
      D.matrix <- DD
    }
    if(sum(a,na.rm = T) == 0){
      DD <- matrix(D.matrix[, 2], nrow(D.matrix), 1)
      colnames(DD) <- "d1"
      row.names(DD) <- row.names(D.matrix)
      D.matrix <- DD
    }
  }

  return(D.matrix)
}
