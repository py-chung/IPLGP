#' LRT Threshold
#'
#' The LRT threshold for QTL interval mapping based on the
#' Gaussian stochastic process (Kao and Ho 2012).
#'
#' @param marker matrix. A k*2 matrix contains the marker information,
#' where the row dimension 'k' represents the number of markers in the
#' chromosomes. The first column labels the chromosomes where the markers
#' are located, and the second column labels the positions of markers (in
#' morgan (M) or centimorgan (cM)). It's important to note that chromosomes
#' and positions must be sorted in order.
#' @param type character. The population type of the dataset. Includes
#' backcross (type="BC"), advanced intercross population (type="AI"), and
#' recombinant inbred population (type="RI"). The default value is "RI".
#' @param ng integer. The generation number of the population type. For
#' instance, in a BC1 population where type="BC", ng=1; in an AI F3
#' population where type="AI", ng=3.
#' @param cM logical. Specify the unit of marker position. If cM=TRUE, it
#' denotes centimorgan; if cM=FALSE, it denotes morgan.
#' @param ns integer. The number of individuals for generating the
#' individual trait values. Changes in this value do not significantly
#' affect the outcome of the LRT threshold value.
#' @param gv numeric. The genetic variance for generating the
#' individual trait values. Changes in this value do not significantly
#' affect the outcome of the LRT threshold value.
#' @param speed numeric. The walking speed of the QTL analysis (in cM).
#' @param simu integer. Determines the number of simulation samples that
#' will be used to compute the LRT threshold using the Gaussian process.
#' @param d.eff logical. Specifies whether the dominant effect will be
#' considered in the parameter estimation for AI or RI population.
#' @param alpha numeric. The type I error rate for the LRT threshold.
#' @param console logical. Determines whether the process of the algorithm
#' will be displayed in the R console or not.
#'
#' @return
#'
#' The LRT threshold for QTL interval mapping.
#'
#' @export
#'
#' @references
#'
#' KAO, C.-H. and H.-A. Ho 2012 A score-statistic approach for determining
#' threshold values in QTL mapping. Frontiers in Bioscience. E4, 2670-2682. <doi: 10.2741/e582>
#'
#' @seealso
#' \code{\link[mvtnorm]{rmvnorm}}
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' LRTthre(marker, type = "RI", ng = 2, speed = 2, simu = 60)
LRTthre <- function(marker, type = "RI", ng = 2, cM = TRUE, ns = 200, gv = 25,
                    speed = 1, simu = 1000, d.eff = FALSE, alpha = 0.05, console = TRUE){

  if(is.null(marker)){
    stop("Input data is missing, please cheak and fix.", call. = FALSE)
  }

  markertest <- c(ncol(marker) != 2, NA %in% marker, marker[,1] != sort(marker[,1]))
  datatry <- try(marker*marker, silent=TRUE)
  if(class(datatry)[1] == "try-error" | TRUE %in% markertest){
    stop("Marker data error, or the number of marker does not match the genetype data.", call. = FALSE)
  }

  if(!cM[1] %in% c(0,1) | length(cM > 1)){cM <- TRUE}

  if(!type[1] %in% c("AI","RI","BC") | length(type) > 1){
    stop("Parameter type error, please input AI, RI, or BC.", call. = FALSE)
  }

  if(!is.numeric(ng) | length(ng) > 1 | min(ng) < 1){
    stop("Parameter ng error, please input a positive integer.", call. = FALSE)
  }
  ng <- round(ng)

  datatry <- try(ns%*%gv, silent=TRUE)
  if(class(datatry)[1] == "try-error" | length(ns) > 1 | length(gv) > 1){
    stop("Parameter ns or gv error, please input a positive integer.", call. = FALSE)
  }

  if(!d.eff[1] %in% c(0,1) | length(d.eff) > 1){d.eff <- TRUE}

  if(!is.numeric(simu) | length(simu) > 1 | min(simu) < 20 | max(simu) > 10^8){
    simu = 1000
  }
  cs <- 100
  if(simu < 200){cs <- 20}

  if(!is.numeric(speed) | length(speed) > 1 | min(speed) < 0){
    stop("Parameter speed error, please input a positive number.", call. = FALSE)
  }

  if(!is.numeric(alpha) | length(alpha) > 1 | min(alpha) < 0 | max(alpha) > 1){
    stop("Parameter alpha error, please input a positive number between 0 and 1.", call. = FALSE)
  }

  if(!console[1] %in% c(0,1) | length(console) > 1){console <- TRUE}

  if(cM){
    marker[, 2] <- marker[, 2]/100
  }

  Diff.P <- function (Freq){
    Diff22 <- (Freq[1]-Freq[3])/sum(Freq[1:3])
    Diff21 <- (Freq[4]-Freq[7])/sum(Freq[4:7])
    Diff20 <- (Freq[8]-Freq[10])/sum(Freq[8:10])
    Diff12 <- (Freq[11]-Freq[14])/sum(Freq[11:14])
    Diff <- c(Diff22, Diff21, Diff20, Diff12, 0, -Diff12, -Diff20, -Diff21, -Diff22)
    return(Diff)
  }

  Diff.k <- function (Freq){
    Diff22 <- (Freq[1]-Freq[3])/sum(Freq[1:3])
    Diff21 <- (Freq[4]-Freq[7])/sum(Freq[4:7])
    Diff20 <- (Freq[8]-Freq[10])/sum(Freq[8:10])
    Diff12 <- (Freq[11]-Freq[14])/sum(Freq[11:14])
    Diff <- c(Diff22, Diff21, Diff20, Diff12, 0, -Diff12, -Diff20, -Diff21, -Diff22)
    return(Diff)
  }

  Diff.l <- function (Freq){
    Diff22 <- (Freq[2]-Freq[1]-Freq[3])/sum(Freq[1:3])
    Diff21 <- (Freq[5]+Freq[6]-Freq[4]-Freq[7])/sum(Freq[4:7])
    Diff20 <- (Freq[9]-Freq[8]-Freq[10])/sum(Freq[8:10])
    Diff12 <- (Freq[12]+Freq[13]-Freq[11]-Freq[14])/sum(Freq[11:14])
    Diff11 <- (sum(Freq[17:20])-2*(Freq[15]+Freq[16]))/(sum(Freq[15:20])+Freq[15]+Freq[16])
    Diff <- c(Diff22, Diff21, Diff20, Diff12, Diff11, Diff12, Diff20, Diff21, Diff22)
    return(Diff)
  }


  cr <- unique(marker[, 1])
  ncr <- length(cr)
  thre <- c()
  if(type == "BC"){
    thre.fun <- function(d, nG, simu, nS, ss, console = FALSE, cr = 0, speed = 1){
      d <- d*100
      marker <- length(d)+1
      D <- sum(d)
      cov <- matrix(0, ncol = 2*marker, nrow = 2*marker)

      bcp3M <- function(d1, d2, nG){
        r1 <- (1-exp(-2*d1/100))/2
        r2 <- (1-exp(-2*d2/100))/2
        pABC <- pABc <- pAbC <- paBC <- pabC <- paBc <- pAbc <- pabc <- c()

        pABC[1] <- (1-r1)*(1-r2)/2
        pABc[1] <- (1-r1)*r2/2
        pAbC[1] <- r1*r2/2
        paBC[1] <- r1*(1-r2)/2
        pabC[1] <- (1-r1)*r2/2
        paBc[1] <- r1*r2/2
        pAbc[1] <- r1*(1-r2)/2
        pabc[1] <- (1-r1)*(1-r2)/2

        if(nG > 1){
          for(i in 2:nG){
            pABC[i] <- pABC[i-1]+1/2*pABc[i-1]+1/2*pAbC[i-1]+1/2*paBC[i-1]+(1-r1)/2*pabC[i-1]+
              ((1-r1)*(1-r2)+r1*r2)/2*paBc[i-1]+(1-r2)/2*pAbc[i-1]+(1-r1)*(1-r2)/2*pabc[i-1]
            pABc[i] <- 1/2*pABc[i-1]+(r1*(1-r2)+(1-r1)*r2)/2*paBc[i-1]+r2/2*pAbc[i-1]+(1-r1)*r2/2*pabc[i-1]
            pAbC[i] <- 1/2*pAbC[i-1]+r1/2*pabC[i-1]+r2/2*pAbc[i-1]+r1*r2/2*pabc[i-1]
            paBC[i] <- 1/2*paBC[i-1]+r1/2*pabC[i-1]+((1-r1)*r2+r1*(1-r2))/2*paBc[i-1]+r1*(1-r2)/2*pabc[i-1]
            pabC[i] <- (1-r1)/2*pabC[i-1]+(1-r1)*r2/2*pabc[i-1]
            paBc[i] <- ((1-r1)*(1-r2)+r1*r2)/2*paBc[i-1]+r1*r2/2*pabc[i-1]
            pAbc[i] <- (1-r2)/2*pAbc[i-1]+r1*(1-r2)/2*pabc[i-1]
            pabc[i] <- (1-r1)*(1-r2)/2*pabc[i-1]
          }
        }
        re <- c(pABC[nG], pABc[nG], pAbC[nG], paBC[nG], pabC[nG], paBc[nG], pAbc[nG], pabc[nG])
        return(re)
      }

      bcp2M <- function(d, nG){
        r <- (1-exp(-2*d/100))/2
        pAB <- pAb <- paB <- pab <- c()

        pAB[1] <- (1-r)/2
        pAb[1] <- r/2
        paB[1] <- r/2
        pab[1] <- (1-r)/2
        if(nG > 1){
          for(i in 2:nG){
            pAB[i] <- pAB[i-1]+1/2*pAb[i-1]+1/2*paB[i-1]+(1-r)/2*pab[i-1]
            pAb[i] <- 1/2*pAb[i-1]+r/2*pab[i-1]
            paB[i] <- 1/2*paB[i-1]+r/2*pab[i-1]
            pab[i] <- (1-r)/2*pab[i-1]
          }
        }

        re <- c(pAB[nG], pAb[nG], paB[nG], pab[nG])
        return(re)
      }

      b3.12 <- bcp3M(d[1], d[2], nG)
      b2.1 <- bcp2M(d[1], nG)
      b2.2 <- bcp2M(d[2], nG)
      if(marker > 2){
        cov[c(1:4), c(5, 6)] <- c(b3.12[1]/sqrt(b2.1[1]*b2.2[1]), 0,
                                  b3.12[4]/sqrt(b2.1[3]*b2.2[1]), 0,
                                  0, b3.12[3]/sqrt(b2.1[2]*b2.2[3]),
                                  0, b3.12[5]/sqrt(b2.1[4]*b2.2[3]))
      }

      if(marker > 3){
        for(k in 4:marker){
          bcp4M <- function(d1, d2, d3, nG){
            r1 <- (1-exp(-2*d1/100))/2
            r2 <- (1-exp(-2*d2/100))/2
            r3 <- (1-exp(-2*d3/100))/2
            pABCD <- pABCd <- pABcD <- pAbCD <- paBCD <- pABcd <- c()
            pAbCd <- paBCd <- pAbcD <- paBcD <- pabCD <- pabcD <- c()
            pabCd <- paBcd <- pAbcd <- pabcd <- c()

            pABCD[1] <- (1-r1)*(1-r2)*(1-r3)/2
            pABCd[1] <- (1-r1)*(1-r2)*r3/2
            pABcD[1] <- (1-r1)*r2*r3/2
            pAbCD[1] <- r1*r2*(1-r3)/2
            paBCD[1] <- r1*(1-r2)*(1-r3)/2
            pABcd[1] <- (1-r1)*r2*(1-r3)/2
            pAbCd[1] <- r1*r2*r3/2
            paBCd[1] <- r1*(1-r2)*r3/2
            pAbcD[1] <- r1*(1-r2)*r3/2
            paBcD[1] <- r1*r2*r3/2
            pabCD[1] <- (1-r1)*r2*(1-r3)/2
            pabcD[1] <- (1-r1)*(1-r2)*r3/2
            pabCd[1] <- (1-r1)*r2*r3/2
            paBcd[1] <- r1*r2*(1-r3)/2
            pAbcd[1] <- r1*(1-r2)*(1-r3)/2
            pabcd[1] <- (1-r1)*(1-r2)*(1-r3)/2

            if(nG > 1){
              for(i in 2:nG){
                pABCD[i]  <- pABCD[i-1]+1/2*pABCd[i-1]+1/2*pABcD[i-1]+1/2*pAbCD[i-1]+1/2*paBCD[i-1]+(1-r3)/2*pABcd[i-1]+
                  ((1-r2)*(1-r3)+r2*r3)/2*pAbCd[i-1]+
                  ((1-r1)*(1-r2)*(1-r3)+r1*r2*(1-r3)+(1-r1)*r2*r3+r1*(1-r2)*r3)/2*paBCd[i-1]+
                  (1-r2)/2*pAbcD[i-1]+((1-r1)*(1-r2)+r1*r2)/2*paBcD[i-1]+ (1-r1)/2*pabCD[i-1]+(1-r1)*(1-r2)/2*pabcD[i-1]+
                  ((1-r1)*(1-r2)*(1-r3)+(1-r1)*r2*r3)/2*pabCd[i-1]+ ((1-r1)*(1-r2)*(1-r3)+r1*r2*(1-r3))/2*paBcd[i-1]+
                  (1-r2)*(1-r3)/2*pAbcd[i-1]+(1-r1)*(1-r2)*(1-r3)/2*pabcd[i-1]

                pABCd[i] <- 1/2*pABCd[i-1]+((1-r2)*r3+r2*r3)/2*pABcd[i-1]+ ((1-r2)*r3+r2*(1-r3))/2*pAbCd[i-1]+
                  (r1*(1-r2)*(1-r3)+r1*r2*r3+(1-r1)*r2*(1-r3)+ (1-r1)*(1-r2)*r3)/2*paBCd[i-1]+
                  (1-r1)*((1-r2)*r3+r2*(1-r3))/2*pabCd[i-1]+((1-r1)*(1-r2)+r1*r2)*r3/2*paBcd[i-1]+
                  (1-r2)*r3/2*pAbcd[i-1]+ (1-r1)*(1-r2)*r3/2*pabcd[i-1]

                pABcD[i] <- 1/2*pABcD[i-1]+((1-r2)*r3+r2*r3)/2*pABcd[i-1]+ ((1-r1)*r2+r1*r2)/2*pAbcD[i-1]+
                  (r1*(1-r2)+(1-r1)*r2)/2*paBcD[i-1]+ (1-r1)*r2/2*pabcD[i-1]+
                  ((1-r1)*r2+r1*(1-r2))*r3/2*paBcd[i-1]+ r2*r3/2*pAbcd[i-1]+(1-r1)*r2*r3/2*pabcd[i-1]

                pAbCD[i] <- 1/2*pAbCD[i-1]+((1-r2)*r3+r2*(1-r3))/2*pAbCd[i-1]+r2/2*pAbcD[i-1]+r1/2*pabCD[i-1]+r1*r2/2*pabcD[i-1]+
                  r1*((1-r2)*r3+r2*(1-r3))/2*pabCd[i-1]+r2*(1-r3)/2*pAbcd[i-1]+ r1*r2*(1-r3)/2*pabcd[i-1]

                paBCD[i] <- 1/2*paBCD[i-1]+(r1*(1-r2)*(1-r3)+r1*r2*r3+(1-r1)*r2*(1-r3)+ (1-r1)*(1-r2)*r3)/2*paBCd[i-1]+
                  (r1*(1-r2)+(1-r1)*r2)/2*paBcD[i-1]+ r1/2*pabCD[i-1]+r1*(1-r2)/2*pabcD[i-1]+ r1*((1-r2)*(1-r3)+r2*r3)/2*pabCd[i-1]+
                  ((1-r1)*r2+r1*(1-r2))*(1-r3)/2*paBcd[i-1]+ r1*(1-r2)*(1-r3)/2*pabcd[i-1]

                pABcd[i] <- (1-r3)/2*pABcd[i-1]+((1-r1)*r2+r1*(1-r2))*(1-r3)/2*paBcd[i-1]+
                  r2*(1-r3)/2*pAbcd[i-1]+(1-r1)*r2*(1-r3)/2*pabcd[i-1]

                pAbCd[i] <- ((1-r2)*(1-r3)+r2*r3)/2*pAbCd[i-1]+ r1*((1-r2)*(1-r3)+r2*r3)/2*pabCd[i-1]+
                  r2*r3/2*pAbcd[i-1]+r1*r2*r3/2*pabcd[i-1]

                paBCd[i] <- ((1-r1)*(1-r2)*(1-r3)+r1*r2*(1-r3)+(1-r1)*r2*r3+ r1*(1-r2)*r3)/2*paBCd[i-1]+
                  r1*((1-r2)*r3+r2*(1-r3))/2*pabCd[i-1]+ ((1-r1)*r2+r1*(1-r2))*r3/2*paBcd[i-1]+r1*(1-r2)*r3/2*pabcd[i-1]

                pAbcD[i] <- (1-r2)/2*pAbcD[i-1]+r1*(1-r2)/2*pabcD[i-1]+(1-r2)*r3/2*pAbcd[i-1]+ r1*(1-r2)*r3/2*pabcd[i-1]

                paBcD[i] <- ((1-r1)*(1-r2)+r1*r2)/2*paBcD[i-1]+r1*r2/2*pabcD[i-1]+
                  ((1-r1)*(1-r2)+r1*r2)*r3/2*paBcd[i-1]+r1*r2*r3/2*pabcd[i-1]

                pabCD[i] <- (1-r1)/2*pabCD[i-1]+(1-r1)*r2/2*pabcD[i-1]+
                  (1-r1)*((1-r2)*r3+r2*(1-r3))/2*pabCd[i-1]+(1-r1)*r2*(1-r3)/2*pabcd[i-1]

                pabcD[i] <- (1-r1)*(1-r2)/2*pabcD[i-1]+(1-r1)*(1-r2)*r3/2*pabcd[i-1]
                pabCd[i] <- (1-r1)*((1-r2)*(1-r3)+r2*r3)/2*pabCd[i-1]+(1-r1)*r2*r3/2*pabcd[i-1]
                paBcd[i] <- ((1-r1)*(1-r2)+r1*r2)*(1-r3)/2*paBcd[i-1]+r1*r2*(1-r3)/2*pabcd[i-1]
                pAbcd[i] <- (1-r2)*(1-r3)/2*pAbcd[i-1]+r1*(1-r2)*(1-r3)/2*pabcd[i-1]
                pabcd[i] <- (1-r1)*(1-r2)*(1-r3)/2*pabcd[i-1]
              }
            }
            re <- c(pABCD[nG], pABCd[nG], pABcD[nG], pAbCD[nG],
                    paBCD[nG], pABcd[nG], pAbCd[nG], paBCd[nG],
                    pAbcD[nG], paBcD[nG], pabCD[nG], pabcD[nG],
                    pabCd[nG], paBcd[nG], pAbcd[nG], pabcd[nG])
            return(re)
          }


          b4.1k1 <- bcp4M(d[1], sum(d[2:(k-2)]), d[k-1], nG)
          b2.k1 <- bcp2M(d[k-1], nG)
          b2.k2 <- bcp2M(d[k-2], nG)
          b3.k2k1 <- bcp3M(d[k-2], d[k-1], nG)
          cov[1:4, c(2*k-1, 2*k)] <- c(b4.1k1[1]/sqrt(b2.1[1]*b2.k1[1]),
                                       b4.1k1[4]/sqrt(b2.1[2]*b2.k1[1]),
                                       b4.1k1[5]/sqrt(b2.1[3]*b2.k1[1]),
                                       b4.1k1[11]/sqrt(b2.1[4]*b2.k1[1]),
                                       b4.1k1[3]/sqrt(b2.1[1]*b2.k1[3]),
                                       b4.1k1[9]/sqrt(b2.1[2]*b2.k1[3]),
                                       b4.1k1[10]/sqrt(b2.1[3]*b2.k1[3]),
                                       b4.1k1[12]/sqrt(b2.1[4]*b2.k1[3]))
          cov[c(2*k-3, 2*k-2), 2*k-1] <- c(b3.k2k1[1]/sqrt(b2.k2[1]*b2.k1[1]),
                                           b3.k2k1[4]/sqrt(b2.k2[3]*b2.k1[1]))
        }
        for(k in 4:(marker-1)){
          for(l in 1:(marker-k)){
            b4.k2kl <- bcp4M(d[k-2], sum(d[(k-1):(k+l-2)]), d[k+l-1], nG)
            b2.k2 <- bcp2M(d[k-2], nG)
            b2.kl <- bcp2M(d[k+l-1], nG)
            cov[c(2*k-3, 2*k-2), c(2*k+2*l-1, 2*k+2*l)] <- c(b4.k2kl[1]/sqrt(b2.k2[1]*b2.kl[1]),
                                                             b4.k2kl[5]/sqrt(b2.k2[3]*b2.kl[1]),
                                                             b4.k2kl[3]/sqrt(b2.k2[1]*b2.kl[3]),
                                                             b4.k2kl[10]/sqrt(b2.k2[3]*b2.kl[3]))
          }
        }
      }
      cov <- cov+t(cov)
      diag(cov) <- 1

      a <- mvtnorm::rmvnorm(simu, rep(0, 2*marker), cov)
      y <- array(NA, dim = c(simu, 4, marker-1))
      for(i in 1:simu){
        y[i, 1, 1] <- a[i, 1]
        y[i, 2, 1] <- a[i, 2]
        y[i, 3, 1] <- a[i, 3]
        y[i, 4, 1] <- a[i, 4]
        if(marker > 2){
          for(j in 2:(marker-1)){
            y[i, 1, j] <- a[i, 2*j+1]
            y[i, 3, j] <- a[i, 2*(j+1)]
            y[i, 2, j] <- 1/sqrt(nS*bcp2M(d[j], nG)[2])*
              (sqrt(nS*bcp2M(d[j-1], nG)[1])*y[i, 1, j-1]+sqrt(nS*bcp2M(d[j-1], nG)[3])*
                 y[i, 3, j-1]-sqrt(nS*bcp2M(d[j], nG)[1])*y[i, 1, j])
            y[i, 4, j] <- 1/sqrt(nS*bcp2M(d[j], nG)[4])*
              (sqrt(nS*bcp2M(d[j-1], nG)[2])*y[i, 2, j-1]+sqrt(nS*bcp2M(d[j-1], nG)[4])*
                 y[i, 4, j-1]-sqrt(nS*bcp2M(d[j], nG)[3])*y[i, 3, j])
          }
        }
        if(i%%cs == 0){
          cat(paste(paste("BC", nG, "\t"), cr, i, "\n", sep = "\t")[console])
        }
      }
      cat(paste(paste("BC", nG, "\t"), cr, "done", "\n", sep = "\t")[console])

      Ds <- seq(0, D, speed)
      Dsc0 <- cumsum(d)
      Dsc1 <- c(0, Dsc0[-length(Dsc0)])
      Dsl <- length(Ds)
      loci0 <- c(Dsc0[1])
      loci1 <- c(0)
      k0 <- c(1)
      for(i in 2:Dsl){
        loci0[i] <- min(Dsc0[Dsc0 >= Ds[i]])
        loci1[i] <- max(Dsc1[Dsc1 < Ds[i]])
        k0[i] <- which(Dsc0 == min(Dsc0[Dsc0 >= Ds[i]]))
      }

      AA <- BB <- CC <- DD <- c()
      for(j in 1:(min(which(Ds >= Dsc0[1]))-1)){
        b3.j1 <- bcp3M(Ds[j], loci0[j]-Ds[j], nG)
        AA[j] <- nS*b2.1[2]*(-b3.j1[1]/b2.1[1]+b3.j1[3]/b2.1[1]+b3.j1[2]/b2.1[2]-b3.j1[7]/b2.1[2])+
          nS*b2.1[3]*(-b3.j1[1]/b2.1[1]+b3.j1[3]/b2.1[1]+b3.j1[4]/b2.1[3]-b3.j1[5]/b2.1[3])+
          nS*b2.1[4]*(-b3.j1[1]/b2.1[1]+b3.j1[3]/b2.1[1]+b3.j1[6]/b2.1[4]-b3.j1[8]/b2.1[4])
        BB[j] <- nS*b2.1[1]*( b3.j1[1]/b2.1[1]-b3.j1[3]/b2.1[1]-b3.j1[2]/b2.1[2]+b3.j1[7]/b2.1[2])+
          nS*b2.1[3]*(-b3.j1[2]/b2.1[2]+b3.j1[7]/b2.1[2]+b3.j1[4]/b2.1[3]-b3.j1[5]/b2.1[3])+
          nS*b2.1[4]*(-b3.j1[2]/b2.1[2]+b3.j1[7]/b2.1[2]+b3.j1[6]/b2.1[4]-b3.j1[8]/b2.1[4])
        CC[j] <- nS*b2.1[1]*( b3.j1[1]/b2.1[1]-b3.j1[3]/b2.1[1]-b3.j1[4]/b2.1[3]+b3.j1[5]/b2.1[3])+
          nS*b2.1[2]*( b3.j1[2]/b2.1[2]-b3.j1[7]/b2.1[2]-b3.j1[4]/b2.1[3]+b3.j1[5]/b2.1[3])+
          nS*b2.1[4]*(-b3.j1[4]/b2.1[3]+b3.j1[5]/b2.1[3]+b3.j1[6]/b2.1[4]-b3.j1[8]/b2.1[4])
        DD[j] <- nS*b2.1[1]*( b3.j1[1]/b2.1[1]-b3.j1[3]/b2.1[1]-b3.j1[6]/b2.1[4]+b3.j1[8]/b2.1[4])+
          nS*b2.1[2]*( b3.j1[2]/b2.1[2]-b3.j1[7]/b2.1[2]-b3.j1[6]/b2.1[4]+b3.j1[8]/b2.1[4])+
          nS*b2.1[3]*( b3.j1[4]/b2.1[3]-b3.j1[5]/b2.1[3]-b3.j1[6]/b2.1[4]+b3.j1[8]/b2.1[4])
      }
      if(marker > 2){
        for(j2 in (j+1):Dsl){
          b2.k <- bcp2M(loci0[j2],nG)
          loci <- Ds[length(AA)+1]
          b3.kkj <- bcp3M(loci-loci1[j2], loci0[j2]-loci, nG)
          AA[j2] <- nS*b2.k[2]*(-b3.kkj[1]/b2.k[1]+b3.kkj[3]/b2.k[1]+b3.kkj[2]/b2.k[2]-b3.kkj[7]/b2.k[2])+
            nS*b2.k[3]*(-b3.kkj[1]/b2.k[1]+b3.kkj[3]/b2.k[1]+b3.kkj[4]/b2.k[3]-b3.kkj[5]/b2.k[3])+
            nS*b2.k[4]*(-b3.kkj[1]/b2.k[1]+b3.kkj[3]/b2.k[1]+b3.kkj[6]/b2.k[4]-b3.kkj[8]/b2.k[4])
          BB[j2] <- nS*b2.k[1]*( b3.kkj[1]/b2.k[1]-b3.kkj[3]/b2.k[1]-b3.kkj[2]/b2.k[2]+b3.kkj[7]/b2.k[2])+
            nS*b2.k[3]*(-b3.kkj[2]/b2.k[2]+b3.kkj[7]/b2.k[2]+b3.kkj[4]/b2.k[3]-b3.kkj[5]/b2.k[3])+
            nS*b2.k[4]*(-b3.kkj[2]/b2.k[2]+b3.kkj[7]/b2.k[2]+b3.kkj[6]/b2.k[4]-b3.kkj[8]/b2.k[4])
          CC[j2] <- nS*b2.k[1]*( b3.kkj[1]/b2.k[1]-b3.kkj[3]/b2.k[1]-b3.kkj[4]/b2.k[3]+b3.kkj[5]/b2.k[3])+
            nS*b2.k[2]*( b3.kkj[2]/b2.k[2]-b3.kkj[7]/b2.k[2]-b3.kkj[4]/b2.k[3]+b3.kkj[5]/b2.k[3])+
            nS*b2.k[4]*(-b3.kkj[4]/b2.k[3]+b3.kkj[5]/b2.k[3]+b3.kkj[6]/b2.k[4]-b3.kkj[8]/b2.k[4])
          DD[j2] <- nS*b2.k[1]*( b3.kkj[1]/b2.k[1]-b3.kkj[3]/b2.k[1]-b3.kkj[6]/b2.k[4]+b3.kkj[8]/b2.k[4])+
            nS*b2.k[2]*( b3.kkj[2]/b2.k[2]-b3.kkj[7]/b2.k[2]-b3.kkj[6]/b2.k[4]+b3.kkj[8]/b2.k[4])+
            nS*b2.k[3]*( b3.kkj[4]/b2.k[3]-b3.kkj[5]/b2.k[3]-b3.kkj[6]/b2.k[4]+b3.kkj[8]/b2.k[4])
        }
      }

      b2M <- array(NA, dim = c(marker-1, 4))
      for(i in 1:(marker-1)){
        for(j in 1:4){
          b2M[i, j] <- bcp2M(d[i], nG)[j]
        }
      }

      umax <- c()
      u <- array(NA, dim = c(simu, Dsl))
      u2 <- array(NA, dim = c(simu, Dsl))
      for(i in 1:simu){
        for(j in 1:(min(which(Ds >= Dsc0[1]))-1)){
          u[i, j] <- (-sqrt(nS*b2M[1, 1])*AA[j]*y[i, 1, 1]-sqrt(nS*b2M[1, 2])*BB[j]*y[i, 2, 1]-
                        sqrt(nS*b2M[1, 3])*CC[j]*y[i, 3, 1]-sqrt(nS*b2M[1, 4])*DD[j]*y[i, 4, 1])/
            sqrt(nS*b2M[1, 1]*AA[j]^2+nS*b2M[1, 2]*BB[j]^2+nS*b2M[1, 3]*CC[j]^2+nS*b2M[1, 4]*DD[j]^2)
          u2[i, j] <- u[i, j]^2
        }
        if(marker > 2){
          for(j2 in (j+1):Dsl){
            k <- k0[j2]
            u[i, j2] <- (-sqrt(nS*b2M[k, 1])*AA[j2]*y[i, 1, k]-sqrt(nS*b2M[k, 2])*BB[j2]*y[i, 2, k]-
                           sqrt(nS*b2M[k, 3])*CC[j2]*y[i, 3, k]-sqrt(nS*b2M[k, 4])*DD[j2]*y[i, 4, k])/
              sqrt(nS*b2M[k, 1]*AA[j2]^2+nS*b2M[k, 2]*BB[j2]^2+nS*b2M[k, 3]*CC[j2]^2+nS*b2M[k, 4]*DD[j2]^2)
            u2[i, j2] <- u[i, j2]^2
          }
        }
        umax[i] <- max(u2[i, ], na.rm = TRUE)
      }
      return(umax)
    }
  } else if (type == "RI"){
    RI2 <- function (d, n.Generation){
      r <- (1-exp(-2*d))/2
      TransMatrix <- matrix(1:25, nrow = 5)
      TransMatrix[1,] <- c(4, 0, 2, (1-r)^2, r^2)/4
      TransMatrix[2,] <- c(0, 4, 2, r^2, (1-r)^2)/4
      TransMatrix[3,] <- c(0, 0, 1, r*(1-r), r*(1-r))/2
      TransMatrix[4,] <- c(0, 0, 0, (1-r)^2, r^2)/2
      TransMatrix[5,] <- c(0, 0, 0, r^2, (1-r)^2)/2
      Freq <- c(0, 0, 0, 1, 0)
      n.Iteration <- n.Generation-1
      ZygoteFreq <- matrix(rep(0, (9*n.Generation-9)), nrow = 9)
      for (i in 1:n.Iteration) {
        Freq <- TransMatrix%*%Freq
        ZygoteFreq[, i]<-c(Freq[c(1, 3, 2, 3)], sum(Freq[4:5]), Freq[c(3, 2, 3, 1)])
      }
      ZygoteFreq[, n.Iteration]
    }

    RI3 <- function (d1, d2, n.Generation){
      r1 <- (1-exp(-2*d1))/2
      r2 <- (1-exp(-2*d2))/2
      r <- (1-exp(-2*(d1+d2)))/2
      f <- r1*(1-r2)+r2*(1-r1)
      TransMatrix <- matrix(rep(0, 20^2), nrow = 20)
      TransMatrix[1, ] <- c(4, 1, 0, 1, (1-r2)^2, r2^2, 0, 0, 0, 0, 1, (1-r1)^2, r1^2,
                            0, (1-f)^2, f^2, ((1-r1)*(1-r2))^2, (r1*r2)^2, ((1-r1)*r2)^2,
                            (r1*(1-r2))^2)/4
      TransMatrix[2, ] <- c(0, 1, 0, 0, r2*(1-r2), r2*(1-r2), 0, 0, 0, 0, 0, r1*(1-r1),
                            r1*(1-r1), 0, 0, 0, rep(r1*(1-r1)*r2*(1-r2), 4))/2
      TransMatrix[3, ] <- c(0, 1, 4, 0, r2^2, (1-r2)^2, 1, 0, 0, 0, 0, r1^2, (1-r1)^2,
                            1, (1-f)^2, f^2, (r1*r2)^2, ((1-r1)*(1-r2))^2, (r1*(1-r2))^2,
                            ((1-r1)*r2)^2)/4
      TransMatrix[4, ] <- c(0, 0, 0, 1, r2*(1-r2), r2*(1-r2), rep(0, 8), f*(1-f), f*(1-f),
                            (1-r1)^2*r2*(1-r2), r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2),
                            r1^2*r2*(1-r2))/2
      TransMatrix[5, ] <- c(0, 0, 0, 0, (1-r2)^2, r2^2, rep(0, 10), r1*(1-r1)*(1-r2)^2,
                            r1*(1-r1)*r2^2, r1*(1-r1)*r2^2, r1*(1-r1)*(1-r2)^2)/2
      TransMatrix[6, ] <- c(0, 0, 0, 0, r2^2, (1-r2)^2, rep(0, 10), r1*(1-r1)*r2^2,
                            r1*(1-r1)*(1-r2)^2, r1*(1-r1)*(1-r2)^2, r1*(1-r1)*r2^2)/2
      TransMatrix[7, ] <- c(0, 0, 0, 0, r2*(1-r2), r2*(1-r2), 1, rep(0, 7), f*(1-f), f*(1-f),
                            r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2), r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2))/2
      TransMatrix[8, ] <- c(0, 0, 0, 1, r2^2, (1-r2)^2, 0, 4, 1, 0, 0, (1-r1)^2, r1^2, 1, f^2,
                            (1-f)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2, ((1-r1)*(1-r2))^2, (r1*r2)^2)/4
      TransMatrix[9, ] <- c(0, 0, 0, 0, r2*(1-r2), r2*(1-r2), 0, 0, 1, 0, 0, r1*(1-r1), r1*(1-r1), 0,
                            0, 0, rep(r1*(1-r1)*r2*(1-r2), 4))/2
      TransMatrix[10, ] <- c(0, 0, 0, 0, (1-r2)^2, r2^2, 1, 0, 1, 4, 1, r1^2, (1-r1)^2, 0, f^2,
                             (1-f)^2, (r1*(1-r2))^2, ((1-r1)*r2)^2, (r1*r2)^2, ((1-r1)*(1-r2))^2)/4
      TransMatrix[11, ] <- c(rep(0, 10), 1, r1*(1-r1), r1*(1-r1), 0, f*(1-f), f*(1-f),
                             r1*(1-r1)*(1-r2)^2, r1*(1-r1)*r2^2, r1*(1-r1)*r2^2, r1*(1-r1)*(1-r2)^2)/2
      TransMatrix[12, ] <- c(rep(0, 11), (1-r1)^2, r1^2, 0, 0, 0, (1-r1)^2*r2*(1-r2), r1^2*r2*(1-r2),
                             (1-r1)^2*r2*(1-r2), r1^2*r2*(1-r2))/2
      TransMatrix[13, ] <- c(rep(0, 11), r1^2, (1-r1)^2, 0, 0, 0, r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2),
                             r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2))/2
      TransMatrix[14, ] <- c(rep(0, 11), r1*(1-r1), r1*(1-r1), 1, f*(1-f), f*(1-f), r1*(1-r1)*r2^2,
                             r1*(1-r1)*(1-r2)^2, r1*(1-r1)*(1-r2)^2, r1*(1-r1)*r2^2)/2
      TransMatrix[15, ] <- c(rep(0, 14), (1-f)^2, f^2, rep(r1*(1-r1)*r2*(1-r2), 4))/2
      TransMatrix[16, ] <- c(rep(0, 14), f^2, (1-f)^2, rep(r1*(1-r1)*r2*(1-r2), 4))/2
      TransMatrix[17, ] <- c(rep(0, 16), ((1-r1)*(1-r2))^2, (r1*r2)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2)/2
      TransMatrix[18, ] <- c(rep(0, 16), (r1*r2)^2, ((1-r1)*(1-r2))^2, (r1*(1-r2))^2, ((1-r1)*r2)^2)/2
      TransMatrix[19, ] <- c(rep(0, 16), ((1-r1)*r2)^2, (r1*(1-r2))^2, ((1-r1)*(1-r2))^2, (r1*r2)^2)/2
      TransMatrix[20, ] <- c(rep(0, 16), (r1*(1-r2))^2, ((1-r1)*r2)^2, (r1*r2)^2, ((1-r1)*(1-r2))^2)/2
      Freq <- c(rep(0, 16), 1, 0, 0, 0)
      n.Iteration <- n.Generation-1
      ZygoteFreq <- matrix(rep(0, (20*n.Generation-20)), nrow = 20)
      for (i in 1:n.Iteration) {
        Freq <- TransMatrix%*%Freq
        ZygoteFreq[, i] <- Freq
      }
      ZygoteFreq[, n.Iteration]
    }

    RI4TM <- function (r1, r2, r3){
      f1 <- r1*(1-r2)+r2*(1-r1)
      f2 <- r2*(1-r3)+r3*(1-r2)
      tm <- matrix(rep(0, 72^2), nrow = 72)
      tm[1,] <- c(4, 1, 0, 1, (1-f2)^2, f2^2, rep(0, 4), 1, (1-r3)^2, r3^2, 0, (1-r2)^2, r2^2,
                  ((1-r2)*(1-r3))^2, (r2*(1-r3))^2, (r2*r3)^2, ((1-r2)*r3)^2, rep(0, 16), 1,
                  ((1-r1)*(1-f2)+r1*f2)^2, ((1-r1)*f2+r1*(1-f2))^2, 0, (1-r1)^2, r1^2,
                  ((1-r1)*(1-f1))^2, (r1*f1)^2, (r1*(1-f1))^2, ((1-r1)*f1)^2, rep(0, 6),
                  (1-f1)^2, f1^2, ((1-f1)*(1-r3))^2, (f1*(1-r3))^2, ((1-f1)*r3)^2, (f1*r3)^2,
                  0, 0, ((1-r1)*(1-r2))^2, (r1*r2)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2,
                  ((1-r1)*(1-r2)*(1-r3))^2, (r1*(1-r2)*(1-r3))^2, ((1-r1)*r2*(1-r3))^2,
                  (r1*r2*r3)^2, ((1-r1)*(1-r2)*r3)^2, (r1*r2*(1-r3))^2, (r1*(1-r2)*r3)^2,
                  ((1-r1)*r2*r3)^2)/4
      tm[2,] <- c(0, 1, 0, 0, (1-f2)*f2, (1-f2)*f2, rep(0, 5), (1-r3)*r3, (1-r3)*r3, 0, 0, 0,
                  (1-r2)^2*(1-r3)*r3, r2^2*(1-r3)*r3, r2^2*(1-r3)*r3, (1-r2)^2*(1-r3)*r3,
                  rep(0, 17), ((1-r1)*(1-f2)+r1*f2)*((1-r1)*f2+r1*(1-f2)),
                  ((1-r1)*(1-f2)+r1*f2)*((1-r1)*f2+r1*(1-f2)), 0, 0, 0, (1-r1)^2*f1*(1-f1),
                  r1^2*f1*(1-f1), r1^2*f1*(1-f1), (1-r1)^2*f1*(1-f1), rep(0, 8),
                  (1-f1)^2*r3*(1-r3), f1^2*r3*(1-r3), (1-f1)^2*r3*(1-r3), f1^2*r3*(1-r3),
                  rep(0, 6), (1-r1)^2*(1-r2)^2*r3*(1-r3), r1^2*(1-r2)^2*r3*(1-r3),
                  (1-r1)^2*r2^2*r3*(1-r3), r1^2*r2^2*r3*(1-r3), (1-r1)^2*(1-r2)^2*r3*(1-r3),
                  r1^2*r2^2*r3*(1-r3), r1^2*(1-r2)^2*r3*(1-r3), (1-r1)^2*r2^2*r3*(1-r3))/2
      tm[3,] <- c(0, 1, 4, 0, f2^2, (1-f2)^2, 1, 0, 0, 0, 0, r3^2, (1-r3)^2, 1, 0, 0,
                  ((1-r2)*r3)^2, (r2*r3)^2, (r2*(1-r3))^2, ((1-r2)*(1-r3))^2, (1-r2)^2, r2^2,
                  rep(0, 15), ((1-r1)*f2+r1*(1-f2))^2, ((1-r1)*(1-f2)+r1*f2)^2, 1, 0, 0,
                  ((1-r1)*f1)^2, (r1*(1-f1))^2, (r1*f1)^2, ((1-r1)*(1-f1))^2, (1-r1)^2, r1^2,
                  rep(0, 6), (r3*(1-f1))^2, (r3*f1)^2, ((1-r3)*(1-f1))^2, ((1-r3)*f1)^2,
                  (1-f1)^2, f1^2, ((1-r1)*(1-r2))^2, (r1*r2)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2,
                  ((1-r1)*(1-r2)*r3)^2, (r1*(1-r2)*r3)^2, ((1-r1)*r2*r3)^2, (r1*r2*(1-r3))^2,
                  ((1-r1)*(1-r2)*(1-r3))^2, (r1*r2*r3)^2, (r1*(1-r2)*(1-r3))^2,
                  ((1-r1)*r2*(1-r3))^2)/4
      tm[4,] <- c(0, 0, 0, 1, rep((1-f2)*f2, 2), rep(0, 8), r2*(1-r2), r2*(1-r2),
                  r2*(1-r2)*c((1-r3)^2, (1-r3)^2, r3^2, r3^2), rep(0, 20), r1*(1-r1), r1*(1-r1),
                  rep(r1*(1-r1)*f1*(1-f1), 4), rep(0, 14), rep(r1*(1-r1)*r2*(1-r2), 4),
                  r1*(1-r1)*r2*(1-r2)*c((1-r3)^2, (1-r3)^2, (1-r3)^2, r3^2, r3^2, (1-r3)^2, r3^2, r3^2))/2
      tm[5,] <- c(0, 0, 0, 0, (1-f2)^2, f2^2, rep(0, 10), rep(r2*(1-r2)*r3*(1-r3), 4), rep(0, 22),
                  r1*(1-r1)*(1-f1)^2, r1*(1-r1)*f1^2, r1*(1-r1)*(1-f1)^2, r1*(1-r1)*f1^2, rep(0, 18),
                  rep(r1*(1-r1)*r2*(1-r2)*r3*(1-r3), 8))/2
      tm[6,] <- c(0, 0, 0, 0, f2^2, (1-f2)^2, rep(0, 10), rep(r2*(1-r2)*r3*(1-r3), 4), rep(0, 22),
                  r1*(1-r1)*f1^2, r1*(1-r1)*(1-f1)^2, r1*(1-r1)*f1^2, r1*(1-r1)*(1-f1)^2,rep(0, 18),
                  rep(r1*(1-r1)*r2*(1-r2)*r3*(1-r3), 8))/2
      tm[7,] <- c(0, 0, 0, 0, f2*(1-f2), f2*(1-f2), 1, rep(0, 9), r2*(1-r2)*r3^2, r2*(1-r2)*r3^2,
                  r2*(1-r2)*(1-r3)^2, r2*(1-r2)*(1-r3)^2, r2*(1-r2), r2*(1-r2), rep(0, 20),
                  rep(r1*(1-r1)*f1*(1-f1), 4), r1*(1-r1), r1*(1-r1), rep(0, 12), rep(r1*(1-r1)*r2*(1-r2),4),
                  r1*(1-r1)*r2*(1-r2)*c(r3^2, r3^2,r3^2, (1-r3)^2, (1-r3)^2, r3^2, (1-r3)^2, (1-r3)^2))/2
      tm[8,] <- c(0, 0, 0, 1, f2^2, (1-f2)^2, 0, 4, 1, rep(0, 5), r2^2, (1-r2)^2, (r2*(1-r3))^2,
                  ((1-r2)*(1-r3))^2, ((1-r2)*r3)^2, (r2*r3)^2, 0, 0, 1, (1-r3)^2, r3^2, rep(0, 15), r1^2,
                  (1-r1)^2, (r1*f1)^2, ((1-r1)*(1-f1))^2, ((1-r1)*f1)^2, (r1*(1-f1))^2, 0, 0, 1,
                  ((1-r1)*(1-f1)+r1*f1)^2, ((1-r1)*f1+r1*(1-f1))^2, 0, 0, 0, ((1-r3)*(1-f1))^2,
                  ((1-r3)*f1)^2, (r3*(1-f1))^2, (r3*f1)^2, (1-f1)^2, f1^2, (r1*r2)^2, ((1-r1)*(1-r2))^2,
                  (r1*(1-r2))^2, ((1-r1)*r2)^2, (r1*r2*(1-r3))^2, ((1-r1)*r2*(1-r3))^2,
                  (r1*(1-r2)*(1-r3))^2, ((1-r1)*(1-r2)*r3)^2, (r1*r2*r3)^2, ((1-r1)*(1-r2)*(1-r3))^2,
                  ((1-r1)*r2*r3)^2, (r1*(1-r2)*r3)^2)/4
      tm[9,] <- c(0, 0, 0, 0, f2*(1-f2), f2*(1-f2), 0, 0, 1, rep(0, 7), r3*(1-r3)*r2^2, r3*(1-r3)*(1-r2)^2,
                  r3*(1-r3)*(1-r2)^2, r3*(1-r3)*r2^2, 0, 0, 0, r3*(1-r3), r3*(1-r3), rep(0, 17),
                  f1*(1-f1)*r1^2, f1*(1-f1)*(1-r1)^2, f1*(1-f1)*(1-r1)^2, f1*(1-f1)*r1^2, 0, 0, 0,
                  rep(((1-r1)*(1-f1)+r1*f1)*((1-r1)*f1+r1*(1-f1)), 2), 0, 0, 0, r3*(1-r3)*(1-f1)^2,
                  r3*(1-r3)*f1^2, r3*(1-r3)*(1-f1)^2, r3*(1-r3)*f1^2, rep(0, 6),
                  r3*(1-r3)*c(r1^2*r2^2, (1-r1)^2*r2^2, r1^2*(1-r2)^2, (1-r1)^2*(1-r2)^2, r1^2*r2^2,
                              (1-r1)^2*(1-r2)^2, (1-r1)^2*r2^2, r1^2*(1-r2)^2))/2
      tm[10,] <- c(0, 0, 0, 0, (1-f2)^2, f2^2, 1, 0, 1, 4, rep(0, 6), r2^2*r3^2, (1-r2)^2*r3^2,
                   (1-r2)^2*(1-r3)^2, r2^2*(1-r3)^2, r2^2, (1-r2)^2, 0, r3^2, (1-r3)^2, 1, rep(0, 16),
                   (r1*(1-f1))^2, ((1-r1)*f1)^2, ((1-r1)*(1-f1))^2, (r1*f1)^2, r1^2, (1-r1)^2, 0,
                   ((1-r1)*f1+r1*(1-f1))^2, ((1-r1)*(1-f1)+r1*f1)^2, 1, (1-f1)^2, f1^2, (r3*(1-f1))^2,
                   (r3*f1)^2, ((1-r3)*(1-f1))^2, ((1-r3)*f1)^2, 0, 0, (r1*r2)^2, ((1-r1)*(1-r2))^2,
                   (r1*(1-r2))^2, ((1-r1)*r2)^2, (r1*r2*r3)^2, ((1-r1)*r2*r3)^2, (r1*(1-r2)*r3)^2,
                   ((1-r1)*(1-r2)*(1-r3))^2, (r1*r2*(1-r3))^2, ((1-r1)*(1-r2)*r3)^2, ((1-r1)*r2*(1-r3))^2,
                   (r1*(1-r2)*(1-r3))^2)/4
      tm[11,] <- c(rep(0, 10), 1, r3*(1-r3), r3*(1-r3), 0, r2*(1-r2), r2*(1-r2), rep(r2*(1-r2)*r3*(1-r3), 4),
                   rep(0, 32), (1-f1)*f1, (1-f1)*f1, rep(r3*(1-r3)*(1-f1)*f1, 4), 0, 0,
                   r2*(1-r2)*c((1-r1)^2, r1^2, (1-r1)^2, r1^2),
                   r2*r3*(1-r2)*(1-r3)*c((1-r1)^2, r1^2, (1-r1)^2, r1^2, (1-r1)^2, r1^2, r1^2, (1-r1)^2))/2
      tm[12,] <- c(rep(0, 11), (1-r3)^2, r3^2, 0, 0, 0, r2*(1-r2)*c((1-r3)^2, (1-r3)^2, r3^2, r3^2),
                   rep(0, 34), f1*(1-f1)*c((1-r3)^2, (1-r3)^2, r3^2, r3^2), rep(0, 6),
                   r2*(1-r2)*c(((1-r1)*(1-r3))^2, (r1*(1-r3))^2, ((1-r1)*(1-r3))^2, (r1*r3)^2, ((1-r1)*r3)^2,
                               (r1*(1-r3))^2, (r1*r3)^2, ((1-r1)*r3)^2))/2
      tm[13,] <- c(rep(0, 11), r3^2, (1-r3)^2, 0, 0, 0, r2*(1-r2)*c(r3^2, r3^2, (1-r3)^2, (1-r3)^2),
                   rep(0, 34), f1*(1-f1)*c(r3^2, r3^2, (1-r3)^2, (1-r3)^2), rep(0, 6),
                   r2*(1-r2)*c(((1-r1)*r3)^2,(r1*r3)^2, ((1-r1)*r3)^2, (r1*(1-r3))^2, ((1-r1)*(1-r3))^2,
                               (r1*r3)^2, (r1*(1-r3))^2, ((1-r1)*(1-r3))^2))/2
      tm[14,] <- c(rep(0, 11), r3*(1-r3), r3*(1-r3), 1, 0, 0, rep(r2*(1-r2)*r3*(1-r3), 4), r2*(1-r2), r2*(1-r2),
                   rep(0, 32), rep(r3*(1-r3)*f1*(1-f1), 4), f1*(1-f1), f1*(1-f1),
                   r2*(1-r2)*c((1-r1)^2, r1^2, (1-r1)^2, r1^2),
                   r2*(1-r2)*r3*(1-r3)*c((1-r1)^2, r1^2, (1-r1)^2, r1^2, (1-r1)^2, r1^2, r1^2, (1-r1)^2))/2
      tm[15,] <- c(rep(0, 14), (1-r2)^2, r2^2, r3*(1-r3)*c((1-r2)^2, r2^2, r2^2, (1-r2)^2), rep(0, 40),
                   r1*(1-r1)*c((1-r2)^2, r2^2, r2^2, (1-r2)^2),
                   r1*(1-r1)*r3*(1-r3)*c((1-r2)^2, (1-r2)^2, r2^2, r2^2, (1-r2)^2, r2^2, (1-r2)^2, r2^2))/2
      tm[16,] <- c(rep(0, 14), r2^2, (1-r2)^2, r3*(1-r3)*c(r2^2, (1-r2)^2, (1-r2)^2, r2^2), rep(0, 40),
                   r1*(1-r1)*c(r2^2, (1-r2)^2, (1-r2)^2, r2^2),
                   r1*(1-r1)*r3*(1-r3)*c(r2^2, r2^2, (1-r2)^2, (1-r2)^2, r2^2, (1-r2)^2, r2^2, (1-r2)^2))/2
      tm[17,] <- c(rep(0, 16), ((1-r2)*(1-r3))^2, (r2*(1-r3))^2, (r2*r3)^2, ((1-r2)*r3)^2, rep(0, 44),
                   r1*(1-r1)*c(((1-r2)*(1-r3))^2, ((1-r2)*(1-r3))^2, (r2*(1-r3))^2, (r2*r3)^2,
                               ((1-r2)*r3)^2, (r2*(1-r3))^2, ((1-r2)*r3)^2, (r2*r3)^2))/2
      tm[18,] <- c(rep(0, 16), (r2*(1-r3))^2, ((1-r2)*(1-r3))^2, ((1-r2)*r3)^2, (r2*r3)^2, rep(0, 44),
                   r1*(1-r1)*c((r2*(1-r3))^2, (r2*(1-r3))^2, ((1-r2)*(1-r3))^2, ((1-r2)*r3)^2, (r2*r3)^2,
                               ((1-r2)*(1-r3))^2, (r2*r3)^2, ((1-r2)*r3)^2))/2
      tm[19,] <- c(rep(0, 16), (r2*r3)^2, ((1-r2)*r3)^2, ((1-r2)*(1-r3))^2, (r2*(1-r3))^2, rep(0, 44),
                   r1*(1-r1)*c((r2*r3)^2, (r2*r3)^2, ((1-r2)*r3)^2, ((1-r2)*(1-r3))^2, (r2*(1-r3))^2,
                               ((1-r2)*r3)^2, (r2*(1-r3))^2, ((1-r2)*(1-r3))^2))/2
      tm[20,] <- c(rep(0, 16), ((1-r2)*r3)^2, (r2*r3)^2, (r2*(1-r3))^2, ((1-r2)*(1-r3))^2, rep(0, 44),
                   r1*(1-r1)*c(((1-r2)*r3)^2, ((1-r2)*r3)^2, (r2*r3)^2, (r2*(1-r3))^2, ((1-r2)*(1-r3))^2,
                               (r2*r3)^2, ((1-r2)*(1-r3))^2, (r2*(1-r3))^2))/2
      tm[21,] <- c(rep(0, 16), r3*(1-r3)*c((1-r2)^2, r2^2, r2^2, (1-r2)^2), (1-r2)^2, r2^2, rep(0, 38),
                   r1*(1-r1)*c((1-r2)^2, r2^2, r2^2, (1-r2)^2),
                   r1*(1-r1)*r3*(1-r3)*c((1-r2)^2, (1-r2)^2, r2^2, r2^2, (1-r2)^2, r2^2, (1-r2)^2, r2^2))/2
      tm[22,] <- c(rep(0, 16), r3*(1-r3)*c(r2^2, (1-r2)^2, (1-r2)^2, r2^2), r2^2, (1-r2)^2, rep(0, 38),
                   r1*(1-r1)*c(r2^2, (1-r2)^2, (1-r2)^2, r2^2),
                   r1*(1-r1)*r3*(1-r3)*c(r2^2, r2^2, (1-r2)^2, (1-r2)^2, r2^2, (1-r2)^2, r2^2, (1-r2)^2))/2
      tm[23,] <- c(rep(0, 14), r2*(1-r2), r2*(1-r2), rep(r2*(1-r2)*r3*(1-r3), 4), 0, 0, 1, r3*(1-r3),
                   r3*(1-r3), rep(0, 29), rep(r3*(1-r3)*f1*(1-f1), 4), f1*(1-f1), f1*(1-f1),
                   r2*(1-r2)*c(r1^2, (1-r1)^2, r1^2, (1-r1)^2),
                   r2*(1-r2)*r3*(1-r3)*c(r1^2, (1-r1)^2, r1^2, (1-r1)^2, r1^2, (1-r1)^2, (1-r1)^2, r1^2))/2
      tm[24,] <- c(rep(0, 16), r2*(1-r2)*c((1-r3)^2, (1-r3)^2, r3^2, r3^2), 0, 0, 0, (1-r3)^2, r3^2,
                   rep(0, 29), f1*(1-f1)*c((1-r3)^2, (1-r3)^2, r3^2, r3^2), rep(0, 6),
                   r2*(1-r2)*c((r1*(1-r3))^2, ((1-r1)*(1-r3))^2, (r1*(1-r3))^2, ((1-r1)*r3)^2,
                               (r1*r3)^2, ((1-r1)*(1-r3))^2, ((1-r1)*r3)^2, (r1*r3)^2))/2
      tm[25,] <- c(rep(0, 16), r2*(1-r2)*c(r3^2, r3^2, (1-r3)^2, (1-r3)^2), 0, 0, 0, r3^2, (1-r3)^2, rep(0, 29),
                   f1*(1-f1)*c(r3^2, r3^2, (1-r3)^2, (1-r3)^2), rep(0, 6),
                   r2*(1-r2)*c((r1*r3)^2, ((1-r1)*r3)^2, (r1*r3)^2, ((1-r1)*(1-r3))^2, (r1*(1-r3))^2,
                               ((1-r1)*r3)^2, ((1-r1)*(1-r3))^2, (r1*(1-r3))^2))/2
      tm[26,] <- c(rep(0, 16), rep(r2*(1-r2)*r3*(1-r3), 4), r2*(1-r2), r2*(1-r2), 0, r3*(1-r3), r3*(1-r3), 1,
                   rep(0, 26), f1*(1-f1), f1*(1-f1), rep(r3*(1-r3)*f1*(1-f1), 4), 0, 0,
                   r2*(1-r2)*c(r1^2, (1-r1)^2, r1^2, (1-r1)^2),
                   r2*(1-r2)*r3*(1-r3)*c(r1^2, (1-r1)^2, r1^2, (1-r1)^2, r1^2, (1-r1)^2, (1-r1)^2, r1^2))/2
      tm[27,] <- c(rep(0, 10), 1, r3^2, (1-r3)^2, 0, r2^2, (1-r2)^2, (r2*r3)^2, ((1-r2)*r3)^2,
                   ((1-r2)* (1-r3))^2, (r2*(1-r3))^2, rep(0, 6), 4, 1, 0, 1, (1-f2)^2, f2^2, rep(0, 10),
                   ((1-r1)*(1-f2))^2, (r1*f2)^2, (r1*(1-f2))^2, ((1-r1)*f2)^2, (1-r1)^2, r1^2, 0,
                   ((1-r1)*(1-f2)+r1*f2)^2, (r1*(1-f2)+(1-r1)*f2)^2, 1, f1^2, (1-f1)^2, (r3*f1)^2, (r3*(1-f1))^2,
                   ((1-r3)*f1)^2, ((1-r3)*(1-f1))^2, 0, 0, ((1-r1)*r2)^2, (r1*(1-r2))^2, ((1-r1)*(1-r2))^2,
                   (r1*r2)^2, ((1-r1)*r2*r3)^2, (r1*r2*r3)^2, ((1-r1)*(1-r2)*r3)^2, (r1*(1-r2)*(1-r3))^2,
                   ((1-r1)*r2*(1-r3))^2, (r1*(1-r2)*r3)^2, (r1*r2*(1-r3))^2, ((1-r1)*(1-r2)*(1-r3))^2)/4
      tm[28,] <- c(rep(0, 11), r3*(1-r3), r3*(1-r3), 0, 0, 0, r3*(1-r3)*c(r2^2, (1-r2)^2, (1-r2)^2, r2^2),
                   rep(0, 7), 1, 0, 0, (1-f2)*f2, (1-f2)*f2, rep(0, 10),
                   f2*(1-f2)*c((1-r1)^2, r1^2, r1^2, (1-r1)^2), 0, 0, 0,
                   rep(((1-r1)*(1-f2)+r1*f2)*(r1*(1-f2)+(1-r1)*f2), 2), 0, 0, 0,
                   r3*(1-r3)*c(f1^2, (1-f1)^2, f1^2, (1-f1)^2), rep(0, 6),
                   r3*(1-r3)*c(((1-r1)*r2)^2, (r1*r2)^2, ((1-r1)*(1-r2))^2, (r1*(1-r2))^2, ((1-r1)*r2)^2,
                               (r1*(1-r2))^2, (r1*r2)^2, ((1-r1)*(1-r2))^2))/2
      tm[29,] <- c(rep(0, 11), (1-r3)^2, r3^2, 1, 0, 0, (r2*(1-r3))^2, ((1-r2)*(1-r3))^2, ((1-r2)*r3)^2,
                   (r2*r3)^2, r2^2, (1-r2)^2, rep(0, 5), 1, 4, 0, f2^2, (1-f2)^2, 1, rep(0, 7), (1-r1)^2,
                   r1^2, ((1-r1)*f2)^2, (r1*(1-f2))^2, (r1*f2)^2, ((1-r1)*(1-f2))^2, 0, 0, 1,
                   (r1*(1-f2)+(1-r1)*f2)^2, ((1-r1)*(1-f2)+r1*f2)^2, 0, 0, 0, ((1-r3)*f1)^2, ((1-r3)*(1-f1))^2,
                   (r3*f1)^2, (r3*(1-f1))^2, f1^2, (1-f1)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2, ((1-r1)*(1-r2))^2,
                   (r1*r2)^2, ((1-r1)*r2*(1-r3))^2, (r1*r2*(1-r3))^2, ((1-r1)*(1-r2)*(1-r3))^2, (r1*(1-r2)*r3)^2,
                   ((1-r1)*r2*r3)^2, (r1*(1-r2)*(1-r3))^2, (r1*r2*r3)^2, ((1-r1)*(1-r2)*r3)^2)/4
      tm[30,] <- c(rep(0, 14), r2*(1-r2), r2*(1-r2), r2*(1-r2)*c(r3^2, r3^2, (1-r3)^2, (1-r3)^2), rep(0, 9), 1,
                   f2*(1-f2), f2*(1-f2), rep(0, 10), rep(r1*(1-r1)*f2*(1-f2), 4), r1*(1-r1), r1*(1-r1),
                   rep(0, 12), rep(r1*(1-r1)*r2*(1-r2), 4),
                   r1*(1-r1)*r2*(1-r2)*c(r3^2, r3^2, r3^2, (1-r3)^2, (1-r3)^2, r3^2, (1-r3)^2, (1-r3)^2))/2
      tm[31,] <- c(rep(0, 16), rep(r2*(1-r2)*r3*(1-r3), 4), rep(0, 10), (1-f2)^2, f2^2, rep(0, 10),
                   r1*(1-r1)*c((1-f2)^2, f2^2, (1-f2)^2, f2^2), rep(0, 18),
                   rep(r1*(1-r1)*r2*(1-r2)*r3*(1-r3), 8))/2
      tm[32,] <- c(rep(0, 16), rep(r2*(1-r2)*r3*(1-r3), 4), rep(0, 10), f2^2, (1-f2)^2, rep(0, 10),
                   r1*(1-r1)*c(f2^2, (1-f2)^2, f2^2, (1-f2)^2), rep(0,18),
                   rep(r1*(1-r1)*r2*(1-r2)*r3*(1-r3),8))/2
      tm[33,] <- c(rep(0, 16), r2*(1-r2)*c((1-r3)^2, (1-r3)^2, r3^2, r3^2), r2*(1-r2), r2*(1-r2), rep(0, 8),
                   f2*(1-f2), f2*(1-f2), 1, rep(0, 7), r1*(1-r1), r1*(1-r1), rep(r1*(1-r1)*f2*(1-f2), 4),
                   rep(0, 14), rep(r1*(1-r1)*r2*(1-r2), 4),
                   r1*(1-r1)*r2*(1-r2)*c((1-r3)^2, (1-r3)^2, (1-r3)^2, r3^2, r3^2, (1-r3)^2, r3^2, r3^2))/2
      tm[34,] <- c(rep(0, 14), (1-r2)^2, r2^2, ((1-r2)*r3)^2, (r2*r3)^2, ((1-r3)*r2)^2, ((1-r2)*(1-r3))^2, 0, 0, 1,
                   r3^2, (1-r3)^2, 0, 0, 0, 0, 1, f2^2, (1-f2)^2, 0, 4, 1, 0, 0, ((1-r1)*(1-f2)+r1*f2)^2,
                   ((1-r1)*f2+r1*(1-f2))^2, 1, 0, 0, (r1*f2)^2, ((1-r1)*(1-f2))^2, ((1-r1)*f2)^2, (r1*(1-f2))^2,
                   r1^2, (1-r1)^2, rep(0, 6), (r3*f1)^2, (r3*(1-f1))^2, ((1-r3)*f1)^2, ((1-r3)*(1-f1))^2, f1^2,
                   (1-f1)^2, (r1*(1-r2))^2, ((1-r1)*r2)^2, (r1*r2)^2, ((1-r1)*(1-r2))^2, (r1*(1-r2)*r3)^2,
                   ((1-r1)*(1-r2)*r3)^2, (r1*r2*r3)^2, ((1-r1)*r2*(1-r3))^2, (r1*(1-r2)*(1-r3))^2, ((1-r1)*r2*r3)^2,
                   ((1-r1)*(1-r2)*(1-r3))^2, (r1*r2*(1-r3))^2)/4
      tm[35,] <- c(rep(0, 16), r3*(1-r3)*c((1-r2)^2, r2^2, r2^2, (1-r2)^2, 0, 0, 0, 1, 1), rep(0, 5), f2*(1-f2),
                   f2*(1-f2), 0, 0, 1, 0, 0, rep(((1-r1)*(1-f2)+r1*f2)*((1-r1)*f2+r1*(1-f2)), 2), 0, 0, 0,
                   f2*(1-f2)*c(r1^2, (1-r1)^2, (1-r1)^2, r1^2), rep(0, 8),
                   r3*(1-r3)*c(f1^2, (1-f1)^2, f1^2, (1-f1)^2), rep(0, 6),
                   r3*(1-r3)*c((r1*(1-r2))^2, ((1-r1)*(1-r2))^2, (r1*r2)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2,
                               ((1-r1)*r2)^2, ((1-r1)*(1-r2))^2, (r1*r2)^2))/2
      tm[36,] <- c(rep(0, 16), ((1-r2)*(1-r3))^2, (r2*(1-r3))^2, (r2*r3)^2, ((1-r2)*r3)^2, (1-r2)^2, r2^2, 0,
                   (1-r3)^2, r3^2, 1, 0, 0, 0, 0, (1-f2)^2, f2^2, 1, 0, 1, 4, 1, ((1-r1)*f2+r1*(1-f2))^2,
                   ((1-r1)*(1-f2)+r1*f2)^2, 0, r1^2, (1-r1)^2, (r1*(1-f2))^2, ((1-r1)*f2)^2, ((1-r1)*(1-f2))^2,
                   (r1*f2)^2, rep(0, 6), f1^2, (1-f1)^2, ((1-r3)*f1)^2, ((1-r3)*(1-f1))^2, (r3*f1)^2, (r3*(1-f1))^2,
                   0, 0, (r1*(1-r2))^2, (r2*(1-r1))^2, (r1*r2)^2, ((1-r1)*(1-r2))^2, (r1*(1-r2)*(1-r3))^2,
                   ((1-r1)*(1-r2)*(1-r3))^2, (r1*r2*(1-r3))^2, ((1-r1)*r2*r3)^2, (r1*(1-r2)*r3)^2,
                   ((1-r1)*r2*(1-r3))^2, ((1-r1)*(1-r2)*r3)^2, (r1*r2*r3)^2)/4
      tm[37,] <- c(rep(0, 36), 1, rep(((1-r1)*(1-f2)+r1*f2)*((1-r1)*f2+r1*(1-f2)), 2), 0,
                   r1*(1-r1)*c(1, 1, (1-f1)^2, f1^2, (1-f1)^2, f1^2), rep(0, 6),
                   f1*(1-f1)*c(1, 1, (1-r3)^2, (1-r3)^2, r3^2, r3^2), 0, 0,
                   r1*(1-r1)*c((1-r2)^2, r2^2, r2^2, (1-r2)^2),
                   r1*(1-r1)*c(((1-r2)*(1-r3))^2, ((1-r2)*(1-r3))^2, (r2*(1-r3))^2, (r2*r3)^2, ((1-r2)*r3)^2,
                               (r2*(1-r3))^2, ((1-r2)*r3)^2, (r2*r3)^2))/2
      tm[38,] <- c(rep(0, 37), ((1-r1)*(1-f2)+r1*f2)^2, ((1-r1)*f2+r1*(1-f2))^2, 0, 0, 0, rep(r1*(1-r1)*f1*(1-f1), 4),
                   rep(0, 8), rep(r3*(1-r3)*f1*(1-f1), 4), rep(0, 6),
                   r1*(1-r1)*r3*(1-r3)*c((1-r2)^2, (1-r2)^2, r2^2, r2^2, (1-r2)^2, r2^2, (1-r2)^2, r2^2))/2
      tm[39,] <- c(rep(0, 37), ((1-r1)*f2+r1*(1-f2))^2, ((1-r1)*(1-f2)+r1*f2)^2, 0, 0, 0, rep(r1*(1-r1)*f1*(1-f1), 4),
                   rep(0, 8), rep(r3*(1-r3)*f1*(1-f1), 4), rep(0, 6),
                   r1*(1-r1)*r3*(1-r3)*c((1-r2)^2, (1-r2)^2, r2^2, r2^2, (1-r2)^2, r2^2, (1-r2)^2, r2^2))/2
      tm[40,] <- c(rep(0, 37), rep(((1-r1)*(1-f2)+r1*f2)*((1-r1)*f2+r1*(1-f2)), 2), 1, 0, 0,
                   r1*(1-r1)*c(f1^2, (1-f1)^2, f1^2, (1-f1)^2, 1, 1), rep(0, 6),
                   f1*(1-f1)*c(r3^2, r3^2, (1-r3)^2, (1-r3)^2, 1, 1), r1*(1-r1)*c((1-r2)^2, r2^2, r2^2, (1-r2)^2),
                   r1*(1-r1)*c(((1-r2)*r3)^2, ((1-r2)*r3)^2, (r2*r3)^2, (r2*(1-r3))^2, ((1-r2)*(1-r3))^2,
                               (r2*r3)^2, ((1-r2)*(1-r3))^2, (r2*(1-r3))^2))/2
      tm[41,] <- c(rep(0, 40), (1-r1)^2, r1^2, f1*(1-f1)*c((1-r1)^2, r1^2, r1^2, (1-r1)^2), rep(0, 14),
                   r2*(1-r2)*c((1-r1)^2, r1^2, (1-r1)^2, r1^2),
                   r2*(1-r2)*c(((1-r1)*(1-r3))^2, (r1*(1-r3))^2, ((1-r1)*(1-r3))^2, (r1*r3)^2,
                               ((1-r1)*r3)^2, (r1*(1-r3))^2, (r1*r3)^2, ((1-r1)*r3)^2))/2
      tm[42,] <- c(rep(0, 40), r1^2, (1-r1)^2, f1*(1-f1)*c(r1^2, (1-r1)^2, (1-r1)^2, r1^2), rep(0, 14),
                   r2*(1-r2)*c(r1^2, (1-r1)^2, r1^2, (1-r1)^2),
                   r2*(1-r2)*c((r1*(1-r3))^2, ((1-r1)*(1-r3))^2, (r1*(1-r3))^2, ((1-r1)*r3)^2, (r1*r3)^2,
                               ((1-r1)*(1-r3))^2, ((1-r1)*r3)^2, (r1*r3)^2))/2
      tm[43,] <- c(rep(0, 42), ((1-r1)*(1-f1))^2, (r1*f1)^2, (r1*(1-f1))^2, ((1-r1)*f1)^2, rep(0, 18),
                   r2*(1-r2)*r3*(1-r3)*((1-r1)^2*c(1, 0, 1, 0, 1, 0, 0, 1)+r1^2*c(0, 1, 0, 1, 0, 1, 1, 0)))/2
      tm[44,] <- c(rep(0, 42), (r1*f1)^2, ((1-r1)*(1-f1))^2, ((1-r1)*f1)^2, (r1*(1-f1))^2, rep(0, 18),
                   r2*(1-r2)*r3*(1-r3)*(r1^2*c(1, 0, 1, 0, 1, 0, 0, 1)+(1-r1)^2*c(0, 1, 0, 1, 0, 1, 1, 0)))/2
      tm[45,] <- c(rep(0, 42), (r1*(1-f1))^2, ((1-r1)*f1)^2, ((1-r1)*(1-f1))^2, (r1*f1)^2, rep(0, 18),
                   r2*(1-r2)*r3*(1-r3)*(r1^2*c(1, 0, 1, 0, 1, 0, 0, 1)+(1-r1)^2*c(0, 1, 0, 1, 0, 1, 1, 0)))/2
      tm[46,] <- c(rep(0, 42), ((1-r1)*f1)^2, (r1*(1-f1))^2, (r1*f1)^2, ((1-r1)*(1-f1))^2, rep(0, 18),
                   r2*(1-r2)*r3*(1-r3)*((1-r1)^2*c(1, 0, 1, 0, 1, 0, 0, 1)+r1^2*c(0, 1, 0, 1, 0, 1, 1, 0)))/2
      tm[47,] <- c(rep(0, 42), f1*(1-f1)*c((1-r1)^2, r1^2, r1^2, (1-r1)^2), (1-r1)^2, r1^2, rep(0, 12),
                   r2*(1-r2)*c((1-r1)^2, r1^2, (1-r1)^2, r1^2),
                   r2*(1-r2)*c(((1-r1)*r3)^2, (r1*r3)^2, ((1-r1)*r3)^2, (r1*(1-r3))^2, ((1-r1)*(1-r3))^2,
                               (r1*r3)^2, (r1*(1-r3))^2, ((1-r1)*(1-r3))^2))/2
      tm[48,] <- c(rep(0, 42), f1*(1-f1)*c(r1^2, (1-r1)^2, (1-r1)^2, r1^2), r1^2, (1-r1)^2, rep(0, 12),
                   r2*(1-r2)*c(r1^2, (1-r1)^2, r1^2, (1-r1)^2),
                   r2*(1-r2)*c((r1*r3)^2, ((1-r1)*r3)^2, (r1*r3)^2, ((1-r1)*(1-r3))^2, (r1*(1-r3))^2,
                               ((1-r1)*r3)^2, ((1-r1)*(1-r3))^2, (r1*(1-r3))^2))/2
      tm[49,] <- c(rep(0, 40), r1*(1-r1)*c(1, 1, f1^2, (1-f1)^2, f1^2, (1-f1)^2), 0, 0, 1,
                   rep(((1-r1)*(1-f1)+r1*f1)*((1-r1)*f1+r1*(1-f1)), 2), rep(0, 3),
                   f1*(1-f1)*c((1-r3)^2, (1-r3)^2, r3^2, r3^2, 1, 1),
                   r1*(1-r1)*c(r2^2, (1-r2)^2, (1-r2)^2, r2^2),
                   r1*(1-r1)*c((r2*(1-r3))^2, (r2*(1-r3))^2, ((1-r2)*(1-r3))^2, ((1-r2)*r3)^2, (r2*r3)^2,
                               ((1-r2)*(1-r3))^2, (r2*r3)^2, ((1-r2)*r3)^2))/2
      tm[50,] <- c(rep(0, 42), rep(r1*(1-r1)*f1*(1-f1), 4), 0, 0, 0, ((1-r1)*(1-f1)+r1*f1)^2,
                   ((1-r1)*f1+r1*(1-f1))^2, 0, 0, 0, rep(r3*(1-r3)*f1*(1-f1), 4), rep(0, 6),
                   r1*(1-r1)*r3*(1-r3)*(r2^2*c(1, 1, 0, 0, 1, 0, 1, 0)+(1-r2)^2*c(0, 0, 1, 1, 0, 1, 0, 1)))/2
      tm[51,] <- c(rep(0, 42), rep(r1*(1-r1)*f1*(1-f1), 4), 0, 0, 0, ((1-r1)*f1+r1*(1-f1))^2,
                   ((1-r1)*(1-f1)+r1*f1)^2, 0, 0, 0, rep(r3*(1-r3)*f1*(1-f1), 4), rep(0, 6),
                   r1*(1-r1)*r3*(1-r3)*(r2^2*c(1, 1, 0, 0, 1, 0, 1, 0)+(1-r2)^2*c(0, 0, 1, 1, 0, 1, 0, 1)))/2
      tm[52,] <- c(rep(0, 42), r1*(1-r1)*c((1-f1)^2, f1^2, (1-f1)^2, f1^2, 1, 1), 0,
                   rep(((1-r1)*(1-f1)+r1*f1)*((1-r1)*f1+r1*(1-f1)), 2), 1,
                   f1*(1-f1)*c(1, 1, r3^2, r3^2, (1-r3)^2, (1-r3)^2), 0, 0,
                   r1*(1-r1)*c(r2^2, (1-r2)^2, (1-r2)^2, r2^2),
                   r1*(1-r1)*c((r2*r3)^2, (r2*r3)^2, ((1-r2)*r3)^2, ((1-r2)*(1-r3))^2, (r2*(1-r3))^2,
                               ((1-r2)*r3)^2, (r2*(1-r3))^2, ((1-r2)*(1-r3))^2))/2
      tm[53,] <- c(rep(0, 52), (1-f1)^2, f1^2, r3*(1-r3)*c((1-f1)^2, f1^2, (1-f1)^2, f1^2), 0, 0,
                   rep(r1*(1-r1)*r2*(1-r2), 4), rep(r1*(1-r1)*r2*(1-r2)*r3*(1-r3), 8))/2
      tm[54,] <- c(rep(0, 52), f1^2, (1-f1)^2, r3*(1-r3)*c(f1^2, (1-f1)^2, f1^2, (1-f1)^2), 0, 0,
                   rep(r1*(1-r1)*r2*(1-r2), 4), rep(r1*(1-r1)*r2*(1-r2)*r3*(1-r3), 8))/2
      tm[55,] <- c(rep(0, 54), ((1-r3)*(1-f1))^2, ((1-r3)*f1)^2, (r3*(1-f1))^2, (r3*f1)^2, rep(0, 6),
                   r1*(1-r1)*r2*(1-r2)*((1-r3)^2*c(1, 1, 1, 0, 0, 1, 0, 0)+r3^2*c(0, 0, 0, 1, 1, 0, 1, 1)))/2
      tm[56,] <- c(rep(0, 54), ((1-r3)*f1)^2, ((1-r3)*(1-f1))^2, (r3*f1)^2, (r3*(1-f1))^2, rep(0, 6),
                   r1*(1-r1)*r2*(1-r2)*((1-r3)^2*c(1, 1, 1, 0, 0, 1, 0, 0)+r3^2*c(0, 0, 0, 1, 1, 0, 1, 1)))/2
      tm[57,] <- c(rep(0, 54), (r3*(1-f1))^2, (r3*f1)^2, ((1-r3)*(1-f1))^2, ((1-r3)*f1)^2, rep(0, 6),
                   r1*(1-r1)*r2*(1-r2)*(r3^2*c(1, 1, 1, 0, 0, 1, 0, 0)+(1-r3)^2*c(0, 0, 0, 1, 1, 0, 1, 1)))/2
      tm[58,] <- c(rep(0, 54), (r3*f1)^2, (r3*(1-f1))^2, ((1-r3)*f1)^2, ((1-r3)*(1-f1))^2, rep(0, 6),
                   r1*(1-r1)*r2*(1-r2)*(r3^2*c(1, 1, 1, 0, 0, 1, 0, 0)+(1-r3)^2*c(0, 0, 0, 1, 1, 0, 1, 1)))/2
      tm[59,] <- c(rep(0, 54), r3*(1-r3)*c((1-f1)^2, f1^2, (1-f1)^2, f1^2), (1-f1)^2, f1^2,
                   rep(r1*(1-r1)*r2*(1-r2), 4), rep(r1*(1-r1)*r2*(1-r2)*r3*(1-r3), 8))/2
      tm[60,] <- c(rep(0, 54), r3*(1-r3)*c(f1^2, (1-f1)^2, f1^2, (1-f1)^2), f1^2, (1-f1)^2,
                   rep(r1*(1-r1)*r2*(1-r2), 4), rep(r1*(1-r1)*r2*(1-r2)*r3*(1-r3), 8))/2
      tm[61,] <- c(rep(0, 60), ((1-r1)*(1-r2))^2, (r1*r2)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2,
                   r3*(1-r3)*c(((1-r1)*(1-r2))^2, (r1*(1-r2))^2, ((1-r1)*r2)^2, (r1*r2)^2,
                               ((1-r1)*(1-r2))^2, (r1*r2)^2, (r1*(1-r2))^2, ((1-r1)*r2)^2))/2
      tm[62,] <- c(rep(0, 60), (r1*r2)^2, ((1-r1)*(1-r2))^2, (r1*(1-r2))^2, ((1-r1)*r2)^2,
                   r3*(1-r3)*c((r1*r2)^2, (r2*(1-r1))^2, ((1-r2)*r1)^2, ((1-r1)*(1-r2))^2, (r1*r2)^2,
                               ((1-r1)*(1-r2))^2, ((1-r1)*r2)^2, (r1*(1-r2))^2))/2
      tm[63,] <- c(rep(0, 60), ((1-r1)*r2)^2, (r1*(1-r2))^2, ((1-r1)*(1-r2))^2, (r1*r2)^2,
                   r3*(1-r3)*c(((1-r1)*r2)^2, (r1*r2)^2, ((1-r1)*(1-r2))^2, (r1*(1-r2))^2,
                               ((1-r1)*r2)^2, (r1*(1-r2))^2, (r1*r2)^2, ((1-r1)*(1-r2))^2))/2
      tm[64,] <- c(rep(0, 60), (r1*(1-r2))^2, ((1-r1)*r2)^2, (r1*r2)^2, ((1-r1)*(1-r2))^2,
                   r3*(1-r3)*c((r1*(1-r2))^2, ((1-r1)*(1-r2))^2, (r1*r2)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2,
                               ((1-r1)*r2)^2, ((1-r1)*(1-r2))^2, (r1*r2)^2))/2
      tm[65,] <- c(rep(0, 64), ((1-r1)*(1-r2)*(1-r3))^2, (r1*(1-r2)*(1-r3))^2, ((1-r1)*r2*(1-r3))^2,
                   (r1*r2*r3)^2, ((1-r1)*(1-r2)*r3)^2, (r1*r2*(1-r3))^2, (r1*(1-r2)*r3)^2,
                   ((1-r1)*r2*r3)^2)/2
      tm[66,] <- c(rep(0, 64), (r1*(1-r2)*(1-r3))^2, ((1-r1)*(1-r2)*(1-r3))^2, (r1*r2*(1-r3))^2, ((1-r1)*r2*r3)^2,
                   (r1*(1-r2)*r3)^2, ((1-r1)*r2*(1-r3))^2, ((1-r1)*(1-r2)*r3)^2, (r1*r2*r3)^2)/2
      tm[67,] <- c(rep(0, 64), ((1-r1)*r2*(1-r3))^2, (r1*r2*(1-r3))^2, ((1-r1)*(1-r2)*(1-r3))^2, (r1*(1-r2)*r3)^2,
                   ((1-r1)*r2*r3)^2, (r1*(1-r2)*(1-r3))^2, (r1*r2*r3)^2, ((1-r1)*(1-r2)*r3)^2)/2
      tm[68,] <- c(rep(0, 64), (r1*r2*r3)^2, ((1-r1)*r2*r3)^2, (r1*(1-r2)*r3)^2, ((1-r1)*(1-r2)*(1-r3))^2,
                   (r1*r2*(1-r3))^2, ((1-r1)*(1-r2)*r3)^2, ((1-r1)*r2*(1-r3))^2, (r1*(1-r2)*(1-r3))^2)/2
      tm[69,] <- c(rep(0, 64), ((1-r1)*(1-r2)*r3)^2, (r1*(1-r2)*r3)^2, ((1-r1)*r2*r3)^2, (r1*r2*(1-r3))^2,
                   ((1-r1)*(1-r2)*(1-r3))^2, (r1*r2*r3)^2, (r1*(1-r2)*(1-r3))^2, ((1-r1)*r2*(1-r3))^2)/2
      tm[70,] <- c(rep(0, 64), (r1*r2*(1-r3))^2, ((1-r1)*r2*(1-r3))^2, (r1*(1-r2)*(1-r3))^2, ((1-r1)*(1-r2)*r3)^2,
                   (r1*r2*r3)^2, ((1-r1)*(1-r2)*(1-r3))^2, ((1-r1)*r2*r3)^2, (r1*(1-r2)*r3)^2)/2
      tm[71,] <- c(rep(0, 64), (r1*(1-r2)*r3)^2, ((1-r1)*(1-r2)*r3)^2, (r1*r2*r3)^2, ((1-r1)*r2*(1-r3))^2,
                   (r1*(1-r2)*(1-r3))^2, ((1-r1)*r2*r3)^2, ((1-r1)*(1-r2)*(1-r3))^2, (r1*r2*(1-r3))^2)/2
      tm[72,] <- c(rep(0, 64), ((1-r1)*r2*r3)^2, (r1*r2*r3)^2, ((1-r1)*(1-r2)*r3)^2, (r1*(1-r2)*(1-r3))^2,
                   ((1-r1)*r2*(1-r3))^2, (r1*(1-r2)*r3)^2, (r1*r2*(1-r3))^2, ((1-r1)*(1-r2)*(1-r3))^2)/2
      return(tm)
    }

    RI4 <- function(d1, d2, d3, n.Generation){
      r1 <- (1-exp(-2*d1))/2
      r2 <- (1-exp(-2*d2))/2
      r3 <- (1-exp(-2*d3))/2
      n.Iteration <- n.Generation-1
      Freq <- rep(0, 72)
      Freq[65] <- 1
      ZygoteFreq <- matrix(rep(0, (72*n.Generation-72)), nrow = 72)
      for (i in 1:n.Iteration) {
        Freq <- RI4TM(r1, r2, r3)%*%Freq
        ZygoteFreq[, i] <- Freq
      }
      re <- ZygoteFreq[, n.Iteration]
      return(re)
    }

    V963R <- function(d1, d2, nG){
      f.left <- RI2(d1, nG)
      f.right <- RI2(d2, nG)
      f <- RI3(d1, d2, nG)
      f.right <- f.right[c(1, 3, 4, 6, 7, 9)]
      floor <- matrix(1/f.left, ncol = 1)%*%matrix(1/f.right, nrow = 1)
      ceiling <- matrix(rep(0, length(f.left)*length(f.right)), nrow = length(f.left))
      ceiling[1, 1:2] <- c(f[1], f[8])
      ceiling[2, 3:4] <- c(f[2], f[9])
      ceiling[3, 5:6] <- c(f[3], f[10])
      ceiling[4, 1:2] <- c(f[11], f[14])
      ceiling[5, 3:4] <- rep((f[12]+f[13]), 2)
      ceiling[6, 5:6] <- c(f[14], f[11])
      ceiling[7, 1:2] <- c(f[10], f[3])
      ceiling[8, 3:4] <- c(f[9], f[2])
      ceiling[9, 5:6] <- c(f[8], f[1])
      block <- floor*ceiling
      return(block)
    }

    V663R <- function(d1, d2, nG){
      f.left <- RI2(d1, nG)
      f.right <- RI2(d2, nG)
      f <- RI3(d1, d2, nG)
      f.left <- f.left[c(1, 3, 4, 6, 7, 9)]
      f.right <- f.right[c(1, 3, 4, 6, 7, 9)]
      floor <- matrix(1/f.left, ncol = 1)%*%matrix(1/f.right, nrow = 1)
      ceiling <- matrix(rep(0, length(f.left)^2), ncol = length(f.left))
      ceiling[1, 1:2] <- c(f[1], f[8])
      ceiling[2, 5:6] <- c(f[3], f[10])
      ceiling[3, 1:2] <- c(f[11], f[14])
      ceiling[4, 5:6] <- c(f[14], f[11])
      ceiling[5, 1:2] <- c(f[10], f[3])
      ceiling[6, 5:6] <- c(f[8], f[1])
      block <- floor*ceiling
      return(block)
    }

    V964R <- function(d1, d2, d3, nG){
      f.left <- RI2(d1, nG)
      f.right <- RI2(d3, nG)
      f <- RI4(d1, d2, d3, nG)
      f.right <- f.right[c(1, 3, 4, 6, 7, 9)]
      floor <- matrix(1/f.left, ncol = 1)%*%matrix(1/f.right, nrow = 1)
      ceiling <- matrix(rep(0, length(f.left)*length(f.right)), nrow = length(f.left))
      ceiling[, 1] <- c(f[1], f[4], f[8], f[37], f[41]+f[42], f[49], f[36], f[33], f[29])
      ceiling[, 2] <- c(f[3], f[7], f[10], f[40], f[47]+f[48], f[52], f[34], f[30], f[27])
      ceiling[, 3] <- c(f[11], f[15]+f[16], f[23], f[53]+f[54], sum(f[61:64]), f[59]+f[60],
                        f[26], f[21]+f[22], f[14])
      ceiling[, 4] <- rev(ceiling[, 3])
      ceiling[, 5] <- rev(ceiling[, 2])
      ceiling[, 6] <- rev(ceiling[, 1])
      block <- floor*ceiling
      return(block)
    }

    V664R <- function(d1, d2, d3, nG){
      f.left <- RI2(d1, nG)
      f.right <- RI2(d3, nG)
      f <- RI4(d1, d2, d3, nG)
      f.left <- f.left[c(1, 3, 4, 6, 7, 9)]
      f.right <- f.right[c(1, 3, 4, 6, 7, 9)]
      floor <- matrix(1/f.left, ncol = 1)%*%matrix(1/f.right, nrow = 1)
      ceiling <- matrix(rep(0, length(f.left)^2), ncol = length(f.left))
      ceiling[, 1] <- c(f[1], f[8], f[37], f[49], f[36], f[29])
      ceiling[, 2] <- c(f[3], f[10], f[40], f[52], f[34], f[27])
      ceiling[, 3] <- c(f[11], f[23], f[53]+f[54], f[59]+f[60], f[26], f[14])
      ceiling[, 4] <- rev(ceiling[, 3])
      ceiling[, 5] <- rev(ceiling[, 2])
      ceiling[, 6] <- rev(ceiling[, 1])
      block <- floor*ceiling
      return(block)
    }

    cov4R=function (d, nG, nS, ss){
      nI <- length(d)
      dim <- 9+6*(nI-1)
      cov <- matrix(rep(0, (dim^2)), ncol = dim)
      cov[(1:9), (10:15)] <- V963R(d[1], d[2], nG)
      for(i in 2:(nI-1)){
        cov[(9+6*(i-2)+(1:6)), (15+6*(i-2)+(1:6))] <- V663R(d[i], d[i+1], nG)
      }
      for(i in 2:(nI-1)){
        cov[(1:9), (15+6*(i-2)+(1:6))] <- V964R(d[1], sum(d[2:i]), d[i+1], nG)
      }
      for(i in 2:(nI-2)){
        for (j in (i+1):(nI-1)) {
          cov[(9+6*(i-2)+1:6), (15+6*(i-1)+6*(j-i-1)+1:6)] <- V664R(d[i], sum(d[(i+1):j]), d[j+1], nG)
        }
      }
      cov <- (cov+t(cov))*ss/nS
      diagnal <- rep(0, dim)
      diagnal[1:9] <- RI2(d[1], nG)
      for(i in 2:nI){
        diagnal[(9+6*(i-2)+(1:6))] <- RI2(d[i], nG)[c(1, 3, 4, 6, 7, 9)]
      }
      diagnal <- 1/diagnal
      diag(cov) <- diagnal*ss/nS
      return(cov)
    }

    cov1R <- function(d, nG, nS, ss){
      cov <- matrix(rep(0, 81), ncol = 9)
      diagnal <- RI2(d, nG)
      diag(cov) <- (1/diagnal)*ss/nS
      return(cov)
    }

    cov2R <- function(d, nG, nS, ss){
      cov <- matrix(rep(0, 225), ncol = 15)
      cov[(1:9), (10:15)] <- V963R(d[1], d[2], nG)
      cov <- (cov+t(cov))*ss/nS
      diagnal <- rep(0, 15)
      diagnal[1:9] <- RI2(d[1], nG)
      diagnal[10:15] <- RI2(d[2], nG)[c(1, 3, 4, 6, 7, 9)]
      diagnal <- 1/diagnal
      diag(cov) <- diagnal*ss/nS
      return(cov)
    }

    cov3R <- function(d, nG, nS, ss){
      cov <- matrix(rep(0, 441), ncol = 21)
      cov[(1:9), (10:15)] <- V963R(d[1], d[2], nG)
      cov[(10:15), (16:21)] <- V663R(d[2], d[3], nG)
      cov[(1:9), (16:21)] <- V964R(d[1], d[2], d[3], nG)
      cov <- (cov+t(cov))*ss/nS
      diagnal <- rep(0, 21)
      diagnal[1:9] <- RI2(d[1], nG)
      diagnal[10:15] <- RI2(d[2], nG)[c(1, 3, 4, 6, 7, 9)]
      diagnal[16:21] <- RI2(d[3], nG)[c(1, 3, 4, 6, 7, 9)]
      diagnal <- 1/diagnal
      diag(cov) <- diagnal*ss/nS
      return(cov)
    }

    sigmaR <- function(d, nG, nS, ss){
      x <- length(d)
      if(x > 3){
        sigma <- cov4R(d, nG, nS, ss)
      } else {
        sigma <- switch(x, cov1R(d, nG, nS, ss), cov2R(d, nG, nS, ss), cov3R(d, nG, nS, ss))
      }
      return(sigma)
    }

    ZR <- function(y, Ls, Ns, Vr, Go, speed){
      Ni <- length(Ls)
      Lsr0 <- cumsum(Ls)
      sl <- sum(Ls)
      Lsp <- seq(0, sl, speed/100)
      Lspl <- length(Lsp)
      Lsr1 <- c(0, Lsr0[-length(Lsr0)])
      loci0 <- c(Lsr0[1])
      loci1 <- c(0)
      k0 <- c(1)
      for(i in 2:Lspl){
        loci0[i] <- min(Lsr0[Lsr0 >= Lsp[i]])
        loci1[i] <- max(Lsr1[Lsr1 < Lsp[i]])
        k0[i] <- which(Lsr0 == min(Lsr0[Lsr0 >= Lsp[i]]))
      }

      TSC <- NULL
      ff0 <- matrix(0, Ni, 9)
      for(i in 1:Ni){
        ff0[i,] <- RI2(Ls[Ni], Go)
      }
      ff <- ff0[1,]

      for(x in 1:(min(which(Lsp>=Ls[1]))-1)){
        loci <- Lsp[x]
        mar <- Lsr0[1]
        K <- Diff.k(RI3(loci, mar-loci, Go))
        L <- Diff.l(RI3(loci, mar-loci, Go))
        u1 <- sum((K-sum(K*ff))*ff*y[, 1])
        u2 <- sum((L-0.5^(Go-2)+1)*ff*y[, 1])
        a <- sum((K-sum(K*ff))^2*ff)
        b <- sum((L-0.5^(Go-2)+1)^2*ff)
        c <- sum((K-sum(K*ff))*(L-0.5^(Go-2)+1)*ff)
        num <- a*u2^2+b*u1^2-2*c*u1*u2
        den <- a*b-c^2
        Usq <- (Ns/Vr)*num/den
        TSC[x] <- Usq
      }
      if(Ni > 1){
        for (x2 in (x+1):Lspl) {
          ff <- ff0[k0[x2],]
          loci <- Lsp[length(TSC)+1]
          K <- Diff.k(RI3(loci-loci1[x2], loci0[x2]-loci, Go))
          L <- Diff.l(RI3(loci-loci1[x2], loci0[x2]-loci, Go))
          u1 <- sum((K-sum(K*ff))*ff*y[, k0[x2]])
          u2 <- sum((L-0.5^(Go-2)+1)*ff*y[, k0[x2]])
          a <- sum((K-sum(K*ff))^2*ff)
          b <- sum((L-0.5^(Go-2)+1)^2*ff)
          c <- sum((K-sum(K*ff))*(L-0.5^(Go-2)+1)*ff)
          num <- a*u2^2+b*u1^2-2*c*u1*u2
          den <- a*b-c^2
          Usq <- (Ns/Vr)*num/den
          TSC[x2] <- Usq
        }
      }
      re <- max(TSC)
      return(re)
    }

    thre.fun <- function(d, nG, simu, nS, ss, console = FALSE, cr = 0, speed = 1){
      cov <- sigmaR(d, nG, nS, ss)
      rawy <- mvtnorm::rmvnorm(simu, rep(0, ncol(cov)), cov)
      TS <- rep(0, simu)

      for(j in 1:simu){
        if(length(d) > 1){
          y <- matrix(rep(1, (9*length(d))), nrow = 9)
          y[, 1] <- rawy[j, 1:9]
          for (i in 2:length(d)) {
            y[c(1, 3), i] <- rawy[j, 9+6*(i-2)+1:2]
            y[c(4, 6), i] <- rawy[j, 9+6*(i-2)+3:4]
            y[c(7, 9), i] <- rawy[j, 9+6*(i-2)+5:6]
          }
          for(i in 1:(length(d)-1)){
            minuend <- RI2(d[i], nG)*y[, i]
            subtrahend <- RI2(d[i+1], nG)*y[, i+1]
            y[2, i+1] <- (sum(minuend[c(1, 4, 7)])-sum(subtrahend[c(1, 3)]))/subtrahend[2]
            y[5, i+1] <- (sum(minuend[c(2, 5, 8)])-sum(subtrahend[c(4, 6)]))/subtrahend[5]
            y[8, i+1] <- (sum(minuend[c(3, 6, 9)])-sum(subtrahend[c(7, 9)]))/subtrahend[8]
          }
        } else {
          y <- matrix(rawy[j, ], nrow = 9)
        }
        TS[j] <- ZR(y, d, nS, ss, nG, speed)
        if(j%%cs == 0){
          cat(paste(paste("RI F", nG, " ad   ", sep = ""), cr, j, "\n", sep = "\t")[console])
        }
      }
      return(TS)
      cat(paste(paste("RI F", nG, " ad   ", sep = ""), cr, "done", "\n", sep = "\t")[console])
    }

    if(!d.eff){
      thre.fun <- function(d, nG, simu, nS, ss, console = FALSE, cr = 0, speed = 1){
        U.RI <- function(y, Ls, Ns, Vr, Go, sp){
          Ni <- length(Ls)
          Lsr0 <- cumsum(Ls)
          sl <- sum(Ls)
          Lsp <- seq(0, sl, sp/100)
          Lspl <- length(Lsp)
          Lsr1 <- c(0, Lsr0[-length(Lsr0)])
          loci0 <- c(Lsr0[1])
          loci1 <- c(0)
          k0 <- c(1)
          for(i in 2:Lspl){
            loci0[i] <- min(Lsr0[Lsr0 >= Lsp[i]])
            loci1[i] <- max(Lsr1[Lsr1 < Lsp[i]])
            k0[i] <- which(Lsr0 == min(Lsr0[Lsr0 >= Lsp[i]]))
          }

          TSC <- NULL
          ff0 <- matrix(0, Ni, 9)
          for(i in 1:Ni){
            ff0[i,] <- RI2(Ls[Ni], Go)
          }
          ff <- ff0[1,]

          for(x in 1:(min(which(Lsp >= Ls[1]))-1)){
            loci <- Lsp[x]
            mar <- Lsr0[1]
            Diff <- Diff.P(RI3(loci, mar-loci, Go))
            num <- sum(Diff*ff*y[, 1])
            den <- sum(Diff^2*ff)
            Usq <- (Ns/Vr)*num^2/den
            TSC[x] <- Usq
          }
          if(Ni > 1){
            for(x2 in (x+1):Lspl){
              ff <- ff0[k0[x2],]
              loci <- Lsp[length(TSC)+1]
              Diff <- Diff.P(RI3(loci-loci1[x2], loci0[x2]-loci, Go))
              num <- sum(Diff*ff*y[, k0[x2]])
              den <- sum(Diff^2*ff)
              Usq <- (Ns/Vr)*num^2/den
              TSC[x2] <- Usq
            }
          }
          max(TSC, na.rm = TRUE)
        }

        cov <- sigmaR(d, nG, nS, ss)
        rawy <- mvtnorm::rmvnorm(simu, rep(0, ncol(cov)), cov)
        TS <- rep(0, simu)


        for (j in 1:simu){
          if(length(d) > 1){
            y <- matrix(rep(1, (9*length(d))), nrow = 9)
            y[, 1] <- rawy[j, 1:9]
            for (i in 2:length(d)) {
              y[c(1, 3), i] <- rawy[j, 9+6*(i-2)+1:2]
              y[c(4, 6), i] <- rawy[j, 9+6*(i-2)+3:4]
              y[c(7, 9), i] <- rawy[j, 9+6*(i-2)+5:6]
            }
            for(i in 1:(length(d)-1)){
              minuend <- RI2(d[i], nG)*y[, i]
              subtrahend <- RI2(d[i+1], nG)*y[, i+1]
              y[2, i+1] <- (sum(minuend[c(1, 4, 7)])-sum(subtrahend[c(1, 3)]))/subtrahend[2]
              y[5, i+1] <- (sum(minuend[c(2, 5, 8)])-sum(subtrahend[c(4, 6)]))/subtrahend[5]
              y[8, i+1] <- (sum(minuend[c(3, 6, 9)])-sum(subtrahend[c(7, 9)]))/subtrahend[8]
            }
          } else {
            y <- matrix(rawy[j, ], nrow = 9)
          }
          TS[j] <- U.RI(y, d, nS, ss, nG, speed)
          if(j%%cs == 0){
            cat(paste(paste("RI F", nG, " a   ", sep = ""), cr, j, "\n", sep = "\t")[console])
          }
        }
        cat(paste(paste("RI F", nG, " a   ", sep = ""), cr, "done", "\n", sep = "\t")[console])
        return(TS)
      }
    }
  } else if (type == "AI"){
    AI2 <- function(d, n.Generation){
      r <- (1-exp(-2*d))/2
      ZygoteFreq <- matrix(rep(0, (9*n.Generation-9)), nrow = 9)
      GameteFreq <- c((1-r)/2, r/2)
      C <- GameteFreq[1]^2
      D <- GameteFreq[2]^2
      E <- 2*GameteFreq[1]*GameteFreq[2]
      F <- 2*GameteFreq[1]^2
      G <- 2*GameteFreq[2]^2
      ZygoteFreq[, 1] <- c(C, E, D, E, F+G, E, D, E, C)
      if (n.Generation > 2) {
        for(i in 3:(n.Generation)){
          GameteFreq <- (1-r)*GameteFreq+r/4
          C <- GameteFreq[1]^2
          D <- GameteFreq[2]^2
          E <- 2*GameteFreq[1]*GameteFreq[2]
          F <- 2*GameteFreq[1]^2
          G <- 2*GameteFreq[2]^2
          ZygoteFreq[, (i-1)] <- c(C, E, D, E, F+G, E, D, E, C)
        }
      }
      ZygoteFreq[, (n.Generation-1)]
    }

    AI3 <- function(d1, d2, n.Generation){
      d <- d1+d2
      r1 <- (1-exp(-2*d1))/2
      r2 <- (1-exp(-2*d2))/2
      r <- (1-exp(-2*d))/2
      ZygoteFreq <- matrix(rep(0, (20*n.Generation-20)), nrow = 20)
      N <- (1-r1)*(1-r2)/2
      R <- (1-r1)*r2/2
      B <- r1*r2/2
      L <- r1*(1-r2)/2
      COEF <- rep(2, 20)
      COEF[c(1, 3, 8, 10)] <- rep(1, 4)
      GS <- c(N, N, B, N, N, R, B, R, R, L, N, N, B, R, N, R, N, B, R, L)
      GE <- c(N, B, B, R, L, B, L, R, L, L, L, R, L, B, B, L, N, B, R, L)
      ZygoteFreq[, 1] <- COEF*GS*GE
      if(n.Generation > 2){
        for(i in 3:n.Generation){
          F23 <- c(N+L, B+R, B+R, N+L)
          F13 <- rep(c(N+B, L+R), 2)
          F12 <- rep(c(N+R, L+B), each = 2)
          NewGameteFreq <- (1-r1)*(1-r2)*c(N, R, B, L)+0.5*r1*(1-r2)*F23+0.5*r1*r2*F13+0.5*(1-r1)*r2*F12
          N <- NewGameteFreq[1]
          R <- NewGameteFreq[2]
          B <- NewGameteFreq[3]
          L <- NewGameteFreq[4]
          GS <- c(N, N, B, N, N, R, B, R, R, L, N, N, B, R, N, R, N, B, R, L)
          GE <- c(N, B, B, R, L, B, L, R, L, L, L, R, L, B, B, L, N, B, R, L)
          ZygoteFreq[, (i-1)] <- COEF*GS*GE
        }
      }
      re <- ZygoteFreq[, (n.Generation-1)]
      return(re)
    }

    AI4 <- function(d1, d2, d3, n.Generation){
      d <- d1+d2+d3
      da <- d1+d2
      db <- d2+d3
      r <- (1-exp(-2*d))/2
      r1 <- (1-exp(-2*d1))/2
      r2 <- (1-exp(-2*d2))/2
      r3 <- (1-exp(-2*d3))/2
      ra <- (1-exp(-2*da))/2
      rb <- (1-exp(-2*db))/2
      ZygoteFreq <- matrix(rep(0, (72*n.Generation-72)), nrow = 72)
      GG1 <- (1-r1)*(1-r2)*(1-r3)/2
      GG2 <- (1-r1)*(1-r2)*r3/2
      GG3 <- (1-r1)*r2*r3/2
      GG4 <- (1-r1)*r2*(1-r3)/2
      GG5 <- r1*r2*(1-r3)/2
      GG6 <- r1*r2*r3/2
      GG7 <- r1*(1-r2)*r3/2
      GG8 <- r1*(1-r2)*(1-r3)/2
      COEF <- rep(2, 72)
      COEF[c(1, 3, 8, 10, 27, 29, 34, 36)] <- rep(1, 8)
      GS <- c(GG1, GG1, GG2, GG1, GG1, GG2, GG2, GG5, GG5, GG6, GG1, GG1, GG2, GG2, GG1, GG3, GG1, GG4, GG6,
              GG7, GG2, GG4, GG5, GG5, GG6, GG6, GG3, GG3, GG4, GG3, GG3, GG4, GG4, GG7, GG7, GG8, GG1, GG1,
              GG2, GG2, GG1, GG5, GG1, GG5, GG6, GG2, GG2, GG6, GG5, GG5, GG6, GG6, GG1, GG3, GG1, GG4, GG2,
              GG3, GG2, GG4, GG1, GG5, GG3, GG7, GG1, GG8, GG4, GG6, GG2, GG5, GG7, GG3)
      GE <- c(GG1, GG2, GG2, GG5, GG6, GG5, GG6, GG5, GG6, GG6, GG3, GG4, GG3, GG4, GG7, GG5, GG8, GG5, GG3,
              GG2, GG8, GG6, GG7, GG8, GG7, GG8, GG3, GG4, GG4, GG7, GG8, GG7, GG8, GG7, GG8, GG8, GG8, GG7,
              GG8, GG7, GG4, GG8, GG3, GG7, GG8, GG4, GG3, GG7, GG4, GG3, GG4, GG3, GG6, GG8, GG5, GG8, GG6,
              GG7, GG5, GG7, GG2, GG6, GG4, GG8, GG1, GG8, GG4, GG6, GG2, GG5, GG7, GG3)
      ZygoteFreq[, 1] <- COEF*GS*GE
      if(n.Generation > 2){
        for(i in 3:n.Generation){
          x <- c(GG1, GG2, GG3, GG4, GG5, GG6, GG7, GG8)
          x1 <- c(GG1+GG8, GG2+GG7, GG3+GG6, GG4+GG5, GG4+GG5, GG3+GG6, GG2+GG7, GG1+GG8)
          x2 <- rep(c(GG1+GG5, GG2+GG6, GG3+GG7, GG4+ GG8), 2)
          x3 <- c(GG1+GG3, GG2+GG4, GG1+GG3, GG2+GG4, GG5+GG7, GG6+GG8, GG5+GG7, GG6+GG8)
          x4 <- c(GG1+GG2, GG1+GG2, GG3+GG4, GG3+GG4, GG5+GG6, GG5+GG6, GG7+GG8, GG7+GG8)
          x12 <- rep(c(GG1+GG2+GG3+GG4, GG5+GG6+ GG7+GG8), each = 4)
          x34 <- rep(c((GG1+GG4+GG5+GG8), (GG2+GG3+GG6+GG7), (GG2+GG3+GG6+GG7), (GG1+GG4+GG5+GG8)), 2)
          x13 <- rep(c(GG1+GG2+GG5+GG6, GG1+GG2+GG5+GG6, GG3+GG4+GG7+GG8, GG3+GG4+GG7+GG8), 2)
          x24 <- c(rep(c(GG1+GG3+GG6+GG8, GG2+GG4+GG5+GG7), 2), rep(c(GG2+GG4+GG5+GG7, GG1+GG3+GG6+GG8), 2))
          x14 <- rep(c(GG1+GG3+GG5+GG7, GG2+GG4+GG6+GG8), 4)
          x23 <- c(rep(GG1+GG2+GG7+GG8, 2), rep(GG3+GG4+GG5+GG6, 4), rep(GG1+GG2+GG7+GG8, 2))
          NewGameteFreq <- (1-r1)*(1-r2)*(1-r3)*x+0.5*r1*(1-r2)*(1-r3)*x1+0.5*r1*r2*(1-r3)*x2+
            0.5*(1-r1)*r2*r3*x3+0.5*(1-r1)*(1-r2)*r3*x4+(1-r1)*r2*(1-r3)*x12*x34+r1*r2*r3*x13*x24+
            r1*(1-r2)*r3*x14*x23
          GG1 <- NewGameteFreq[1]
          GG2 <- NewGameteFreq[2]
          GG3 <- NewGameteFreq[3]
          GG4 <- NewGameteFreq[4]
          GG5 <- NewGameteFreq[5]
          GG6 <- NewGameteFreq[6]
          GG7 <- NewGameteFreq[7]
          GG8 <- NewGameteFreq[8]
          GS <- c(GG1, GG1, GG2, GG1, GG1, GG2, GG2, GG5, GG5, GG6, GG1, GG1, GG2, GG2, GG1, GG3, GG1, GG4,
                  GG6, GG7, GG2, GG4, GG5, GG5, GG6, GG6, GG3, GG3, GG4, GG3, GG3, GG4, GG4, GG7, GG7, GG8,
                  GG1, GG1, GG2, GG2, GG1, GG5, GG1, GG5, GG6, GG2, GG2, GG6, GG5, GG5, GG6, GG6, GG1, GG3,
                  GG1, GG4, GG2, GG3, GG2, GG4, GG1, GG5, GG3, GG7, GG1, GG8, GG4, GG6, GG2, GG5, GG7, GG3)
          GE <- c(GG1, GG2, GG2, GG5, GG6, GG5, GG6, GG5, GG6, GG6, GG3, GG4, GG3, GG4, GG7, GG5, GG8, GG5,
                  GG3, GG2, GG8, GG6, GG7, GG8, GG7, GG8, GG3, GG4, GG4, GG7, GG8, GG7, GG8, GG7, GG8, GG8,
                  GG8, GG7, GG8, GG7, GG4, GG8, GG3, GG7, GG8, GG4, GG3, GG7, GG4, GG3, GG4, GG3, GG6, GG8,
                  GG5, GG8, GG6, GG7, GG5, GG7, GG2, GG6, GG4, GG8, GG1, GG8, GG4, GG6, GG2, GG5, GG7, GG3)
          ZygoteFreq[, (i-1)] <- COEF*GS*GE
        }
      }
      ZygoteFreq[, (n.Generation-1)]
    }

    V963A <- function(d1, d2, nG){
      f.left <- AI2(d1, nG)
      f.right <- AI2(d2, nG)
      f <- AI3(d1, d2, nG)
      f.right <- f.right[c(1, 3, 4, 6, 7, 9)]
      floor <- matrix(1/f.left, ncol = 1)%*%matrix(1/f.right, nrow = 1)
      ceiling <- matrix(rep(0, length(f.left)*length(f.right)), nrow = length(f.left))
      ceiling[1, 1:2] <- c(f[1], f[8])
      ceiling[2, 3:4] <- c(f[2], f[9])
      ceiling[3, 5:6] <- c(f[3], f[10])
      ceiling[4, 1:2] <- c(f[11], f[14])
      ceiling[5, 3:4] <- rep((f[12]+f[13]), 2)
      ceiling[6, 5:6] <- c(f[14], f[11])
      ceiling[7, 1:2] <- c(f[10], f[3])
      ceiling[8, 3:4] <- c(f[9], f[2])
      ceiling[9, 5:6] <- c(f[8], f[1])
      block <- floor*ceiling
      return(block)
    }

    V663A <- function(d1, d2, nG){
      f.left <- AI2(d1, nG)
      f.right <- AI2(d2, nG)
      f <- AI3(d1, d2, nG)
      f.left <- f.left[c(1, 3, 4, 6, 7, 9)]
      f.right <- f.right[c(1, 3, 4, 6, 7, 9)]
      floor <- matrix(1/f.left, ncol = 1)%*%matrix(1/f.right, nrow = 1)
      ceiling <- matrix(rep(0, length(f.left)^2), ncol = length(f.left))
      ceiling[1, 1:2] <- c(f[1], f[8])
      ceiling[2, 5:6] <- c(f[3], f[10])
      ceiling[3, 1:2] <- c(f[11], f[14])
      ceiling[4, 5:6] <- c(f[14], f[11])
      ceiling[5, 1:2] <- c(f[10], f[3])
      ceiling[6, 5:6] <- c(f[8], f[1])
      block <- floor*ceiling
      return(block)
    }

    V964A <- function(d1, d2, d3, nG){
      f.left <- AI2(d1, nG)
      f.right <- AI2(d3, nG)
      f <- AI4(d1, d2, d3, nG)
      f.right <- f.right[c(1, 3, 4, 6, 7, 9)]
      floor <- matrix(1/f.left, ncol = 1) %*% matrix(1/f.right, nrow = 1)
      ceiling <- matrix(rep(0, length(f.left)*length(f.right)), nrow = length(f.left))
      ceiling[, 1] <- c(f[1], f[4], f[8], f[37], f[41]+f[42], f[49], f[36], f[33], f[29])
      ceiling[, 2] <- c(f[3], f[7], f[10], f[40], f[47]+f[48], f[52], f[34], f[30], f[27])
      ceiling[, 3] <- c(f[11], f[15]+f[16], f[23], f[53]+f[54], sum(f[61:64]), f[59]+f[60], f[26],
                        f[21]+f[22], f[14])
      ceiling[, 4] <- rev(ceiling[, 3])
      ceiling[, 5] <- rev(ceiling[, 2])
      ceiling[, 6] <- rev(ceiling[, 1])
      block <- floor*ceiling
      return(block)
    }

    V664A <- function(d1, d2, d3, nG){
      f.left <- AI2(d1, nG)
      f.right <- AI2(d3, nG)
      f <- AI4(d1, d2, d3, nG)
      f.left <- f.left[c(1, 3, 4, 6, 7, 9)]
      f.right <- f.right[c(1, 3, 4, 6, 7, 9)]
      floor <- matrix(1/f.left, ncol = 1) %*% matrix(1/f.right, nrow = 1)
      ceiling <- matrix(rep(0, length(f.left)^2), ncol = length(f.left))
      ceiling[, 1] <- c(f[1], f[8], f[37], f[49], f[36], f[29])
      ceiling[, 2] <- c(f[3], f[10], f[40], f[52], f[34], f[27])
      ceiling[, 3] <- c(f[11], f[23], f[53]+f[54], f[59]+f[60], f[26], f[14])
      ceiling[, 4] <- rev(ceiling[, 3])
      ceiling[, 5] <- rev(ceiling[, 2])
      ceiling[, 6] <- rev(ceiling[, 1])
      block <- floor*ceiling
      return(block)
    }

    cov4A <- function(d, nG, nS, ss){
      nI <- length(d)
      dim <- 9+6*(nI-1)
      cov <- matrix(rep(0, (dim^2)), ncol = dim)
      cov[(1:9), (10:15)] <- V963A(d[1], d[2], nG)
      for(i in 2:(nI-1)){
        cov[(9+6*(i-2)+(1:6)), (15+6*(i-2)+(1:6))] <- V663A(d[i], d[i+1], nG)
      }
      for(i in 2:(nI-1)){
        cov[(1:9), (15+6*(i-2)+(1:6))] <- V964A(d[1], sum(d[2:i]), d[i+1], nG)
      }
      for(i in 2:(nI-2)){
        for(j in (i+1):(nI-1)){
          cov[(9+6*(i-2)+1:6), (15+6*(i-1)+6*(j-i-1)+1:6)] <- V664A(d[i], sum(d[(i+1):j]), d[j+1], nG)
        }
      }
      cov <- (cov+t(cov))*ss/nS
      diagnal <- rep(0, dim)
      diagnal[1:9] <- AI2(d[1], nG)
      for (i in 2:nI) {
        diagnal[(9+6*(i-2)+(1:6))] <- AI2(d[i], nG)[c(1, 3, 4, 6, 7, 9)]
      }
      diagnal <- 1/diagnal
      diag(cov) <- diagnal*ss/nS
      return(cov)
    }

    cov1A <- function(d, nG, nS, ss){
      cov <- matrix(rep(0, 81), ncol = 9)
      diagnal <- AI2(d, nG)
      diag(cov) <- (1/diagnal)*ss/nS
      return(cov)
    }

    cov2A <- function(d, nG, nS, ss){
      cov <- matrix(rep(0, 225), ncol = 15)
      cov[(1:9), (10:15)] <- V963A(d[1], d[2], nG)
      cov <- (cov+t(cov))*ss/nS
      diagnal <- rep(0, 15)
      diagnal[1:9] <- AI2(d[1], nG)
      diagnal[10:15] <- AI2(d[2], nG)[c(1, 3, 4, 6, 7, 9)]
      diagnal <- 1/diagnal
      diag(cov) <- diagnal*ss/nS
      return(cov)
    }

    cov3A <- function(d, nG, nS, ss){
      cov <- matrix(rep(0, 441), ncol = 21)
      cov[(1:9), (10:15)] <- V963A(d[1], d[2], nG)
      cov[(10:15), (16:21)] <- V663A(d[2], d[3], nG)
      cov[(1:9), (16:21)] <- V964A(d[1], d[2], d[3], nG)
      cov <- (cov+t(cov))*ss/nS
      diagnal <- rep(0, 21)
      diagnal[1:9] <- AI2(d[1], nG)
      diagnal[10:15] <- AI2(d[2], nG)[c(1, 3, 4, 6, 7, 9)]
      diagnal[16:21] <- AI2(d[3], nG)[c(1, 3, 4, 6, 7, 9)]
      diagnal <- 1/diagnal
      diag(cov) <- diagnal*ss/nS
      return(cov)
    }

    sigmaA <- function(d, nG, nS, ss){
      x <- length(d)
      if(x > 3){
        sigma <- cov4A(d, nG, nS, ss)
      } else {
        sigma <- switch(x, cov1A(d, nG, nS, ss), cov2A(d, nG, nS, ss), cov3A(d, nG, nS, ss))
      }
      return(sigma)
    }

    ZA <- function(y, Ls, Ns, Vr, Go, sp){
      Ni <- length(Ls)
      Lsr0 <- cumsum(Ls)
      sl <- sum(Ls)
      Lsp <- seq(0, sl, sp/100)
      Lspl <- length(Lsp)
      Lsr1 <- c(0,Lsr0[-length(Lsr0)])
      loci0 <- c(Lsr0[1])
      loci1 <- c(0)
      k0 <- c(1)
      for(i in 2:Lspl){
        loci0[i] <- min(Lsr0[Lsr0 >= Lsp[i]])
        loci1[i] <- max(Lsr1[Lsr1 < Lsp[i]])
        k0[i] <- which(Lsr0 == min(Lsr0[Lsr0 >= Lsp[i]]))
      }

      TSC <- NULL
      ff0 <- matrix(0, Ni, 9)
      for(i in 1:Ni){
        ff0[i,] <- AI2(Ls[Ni], Go)
      }
      ff <- ff0[1,]

      for(x in 1:(min(which(Lsp >= Ls[1]))-1)){
        loci <- Lsp[x]
        mar <- Lsr0[1]
        K <- Diff.k(AI3(loci, mar-loci, Go))
        L <- Diff.l(AI3(loci, mar-loci, Go))
        u1 <- sum((K-sum(K*ff))*ff*y[, 1])
        u2 <- sum((L-sum(L*ff))*ff*y[, 1])
        a <- sum((K-sum(K*ff))^2*ff)
        b <- sum((L-sum(L*ff))^2*ff)
        c <- sum((K-sum(K*ff))*(L-sum(L*ff))*ff)
        num <- a*u2^2+b*u1^2-2*c*u1*u2
        den <- a*b-c^2
        Usq <- (Ns/Vr)*num/den
        TSC[x] <- Usq
      }
      if(Ni > 1){
        for(x2 in (x+1):Lspl){
          ff <- ff0[k0[x2],]
          loci <- Lsp[length(TSC)+1]
          K <- Diff.k(AI3(loci-loci1[x2], loci0[x2]-loci, Go))
          L <- Diff.l(AI3(loci-loci1[x2], loci0[x2]-loci, Go))
          u1 <- sum((K-sum(K*ff))*ff*y[, k0[x2]])
          u2 <- sum((L-sum(L*ff))*ff*y[, k0[x2]])
          a <- sum((K-sum(K*ff))^2*ff)
          b <- sum((L-sum(L*ff))^2*ff)
          c <- sum((K-sum(K*ff))*(L-sum(L*ff))*ff)
          num <- a*u2^2+b*u1^2-2*c*u1*u2
          den <- a*b-c^2
          Usq <- (Ns/Vr)*num/den
          TSC[x2] <- Usq
        }
      }
      re <- max(TSC)
      return(re)
    }

    thre.fun <- function(d, nG, simu, nS, ss, console = FALSE, cr = 0, speed = 1){
      cov <- sigmaA(d, nG, nS, ss)
      rawy <- mvtnorm::rmvnorm(simu, rep(0, ncol(cov)), cov)
      TS <- rep(0, simu)

      for(j in 1:simu){
        if(length(d) > 1){
          y <- matrix(rep(1, (9*length(d))), nrow = 9)
          y[, 1] <- rawy[j, 1:9]
          for(i in 2:length(d)){
            y[c(1, 3), i] <- rawy[j, 9+6*(i-2)+1:2]
            y[c(4, 6), i] <- rawy[j, 9+6*(i-2)+3:4]
            y[c(7, 9), i] <- rawy[j, 9+6*(i-2)+5:6]
          }
          for(i in 1:(length(d)-1)){
            minuend <- AI2(d[i], nG)*y[, i]
            subtrahend <- AI2(d[i+1], nG)*y[, i+1]
            y[2, i+1] <- (sum(minuend[c(1, 4, 7)])-sum(subtrahend[c(1, 3)]))/subtrahend[2]
            y[5, i+1] <- (sum(minuend[c(2, 5, 8)])-sum(subtrahend[c(4, 6)]))/subtrahend[5]
            y[8, i+1] <- (sum(minuend[c(3, 6, 9)])-sum(subtrahend[c(7, 9)]))/subtrahend[8]
          }
        } else {
          y <- matrix(rawy[j, ], nrow = 9)
        }
        TS[j] <- ZA(y, d, nS, ss, nG, speed)
        if(j%%cs == 0){
          cat(paste(paste("AI F", nG, " ad   ", sep = ""), cr, j, "\n", sep = "\t")[console])
        }
      }
      cat(paste(paste("AI F", nG, " ad   ", sep = ""), cr, "done", "\n", sep = "\t")[console])
      return(TS)
    }

    if(!d.eff){
      thre.fun <- function(d, nG, simu, nS, ss, console = FALSE, cr = 0, speed = 1){
        U.AI <- function(y, Ls, Ns, Vr, Go, sp){
          Ni <- length(Ls)
          Lsr0 <- cumsum(Ls)
          sl <- sum(Ls)
          Lsp <- seq(0,sl,sp/100)
          Lspl <- length(Lsp)
          Lsr1 <- c(0,Lsr0[-length(Lsr0)])
          loci0 <- c(Lsr0[1])
          loci1 <- c(0)
          k0 <- c(1)
          for(i in 2:Lspl){
            loci0[i] <- min(Lsr0[Lsr0 >= Lsp[i]])
            loci1[i] <- max(Lsr1[Lsr1 < Lsp[i]])
            k0[i] <- which(Lsr0 == min(Lsr0[Lsr0 >= Lsp[i]]))
          }

          TSC <- NULL
          ff0 <- matrix(0, Ni, 9)
          for(i in 1:Ni){
            ff0[i,] <- AI2(Ls[Ni], Go)
          }
          ff <- ff0[1,]

          for(x in 1:(min(which(Lsp >= Ls[1]))-1)){
            loci <- Lsp[x]
            mar <- Lsr0[1]
            Diff <- Diff.P(AI3(loci, mar-loci, Go))
            num <- sum(Diff*ff*y[, 1])
            den <- sum(Diff^2*ff)
            Usq <- (Ns/Vr)*num^2/den
            TSC[x] <- Usq
          }
          if(Ni > 1){
            for(x2 in (x+1):Lspl){
              ff <- ff0[k0[x2],]
              loci <- Lsp[length(TSC)+1]
              Diff <- Diff.P(AI3(loci-loci1[x2], loci0[x2]-loci, Go))
              num <- sum(Diff*ff*y[, k0[x2]])
              den <- sum(Diff^2*ff)
              Usq <- (Ns/Vr)*num^2/den
              TSC[x2] <- Usq
            }
          }
          re <- max(TSC, na.rm = TRUE)
          return(re)
        }

        cov <- sigmaA(d, nG, nS, ss)
        rawy <- mvtnorm::rmvnorm(simu, rep(0, ncol(cov)), cov)
        TS <- rep(0, simu)

        for(j in 1:simu){
          if(length(d) > 1){
            y <- matrix(rep(1, (9*length(d))), nrow = 9)
            y[, 1] <- rawy[j, 1:9]
            for(i in 2:length(d)){
              y[c(1, 3), i] <- rawy[j, 9+6*(i-2)+1:2]
              y[c(4, 6), i] <- rawy[j, 9+6*(i-2)+3:4]
              y[c(7, 9), i] <- rawy[j, 9+6*(i-2)+5:6]
            }
            for(i in 1:(length(d)-1)){
              minuend <- AI2(d[i], nG)*y[, i]
              subtrahend <- AI2(d[i+1], nG)*y[, i+1]
              y[2, i+1] <- (sum(minuend[c(1, 4, 7)])-sum(subtrahend[c(1, 3)]))/subtrahend[2]
              y[5, i+1] <- (sum(minuend[c(2, 5, 8)])-sum(subtrahend[c(4, 6)]))/subtrahend[5]
              y[8, i+1] <- (sum(minuend[c(3, 6, 9)])-sum(subtrahend[c(7, 9)]))/subtrahend[8]
            }
          } else {
            y <- matrix(rawy[j, ], nrow = 9)
          }
          TS[j] <- U.AI(y, d, nS, ss, nG, speed)
          if(j%%cs == 0){
            cat(paste(paste("AI F", nG, " a   ", sep = ""), cr, j, "\n", sep = "\t")[console])
          }
        }
        cat(paste(paste("AI F", nG, " a   ", sep = ""), cr, "done", "\n", sep = "\t")[console])
        return(TS)
      }
    }
  }

  cat(paste("Generation", "chr", "simulation times", "\n", sep = "\t")[console])
  for(i in 1:ncr){
    d <- marker[marker[, 1] == cr[i], 2]
    d <- d[-1]-d[-length(d)]
    thre0 <- thre.fun(d, ng, simu, ns, gv, console, cr[i], speed)
    thre <- cbind(thre, thre0)
  }

  max0 <- function(x){max(x, na.rm = TRUE)}
  permu.max <- apply(thre,1,max0)

  thre.all <- stats::quantile(permu.max, (1-alpha), type = 3)
  return(thre.all)
}
