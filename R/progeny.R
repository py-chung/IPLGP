#' Progeny Simulation
#'
#' Generate simulated phenotype and genotype data for a specified
#' generation from various breeding schemes.
#'
#' @param QTL matrix. A q*2 matrix contains the QTL information, where the
#' row dimension 'q' represents the number of QTLs in the chromosomes. The
#' first column labels the chromosomes where the QTLs are located, and the
#' second column labels the positions of QTLs (in morgan (M) or centimorgan
#' (cM)).
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
#' @param  E.vector vector. Set the effect of QTLs. It should be a named
#' vector, where the names of elements represent the effects of QTLs and
#' their interactions. For example: the additive effect of QTL1 is coded as
#' "a1"; the dominant effect of QTL2 is coded as "d2"; the interaction of
#' the additive effect of QTL2 and the dominant effect of QTL1 is coded as
#' "a2:d1". So, if the additive effect of QTL1 is 2, the dominant effect of
#' QTL2 is 5, and the interaction of the additive effect of QTL2 and the
#' dominant effect of QTL1 is 3, the user should input: E.vector = c("a1"=2,
#' "d2"=5, "a2:d1"=3). If E.vector=NULL, the phenotypic value will not be
#' simulated.
#' @param h2 numeric. Set the heritability for simulated phenotypes. It
#' should be a number between 0 and 1.
#' @param size numeric. The population size of simulated progeny.
#'
#' @return
#' \item{phe}{The phenotypic value of each simulated progeny.}
#' \item{E.vector}{The effect vector used in this simulation.}
#' \item{marker.prog}{The marker genotype of each simulated progeny.}
#' \item{QTL.prog}{The QTL genotype of each simulated progeny.}
#' \item{VG}{The genetic variance of this population.}
#' \item{VE}{The environmental variance of this population.}
#' \item{genetic.value}{The genetic value of each simulated progeny.}
#'
#' @export
#'
#' @references
#'
#' Haldane J.B.S. 1919. The combination of linkage values and the calculation
#' of distance between the loci for linked factors. Genetics 8: 299–309.
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' result <- progeny(QTL, marker, type = "RI", ng = 5, E.vector = c("a1" = 2, "d2" = 5, "a2:d1" = 3),
#' h2 = 0.5, size = 200)
#' result$phe
progeny <- function(QTL, marker, type = "RI", ng = 2, cM = TRUE, E.vector = NULL, h2 = 0.5, size = 200){

  if(is.null(QTL) | is.null(marker)){
    stop("Input data is missing, please cheak and fix.", call. = FALSE)
  }

  marker <- as.matrix(marker)
  markertest <- c(ncol(marker) != 2, NA %in% marker, marker[,1] != sort(marker[,1]))

  datatry <- try(marker*marker, silent=TRUE)
  if(class(datatry)[1] == "try-error" | TRUE %in% markertest){
    stop("Marker data error, or the number of marker does not match the genetype data.", call. = FALSE)
  }

  datatry <- try(QTL*QTL, silent=TRUE)
  if(class(datatry)[1] == "try-error" | ncol(QTL) != 2 | NA %in% QTL | max(QTL[, 1]) > max(marker[, 1])){
    stop("QTL data error, please cheak your QTL data.", call. = FALSE)
  }

  if(!type[1] %in% c("AI","RI","BC") | length(type) > 1){
    stop("Parameter type error, please input AI, RI, or BC.", call. = FALSE)
  }

  if(!is.numeric(ng) | length(ng) > 1 | min(ng) < 1){
    stop("Parameter ng error, please input a positive integer.", call. = FALSE)
  }
  ng <- round(ng)

  if(!cM[1] %in% c(0,1) | length(cM > 1)){cM <- TRUE}

  if(!is.null(E.vector)){
    E.vector <- c(E.vector)
    if(!is.numeric(E.vector) | is.null(names(E.vector))){
      stop("Parameter E.vector error, please check and fix.", call. = FALSE)
    }
  }

  if(!is.numeric(h2) | length(h2) > 1 | min(h2) < 0 | max(h2) > 1){
    stop("Parameter h2 error, please input a positive number between 0 and 1.", call. = FALSE)
  }

  if(!is.numeric(size) | length(size) > 1 | min(size) < 1){
    stop("Parameter size error, please input a positive integer.", call. = FALSE)
  }

  if(!cM){
    QTL <- QTL*100
    marker[, 2] <- marker[, 2]*100
  }

  n <- as.numeric(nrow(QTL))
  QTLs <- cbind(QTL, 0)
  marker <- cbind(marker, 1:nrow(marker))
  marker0 <- rbind(as.matrix(marker), as.matrix(QTLs))
  marker0 <- marker0[order(marker0[, 1], marker0[, 2]), ]
  dis <- marker0[, 2]
  p1 <- rep(2, nrow(marker0))
  p2 <- rep(0, nrow(marker0))
  F1 <- matrix(c(p1, p2), nrow(marker0), 2)

  D.matrix <- D.make(n, type = type, aa = TRUE, dd = TRUE, ad = TRUE)
  if(is.null(E.vector)){E.vector <- rep(0, ncol(D.matrix))
  } else if (length(E.vector) != ncol(D.matrix)){
    E0 <- rep(0, ncol(D.matrix))
    for(i in 1:length(E.vector)){
      E0[colnames(D.matrix) %in% (names(E.vector)[i])] <- E.vector[i]
      if(!((names(E.vector)[i]) %in% colnames(D.matrix))){
        name0 <- strsplit(names(E.vector)[i], split = "")
        name1 <- c()
        for(k in c(4, 5, 3, 1, 2)){
          name1 <- paste(name1, name0[[1]][k], sep = "")
        }
        names(E.vector)[i] <- name1
        E0[colnames(D.matrix) %in% (names(E.vector)[i])] <- E.vector[i]
      }
    }
    E.vector <- E0
  }
  names(E.vector) <- colnames(D.matrix)

  simulate.progeny <- function(DNA0, distance){
    distance.c <- distance[-1]-distance[-length(distance)]
    distance.r <- (1-exp(-2*distance.c/100))/2
    distance.r[distance.r < 0] <- 0.5

    index <- sample(1:2, 1)
    DNA1 <- DNA0[, index]
    DNA2 <- DNA0[, -index]

    DNA.new <- DNA1[1]
    for(i in 2:nrow(DNA0)){
      x <- stats::runif(1)
      if(x < distance.r[i-1]){
        DNA.new[i] <- DNA2[i]
        D <- DNA1
        DNA1 <- DNA2
        DNA2 <- D
      } else { DNA.new[i] <- DNA1[i] }
    }

    return(DNA.new)
  }

  VG.RI <- function(QTL, E.vector, ng){
    n <- as.numeric(nrow(QTL))

    r.matrix <- matrix(0, n, n)
    for(i in 1:n){
      for(j in 1:i){
        Qij <- QTL[c(i, j),]
        if(Qij[1, 1] != Qij[2, 1]){
          r.matrix[i, j] <- 0.5
        } else {
          r.matrix[i, j] <- (1-exp(-2*(abs(Qij[1, 2]-Qij[2, 2]))/100))/2
        }
      }
    }
    r.matrix <- r.matrix+t(r.matrix)
    lamda <- (1-2*r.matrix)*(1-r.matrix)^(ng-2)
    lamda.aadd <- c()
    for(i in 1:(n-1)){
      lamda.aadd <- c(lamda.aadd, lamda[i, (i+1):n])
    }
    lamda.adda <- c()
    for(i in 1:n){
      lamda.adda <- c(lamda.adda, lamda[i, -i])
    }

    E.a <- E.vector[(1:n)*2-1]
    E.d <- E.vector[(1:n)*2]
    E.iaa <- E.vector[(n*2+1):(n*2+choose(n, 2))]
    E.iad <- E.vector[(n*2+choose(n, 2)*2+1):length(E.vector)]
    E.idd <- E.vector[(n*2+choose(n, 2)+1):(n*2+choose(n, 2)*2)]
    E.iad2 <- c()
    for(i in 1:n){
      s0 <- seq(i, n*(n-1), (n-1))
      s0[1:i] <- s0[1:i]-1
      s0 <- s0[-i]
      E.iad2 <- c(E.iad2, E.iad[s0])
    }
    E.iaa2 <- matrix(0, n, n)
    k0 <- 1
    for(i in 1:n){
      E.iaa2[i, -(1:i)] <- E.iaa[k0:(k0+n-i-1)]
      k0 <- k0+n-i
    }
    E.iaa2 <- E.iaa2+t(E.iaa2)
    diag(E.iaa2) <- NA
    E.iaa2 <- c(E.iaa2)
    E.iaa2 <- E.iaa2[!is.na(E.iaa2)]

    va <- sum(E.a^2/2)
    vd <- sum(E.d^2/4)
    viaa <- sum(E.iaa^2/4)
    viad <- sum(E.iad^2/8)
    vidd <- sum((1-lamda.aadd^4)*E.idd^2/16)
    vaa <- c()
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        vaa <- c(vaa, lamda.aadd[length(vaa)+1]*E.a[i]*E.a[j])
      }
    }
    vaa <- sum(vaa)
    vdd <- c()
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        vdd <- c(vdd, (lamda.aadd[length(vdd)+1]^2/2)*E.d[i]*E.d[j])
      }
    }
    vdd <- sum(vdd)
    vaiad <- c()
    for(i in 1:n){
      vaiad <- c(vaiad, (E.a[i]*E.iad[((i-1)*(n-1)+1):(i*(n-1))]*lamda.adda[((i-1)*(n-1)+1):(i*(n-1))]^2)/2)
    }
    vaiad <- -sum(vaiad)
    vaida <- c()
    for(i in 1:n){
      vaida <- c(vaida, (E.a[i]*E.iad2[((i-1)*(n-1)+1):(i*(n-1))]*lamda.adda[((i-1)*(n-1)+1):(i*(n-1))])/2)
    }
    vaida <- -sum(vaida)
    vdiaa <- c()
    for(i in 1:n){
      vdiaa <- c(vdiaa, (E.d[i]*E.iad2[((i-1)*(n-1)+1):(i*(n-1))]*lamda.adda[((i-1)*(n-1)+1):(i*(n-1))])/2)
    }
    vdiaa <- -sum(vdiaa)
    viaaidd <- sum(lamda.aadd*(1-lamda.aadd^2)*E.iaa*E.idd/2)
    viadida <- sum(lamda.adda*E.iad*E.iad2/4)
    VG <- va+vd+viaa+viad+vidd+vaa+vdd+vaiad+vaida+vdiaa+viaaidd+viadida
    return(VG)
  }

  VG.BC <- function(QTL, E.vector, ng){
    n <- nrow(QTL)

    r.matrix <- matrix(0, n, n)
    for(i in 1:n){
      for(j in 1:i){
        Qij <- QTL[c(i, j), ]
        if(Qij[1, 1] != Qij[2, 1]){
          r.matrix[i, j] <- r.matrix[j, i] <- 0.5
        } else {
          r.matrix[i, j] <- r.matrix[j, i] <- (1-exp(-2*(abs(Qij[1, 2]-Qij[2, 2]))/100))/2
        }
      }
    }
    lamda <- (1-2*r.matrix)*(1-r.matrix)^(ng-1)
    lamda.aadd <- c()
    for(i in 1:(n-1)){
      lamda.aadd <- c(lamda.aadd, lamda[i, (i+1):n])
    }

    E.a <- E.vector[1:n]
    E.iaa <- E.vector[-(1:n)]

    va <- sum(E.a^2/4)
    viaa <- sum(E.iaa^2/16)
    vaa <- c()
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        vaa <- c(vaa, lamda.aadd[length(vaa)+1]*E.a[i]*E.a[j]/2)
      }
    }
    vaa <- sum(vaa)

    VG <- va+viaa+vaa
    return(VG)
  }

  prog.g <- list()
  if (type == "BC"){
    p0 <- cbind(p1, p1)
    for(i in 1:size){
      pr1 <- simulate.progeny(F1, dis)
      pr2 <- simulate.progeny(p0, dis)
      prog.g[[i]] <- cbind(pr1, pr2)
    }
    if(ng > 1){
      for(g in 2:ng){
        prog.g2 <- list()
        for(i in 1:size){
          sele <- sample(1:size, 1)
          pr1 <- simulate.progeny(prog.g[[sele]], dis)
          pr2 <- simulate.progeny(p0, dis)
          prog.g2[[i]] <- cbind(pr1, pr2)
        }
        prog.g <- prog.g2
      }
    }
    if(n == 1){
      VG <- sum((E.vector[1:n])^2/4)
    } else {
      VG <- VG.BC(QTL, E.vector, ng)
    }
  } else {
    for(i in 1:size){
      pr1 <- simulate.progeny(F1, dis)
      pr2 <- simulate.progeny(F1, dis)
      prog.g[[i]] <- cbind(pr1, pr2)
    }
    if(ng > 2){
      if(type == "RI"){
        for(g in 3:ng){
          prog.g2 <- list()
          for(i in 1:size){
            sele <- sample(1:size, 1)
            pr1 <- simulate.progeny(prog.g[[sele]], dis)
            pr2 <- simulate.progeny(prog.g[[sele]], dis)
            prog.g2[[i]] <- cbind(pr1, pr2)
          }
          prog.g <- prog.g2
        }
      } else if (type == "AI"){
        for(g in 3:ng){
          prog.g2 <- list()
          for(i in 1:size){
            sele <- sample(1:size, 2)
            pr1 <- simulate.progeny(prog.g[[sele[1]]], dis)
            pr2 <- simulate.progeny(prog.g[[sele[2]]], dis)
            prog.g2[[i]] <- cbind(pr1, pr2)
          }
          prog.g <- prog.g2
        }
      }
    } else if (ng == 1){
      for(i in 1:size){
        prog.g[[i]] <- F1
      }
    }

    if(n == 1){
      VG <- sum((E.vector[(1:n)*2-1])^2/2)+sum((E.vector[(1:n)*2])^2/4)
    } else {
      VG <- VG.RI(QTL, E.vector, ng)
    }
  }
  prog <- matrix(0, size, length(dis))
  for(i in 1:size){
    prog[i,] <- apply(prog.g[[i]], 1, mean)
  }

  VE <- VG*(1-h2)/h2
  phe <- c()
  genetic.value <- c()
  for(i in 1:size){
    genoQ <- prog[i, marker0[, 3] == 0]
    genoQ1 <- genoQ[1]
    if(n>1){
      for(j in 2:n){
        genoQ1 <- paste(genoQ1, genoQ[j], sep = "")
      }
    }
    eff <- D.matrix%*%E.vector
    eff <- eff[row.names(eff) == genoQ1,]
    genetic.value[i] <- eff
    phe[i] <- eff+stats::rnorm(1, 0, sqrt(VE))
  }

  marker.prog <- prog[, marker0[, 3] != 0]
  QTL.prog <- prog[, marker0[, 3] == 0]

  output <- list(phe = phe, E.vector = E.vector, marker.prog = marker.prog,
                 QTL.prog = QTL.prog, genetic.value = genetic.value, VG = VG, VE = VE)
  return(output)
}
