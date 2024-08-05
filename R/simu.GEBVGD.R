#' Simulate Progeny with GEBV-GD Strategy
#'
#' Identify parental lines based on GEBV-GD strategy and simulate their offsprings.
#'
#' @param fittedA.t matrix. An n*t matrix denotes the fitted values of each traits
#' of the training population. The missing value must have been already imputed.
#' If outcross is set to be TRUE, this argument must be the additive effect part
#' of fitted values.
#' @param fittedD.t matrix. An n*t matrix denotes the dominance effect part of
#' fitted values when outcross is set to be TRUE. The missing value must have been
#' already imputed.
#' @param fittedmu.t numeric or vector. A p*1 vector denote the average value of
#' fitted values when outcross is set to be TRUE. The length must be the same as
#' the number of traits.
#' @param geno.t matrix. An n*p matrix denotes the marker score matrix of the
#' training population. The markers must be coded as 1, 0, or -1 for alleles
#' AA, Aa, or aa. The missing value must have been already imputed.
#' @param marker matrix. A p*2 matrix whose first column indicates the chromosome
#' number to which a marker belongs; and second column indicates the position of
#' the marker in centi-Morgan (cM).
#' @param geno.c matrix. An nc*p matrix denotes the marker score matrix of the
#' candidate population with nc individuals and p markers. It should be pure lines
#' and markers must be coded as 1, or -1 for alleles AA, or aa. The missing value
#' must have been already imputed. If geno.c is set to be NULL, the candidate
#' population is exactly the training population.
#' @param npl integer. An integer indicates the number of individuals who will
#' be chosen as the parental lines. If npl = NULL, it will be 4 times the number
#' of traits.
#' @param better.c logical. A logical variable, if better.c is set to be TRUE,
#' the candidate individuals with GEBVs better than average for all the target
#' traits will comprise the candidate set. Otherwise, all the candidate
#' individuals will comprise the candidate set.
#' @param npl.best integer. A integer indicates the numbers of the candidate
#' individuals with the top GEBV index will be retained. If npl.best is set to
#' be NULL, it will be 2 times the number of traits.
#' @param weight vector. A vector with length t indicates the weights of target
#' traits in selection index. If weight is set to be NULL, the equal weight will
#' be assigned to all the target traits. The weights should be a positive number.
#' @param direction vector. A vector with length t indicates the selecting
#' directions for target traits. The elements of direction are Inf, or -Inf
#' representing the rule that the larger the better; or the smaller the better.
#' Or if the element is a number, it will select the individuals with the trait
#' value close to the number. If direction is set to be NULL, the selecting
#' direction will be the larger the better for all trait.
#' @param outcross logical. A logical variable, if outcross is set to be TRUE,
#' the crop is regarded as an outcross crop. The kinship matrix of dominance
#' effects are also considered in the model, and crossing and selection will be
#' performed in F1 generation. The detail can be seen in the references.
#' @param nprog integer. An integer indicates the number of progenies which
#' will be produced for each of the best individuals at every generation.
#' @param nsele integer. An integer indicates the number of the best individuals
#' which will be selected at each generation. If nsele is set to be NULL, the
#' number will be the same as the number of F1 individuals.
#' @param ngen integer. An integer indicates the number of generations in the
#' simulation process.
#' @param nrep integer. An integer indicates the number of repetitions in the
#' simulation process.
#' @param cri integer. An integer indicates the stopping criterion, note that
#' cri < 1e+06. The genetic algorithm will stop if the number of iterations
#' reaches cri.
#' @param console logical. A logical variable, if console is set to be TRUE,
#' the simulation process will be shown in the R console.
#'
#' @return
#' \item{method}{The GEBV-GD strategy.}
#' \item{weight}{The weights of target traits in selection index.}
#' \item{direction}{The selecting directions of target traits in selection index.}
#' \item{mu}{The mean vector of target traits.}
#' \item{sd}{The standard deviation vector of target traits.}
#' \item{GEBV.value}{The GEBVs of target traits in each generation and each
#' repetition.}
#' \item{parental.lines}{The IDs and D-score of parental lines selected in
#' each repetition.}
#' \item{suggested.subset}{The most frequently selected parental lines by this
#' strategy.}
#'
#' @note
#' The function output.best and output.gain can be used to summarize the result.
#'
#' The fitted value data in the input data can be obtained by the function
#' GBLUP.fit and mmer, that can be seen in the Examples shown below.
#'
#' @export
#'
#' @references
#'
#' Chung PY, Liao CT. 2020. Identification of superior parental lines for
#' biparental crossing via genomic prediction. PLoS ONE 15(12):e0243159.
#'
#' @seealso
#' \code{\link[sommer]{mmer}}
#' \code{\link[IPLGP]{GBLUP.fit}}
#' \code{\link[IPLGP]{GA.Dscore}}
#' \code{\link[IPLGP]{simu.gamete}}
#' \code{\link[IPLGP]{simu.GEBVO}}
#' \code{\link[IPLGP]{simu.GEBVGD}}
#' \code{\link[IPLGP]{output.best}}
#' \code{\link[IPLGP]{output.gain}}
#'
#' @examples
#' # generate simulated data
#' set.seed(6000)
#' geno.test <- matrix(sample(c(1, -1), 200, replace = TRUE), 10, 20)
#' t1 <- 5*geno.test[,3]+3*geno.test[,7]-geno.test[,11]+rnorm(10,30,10)
#' t2 <- 3*geno.test[,3]+geno.test[,12]-2*geno.test[,18]+rnorm(10,10,5)
#' t3 <- NULL
#' t4 <- NULL
#' t5 <- NULL
#' marker.test <- cbind(rep(1:2, each=10), rep(seq(0, 90, 10), 2))
#' fit <- GBLUP.fit(t1, t2, t3, t4, t5, geno = geno.test)
#' fitvalue <- fit$fitted.value
#'
#' geno.candidate <- matrix(sample(c(1,-1), 300, replace = TRUE), 15, 20)
#'
#' # run and output
#' result <- simu.GEBVGD(fitvalue, geno.t = geno.test, marker = marker.test,
#' geno.c = geno.candidate, nprog = 5, nsele = 10, ngen = 5, nrep = 5, cri = 250)
#' result$suggested.subset
#'
#'
#'
#' # other method: use mmer to obtain the fitted value
#' \dontrun{
#' set.seed(6000)
#' geno.test <- matrix(sample(c(1, -1), 200, replace = TRUE), 10, 20)
#' t1 <- 5*geno.test[,3]+3*geno.test[,7]-geno.test[,11]+rnorm(10,30,10)
#' t2 <- 3*geno.test[,3]+geno.test[,12]-2*geno.test[,18]+rnorm(10,10,5)
#' phe <- cbind(t1, t2)
#' nt <- ncol(phe)
#' marker.test <- cbind(rep(1:2, each=10), rep(seq(0, 90, 10), 2))
#' rownames(geno.test) <- 1:nrow(geno.test)
#' id <- rownames(geno.test)
#' K0 <- geno.test%*%t(geno.test)/ncol(geno.test)
#'
#' dat <- data.frame(id, phe)
#' fit0 <- sommer::mmer(cbind(t1, t2)~1,
#'       random = ~sommer::vsr(id, Gu = K0, Gtc = sommer::unsm(nt)),
#'       rcov = ~sommer::vsr(units, Gtc = sommer::unsm(nt)),
#'       data = dat,
#'       tolParInv = 0.1)
#'
#' u0 <- fit0$U$`u:id`
#' fit <- matrix(unlist(u0), ncol = nt)
#' colnames(fit) <- names(u0)
#'
#' fit <- fit+matrix(fit0$fitted[1,], nrow(fit), nt, byrow = TRUE)
#' fitvalue <- fit[order(as.numeric(names((u0[[1]])))),]
#' }
simu.GEBVGD <- function(fittedA.t, fittedD.t = NULL, fittedmu.t =NULL, geno.t, marker, geno.c = NULL, npl = NULL,
                        better.c = FALSE, npl.best = NULL, weight = NULL, direction = NULL, outcross = FALSE,
                        nprog = 50, nsele = NULL, ngen = 10, nrep = 30, cri = 10000, console = TRUE){
  phe.t <- fittedA.t
  datatry <- try(phe.t*phe.t, silent=TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% phe.t){
    stop("fittedA.t error, please cheak your data and impute the missing value.", call. = FALSE)
  }

  datatry <- try(geno.t%*%t(geno.t), silent = TRUE)
  if(class(datatry)[1] == "try-error" | length(geno.t[geno.t != 1 & geno.t != 0 & geno.t != -1]) > 0){
    stop("Genotype data of training population error or have not been imputed.", call. = FALSE)
  }

  nt <- ncol(phe.t)
  ind.t <- nrow(phe.t)
  if(is.null(npl) | !is.numeric(npl)){
    npl <- 4*nt
  } else {npl <- npl[1]}
  if(is.null(npl.best) | !is.numeric(npl.best)){
    npl.best <- 2*nt
  } else {npl.best <- npl.best[1]}
  if(!is.numeric(weight) | length(weight) < nt){
    weight <- rep(1, nt)/nt
  } else {
    weight <- weight[1:nt]
    weight <- abs(weight)
  }
  if(!better.c[1] %in% c(0,1)){
    better.c <- FALSE
  } else {better.c <- better.c[1]}
  if(!is.numeric(direction) | length(direction) < nt){
    direction <- rep(Inf, nt)
  } else {direction <- direction[1:nt]}
  if(!outcross[1] %in% c(0,1)){
    outcross <- FALSE
  } else {outcross <- outcross[1]}
  if(!console[1] %in% c(0,1)){
    console <- TRUE
  } else {console <- console[1]}

  nt <- ncol(phe.t)
  nf1 <- choose(npl,2)
  if(is.null(nsele) | !is.numeric(npl)){
    nsele <- nf1
  } else {nsele <- nsele[1]}

  markertest <- c(nrow(marker) != ncol(geno.t), NA%in%marker[,2], marker[,1] != sort(marker[,1]))
  datatry <- try(marker[,2]%*%marker[,2], silent = TRUE)
  if(class(datatry)[1] == "try-error" | T%in%markertest){
    stop("Marker data error, please cheak your marker data. Or the number of marker does not match the genetype data.", call. = FALSE)
  }

  datatry <- try(npl*weight*direction*nprog*nsele*ngen*nrep, silent = TRUE)
  if(class(datatry)[1] == "try-error" | NA%in%datatry){
    stop("Argument error, please cheak your argument.", call. = FALSE)
  }

  if(nprog*nf1 < nsele){
    stop("Argument error, 'nprog' too small or 'nsele' too big.", call. = FALSE)
  }

  if(outcross){
    if(!is.null(fittedD.t)){
      datatry <- try(fittedD.t*fittedD.t, silent=TRUE)
      if(class(datatry)[1] == "try-error" | NA %in% fittedD.t | nrow(fittedD.t) != nrow(fittedA.t) | ncol(fittedD.t) != ncol(fittedA.t)){
        stop("fittedD.t error, please cheak your data and impute the missing value.", call. = FALSE)
      }
    }

    if(!is.null(fittedmu.t)){
      datatry <- try(fittedmu.t*fittedmu.t, silent=TRUE)
      if(class(datatry)[1] == "try-error" | NA %in% fittedmu.t | length(fittedmu.t) != ncol(fittedA.t)){
        stop("fittedmu.t error, please cheak your data and fix.", call. = FALSE)
      }
    }

    pheA <- phe.sd(fittedA.t)
    pheD <- phe.sd(fittedD.t)
    phe1 <- phe.sd(fittedA.t+fittedD.t)
    sd0 <- pheA[[3]]
    sdD <- pheD[[3]]
    mu0 <- fittedmu.t+pheA[[2]]+pheD[[2]]
    fit <- pheA[[1]]
    fitD <- pheD[[1]]

    gd.t <- geno.d(geno.t)
    K0 <- gd.t$KA
    gd.A <- gd.t$genoA
    KD <- gd.t$KD
    gd.D <- gd.t$genoD

    d0 <- det(K0)
    i0 <- 10^-10
    while(d0 == 0){
      i0 <- i0*10
      d0 <- det(K0+diag(i0, nrow(K0)))
    }
    K00 <- solve(K0+diag(i0, nrow(K0)))

    d0 <- det(KD)
    i0 <- 10^-10
    while(d0 == 0){
      i0 <- i0*10
      d0 <- det(KD+diag(i0, nrow(K0)))
    }
    KD0 <- solve(KD+diag(i0, nrow(K0)))
  } else {
    phe1 <- phe.sd(phe.t)
    mu0 <- phe1[[2]]
    sd0 <- phe1[[3]]
    fit <- phe1[[1]]

    sdD <- 0
    fitD <- fit
    fitD[] <- 0

    K0 <- geno.t%*%t(geno.t)/ncol(geno.t)

    d0 <- det(K0)
    i0 <- 10^-10
    while(d0 == 0){
      i0 <- i0*10
      d0 <- det(K0+diag(i0, nrow(K0)))
    }
    K00 <- solve(K0+diag(i0, nrow(K0)))

    KD0 <- K00
    KD0[] <- 0
  }
  sds <- phe1[[3]]

  row.names(fit) <- 1:nrow(fittedA.t)
  row.names(fitD) <- 1:nrow(fittedA.t)

  if(is.null(geno.c)){
    geno.c <- geno.t
    p.c <- fit
    p.cD <- fitD
  } else {
    datatry <- try(geno.c%*%t(geno.c), silent = TRUE)
    if(class(datatry)[1] == "try-error"){
      stop("Genotype data of candidate population error or have not been imputed.", call. = FALSE)
    }

    datatry <- try(geno.t%*%t(geno.c), silent=TRUE)
    if(class(datatry)[1] == "try-error"){
      stop("Candidate set genotype data error, please cheak your candidate set genotype data.", call. = FALSE)
    }

    geno.c2 <- rbind(geno.t,geno.c)

    if(outcross){
      geno.c20 <- geno.d(geno.c2)

      gd.A2 <- geno.c20$genoA
      KptpA <- t(gd.A%*%t(gd.A2))/ncol(gd.A2)
      gd.D2 <- geno.c20$genoD
      KptpD <- t(gd.D%*%t(gd.D2))/ncol(gd.D2)
    } else {
      KptpA <- t(geno.t%*%t(geno.c2))/ncol(geno.c2)
      KptpD <- KptpA
      KptpD[] <- 0
    }
    p.c <- KptpA%*%K00%*%fit
    p.cD <- KptpD%*%KD0%*%fitD
    npc <- nrow(p.c)-ind.t
    p.c <- cbind(matrix(0, npc, ind.t), diag(npc))%*%p.c
    p.cD <- cbind(matrix(0, npc, ind.t), diag(npc))%*%p.cD
  }

  if(length(geno.c[geno.c != 1 & geno.c != -1]) > 0){
    stop("Genotype data of candidate population error, the candidate population should be pure lines.", call. = FALSE)
  }

  if(is.null(row.names(geno.c))){
    row.names(p.c) <- 1:nrow(p.c)
    row.names(p.cD) <- 1:nrow(p.cD)
    row.names(geno.c) <- row.names(p.c)
  } else {
    row.names(p.c) <- row.names(geno.c)
    row.names(p.cD) <- row.names(geno.c)
  }

  sele <- c()
  sele.c <- c()
  for(i in 1:nt){
    if(direction[i] == Inf){
      sele0 <- (p.c[,i]+p.cD[,i])*weight[i]
    } else if(direction[i] == -Inf){
      sele0 <- (p.c[,i]+p.cD[,i])*weight[i]
      sele0 <- -sele0
    } else if(abs(direction[i]) != Inf){
      sele0 <- (p.c[,i]+p.cD[,i]+(mu0[i]-direction[i])/sds[i])*weight[i]
      sele0 <- -abs(sele0)
    }
    sele1 <- rep(1, nrow(p.c))
    sele1[sele0 < mean(sele0)] <- 0
    sele.c <- cbind(sele.c, sele1)
    sele <- cbind(sele, sele0)
  }
  sele <- apply(sele, 1, sum)

  if(better.c){
    sele.c <- apply(sele.c, 1, prod)
    fit.c.better <- p.c
    fit.c.better <- row.names(fit.c.better)[sele.c == 1]
    if(length(fit.c.better) >= 2*npl){
      geno.c <- geno.c[row.names(geno.c)%in%fit.c.better,]
      better.m0 <- matrix(0, length(fit.c.better), nrow(p.c))
      for(i in 1:length(fit.c.better)){
        better.m0[i, fit.c.better[i]==rownames(p.c)] <- 1
      }
      p.c <- better.m0%*%p.c
      p.cD <- better.m0%*%p.cD
      rownames(p.c) <- fit.c.better
      rownames(p.cD) <- fit.c.better
      sele.k <- sele[names(sele) %in% fit.c.better]
    } else {
      warning("The intersection of top individuals from each trait is too small. The function takes all candidate set to carry out the simulation.")
      sele.k <- sele
    }
  } else {sele.k <- sele}

  if(outcross){
    gd.c0 <- geno.d(geno.c)
    kp0 <- gd.c0$KA
    gd.c0A <- gd.c0$genoA
    kpD <- gd.c0$KD
    gd.c0D <- gd.c0$genoD
  } else {
    kp0 <- geno.c%*%t(geno.c)/ncol(geno.c)

    kpD <- kp0
    kpD[] <- 0
  }

  nGA <- npl
  if(npl > 20){nGA <- 20}
  mut <- 3
  if(npl <= 5){mut <- 1}

  pk <- order(sele.k, decreasing = TRUE)[1:npl.best]

  gvalue.result <- list()
  p.result <- matrix(0, nrep, npl)
  Dscore <- c()

  for(m in 1:nrep){
    nf1 <- choose(npl, 2)
    GA0 <- GA.Dscore(kp0, npl, keep = pk, n0 = nGA, mut = mut ,cri = cri)
    p0 <- GA0[[1]]
    Dscore[m] <- GA0[[2]]

    phe.m0 <- matrix(0, length(p0), nrow(p.c))
    for(i in 1:length(p0)){
      phe.m0[i, p0[i]] <- 1
    }
    phe.p0 <- phe.m0%*%p.c
    phe.p0D <- phe.m0%*%p.cD
    rownames(phe.p0) <- rownames(p.c)[p0]
    rownames(phe.p0D) <- rownames(p.c)[p0]
    geno.p0 <- geno.c[p0,]

    GEBVGD.p0 <- cbind(1:nrow(phe.p0), phe.p0, geno.p0)

    GEBVGD.gvalue <- list()
    GEBVGD.SNP <- list()

    p.result[m,] <- rownames(GEBVGD.p0)

    GEBV0 <- phe.p0
    GEBV0D <- phe.p0D
    if(is.null(colnames(phe.t))){
      colnames(GEBV0)<-paste("t", 1:nt, sep = "")
      colnames(GEBV0D)<-paste("t", 1:nt, sep = "")
    } else {
      colnames(GEBV0) <- colnames(phe.t)
      colnames(GEBV0D) <- colnames(phe.t)
    }
    GEBVGD.gvalue[[1]] <- GEBV0*matrix(sd0, nrow(GEBVGD.p0), nt, byrow = TRUE)+
      GEBV0D*matrix(sdD, nrow(GEBVGD.p0), nt, byrow = TRUE)+matrix(mu0, nrow(GEBVGD.p0), nt, byrow = TRUE)
    GEBVGD.p <- GEBVGD.p0[,(nt+2):ncol(GEBVGD.p0)]

    GEBVGD.SNP[[1]] <- GEBVGD.p
    GEBVGD.F1 <- list()
    k1 <- 1
    for(i in 1:(npl-1)){
      for(j in (i+1):npl){
        GEBVGD.F1[[k1]] <- as.matrix(GEBVGD.p[c(i,j),])
        k1 <- k1+1
      }
    }

    GEBVGD.F1.SNP <- c()
    nf1 <- choose(npl, 2)
    for(i in 1:nf1){
      F1_0 <- (GEBVGD.F1[[i]][1,]+GEBVGD.F1[[i]][2,])/2
      GEBVGD.F1.SNP <- cbind(GEBVGD.F1.SNP,F1_0)
    }
    GEBVGD.SNP[[2]] <- GEBVGD.F1.SNP
    GEBVGD.F1.SNP2 <- rbind(geno.t,t(GEBVGD.F1.SNP))

    if(outcross){
      gd.F1.SNP2 <- geno.d(GEBVGD.F1.SNP2)
      gd.F1.SNP2A <- gd.F1.SNP2$genoA
      Kpt <- t(gd.A%*%t(gd.F1.SNP2A))/ncol(gd.F1.SNP2A)
      gd.F1.SNP2D <- gd.F1.SNP2$genoD
      KptD <- t(gd.D%*%t(gd.F1.SNP2D))/ncol(gd.F1.SNP2D)
    } else {
      Kpt <- t(geno.t%*%t(GEBVGD.F1.SNP2))/ncol(GEBVGD.F1.SNP2)
      KptD <- Kpt
      KptD[] <- 0
    }

    fity.F1 <- Kpt%*%K00%*%fit
    fity.F1D <- KptD%*%KD0%*%fitD
    fity.F1 <- cbind(matrix(0, nf1, ind.t), diag(nf1))%*%fity.F1
    fity.F1D <- cbind(matrix(0, nf1, ind.t), diag(nf1))%*%fity.F1D
    row.names(fity.F1) <- NULL
    row.names(fity.F1D) <- NULL
    if(is.null(colnames(phe.t))){
      colnames(fity.F1) <- paste("t",1:nt,sep="")
      colnames(fity.F1D) <- paste("t",1:nt,sep="")
    } else {
      colnames(fity.F1) <- colnames(phe.t)
      colnames(fity.F1D) <- colnames(phe.t)
    }
    GEBVGD.gvalue[[2]] <- fity.F1*matrix(sd0, nrow(fity.F1), nt, byrow = TRUE)+
      fity.F1D*matrix(sdD, nrow(fity.F1), nt, byrow = TRUE)+matrix(mu0, nrow(fity.F1), nt, byrow = TRUE)

    if(console){
      cat("Method", "Repeat", "Generation", "\n")
      cat("GEBVGD", m, paste("F", 1, sep = ""), "\n", sep = "\t")
    }

    p.snp <- GEBVGD.F1
    g0 <- 3

    if(outcross){
      sele <- c()
      for(i in 1:nt){
        if(direction[i] == Inf){
          sele0 <- (fity.F1[,i]+fity.F1D[,i])*weight[i]
        } else if (direction[i] == -Inf){
          sele0 <- (fity.F1[,i]+fity.F1D[,i])*weight[i]
          sele0 <- -sele0
        } else {
          sele0 <- (fity.F1[,i]+fity.F1D[,i]+(mu0[i]-direction[i])/sds[i])*weight[i]
          sele0 <- -abs(sele0)
        }
        sele <- cbind(sele, sele0)
      }
      sele <- order(apply(sele, 1, sum), decreasing = TRUE)

      F1.snp <- list()
      for(i in 1:npl){
        F1.snp[[i]] <- p.snp[[sele[i]]]
      }

      GEBVGD.F2 <- list()
      k1 <- 1
      for(i in 1:(npl-1)){
        for(j in (i+1):npl){
          for(k in 1:nprog){
            marker.ind1 <- cbind(marker, t(F1.snp[[i]]))
            marker.ind2 <- cbind(marker, t(F1.snp[[j]]))
            F2_1 <- simu.gamete(marker.ind1)
            F2_2 <- simu.gamete(marker.ind2)
            F2 <- cbind(F2_1, F2_2)
            GEBVGD.F2[[k1]] <- t(F2)
            k1 <- k1+1
          }
        }
      }

      GEBVGD.F2.SNP <- c()
      for(i in 1:(nf1*nprog)){
        F2_0 <- (GEBVGD.F2[[i]][1,]+GEBVGD.F2[[i]][2,])/2
        GEBVGD.F2.SNP <- cbind(GEBVGD.F2.SNP, F2_0)
      }

      GEBVGD.SNP[[3]] <- GEBVGD.F2.SNP
      GEBVGD.F2.SNP2 <- rbind(geno.t, t(GEBVGD.F2.SNP))
      gd.F2.SNP2 <- geno.d(GEBVGD.F2.SNP2)
      gd.F2.SNP2A <- gd.F2.SNP2$genoA
      Kpt <- t(gd.A%*%t(gd.F2.SNP2A))/ncol(gd.F2.SNP2A)
      gd.F2.SNP2D <- gd.F2.SNP2$genoD
      KptD <- t(gd.D%*%t(gd.F2.SNP2D))/ncol(gd.F2.SNP2D)

      fity.F2 <- Kpt%*%K00%*%fit
      fity.F2D <- KptD%*%KD0%*%fitD
      nf2 <- nrow(fity.F2)-ind.t
      fity.F2 <- cbind(matrix(0, nf2, ind.t), diag(nf2))%*%fity.F2
      fity.F2D <- cbind(matrix(0, nf2, ind.t), diag(nf2))%*%fity.F2D
      row.names(fity.F2) <- NULL
      row.names(fity.F2D) <- NULL
      if(is.null(colnames(phe.t))){
        colnames(fity.F2) <- paste("t",1:nt,sep="")
        colnames(fity.F2D) <- paste("t",1:nt,sep="")
      } else {
        colnames(fity.F2) <- colnames(phe.t)
        colnames(fity.F2D) <- colnames(phe.t)
      }

      GEBVGD.gvalue[[3]] <- fity.F2*matrix(sd0, nrow(fity.F2), nt, byrow = TRUE)+
        fity.F2D*matrix(sdD, nrow(fity.F2), nt, byrow = TRUE)+matrix(mu0, nrow(fity.F2), nt, byrow = TRUE)

      if(console){cat("GEBVGD", m, paste("F", 2, sep = ""), "\n", sep = "\t")}

      sele <- c()
      for(i in 1:nt){
        if(direction[i] == Inf){
          sele0 <- (fity.F2[,i]+fity.F2D[,i])*weight[i]
        } else if (direction[i] == -Inf){
          sele0 <- (fity.F2[,i]+fity.F2D[,i])*weight[i]
          sele0 <- -sele0
        } else {
          sele0 <- (fity.F2[,i]+fity.F2D[,i]+(mu0[i]-direction[i])/sds[i])*weight[i]
          sele0 <- -abs(sele0)
        }
        sele <- cbind(sele, sele0)
      }
      sele <- order(apply(sele, 1, sum), decreasing = TRUE)
      p.snp <- list()
      nf1 <- nsele
      for(i in 1:nf1){
        p.snp[[i]] <- GEBVGD.F2[[sele[i]]]
      }
      g0 <- 4
    }


    for(g in g0:(ngen+1)){
      GEBVGD.F2 <- list()
      k2 <- 1
      for(i in 1:nf1){
        for(j in 1:nprog){
          marker.ind <- cbind(marker,t(p.snp[[i]]))
          F2_1 <- simu.gamete(marker.ind)
          F2_2 <- simu.gamete(marker.ind)
          F2 <- cbind(F2_1, F2_2)
          GEBVGD.F2[[k2]] <- t(F2)
          k2 <- k2+1
        }
      }
      GEBVGD.F2.SNP <- matrix(0,length(GEBVGD.F2[[1]][1,]), (nf1*nprog))
      for(i in 1:(nf1*nprog)){
        F2_0 <- (GEBVGD.F2[[i]][1,]+GEBVGD.F2[[i]][2,])/2
        GEBVGD.F2.SNP[,i] <- F2_0
      }
      GEBVGD.SNP[[g]] <- GEBVGD.F2.SNP
      GEBVGD.F2.SNP2 <- rbind(geno.t, t(GEBVGD.F2.SNP))

      if(outcross){
        gd.F2.SNP2 <- geno.d(GEBVGD.F2.SNP2)
        gd.F2.SNP2A <- gd.F2.SNP2$genoA
        Kpt2 <- t(gd.A%*%t(gd.F2.SNP2A))/ncol(gd.F2.SNP2A)
        gd.F2.SNP2D <- gd.F2.SNP2$genoD
        Kpt2D <- t(gd.D%*%t(gd.F2.SNP2D))/ncol(gd.F2.SNP2D)
      } else {
        Kpt2 <- t(geno.t%*%t(GEBVGD.F2.SNP2))/ncol(GEBVGD.F2.SNP2)
        Kpt2D <- Kpt2
        Kpt2D[] <- 0
      }

      fity.F2 <- Kpt2%*%K00%*%fit
      fity.F2D <- Kpt2D%*%KD0%*%fitD
      nf2 <- nrow(fity.F2)-ind.t
      fity.F2 <- cbind(matrix(0, nf2, ind.t), diag(nf2))%*%fity.F2
      fity.F2D <- cbind(matrix(0, nf2, ind.t), diag(nf2))%*%fity.F2D
      row.names(fity.F2) <- NULL
      if(is.null(colnames(phe.t))){
        colnames(fity.F2) <- paste("t", 1:nt, sep = "")
        colnames(fity.F2D) <- paste("t", 1:nt, sep = "")
      } else {
        colnames(fity.F2) <- colnames(phe.t)
        colnames(fity.F2D) <- colnames(phe.t)
      }
      GEBVGD.gvalue[[g]] <- fity.F2*matrix(sd0, nrow(fity.F2), nt, byrow = TRUE)+
        fity.F2D*matrix(sdD, nrow(fity.F2), nt, byrow = TRUE)+matrix(mu0, nrow(fity.F2), nt, byrow = TRUE)

      if(console){cat("GEBVGD", m, paste("F", g-1, sep = ""), "\n", sep = "\t")}

      sele <- c()
      for(i in 1:nt){
        if(direction[i] == Inf){
          sele0 <- (fity.F2[,i]+fity.F2D[,i])*weight[i]
        } else if (direction[i] == -Inf){
          sele0 <- (fity.F2[,i]+fity.F2D[,i])*weight[i]
          sele0 <- -sele0
        } else {
          sele0 <- (fity.F2[,i]+fity.F2D[,i]+(mu0[i]-direction[i])/sds[i])*weight[i]
          sele0 <- -abs(sele0)
        }
        sele <- cbind(sele, sele0)
      }
      sele <- order(apply(sele, 1, sum), decreasing = TRUE)
      p.snp <- list()
      nf1 <- nsele
      for(i in 1:nf1){
        p.snp[[i]] <- GEBVGD.F2[[sele[i]]]
      }
    }
    names(GEBVGD.gvalue) <- c("p", paste("F", 1:ngen, sep = ""))
    gvalue.result[[m]] <- GEBVGD.gvalue
  }

  bestsub <- table(p.result)
  bestsub <- data.frame(parental.lines = names(bestsub), chosen.ratio = as.numeric(bestsub)/nrep)
  bestsub <- bestsub[order(bestsub[,2], decreasing = TRUE),]
  bestsub <- bestsub[1:npl,]
  parental.lines <- list(parental.lines = p.result, D.score = Dscore)

  return(list(method = "GEBV-GD", weight = weight, direction = direction, mu = phe1[[2]], sd = sd0,
              GEBV.value = gvalue.result, parental.lines = parental.lines, suggested.subset = bestsub))
}
