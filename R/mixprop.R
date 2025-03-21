mixprop <- function(QTL, nu, marker, geno, model, cM = TRUE, type = "RI", ng = 2, cp.matrix = NULL){
  nQTL <- nrow(QTL)
  if(nQTL > 1){QTL <- QTL[order(QTL[, 1], QTL[, 2]),]}
  marker <- marker[order(marker[, 1], marker[, 2]),]

  Q3 <- 3^nQTL
  Q2 <- 2^nQTL
  if(is.null(cp.matrix)){
    cp.matrix <- Q.make(QTL, marker, geno, cM = cM, type = type, ng = ng)[[(nQTL+1)]]
  }
  QTL.freq <- cp.matrix

  if(cM){
    QTL[, 2] <- QTL[, 2]/100
    marker[, 2] <- marker[, 2]/100
  }

  if(model == 2){
    mix.prop <- QTL.freq
    popu.freq <- NULL
    freq.u <- NULL
  } else {
    N <- nrow(QTL.freq)+nu
    freq.s <- colSums(QTL.freq)/N

    if(nQTL == 1){
      if(type == "BC"){
        popu.freq <- c(1-0.5^ng, 0.5^ng)
        Qn <- 2
      } else if (type == "RI"){
        hetro <- 0.5^(ng-1)
        popu.freq <- c((1-hetro)/2, hetro, (1-hetro)/2)
        Qn <- 3
      } else {
        popu.freq <- c(0.25,0.5,0.25)
        Qn <- 3
      }

    } else {
      M <- c(1,0)
      GAM <- matrix(0, Q2, nQTL)
      Pg <- c(2, 1, 0)
      PG <- matrix(0, Q3, nQTL)
      rmn <- rep(0, nQTL-1)
      for(n0 in 1:nQTL){
        GAM[, n0] <- rep(rep(M, 2^(n0-1)), each = 2^(nQTL-n0))
        PG[, n0] <- rep(rep(Pg, 3^(n0-1)), each = 3^(nQTL-n0))
        if(n0 == nQTL){break}
        if(QTL[n0, 1] == QTL[(n0+1), 1]){
          rmd <- QTL[(n0+1), 2]-QTL[n0, 2]
          rmn[n0] <- (1-exp(-2*rmd))/2
        }else{
          rmn[n0] <- 0.5
        }
      }

      if(type == "AI" & ng > 2){
        k0 <- seq(0, (ng-1), 2)
        rmn0 <- matrix(0, (nQTL-1), length(k0))
        for(i in 1:length(k0)){
          rmn0[,i] <- choose((ng-1), k0[i])*rmn^k0[i]*(1-rmn)^(ng-1-k0[i])
        }
        rmn <- 1-apply(rmn0, 1, sum)
      }

      gamf <- matrix(0, Q2, nQTL-1)
      GAM.f <- matrix(0, Q2, nQTL-1)
      GAM.freq <- rep(1, Q2)
      for(nr in 1:(nQTL-1)){
        gamf[, nr] <- GAM[, nr] == GAM[, nr+1]
        for(ngam in 1:Q2){
          GAM.f[ngam,nr] <- ifelse(gamf[ngam, nr] == 1, 1-rmn[nr], rmn[nr])
        }
        GAM.freq <- GAM.freq*GAM.f[, nr]
      }

      if(type == "BC"){
        popu.freq <- GAM.freq/2
        if(ng > 1){
          G.matrix <- matrix(0, Q2, Q2)
          G.matrix[1, 1] <- 1
          G.matrix[Q2, ] <- GAM.freq/2
          for(i in 2:(Q2-1)){
            p0 <- apply(matrix(t(GAM[, which(GAM[i,] != 0)]) == GAM[i, which(GAM[i,]!=0)],
                               nrow = Q2, byrow = TRUE), 1, sum) == sum(GAM[i,] != 0)
            p1 <- popu.freq[p0]/sum(popu.freq[p0])
            G.matrix[i, p0] <- p1
          }
          for(i in 1:(ng-1)){
            popu.freq <- crossprod(G.matrix, popu.freq)
          }
        }
        Qn <- Q2
      }else{
        Gfreq <- rep(0, (Q2)*(Q2))
        Gtype <- matrix(0, (Q2)*(Q2), nQTL)
        wg <- 0
        for(ga1 in 1:Q2){
          for(ga2 in 1:Q2){
            wg <- wg+1
            Gfreq[wg] <- (GAM.freq[ga1]/2)*(GAM.freq[ga2]/2)
            Gtype[wg,] <- GAM[ga1,]+GAM[ga2,]
          }
        }
        popu.freq <- rep(0, Q3)
        ind <- rep(0, (Q2)*(Q2))
        for(pfi in 1:Q3){
          for(gti in 1:((Q2)*(Q2))){
            if(sum(Gtype[gti,] == PG[pfi,]) == nQTL){ind[gti] <- pfi}
          }
          popu.freq[pfi] <- sum(Gfreq[ind == pfi])
        }

        if(type == "RI" & ng > 2){
          G.matrix <- matrix(0, Q3, Q3)
          G.matrix[1, 1] <- 1
          G.matrix[Q3, Q3] <- 1
          for(i in 2:(Q3-1)){
            p0 <- apply(matrix(t(PG[, which(PG[i,] != 1)]) == PG[i, which(PG[i,] != 1)],
                               nrow = Q3, byrow = TRUE), 1, sum) == sum(PG[i,] != 1)
            p1 <- popu.freq[p0]/sum(popu.freq[p0])
            G.matrix[i, p0] <- p1
          }
          for(i in 1:(ng-1)){
            popu.freq <- crossprod(G.matrix, popu.freq)
          }
        }
        Qn <- Q3
      }
    }

    freq.u <- as.vector(popu.freq)-freq.s
    if(model == 3){
      freq.u <- as.vector(popu.freq)
      names(freq.u) <- names(freq.s)
    }
    freq.u[freq.u < 0] <- 0
    Freq.u <- matrix(rep(freq.u/sum(freq.u), each = nu), nu, Qn)
    mix.prop <- rbind(QTL.freq, Freq.u)
  }
  re <- list(popu.freq = popu.freq, mix.prop = mix.prop, QTL = QTL, marker = marker, freq.u = freq.u)
  return(re)
}
