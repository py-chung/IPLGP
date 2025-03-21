methSG <- function(QTL, marker, geno, D.matrix, y, yu, tL, tR, type, ng,sele.g, crit, cM, ya){
  mp <- switch(sele.g,
               p = mixprop(QTL, length(yu), marker, geno, model = 1, cM = cM, type = type, ng = ng)[[2]],
               t = mixprop(QTL, length(yu), marker, geno, model = 2, cM = cM, type = type, ng = ng)[[2]],
               f = mixprop(QTL, length(yu), marker, geno, model = 3, cM = cM, type = type, ng = ng)[[2]],
               n = Q.make(QTL, marker, geno, type = type, ng = ng)$cp.matrix)
  X <- mp%*%D.matrix
  fit <- stats::lm(ya~X)
  eff <- as.numeric(fit$coefficients[-1])
  mu0 <- as.numeric(fit$coefficients[1])
  ms <- stats::anova(fit)$`Mean Sq`
  sigma <- ms[2]^0.5
  R2 <- summary(fit)$r.squared

  L0 <- c()
  L1 <- c()
  for(k in 1:nrow(mp)){
    L00 <- c()
    L01 <- c()
    for(m in 1:nrow(D.matrix)){
      L00[m] <- mp[k, m]*stats::dnorm((ya[k]-mu0)/sigma)
      L01[m] <- mp[k, m]*stats::dnorm((ya[k]-(mu0+D.matrix[m,]%*%eff))/sigma)
    }
    L0[k] <- sum(L00)
    L1[k] <- sum(L01)
  }
  LRT <- -2*sum(log(L0[!is.na(L0) & !is.na(L1)]/L1[!is.na(L0) & !is.na(L1)]))
  like <- sum(log(L1[!is.na(L0) & !is.na(L1)]))

  result <- list(eff, mu0, sigma, LRT, like, R2, model = "regression interval mapping model")
  return(result)
}
