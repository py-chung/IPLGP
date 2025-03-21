methEM <- function(QTL, marker, geno, D.matrix, y, yu, tL, tR, type, ng, sele.g, crit, cM, ya){
  EM <- EM.MIM(QTL, marker, geno, D.matrix, y = y, yu = yu, tL = tL, tR = tR, type = type,
                ng = ng, sele.g = sele.g, crit = crit, console = FALSE)
  eff <- as.numeric(EM$E.vector)
  mu0 <- as.numeric(EM$beta)
  sigma <- sqrt(as.numeric(EM$variance))
  LRT <- EM$LRT
  model <- EM$model
  R2 <- EM$R2
  like <- EM$log.likelihood
  result <- list(eff, mu0, sigma, LRT, like, R2, model)
  return(result)
}
