EM.MIM0 <- function(D.matrix, cp.matrix, y, E.vector0 = NULL, X = NULL, beta0 = NULL,
                   variance0 = NULL, crit = 10^-5, stop = 1000, conv = TRUE, console = TRUE){

  Y <- c(y)
  ind <- length(Y)
  g <- nrow(D.matrix)
  eff <- ncol(D.matrix)

  Y[is.na(Y)] <- mean(Y,na.rm = TRUE)

  E.vector <- E.vector0
  beta <- beta0
  variance <- variance0

  sigma <- sqrt(variance)

  Delta <- 1
  number <- 0

  if(length(colnames(D.matrix)) == ncol(D.matrix)){
    effectname <- colnames(D.matrix)
  }

  cat(paste("number", "var", paste(effectname, collapse = "\t"), "\n", sep = "\t")[console])

  Yt <- as.matrix(Y)
  indvec <- matrix(1, 1, ind)
  gvec <- matrix(1, 1, g)

  while (max(abs(Delta)) > crit & number < stop) {

    Et <- as.matrix(E.vector)
    bt<- as.matrix(beta)

    muji.matrix <- t(D.matrix%*%E.vector%*%indvec)+X%*%beta%*%gvec

    P0.matrix <- cp.matrix*stats::dnorm(y, muji.matrix, c(sigma))
    P0s <- apply(P0.matrix, 1, sum)
    P0.matrix[P0s == 0, ] <- rep(1, g)
    PIt <- P0.matrix*(1/P0s%*%gvec)

    PD <- t(Yt-X%*%bt)%*%PIt%*%D.matrix
    PDD <- indvec%*%PIt%*%(D.matrix^2)
    r.vector <- c(PD/PDD)

    M.matrix <- matrix(0, eff, eff)
    V.matrix <- matrix(0, eff, eff)
    for(i in 1:eff){
      Di <- D.matrix[, i]
      Dij <- D.matrix
      for(j in 1:eff){
        Dij[, j] <- Di*D.matrix[, j]
      }
      iPI <- indvec%*%PIt
      M01 <- c(iPI)%*%Dij
      V.matrix[i, ] <- c(M01)
      M.matrix[, i] <- c(M01)/c(PDD)
    }
    diag(M.matrix) <- 0

    E.t <- r.vector-M.matrix%*%E.vector
    beta.t <- solve(t(X)%*%X)%*%t(X)%*%(c(Yt)-PIt%*%D.matrix%*%E.t)

    YXb <- Y-X%*%c(beta.t)
    EVE <- t(E.t)%*%V.matrix%*%E.t
    sigma.t <- sqrt((t(YXb)%*%(YXb)-t(YXb)%*%PIt%*%D.matrix%*%E.t*2+EVE)/ind)

    Delta <- c(E.t-E.vector)
    if(NaN %in% Delta){
      break()
    }
    number <- number+1
    Ep <- round(c(E.t), 3)
    sp <- round(c(sigma.t)^2, 3)
    cat(paste(number, sp, paste(Ep, collapse = "\t"), "\n", sep = "\t")[console])

    E.vector <- c(E.t)
    beta <- c(beta.t)
    sigma <- c(sigma.t)

  }

  muji.matrix <- t(D.matrix%*%E.vector%*%indvec)+X%*%beta%*%gvec
  P0.matrix <- cp.matrix*stats::dnorm(y, muji.matrix, c(sigma))
  P0s <- apply(P0.matrix, 1, sum)
  P0.matrix[P0s == 0, ] <- rep(1, g)
  PI.matrix <- P0.matrix*(1/P0s%*%gvec)

  names(E.vector) <- effectname
  colnames(PI.matrix) <- colnames(cp.matrix)

  variance <- sigma^2

  L0 <- rep(0, ind)
  L1 <- rep(0, ind)
  Xb <- X%*%beta
  for(m in 1:g){
    L0 <- L0+(cp.matrix[, m]*stats::dnorm(Y, mean(Xb), sigma))
    L1 <- L1+(cp.matrix[, m]*stats::dnorm(Y, mean(Xb)+D.matrix[m,]%*%E.vector, sigma))
  }
  like0 <- sum(log(L0))
  like1 <- sum(log(L1))
  LRT <- 2*(like1-like0)

  y.hat <- PI.matrix%*%D.matrix%*%E.vector+Xb
  r2 <- c(stats::var(y.hat)/stats::var(y))

  if(number == stop){
    if(conv){
      E.vector <- rep(0, length(E.vector))
      beta <- 0
      variance <- 0
      PI.matrix <- matrix(0, nrow(PI.matrix), ncol(PI.matrix))
      like1 <- -Inf
      LRT <- 0
      r2 <- 0
    }
    warning("EM algorithm fails to converge, please check the input data or adjust
            the convergence criterion and stopping criterion.")
  }

  result <- list(E.vector = E.vector, beta = as.numeric(beta), variance = as.numeric(variance),
                 PI.matrix = PI.matrix, log.likelihood = like1, LRT = LRT, R2 = r2,
                 y.hat = y.hat, iteration.number = number)

  return(result)
}
