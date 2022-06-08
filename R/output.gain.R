#' Summary For Genetic Gain
#'
#' Output the GEBV average of parental lines, the GEBV average of the last
#' generation in simulation process, and the genetic gain average over repetitions
#' for each target trait.
#'
#' @param result list. The data list of the output from simu.GEBVO, simu.GDO,
#' or simu.GEBVGD.
#'
#' @return
#' The output contains the table of the GEBV average of parental lines, the GEBV
#' average of the last generation in simulation process, and the genetic gain
#' average over repetitions for each target trait.
#'
#' @export
#'
#' @references
#'
#' Chung PY, Liao CT. 2020. Identification of superior parental lines for
#' biparental crossing via genomic prediction. PLoS ONE 15(12):e0243159.
#'
#' @seealso
#' \code{\link[IPLGP]{simu.GEBVO}}
#' \code{\link[IPLGP]{simu.GDO}}
#' \code{\link[IPLGP]{simu.GEBVGD}}
#'
#' @examples
#' # generate simulated data
#' set.seed(2000)
#' t1 <- rnorm(10,30,10)
#' t2 <- rnorm(10,10,5)
#' t3 <- NULL
#' t4 <- NULL
#' t5 <- NULL
#' geno.test <- matrix(sample(c(1, -1), 200, replace = TRUE), 10, 20)
#' marker.test <- cbind(rep(1:2, each=10), rep(seq(0, 90, 10), 2))
#' fit <- GBLUP.fit(t1, t2, t3, t4, t5, geno = geno.test)
#' fitvalue <- fit$fitted.value
#'
#' geno.candidate <- matrix(sample(c(1,-1), 300, replace = TRUE), 15, 20)
#'
#' # run
#' result <- simu.GEBVO(fitvalue, geno.t = geno.test, marker = marker.test,
#' geno.c = geno.candidate, nprog = 5, nsele = 10, ngen = 5, nrep = 5)
#'
#' # summary for genetic gain
#' output <- output.gain(result)
#' output
output.gain <- function(result){

  datatest <- names(result) != c("method", "weight", "direction", "mu", "sd", "GEBV.value",
                                 "parental.lines", "suggested.subset")
  if(TRUE %in% (datatest) | length(datatest) != 8){
    stop("Input data error, please input the original output data of simu.GEBVO, simu.GDO, or simu.GEBVGD.",
         call. = FALSE)
  }

  GEBV <- result$GEBV.value
  t.n <- colnames(GEBV[[1]][[1]])
  method <- result$method
  weight <- result$weight
  direction <- result$direction
  mu <- result$mu
  sd <- result$sd
  nrep <- length(GEBV)
  ngen <- length(GEBV[[1]])-1
  nt <- ncol(GEBV[[1]][[1]])

  datatry <- try(weight*direction*mu*sd, silent = TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% datatry){
    stop("Input data error, please input the original output data of simu.GEBVO, simu.GDO, or simu.GEBVGD.",
         call. = FALSE)
  }

  GEBV.mean <- list()
  for(i in 1:nrep){
    mean.gvalue0 <- matrix(0, ngen+1, nt)
    for(j in 1:(ngen+1)){
      mean.gvalue <- GEBV[[i]][[j]]
      datatry <- try(mean.gvalue*mean.gvalue, silent = TRUE)
      if(class(datatry)[1] == "try-error" | NA %in% mean.gvalue){
        stop("Input data error, please input the original output data of simu.GEBVO, simu.GDO, or simu.GEBVGD.",
             call. = FALSE)
      }

      for(k in 1:nt){
        mean.gvalue0[j, k] <- mean(mean.gvalue[, k])
      }
    }
    GEBV.mean[[i]] <- mean.gvalue0
  }

  GEBV.mean.ave <- list()
  for(i in 1:nt){
    ave0 <- matrix(0, (ngen+1), nrep)
    for(j in 1:nrep){
      ave0[, j] <- GEBV.mean[[j]][, i]
    }
    GEBV.mean.ave[[i]] <- ave0
  }

  GEBV.all <- list()
  gain <- data.frame(matrix(NA, nt, 4))
  for(i in 1:nt){
    GEBV.all0 <- data.frame(
      generation <- 0:ngen,
      best.GEBV.average <- apply(GEBV.mean.ave[[i]], 1, mean),
      GEBV.sd <- apply(GEBV.mean.ave[[i]], 1, sd)
    )
    colnames(GEBV.all0) <- c("generation", "mean", "standard deviation")

    gain[i, c(2, 3)] <- GEBV.all0[c(1, (ngen+1)), 2]
    if(abs(direction[i]) == Inf){
      gain[i, 4] <- gain[i, 3]-gain[i, 2]
    } else {
      dir0 <- direction[i]
      g0 <- c()
      for(j in 1:nrep){
        g0[j] <- abs(dir0-GEBV.mean.ave[[i]][1, j])-abs(dir0-GEBV.mean.ave[[i]][ngen+1, j])
      }
      gain[i, 4] <- mean(g0)
    }

    GEBV.all0[, 1] <- c("P", paste("F", 1:ngen, sep = ""))
    GEBV.all[[i+1]] <- GEBV.all0
  }
  gain[, 1] <- t.n
  colnames(gain) <- c("trait", "parental line", paste("F", ngen, sep = ""), "genetic gain")
  GEBV.all[[1]] = gain

  names(GEBV.all) <- c("genetic gain", t.n)
  return(GEBV.all)
}
