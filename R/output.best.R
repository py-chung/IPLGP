#' Summary For The Best Individuals
#'
#' Output the GEBV average curves and the summary statistics for the best
#' individuals selected over generations.
#'
#' @param result list. The data list of the output from simu.GEBVO, simu.GDO,
#' or simu.GEBVGD.
#' @param save.pdf logical. A logical variable, if save.pdf is set to be TRUE,
#' the pdf file of plots will be saved in the working directory instead of
#' being shown in the console.
#'
#' @return
#' The GEBV averages of the best individuals among the repetitions over
#' generations for each trait.
#'
#' @note
#' The figure output contains the plots of GEBV averages of the best individuals
#' selected over generations for each trait. If save.pdf is set to be TRUE, the
#' pdf file of plots will be saved in the working directory instead of being
#' shown in the console.
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
#' \code{\link[ggplot2]{ggplot}}
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
#' # summary for the best individuals
#' output <- output.best(result)
#' output
output.best <- function(result, save.pdf =  FALSE){

  if(!save.pdf[1] %in% c(0, 1) | length(save.pdf) > 1){save.pdf <- FALSE}
  datatest <- names(result) != c("method", "weight", "direction", "mu", "sd", "GEBV.value",
                                 "parental.lines", "suggested.subset")
  if(TRUE %in% (datatest) | length(datatest)!=8){
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

  datatry <- try(weight*direction*mu*sd, silent =  TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% datatry){
    stop("Input data error, please input the original output data of simu.GEBVO, simu.GDO, or simu.GEBVGD.",
         call. = FALSE)
  }

  GEBV.max <- list()
  for(i in 1:nrep){
    max.gvalue0 <- matrix(0, ngen+1, nt)
    for(j in 1:(ngen+1)){
      max.gvalue <- GEBV[[i]][[j]]
      datatry <- try(max.gvalue*max.gvalue, silent =  TRUE)
      if(class(datatry)[1] == "try-error" | NA %in% max.gvalue){
        stop("Input data error, please input the original output data of simu.GEBVO, simu.GDO, or simu.GEBVGD.",
             call. = FALSE)
      }

      sele <- c()
      for(k in 1:nt){
        if(direction[k] == Inf){
          sele0 <- max.gvalue[,k]/sd[k]*weight[k]
        } else if (direction[k] == -Inf){
          sele0 <- max.gvalue[,k]/sd[k]*weight[k]
          sele0 <- -sele0
        } else {
          sele0 <- ((max.gvalue[,k]-direction[k])/sd[k])*weight[k]
          sele0 <- -abs(sele0)
        }
        sele <- cbind(sele, sele0)
      }
      sele <- apply(sele, 1, sum)

      max.gvalue0[j,] <- max.gvalue[which.max(sele),]
    }
    GEBV.max[[i]] <- max.gvalue0
  }

  GEBV.max.ave <- list()
  for(i in 1:nt){
    ave0 <- matrix(0, (ngen+1), nrep)
    for(j in 1:nrep){
      ave0[, j] <- GEBV.max[[j]][, i]
    }
    GEBV.max.ave[[i]] <- ave0
  }

  GEBV.all <- list()

  if(save.pdf){grDevices::pdf("IPLGP.GEBVplot.pdf", width = 8, height = 5)}
  for(i in 1:nt){
    GEBV.all0 <- data.frame(
      generation <- 0:ngen,
      best.GEBV.average <- apply(GEBV.max.ave[[i]], 1, mean),
      GEBV.sd <- apply(GEBV.max.ave[[i]], 1, sd)
    )
    colnames(GEBV.all0) =  c("generation", "mean", "standard deviation")

    plot0 <- ggplot2::ggplot(GEBV.all0, ggplot2::aes(x =  generation, y =  best.GEBV.average), main =  "GEBV.all") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin =  best.GEBV.average-GEBV.sd, ymax =  best.GEBV.average+GEBV.sd),
                             width =  .1) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size =  1.5)+
      ggplot2::labs(title = paste(t.n[i], " (", method, ")", sep =  ""))+
      ggplot2::scale_x_continuous(breaks =  seq(0, ngen, 1), labels = c("P", paste("F", 1:ngen, sep = "")))+
      ggplot2::theme_bw()

    print(plot0)
    GEBV.all0[, 1] <- c("P", paste("F", 1:ngen, sep = ""))
    GEBV.all[[i]] <- GEBV.all0
  }
  if(save.pdf){grDevices::dev.off()}

  names(GEBV.all) <- t.n
  return(GEBV.all)
}
