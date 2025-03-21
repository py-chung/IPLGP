#' QTL search by IM
#'
#' Expectation-maximization algorithm for QTL interval mapping to search
#' for possible positions of QTL in all chromosomes. It can handle genotype
#' data which is selective genotyping too.
#'
#' @param marker matrix. A k*2 matrix contains the marker information,
#' where the row dimension 'k' represents the number of markers in the
#' chromosomes. The first column labels the chromosomes where the markers
#' are located, and the second column labels the positions of markers (in
#' morgan (M) or centimorgan (cM)). It's important to note that chromosomes
#' and positions must be sorted in order.
#' @param geno matrix. A n*k matrix contains the genotypes of k markers
#' for n individuals. The marker genotypes of P1 homozygote (MM),
#' heterozygote (Mm), and P2 homozygote (mm) are coded as 2, 1, and 0,
#' respectively, with NA indicating missing values.
#' @param y vector. A vector that contains the phenotype values of
#' individuals with genotypes.
#' @param yu vector. A vector that contains the phenotype values of
#' individuals without genotypes.
#' @param sele.g character. Determines the type of data being analyzed:
#' If sele.g="n", it considers the data as complete genotyping data. If
#' sele.g="f", it treats the data as selective genotyping data and utilizes
#' the proposed corrected frequency model (Lee 2014) for analysis; If
#' sele.g="t", it considers the data as selective genotyping data and uses
#' the truncated model (Lee 2014) for analysis; If sele.g="p", it treats
#' the data as selective genotyping data and uses the population
#' frequency-based model (Lee 2014) for analysis. Note that the 'yu'
#' argument must be provided when sele.g="f" or "p".
#' @param tL numeric. The lower truncation point of phenotype value when
#' sele.g="t". When sele.g="t" and tL=NULL, the 'yu' argument must be
#' provided. In this case, the function will consider the minimum of 'yu'
#' as the lower truncation point.
#' @param tR numeric. The upper truncation point of phenotype value when
#' sele.g="t". When sele.g="t" and tR=NULL, the 'yu' argument must be
#' provided. In this case, the function will consider the maximum of 'yu'
#' as the upper truncation point.
#' @param method character. When method="EM", it indicates that the interval
#' mapping method by Lander and Botstein (1989) is used in the analysis.
#' Conversely, when method="REG", it indicates that the approximate regression
#' interval mapping method by Haley and Knott (1992) is used in the analysis.
#' @param type character. The population type of the dataset. Includes
#' backcross (type="BC"), advanced intercross population (type="AI"), and
#' recombinant inbred population (type="RI"). The default value is "RI".
#' @param D.matrix matrix. The design matrix of the IM model. If D.matrix=NULL,
#' the design matrix will be constructed using Cockerham’s model: In BC
#' population, it is a 2*1 matrix with values 0.5 and -0.5 for the additive
#' effect; In RI or AI population, it is a 3*2 matrix. The first column
#' consists of 1, 0, and -1 for the additive effect, and the second column
#' consists of 0.5, -0.5, and 0.5 for the dominant effect.
#' @param ng integer. The generation number of the population type. For
#' instance, in a BC1 population where type="BC", ng=1; in an AI F3
#' population where type="AI", ng=3.
#' @param cM logical. Specify the unit of marker position. If cM=TRUE, it
#' denotes centimorgan; if cM=FALSE, it denotes morgan.
#' @param speed numeric. The walking speed of the QTL search (in cM).
#' @param crit numeric. The convergence criterion of EM algorithm.
#' The E and M steps will iterate until a convergence criterion is met.
#' It must be a value between 0 and 1.
#' @param d.eff logical. Specifies whether the dominant effect will be
#' considered in the parameter estimation for AI or RI population.
#' @param LRT.thre logical or numeric. If set to TRUE, the LRT threshold
#' will be computed based on the Gaussian stochastic process (Kao and Ho 2012).
#' Alternatively, users can input a numerical value as the LRT threshold to
#' evaluate the significance of QTL detection.
#' @param simu integer. Determines the number of simulation samples that
#' will be used to compute the LRT (Likelihood Ratio Test) threshold using
#' the Gaussian process. It must be a value between 50 and 10^8.
#' @param alpha numeric. The type I error rate for the LRT threshold.
#' @param detect logical. Determines whether the significant QTL, whose LRT
#' statistic is larger than the LRT threshold, will be displayed in the
#' output dataset or not.
#' @param QTLdist numeric. The minimum distance (in cM) among different
#' linked significant QTL.
#' @param plot.all logical. When set to TRUE, it directs the function to
#' output the profile of LRT statistics for the genome in one figure.
#' @param plot.chr logical. When set to TRUE, it instructs the function to
#' output the profile of LRT statistics for the chromosomes.
#' @param console logical. Determines whether the process of the algorithm
#' will be displayed in the R console or not.
#'
#' @return
#' \item{effect}{The estimated effects and LRT statistics of all positions.}
#' \item{LRT.threshold}{The LRT threshold value computed for the data using the
#' Gaussian stochastic process (Kuo 2011; Kao and Ho 2012).}
#' \item{detect.QTL}{The positions, effects and LRT statistics of the detected
#' QTL significant using the obtained LRT threshold value.}
#' \item{model}{The model of selective genotyping data in this analyze.}
#' \item{inputdata}{The input data of this analysis. It contains marker, geno,
#' y, yu, sele.g, type, ng, cM, and d.eff. The parameters not provided by the
#' user will be output with default values.}
#'
#' Graphical outputs including LOD value and effect of each position.
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
#' H.-I LEE, H.-A. HO and C.-H. KAO 2014 A new simple method for improving
#' QTL mapping under selective genotyping. Genetics 198: 1685-1698. <doi: 10.1534/genetics.114.168385.>
#'
#' KAO, C.-H. and H.-A. Ho 2012 A score-statistic approach for determining
#' threshold values in QTL mapping. Frontiers in Bioscience. E4, 2670-2682. <doi: 10.2741/e582>
#'
#' @seealso
#' \code{\link[QTLEMM]{EM.MIM}}
#' \code{\link[QTLEMM]{LRTthre}}
#'
#' @examples
#'
#' # load the example data
#' load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' result <- IM.search(marker, geno, y, type = "RI", ng = 2, speed = 7, crit = 10^-3, LRT.thre = 10)
#' result$detect.QTL
#'
#' \dontrun{
#' # Example for selective genotyping data
#' # load the example data
#' load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))
#'
#' # make the seletive genotyping data
#' ys <- y[y > quantile(y)[4] | y < quantile(y)[2]]
#' yu <- y[y >= quantile(y)[2] & y <= quantile(y)[4]]
#' geno.s <- geno[y > quantile(y)[4] | y < quantile(y)[2],]
#'
#' # run and result
#' result <- IM.search(marker, geno.s, ys, yu, sele.g = "f", type = "RI", ng = 2,
#' speed = 7, crit = 10^-3, LRT.thre = 10)
#' result$detect.QTL
#' }
#'
IM.search <- function(marker, geno, y, yu = NULL, sele.g = "n", tL = NULL, tR = NULL, method = "EM",
                      type = "RI", D.matrix = NULL, ng = 2, cM = TRUE, speed = 1, crit = 10^-5,
                      d.eff = FALSE, LRT.thre = TRUE, simu = 1000, alpha = 0.05, detect = TRUE,
                      QTLdist = 15, plot.all = TRUE, plot.chr = TRUE, console = TRUE){

  if(is.null(marker) | is.null(geno) | is.null(y)){
    stop("Input data is missing. The argument marker, geno, and y must be input.", call. = FALSE)
  }

  genotest <- table(c(geno))
  datatry <- try(geno*geno, silent=TRUE)
  if(class(datatry)[1] == "try-error" | FALSE %in% (names(genotest) %in% c(0, 1, 2))  | length(dim(geno)) != 2){
    stop("Genotype data error, please cheak your genotype data.", call. = FALSE)
  }

  marker <- as.matrix(marker)
  markertest <- c(ncol(marker) != 2, NA %in% marker, marker[,1] != sort(marker[,1]), nrow(marker) != ncol(geno))
  datatry <- try(marker*marker, silent=TRUE)
  if(class(datatry)[1] == "try-error" | TRUE %in% markertest){
    stop("Marker data error, or the number of marker does not match the genotype data.", call. = FALSE)
  }

  y[is.na(y)] <- mean(y,na.rm = TRUE)

  if(!is.null(yu)){
    yu <- yu[!is.na(yu)]
    datatry <- try(yu%*%yu, silent=TRUE)
    if(class(datatry)[1] == "try-error"){
      stop("yu data error, please check your yu data.", call. = FALSE)
    }
  }

  y <- c(y)
  datatry <- try(y%*%geno, silent=TRUE)
  if(class(datatry)[1] == "try-error"){
    stop("Phenotype data error, or the number of individual does not match the genetype data.", call. = FALSE)
  }

  if(!sele.g[1] %in% c("n", "t", "p", "f") | length(sele.g)!=1){
    stop("Parameter sele.g error, please check and fix.", call. = FALSE)
  }

  if(sele.g == "t"){
    lrtest <- c(tL[1] < min(c(y,yu)), tR[1] > max(c(y,yu)), tR[1] < tL[1],
                length(tL) > 1, length(tR) > 1)
    datatry <- try(tL[1]*tR[1], silent=TRUE)
    if(class(datatry)[1] == "try-error" | TRUE %in% lrtest){
      stop("Parameter tL or tR error, please check and fix.", call. = FALSE)
    }
    if(is.null(tL) | is.null(tR)){
      if(is.null(yu)){
        stop("yu data error, the yu data must be input for truncated model when parameter tL or tR is not set.", call. = FALSE)
      }
    }
  }

  if(!method[1] %in% c("EM","REG") | length(method) > 1){method <- "EM"}

  if(!type[1] %in% c("AI","RI","BC") | length(type) > 1){
    stop("Parameter type error, please input AI, RI, or BC.", call. = FALSE)
  }

  if(!is.null(D.matrix)){
    if(is.matrix(D.matrix)){
      if(type == "BC"){
        if(FALSE %in% (dim(D.matrix) == c(2,1))){
          stop("Design matrix error, please cheak and fix", call. = FALSE)
        }
      } else {
        if(FALSE %in% (dim(D.matrix) == c(3,2))){
          stop("Design matrix error, please cheak and fix", call. = FALSE)
        }
      }
    } else {
      stop("Design matrix error, please cheak and fix", call. = FALSE)
    }
  }

  if(!is.numeric(ng) | length(ng) > 1 | min(ng) < 1){
    stop("Parameter ng error, please input a positive integer.", call. = FALSE)
  }
  ng <- round(ng)

  if(!cM[1] %in% c(0,1) | length(cM > 1)){cM <- TRUE}

  if(!is.numeric(speed) | length(speed) > 1 | min(speed) < 0){
    stop("Parameter speed error, please input a positive number.", call. = FALSE)
  }

  if(!is.numeric(crit) | length(crit) > 1 | min(crit) <= 0 | max(crit) >= 1){
    stop("Parameter crit error, please input a positive number between 0 and 1.", call. = FALSE)
  }

  if(!d.eff[1] %in% c(0,1) | length(d.eff) > 1){d.eff <- TRUE}

  if(min(LRT.thre) < 0 | length(LRT.thre) > 1 | is.character(LRT.thre)){LRT.thre <- TRUE}

  if(!is.numeric(simu) | length(simu) > 1 | min(simu) < 50 | max(simu) > 10^8){
    simu = 1000
  }

  if(!is.numeric(alpha) | length(alpha) > 1 | min(alpha) <= 0 | max(alpha) >= 1){
    stop("Parameter alpha error, please input a positive number between 0 and 1.", call. = FALSE)
  }

  if(!detect[1] %in% c(0,1) | length(detect > 1)){detect <- TRUE}

  if(!is.numeric(QTLdist) | length(QTLdist) > 1 | min(QTLdist) < speed*2){
    stop("Parameter QTLdist error, please input a more suitable positive number.", call. = FALSE)
  }

  if(!plot.all[1] %in% c(0,1) | length(plot.all) > 1){plot.all <- TRUE}
  if(!plot.chr[1] %in% c(0,1) | length(plot.chr) > 1){plot.chr <- TRUE}
  if(!console[1] %in% c(0,1) | length(console) > 1){console <- TRUE}

  if(!cM){
    marker[, 2] <- marker[, 2]*100
  }

  if(is.null(D.matrix)){
    if(type == "BC"){
      D.matrix <- matrix(c(0.5, -0.5), 2, 1)
      row.names(D.matrix) <- c(2, 1)
      colnames(D.matrix) <- "a1"
    } else {
      a1 <- c(1, 0, -1)
      d1 <- c(-0.5, 0.5, -0.5)
      D.matrix <- cbind(a1, d1)
      row.names(D.matrix) <- c(2, 1, 0)
    }
  }
  if(!d.eff & type != "BC"){
    D.matrix <- matrix(D.matrix[, 1], 3, 1)
    row.names(D.matrix) <- c(2, 1, 0)
    colnames(D.matrix) <- "a1"
  }

  if(method == "EM"){
    meth <- methEM
    if(sele.g == "p" | sele.g == "f"){
      ya <- c(y, yu)
    } else {
      ya <- y
    }
  } else if (method=="REG"){
    if(sele.g == "p" | sele.g == "f"){
      ya <- c(y, yu)
    } else {
      ya <- y
    }

    meth <- methSG
  }

  cat(paste("chr", "cM", "LRT", "\n", sep = "\t")[console])
  effect <- c()
  cr0 <- unique(marker[, 1])
  for(i in cr0){
    cr <- marker[marker[, 1] == i,]
    minpos <- min(cr[, 2])
    if(speed%%1 == 0){minpos = ceiling(min(cr[, 2]))}
    for(j in seq(minpos, (max(cr[, 2])), speed)){
      QTL <- matrix(c(i, j), 1, 2)

      result <- meth(QTL, marker, geno, D.matrix, y, yu, tL, tR, type, ng, sele.g, crit, cM, ya)
      eff <- result[[1]]
      mu0 <- result[[2]]
      sigma <- result[[3]]
      LRT <- result[[4]]
      R2 <- result[[6]]
      model <- result[[7]]

      eff0 <- c(i, j, eff, LRT, R2)
      effect <- rbind(effect, eff0)
      cat(paste(i, j, LRT, "\n", sep = "\t")[console])
    }
  }
  row.names(effect) <- 1:nrow(effect)
  colnames(effect) <- c("chr", "cM", colnames(D.matrix), "LRT", "R2")
  effect <- data.frame(effect)
  effect[effect[, ncol(effect)] == Inf & !is.na(effect[, ncol(effect)]), 3:ncol(effect)] <- 0


  LRT.threshold <- NULL
  if(LRT.thre == TRUE){
    cat("LRT threshold calculation \n"[console])
    LRT.threshold <- LRTthre(marker, type = type, ng = ng, ns = length(y), gv = 25,
                             cM = cM, d.eff = d.eff, simu = simu, speed = speed,
                             alpha = alpha, console = console)
  } else if (is.numeric(LRT.thre)){
    LRT.threshold <- LRT.thre
  }

  detectQTL <- function(effect, LRT.threshold, QTLdist = 10){
    LRT <- effect[, c(1, 2, ncol(effect)-1)]
    LRT[LRT[, 3] < LRT.threshold, 3] <- 0
    det0 <- c()
    for(i in unique(LRT[, 1])){
      LRTcr <- LRT[LRT[, 1] == i, 2:3]
      detcr <- rep(0, nrow(LRTcr))
      if(sum(LRTcr[, 2]) > 0){
        for(j in 1:nrow(LRTcr)){
          LRTdet <- LRTcr[LRTcr[,1 ] >= max(c(min(LRTcr[, 1]), LRTcr[j, 1]-QTLdist)) &
                            LRTcr[, 1] <= min(c(max(LRTcr[, 1]), LRTcr[j, 1]+QTLdist)), 2]
          LRT0 <- LRTcr[j, 2]
          if(LRT0 == max(LRTdet) & LRT0 > 0){detcr[j] <- 1
          } else {detcr[j] <- 0}
        }
      }
      det0 <- c(det0, detcr)
    }
    detect.QTL <- effect[det0 == 1,]
    detect.QTL
  }

  detect.QTL <- NULL
  if(detect == TRUE & is.numeric(LRT.threshold)){
    detect.QTL <- detectQTL(effect, LRT.threshold, QTLdist)
  }

  if(plot.chr | plot.all){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
  }

  nc0 <- length(cr0)
  if(plot.chr | (nc0 == 1 & plot.all)){
    for(k in cr0){
      graphics::par(mar = c(5, 5, 4, 2))
      graphics::par(fig = c(0, 1, 0.35, 1))
      plot(effect[effect$chr == k,]$cM, effect[effect$chr == k,]$LRT, main = "", ylab = "LRT statistic",
           xlab = "", type = "l", ylim = c(0, max(c(effect$LRT, LRT.threshold), na.rm = TRUE)), bty = "l")
      if(!is.null(LRT.threshold)){
        graphics::abline(h = LRT.threshold, col = "red", lwd = 2)
      }
      graphics::par(mar = c(2, 5, 4, 2))
      graphics::par(fig = c(0, 1, 0, 0.45), new = TRUE)
      plot(effect[effect$chr == k,]$cM, effect[effect$chr == k,]$a1, ylab = "effect", main = paste("Chromosome", k),
           xlab = "", type = "l", ylim = c(min(c(effect$a1, 0), na.rm = TRUE)*1.4, max(c(effect$a1, 0), na.rm = TRUE)*1.2), xaxt = "n", bty = "n")
      graphics::abline(h=0,col="black",lwd=2)
    }
  }
  if(plot.all & nc0 > 1){
    ncr <- table(effect$chr)
    graphics::par(mar = c(5, 5, 4, 2))
    graphics::par(fig = c(0, 1, 0.3, 1))

    lcr <- c(0, cumsum(ncr))
    x0 <- c()
    s0 <- c()
    cut0 <- (max(lcr)*speed/length(lcr)/5)*(cr0-1)
    xn <- 0
    for(i in 1:length(ncr)){
      x1 <- seq(speed, ncr[i]*speed, speed)+lcr[i]*speed+cut0[i]
      x0 <- c(x0, x1)
      s0 <- c(s0, seq(min(x1), max(x1), 1))
      xn <- c(xn, length(x0))
    }
    xn1 <- (lcr[-1]-ncr/2)*speed+cut0
    yli <- max(c(effect$LRT, LRT.threshold), na.rm = TRUE)
    plot(x0, effect$LRT, main = "", ylab = "LRT statistic", xlab = "", type = "n",
         ylim = c(-yli/10, yli), xaxt = "n", bty = "l",axes = FALSE)
    for(k in 1:(length(xn)-1)){
      graphics::points(x0[(xn[k]+1):xn[k+1]], effect$LRT[(xn[k]+1):xn[k+1]], main = "", ylab = "LRT statistic",
                       xlab = "", type = "l", ylim = c(-yli/10, yli), xaxt = "n", bty = "l")
    }
    if(!is.null(LRT.threshold)){graphics::abline(h = LRT.threshold, col = "red", lwd = 2)}
    graphics::axis(2, seq(0, yli*1.2, 5))
    graphics::axis(1, cr0, at = xn1, cex = 1.2, tick = FALSE)
    lse <- 4000/max(x0)
    if(speed < 1){
      s0 <- x0
    }
    graphics::segments(s0, rep(-yli/10, length(s0)), s0, rep(-yli/5, length(s0)), lwd = lse)


    graphics::par(mar = c(2, 5, 4, 2))
    graphics::par(fig = c(0, 1, 0, 0.5), new = TRUE)
    plot(x0, effect$a1, main = "", ylab = "effect", xlab = "", type = "n",
         ylim = c(min(c(effect$a1, 0), na.rm = TRUE)*1.2, max(c(effect$a1, 0), na.rm = TRUE)*1.2), xaxt = "n", bty = "n")
    graphics::abline(h = 0, col = "blue", lwd = 2)
    for(k in 1:(length(xn)-1)){
      graphics::points(x0[(xn[k]+1):xn[k+1]], effect$a1[(xn[k]+1):xn[k+1]], main = "", ylab = "effect",
                       xlab = "", type = "l", ylim = c(min(c(effect$a1, 0), na.rm = TRUE)*1.2, max(c(effect$a1, 0))*1.2),
                       xaxt = "n", bty = "n")
    }
  }

  inputdata <- list(marker = marker, geno = geno, y = y, yu = yu, sele.g = sele.g, type = type, ng = ng, cM = cM,
                    d.eff = d.eff)

  result <- list(effect = effect, LRT.threshold = LRT.threshold, detect.QTL = detect.QTL, model = model, inputdata = inputdata)
  return(result)
}
