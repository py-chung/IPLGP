#' EM Algorithm for QTL MIM
#'
#' Expectation-maximization algorithm for QTL multiple interval mapping.
#' It can handle genotype data which is selective genotyping too.
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
#' @param geno matrix. A n*k matrix contains the genotypes of k markers
#' for n individuals. The marker genotypes of P1 homozygote (MM),
#' heterozygote (Mm), and P2 homozygote (mm) are coded as 2, 1, and 0,
#' respectively, with NA indicating missing values.
#' @param D.matrix matrix. The design matrix of QTL effects is a g*p
#' matrix, where g is the number of possible QTL genotypes, and p is the
#' number of effects considered in the MIM model. This design matrix can
#' be conveniently generated using the function D.make().
#' @param cp.matrix matrix. The conditional probability matrix is an
#' n*g matrix, where n is the number of genotyped individuals, and g is
#' the number of possible genotypes of QTLs. If cp.matrix=NULL, the
#' function will calculate the conditional probability matrix for selective
#' genotyping.
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
#' @param type character. The population type of the dataset. Includes
#' backcross (type="BC"), advanced intercross population (type="AI"), and
#' recombinant inbred population (type="RI"). The default value is "RI".
#' @param ng integer. The generation number of the population type. For
#' instance, in a BC1 population where type="BC", ng=1; in an AI F3
#' population where type="AI", ng=3.
#' @param cM logical. Specify the unit of marker position. If cM=TRUE, it
#' denotes centimorgan; if cM=FALSE, it denotes morgan.
#' @param E.vector0 vector. The initial value for QTL effects. The
#' number of elements corresponds to the column dimension of the design
#' matrix. If E.vector0=NULL, the initial value for all effects will be
#' set to 0.
#' @param X matrix. The design matrix of the fixed factors except for
#' QTL effects. It is an n*k matrix, where n is the number of
#' individuals, and k is the number of fixed factors. If X=NULL,
#' the matrix will be an n*1 matrix where all elements are 1.
#' @param beta0 vector. The initial value for effects of the fixed
#' factors. The number of elements corresponds to the column dimension
#' of the fixed factor design matrix.  If beta0=NULL, the initial value
#' will be set to the average of y.
#' @param variance0 numeric. The initial value for variance. If
#' variance0=NULL, the initial value will be set to the variance of
#' phenotype values.
#' @param crit numeric. The convergence criterion of EM algorithm.
#' The E and M steps will iterate until a convergence criterion is met.
#' It must be a value between 0 and 1.
#' @param stop numeric. The stopping criterion of EM algorithm. The E and
#' M steps will halt when the iteration number reaches the stopping
#' criterion, treating the algorithm as having failed to converge.
#' @param conv logical. If set to False, it will disregard the failure to
#' converge and output the last result obtained during the EM algorithm
#' before reaching the stopping criterion.
#' @param console logical. Determines whether the process of the algorithm
#' will be displayed in the R console or not.
#' @param IMresult list. The data list of the output from IM.search(). The
#' required parameters for this function will be extracted from the data list.
#' @param MIMresult list. The data list of the output from MIM.search() or
#' MIM.points(). The required parameters for this function will be extracted
#' from the data list.
#'
#' @return
#' \item{QTL}{The QTL imformation of this analysis.}
#' \item{E.vector}{The QTL effects are calculated by the EM algorithm.}
#' \item{beta}{The effects of the fixed factors are calculated by the EM
#' algorithm.}
#' \item{variance}{The variance is calculated by the EM algorithm.}
#' \item{PI.matrix}{The posterior probabilities matrix after the
#' process of the EM algorithm.}
#' \item{log.likelihood}{The log-likelihood value of this model.}
#' \item{LRT}{The LRT statistic of this model.}
#' \item{R2}{The coefficient of determination of this model. This
#' can be used as an estimate of heritability.}
#' \item{y.hat}{The fitted values of trait values with genotyping are
#' calculated by the estimated values from the EM algorithm.}
#' \item{yu.hat}{The fitted values of trait values without genotyping
#' are calculated by the estimated values from the EM algorithm.}
#' \item{iteration.number}{The iteration number of the EM algorithm.}
#' \item{model}{The model of this analysis, contains complete a
#' genotyping model, a proposed model, a truncated model, and a
#' frequency-based model.}
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
#' @seealso
#' \code{\link[QTLEMM]{D.make}}
#' \code{\link[QTLEMM]{Q.make}}
#' \code{\link[QTLEMM]{IM.search}}
#' \code{\link[QTLEMM]{MIM.search}}
#' \code{\link[QTLEMM]{MIM.points}}
#'
#' @examples
#'
#' # load the example data
#' load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' D.matrix <- D.make(3, type = "RI", aa = c(1, 3, 2, 3), dd = c(1, 2, 1, 3), ad = c(1, 2, 2, 3))
#' cp.matrix <- Q.make(QTL, marker, geno, type = "RI", ng = 2)$cp.matrix
#' result <- EM.MIM(D.matrix = D.matrix, cp.matrix = cp.matrix, y = y)
#' result$E.vector
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
#' D.matrix <- D.make(3, type = "RI", aa = c(1, 3, 2, 3), dd = c(1, 2, 1, 3), ad = c(1, 2, 2, 3))
#' result <- EM.MIM(QTL, marker, geno.s, D.matrix, y = ys, yu = yu, sele.g = "f")
#' result$E.vector
#' }
#'
EM.MIM <- function(QTL = NULL, marker = NULL, geno = NULL, D.matrix = NULL, cp.matrix = NULL, y = NULL,
                   yu = NULL, sele.g = "n", tL = NULL, tR = NULL, type = "RI", ng = 2, cM = TRUE,
                   E.vector0 = NULL, X = NULL, beta0 = NULL, variance0 = NULL, crit = 10^-5,
                   stop = 1000, conv = TRUE, console = TRUE, IMresult = NULL, MIMresult = NULL){

  if(!is.null(IMresult)){
    check1 <- names(IMresult) == c("effect", "LRT.threshold", "detect.QTL", "model", "inputdata")
    if(!is.list(IMresult) | sum(check1) != 5){
      stop("IMresult data error, please input all of the original output data of IM.search().", call. = FALSE)
    } else {
      check2 <- names(IMresult$inputdata) == c("marker", "geno", "y", "yu", "sele.g", "type", "ng", "cM", "d.eff" )
      if(!is.list(IMresult$inputdata) | sum(check2) != 9){
        stop("IMresult data error, please input all of the original output data of IM.search().", call. = FALSE)
      } else {
        if(is.data.frame(IMresult$detect.QTL) | is.matrix(IMresult$detect.QTL)){
          QTL <- as.matrix(IMresult$detect.QTL[,1:2])
        } else {stop("IMresult data error, please input all of the original output data of IM.search().", call. = FALSE)}
        marker <- IMresult$inputdata$marker
        geno <- IMresult$inputdata$geno
        y <- IMresult$inputdata$y
        yu <- IMresult$inputdata$yu
        sele.g <- IMresult$inputdata$sele.g
        type <- IMresult$inputdata$type
        ng <- IMresult$inputdata$ng
        cM <- IMresult$inputdata$cM
        if(is.null(D.matrix)){
          D.matrix <- D.make(nrow(QTL), type = type)
          if(length(IMresult$inputdata$d.eff) != 0){
            if(IMresult$inputdata$d.eff[1] == 0){
              D.matrix <- D.make(nrow(QTL), type = type, d = 0)
            }
          }
        }
        cat("Use the output data of IM.search() for analysis. \n")
        cat("If an error occurs, please check whether the original output data of IM.search() is used. \n")
      }
    }
  } else if(!is.null(MIMresult)){
    check1 <- names(MIMresult) == c("effect", "QTL.best", "effect.best", "model", "inputdata")
    if(!is.list(MIMresult) | sum(check1) != 5){
      stop("MIMresult data error, please input all of the original output data of MIM.search()/MIM.points().", call. = FALSE)
    } else {
      check2 <- names(MIMresult$inputdata) == c("marker", "geno", "y", "yu", "sele.g", "type", "ng", "cM", "D.matrix" )
      if(!is.list(MIMresult$inputdata) | sum(check2) != 9){
        stop("MIMresult data error, please input all of the original output data of MIM.search()/MIM.points().", call. = FALSE)
      } else {
        if(is.data.frame(MIMresult$QTL.best) | is.matrix(MIMresult$QTL.best)){
          QTL <- as.matrix(MIMresult$QTL.best[,1:2])
        } else {stop("MIMresult data error, please input all of the original output data of MIM.search()/MIM.points().", call. = FALSE)}
        marker <- MIMresult$inputdata$marker
        geno <- MIMresult$inputdata$geno
        y <- MIMresult$inputdata$y
        yu <- MIMresult$inputdata$yu
        sele.g <- MIMresult$inputdata$sele.g
        type <- MIMresult$inputdata$type
        ng <- MIMresult$inputdata$ng
        cM <- MIMresult$inputdata$cM
        if(is.null(D.matrix)){
          D.matrix <- MIMresult$inputdata$D.matrix
        }
        cat("Use the output data of MIM.search()/MIM.points() for analysis. \n")
        cat("If an error occurs, please check whether the original output data of MIM.search()/MIM.points() is used. \n")
      }
    }
  }

  if(!sele.g[1] %in% c("n", "t", "p", "f") | length(sele.g)!=1){
    stop("Parameter sele.g error, please check and fix. It must be 'n', 't', 'p', or 'f'.", call. = FALSE)
  }

  if((is.null(cp.matrix) | sele.g != "n") & (is.null(QTL) | is.null(marker) | is.null(geno))){
    stop("Input data is missing. The argument QTL, marker, and geno must be input when analize the selective genotyping data or cp.matrix is set to null.", call. = FALSE)
  }

  if(is.null(D.matrix)){
    stop("Input data is missing. The argument D.matrix must be input.", call. = FALSE)
  }

  if(!is.null(geno)){
    genotest <- table(c(geno))
    datatry <- try(geno*geno, silent=TRUE)
    if(class(datatry)[1] == "try-error" | FALSE %in% (names(genotest) %in% c(0, 1, 2))  | length(dim(geno)) != 2){
      stop("Genotype data error, please cheak your genotype data.", call. = FALSE)
    }
  }

  if(!is.null(marker)){
    marker <- as.matrix(marker)
    markertest <- c(ncol(marker) != 2, NA %in% marker, marker[,1] != sort(marker[,1]), nrow(marker) != ncol(geno))
    datatry <- try(marker*marker, silent=TRUE)
    if(class(datatry)[1] == "try-error" | TRUE %in% markertest){
      stop("Marker data error, or the number of marker does not match the genotype data.", call. = FALSE)
    }
  }

  if(!is.null(QTL)){
    QTL <- as.matrix(QTL)
    datatry <- try(QTL*QTL, silent=TRUE)
    if(class(datatry)[1] == "try-error" | ncol(QTL) != 2 | NA %in% QTL | max(QTL[, 1]) > max(marker[, 1])){
      stop("QTL data error, please cheak your QTL data.", call. = FALSE)
    }

    for(i in 1:nrow(QTL)){
      ch0 <- QTL[i, 1]
      q0 <- QTL[i, 2]
      if(!ch0 %in% marker[, 1]){
        stop("The specified QTL is not in the position that can estimate by the marker data.", call. = FALSE)
      }
      if(q0 > max(marker[marker[, 1] == ch0, 2]) | q0 < min(marker[marker[, 1] == ch0, 2])){
        stop("The specified QTL is not in the position that can estimate by the marker data.", call. = FALSE)
      }
    }

    nQTL <- nrow(QTL)
  }

  if(!is.null(cp.matrix)){
    datatry <- try(y%*%cp.matrix%*%D.matrix, silent=TRUE)
    if(class(datatry)[1] == "try-error" | NA %in% D.matrix | NA %in% cp.matrix){
      stop("Input data error, please check your input data y, D.matrix, and cp.matrix.", call. = FALSE)
    }
  } else {
    datatry <- try(y%*%t(y), silent=TRUE)
    datatry1 <- try(D.matrix*D.matrix, silent=TRUE)
    if(class(datatry)[1] == "try-error" | class(datatry1)[1] == "try-error" | NA %in% D.matrix | length(y) < 2){
      stop("Input data error, please check your input data y and D.matrix.", call. = FALSE)
    }
  }

  ind <- length(y)
  g <- nrow(D.matrix)
  eff <- ncol(D.matrix)

  y[is.na(y)] <- mean(y,na.rm = TRUE)

  E.vector <- E.vector0
  beta <- beta0
  variance <- variance0

  if(!is.null(yu)){
    yu <- yu[!is.na(yu)]
    datatry <- try(yu%*%yu, silent=TRUE)
    if(class(datatry)[1] == "try-error"){
      stop("yu data error, please check your yu data.", call. = FALSE)
    }
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

  if(!type[1] %in% c("AI","RI","BC") | length(type) > 1){
    stop("Parameter type error, please input AI, RI, or BC.", call. = FALSE)
  }

  if(!is.numeric(ng) | length(ng) > 1 | min(ng) < 1){
    stop("Parameter ng error, please input a positive integer.", call. = FALSE)
  }
  ng <- round(ng)

  if(!cM[1] %in% c(0,1) | length(cM) > 1){cM <- TRUE}

  if(sele.g == "p" | sele.g == "f"){
    if(is.null(yu)){
      stop("yu data error, the yu data must be input for proposed model and frequency-based model", call. = FALSE)
    }
    ya <- c(y, yu)

    QTL.p <- function(QTL, marker, geno, D.matrix, ys, yu, cM = TRUE, type = "RI", ng = 2,
                      E.vector = NULL, X = NULL, beta = NULL, variance = NULL, cp.matrix = NULL,
                      crit = 10^-5, stop = 1000, conv = TRUE, console = FALSE, pf = FALSE){
      y <- c(ys,yu)
      N <- length(y)
      nu <- length(yu)
      nQTL <- nrow(QTL)
      model <- 1
      if(pf == T){model <- 3}
      mp <- mixprop(QTL, nu, marker, geno, model = model, cM = cM, type = type, ng = ng, cp.matrix = cp.matrix)
      Freq <- mp[[2]]
      QTL <- mp[[3]]
      marker <- mp[[4]]
      n.para <- ncol(D.matrix)
      Qn <- nrow(D.matrix)

      L0 <- sum(log(stats::dnorm(y,mean = mean(y), sd = stats::var(c(y))^0.5)))
      if(is.null(X)){
        X <- matrix(1, N, 1)
      } else if (is.vector(X)){
        X <- matrix(X, N, 1)
      }
      if(is.null(beta)){
        mut <- matrix(rep(mean(y), ncol(X)), ncol(X), 1)
      } else if (is.numeric(beta)){
        mut <- matrix(rep(beta, ncol(X)), ncol(X), 1)
      } else {mut <- beta}
      if(is.null(E.vector)){
        Et <- rep(0, n.para)
      } else {Et <- E.vector}
      if(is.null(variance)){
        sigt <- stats::var(y)
      } else {sigt <- variance}
      MUt <- t(D.matrix%*%Et%*%matrix(1, 1, N))+X%*%mut%*%matrix(1, 1, Qn)
      PIup <- matrix(0, N, Qn)
      for(j in 1:N){
        P0 <- c()
        for(i.p in 1:Qn){
          P0[i.p] <- Freq[j, i.p]*stats::dnorm(y[j], mean = MUt[j, i.p], sd = sigt^0.5)
        }
        if(sum(P0) != 0){P0 <- P0/sum(P0)
        } else {P0 <- rep(1/g,g)}
        PIup[j,] <- P0
      }
      PI <- PIup

      R <- matrix(0, n.para, 1)
      M <- matrix(0, n.para, n.para)
      V <- matrix(0, n.para, n.para)
      one <- rep(1, N)
      for(m1 in 1:n.para){
        R[m1,] <- (t(y-X%*%mut)%*%PI%*%D.matrix[, m1])/(one%*%PI%*%(D.matrix[, m1]*D.matrix[, m1]))
        for(m2 in 1:n.para){
          V[m1, m2] <- one%*%PI%*%(D.matrix[, m1]*D.matrix[, m2])
          if(m1 != m2){
            M[m1, m2] <- (one%*%PI%*%(D.matrix[, m1]*D.matrix[, m2]))/(one%*%PI%*%(D.matrix[, m1]*D.matrix[, m1]))
          }
        }
      }
      Et1 <- R-M%*%Et
      mut1 <- solve(t(X)%*%X)%*%t(X)%*%(y-PI%*%D.matrix%*%Et1)
      sigt1 <- (t(y-X%*%mut1)%*%(y-X%*%mut1)-2*(t(y-X%*%mut1)%*%PI%*%D.matrix%*%Et1)+(t(Et1)%*%V%*%Et1))/N
      MUt1 <- t(D.matrix%*%Et1%*%matrix(1, 1, N))+X%*%mut1%*%matrix(1, 1, Qn)
      Lih10 <- matrix(0, N, Qn)
      Lih1 <- matrix(0, N, Qn)
      for(i.L in 1:Qn){
        Lih10[, i.L] <- Freq[, i.L]*stats::dnorm(y, mean = MUt[, i.L], sd = sigt^0.5)
        Lih1[, i.L] <- Freq[, i.L]*stats::dnorm(y, mean = MUt1[, i.L], sd = sigt1^0.5)
      }
      L10 <- sum(log(rowSums(Lih10)))
      L1 <- sum(log(rowSums(Lih1)))
      number <- 1

      cat(paste("number", "var", paste(effectname, collapse = "\t"), "\n", sep = "\t")[console])
      while(max(abs(c(Et1-Et,sigt1-sigt))) >= crit & number < stop){
        repeat{
          mut <- mut1
          Et <- Et1
          sigt <- sigt1
          MUt <- MUt1
          L10 <- L1
          PIup <- matrix(0, N, Qn)
          for(j in 1:N){
            P0 <- c()
            for(i.p in 1:Qn){
              P0[i.p] <- Freq[j, i.p]*stats::dnorm(y[j], mean = MUt[j, i.p], sd = sigt^0.5)
            }
            if(sum(P0) != 0){P0 <- P0/sum(P0)
            } else {P0 <- rep(1/g,g)}
            PIup[j,] <- P0
          }
          PI <- PIup
          for(m1 in 1:n.para){
            R[m1,] <- (t(y-X%*%mut)%*%PI%*%D.matrix[, m1])/(one%*%PI%*%(D.matrix[, m1]*D.matrix[, m1]))
            for(m2 in 1:n.para){
              V[m1, m2] <- one%*%PI%*%(D.matrix[, m1]*D.matrix[, m2])
              if(m1 != m2){
                M[m1, m2] <- (one%*%PI%*%(D.matrix[, m1]*D.matrix[, m2]))/(one%*%PI%*%(D.matrix[, m1]*D.matrix[, m1]))
              }
            }
          }
          Et1 <- R-M%*%Et
          mut1 <- solve(t(X)%*%X)%*%t(X)%*%(y-PI%*%D.matrix%*%Et1)
          sigt1 <- (t(y-X%*%mut1)%*%(y-X%*%mut1)-2*(t(y-X%*%mut1)%*%PI%*%D.matrix%*%Et1)+(t(Et1)%*%V%*%Et1))/N
          MUt1 <- t(D.matrix%*%Et1%*%matrix(1, 1, N))+X%*%mut1%*%matrix(1, 1, Qn)
          Lih10 <- matrix(0, N, Qn)
          Lih1 <- matrix(0, N, Qn)
          for(i.L in 1:Qn){
            Lih10[, i.L] <- Freq[, i.L]*stats::dnorm(y, mean = MUt[, i.L], sd = sigt^0.5)
            Lih1[, i.L] <- Freq[, i.L]*stats::dnorm(y, mean = MUt1[, i.L], sd = sigt1^0.5)
          }
          L10 <- sum(log(rowSums(Lih10)))
          L1 <- sum(log(rowSums(Lih1)))
          break(max(abs(c(Et1-Et,sigt1-sigt))) < crit)
        }
        Ep <- round(Et1, 3)
        sp <- round(sigt1, 3)
        cat(paste(number, sp, paste(Ep, collapse = "\t"), "\n", sep = "\t")[console])
        number <- number+1

        if(NaN %in% (Et1-Et)){
          Et1 <- Et
          break()
        }
      }
      number <- number-1

      like0 <- L0
      like1 <- L1
      LRT <- 2*(like1-like0)
      parameter <- c(mut1, Et1, sigt1)

      y.hat <- PI%*%D.matrix%*%Et1+X%*%mut1
      r2 <- stats::var(y.hat)/stats::var(y)

      Et1 <- parameter[2:(1+ncol(D.matrix))]
      mut1 <- parameter[1]
      sigt1 <- parameter[length(parameter)]
      colnames(PI) <- colnames(mp)

      if(number == stop){
        if(conv){
          Et1 <- rep(0, length(E.vector))
          mut1 <- 0
          sigt1 <- 0
          PI <- matrix(0, nrow(PI), ncol(PI))
          like1 <- -Inf
          LRT <- 0
          r2 <- 0
        }
        warning("EM algorithm fails to converge, please check the input data or adjust the convergence criterion.")
      }

      output <- list(QTL = QTL, E.vector = Et1, beta = mut1, variance = sigt1, PI.matrix = PI,
                     log.likelihood = like1, LRT = LRT, R2 = r2, y.hat = y.hat[1:length(ys)],
                     yu.hat = y.hat[-(1:length(ys))], iteration.number = number)
      return(output)
    }
  } else {ya <- y}

  if(is.null(E.vector)){E.vector <- rep(0, eff)}
  if(is.null(X)){
    X <- matrix(1,length(ya),1)
  }else if(is.vector(X)){
    X <- matrix(X, length(ya), 1)
  }
  if(is.null(beta)){
    beta <- matrix(rep(mean(y), ncol(X)), ncol(X), 1)
  }else if(is.numeric(beta)){
    beta <- matrix(rep(beta, ncol(X)), ncol(X), 1)
  }
  if(is.null(variance)){variance <- stats::var(ya)}
  if(!console[1] %in% c(0,1) | length(console) > 1){console <- TRUE}
  if(!conv[1] %in% c(0,1) | length(conv) > 1){conv <- TRUE}

  datatry <- try(D.matrix%*%E.vector, silent=TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% E.vector){
    stop("Parameter E.vector0 error, please check and fix.", call. = FALSE)
  }

  datatry <- try(ya%*%X%*%beta, silent=TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% X | NA %in% beta){
    stop("Parameter X or bata0 error, please check and fix.", call. = FALSE)
  }

  if(!is.numeric(variance) | length(variance) > 1 | min(variance) < 0){
    stop("Parameter variance0 error, please input a positive number.", call. = FALSE)
  }
  sigma <- sqrt(variance)

  if(!is.numeric(crit) | length(crit) > 1 | min(crit) <= 0 | max(crit) >= 1){
    stop("Parameter crit error, please input a positive number between 0 and 1.", call. = FALSE)
  }

  if(!is.numeric(stop) | length(stop) > 1 | min(crit) <= 0){
    stop = 1000
    warning("Parameter stop error, adjust to 1000.")
  }

  if(length(colnames(D.matrix)) == ncol(D.matrix)){
    effectname <- colnames(D.matrix)
  }

  if(sele.g == "f"){
    if(nQTL > 1){QTL <- QTL[order(QTL[, 1], QTL[, 2]),]}
    marker <- marker[order(marker[, 1],marker[, 2]),]
    result <- QTL.p(QTL, marker, geno, D.matrix, ys = y, yu = yu, cM = cM, type = type, ng = ng,
                    E.vector = E.vector, X = X, beta = beta, variance = variance,
                    cp.matrix = cp.matrix, crit = crit, stop = stop, conv = conv, console = console)
    model <- "proposed corrected frequency model (Lee 2014) of selective genotyping"
    result[[12]] <- model
    names(result)[12] <- "model"
  } else if (sele.g == "p"){
    if(nQTL > 1){QTL <- QTL[order(QTL[, 1], QTL[, 2]),]}
    marker <- marker[order(marker[, 1], marker[, 2]),]
    result <- QTL.p(QTL, marker, geno, D.matrix, ys = y, yu = yu, cM = cM, type = type, ng = ng,
                    E.vector = E.vector, X = X, beta = beta, variance = variance,
                    cp.matrix = cp.matrix, crit = crit, stop = stop, conv = conv, console = console, pf = TRUE)
    model <- "population frequency-based model of selective genotyping"
    result[[12]] <- model
    names(result)[12] <- "model"
  } else if (sele.g == "t"){
    QTL.t <- function(QTL, marker, geno, D.matrix, ys, yu, tL, tR, cM = TRUE, type = "RI", ng = 2,
                      E.vector = NULL, beta = NULL,  X=NULL, variance = NULL, cp.matrix = NULL,
                      crit = 10^-5, stop = 1000, conv =TRUE, console = FALSE){
      maxit <- stop
      nS <- length(ys)
      nType <- nrow(D.matrix)
      nE <- ncol(D.matrix)
      nu <- length(yu)
      mp <- mixprop(QTL, nu, marker, geno, model = 2, cM = cM, type = type, ng = ng, cp.matrix = cp.matrix)
      Freq <- mp[[2]]
      QTL <- mp[[3]]
      marker <- mp[[4]]
      y <- c(ys, yu)
      init <- c(mean(y), rep(0, nE), stats::sd(y))
      negLL0 <- function(theta) {
        mu <- theta[1]
        sigma <- theta[2]
        tauL <- (tL - mu)/sigma
        tauR <- (tR - mu)/sigma
        U <- 1 + stats::pnorm(tauL) - stats::pnorm(tauR)
        LL <- sum(log(stats::dnorm(ys, mu, sigma)/U))
        return(-LL)
      }
      temp <- stats::optim(c(mean(y), stats::sd(y)), negLL0, method = "Nelder-Mead", control = list(maxit = 3000))
      LL0 <- -1 * (temp$value)

      err <- 1
      if(is.null(X)){
        X <- matrix(1, nS, 1)
        if(length(beta > 1)){
          beta <- mean(beta)
        }
      } else if (is.vector(X)){
        X <- matrix(1, nS, 1)
      }
      if(is.null(beta)){
        mu0 <- matrix(rep(mean(y), ncol(X)), ncol(X), 1)
      } else if (is.numeric(beta)){
        mu0 <- matrix(rep(beta, ncol(X)), ncol(X), 1)
      } else {mut <- beta}
      if(is.null(E.vector)){
        E <- rep(0, nE)
      } else {E <- E.vector}
      if(is.null(variance)){
        sigma <- stats::sd(y)
      } else {sigma <- sqrt(variance)}
      mu <- t(D.matrix%*%E%*%matrix(1, 1, nS))+X%*%mu0%*%matrix(1, 1, nType)

      tauL <- (tL - X%*%mu0%*%matrix(1, 1, nType) - t(D.matrix%*%E%*%matrix(1, 1, nS)))/sigma
      tauR <- (tR - X%*%mu0%*%matrix(1, 1, nType) - t(D.matrix%*%E%*%matrix(1, 1, nS)))/sigma
      U <- 1+stats::pnorm(tauL)-stats::pnorm(tauR)
      LL <- sum(log(rowSums(Freq*stats::dnorm(ys, mu, sigma)/U)))
      Et <- matrix(NA, 3000, nE)
      Mt <- c()
      St <- c()

      cat(paste("number", "var", paste(effectname, collapse = "\t"), "\n", sep = "\t")[console])
      for (iter in 1:maxit) {
        E0 <- E
        sigma0 <- sigma
        mu1 <- mu0
        Pi1 <- Freq*stats::dnorm(ys, mu, sigma)/U
        if (NaN %in% Pi1){
          break()
        }
        Pi <- Pi1
        Pi[rowSums(Pi) == 0,] <- 1
        Pi <- Pi/rowSums(Pi)
        Amy <- (stats::dnorm(tauL)-stats::dnorm(tauR))/U
        Bob <- (tauL*stats::dnorm(tauL)-tauR*stats::dnorm(tauR))/U
        V <- matrix(0, nE, nE)
        for (i in 1:nE) {
          V[i, ] <- colSums(Pi)%*%(D.matrix[, i]*D.matrix)
        }
        M <- V/diag(V)-diag(1, nE)
        r <- (t(ys-X%*%mu0)%*%Pi+sigma*colSums(Pi*Amy))%*%D.matrix/colSums(Pi%*%D.matrix^2)
        E <- t(r)-M%*%E
        mu0 <- solve(t(X)%*%X)%*%t(X)%*%(ys-Pi%*%D.matrix%*%E+sigma*rowSums(Pi*Amy))
        sigma <- (t(ys-X%*%mu0)%*%(ys-X%*%mu0)-2*t(ys-X%*%mu0)%*%Pi%*%D.matrix%*%E+t(E)%*%V%*%E)/(nS-sum(Pi*Bob))
        sigma <- sqrt(as.vector(sigma))
        mu <- t(D.matrix%*%E%*%matrix(1, 1, nS))+X%*%mu0%*%matrix(1, 1, nType)
        tauL <- (tL-X%*%mu0%*%matrix(1, 1, nType)-t(D.matrix%*%E%*%matrix(1, 1, nS)))/sigma
        tauR <- (tR-X%*%mu0%*%matrix(1, 1, nType)-t(D.matrix%*%E%*%matrix(1, 1, nS)))/sigma
        U <- 1+stats::pnorm(tauL)-stats::pnorm(tauR)
        LL1 <- sum(log(rowSums(Freq*stats::dnorm(ys, mu, sigma)/U)))
        err <- abs(LL1 - LL)
        LL <- LL1
        Ep <- round(E, 3)
        sp <- round(sigma^2, 3)
        cat(paste(iter, sp, paste(Ep, collapse = "\t"), "\n", sep = "\t")[console])

        Et[iter,] <- E
        Mt[iter] <- mu0
        St[iter] <- sigma
        if (max(abs(c(E-E0, mu0-mu1, sigma-sigma0))) < crit){
          break()
        }
        if(iter == maxit){
          if(ncol(Et) > 1){
            E <- apply(Et[(maxit-119):maxit,], 2, mean)
          } else {E = mean(Et)}
          mu0 <- mean(Mt[(maxit-119):maxit])
          sigma <- mean(St[(maxit-119):maxit])
        }
      }
      parameter <- c(mu0, E, sigma^2)
      number <- iter

      like0 <- 2*LL0
      like1 <- 2*LL1
      LRT <- like1 - like0

      Et1 <- parameter[2:(1+ncol(D.matrix))]
      mut1 <- parameter[1]
      sigt1 <- parameter[length(parameter)]
      colnames(Pi) <- colnames(mp)

      y.hat <- Pi%*%D.matrix%*%Et1+X%*%mut1
      r2 <- stats::var(y.hat)/stats::var(ys)

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
        warning("EM algorithm fails to converge, please check the input data or adjust the convergence criterion and stopping criterion.")
      }

      output <- list(QTL = QTL, E.vector = Et1, beta = mut1, variance = sigt1, PI.matrix = Pi,
                     log.likelihood = like1, LRT = LRT, R2 = r2, y.hat = y.hat, yu.hat = NULL,
                     iteration.number = number)
      return(output)
    }

    if(nQTL>1){QTL <- QTL[order(QTL[, 1], QTL[, 2]),]}
    marker <- marker[order(marker[, 1], marker[, 2]),]
    if(is.null(tL)){tL <- min(yu)}
    if(is.null(tR)){tR <- max(yu)}
    result <- QTL.t(QTL, marker, geno, D.matrix, ys = y, yu = yu, tL, tR, cM = cM, type = type, ng = ng,
                    E.vector = E.vector, beta = beta, X = X, variance = variance,
                    cp.matrix = cp.matrix, crit = crit, stop = stop, conv = conv, console = console)
    model <- "truncated model of selective genotyping"
    result[[12]] <- model
    names(result)[12] <- "model"
  } else {
    if(is.null(cp.matrix)){
      cp.matrix <- Q.make(QTL, marker, geno, cM = cM, type = type, ng = ng)[[(nrow(QTL)+1)]]
    }
    re <- EM.MIM0(D.matrix, cp.matrix, y, E.vector0 = E.vector, X = X, beta0 = beta, variance0 = variance,
                  crit = crit, stop = stop, conv = conv, console = console)
    result <- list(QTL = QTL, E.vector = re[[1]], beta = re[[2]], variance = re[[3]], PI.matrix = re[[4]],
                   log.likelihood = re[[5]], LRT = re[[6]], R2 = re[[7]], y.hat = re[[8]], yu.hat = NULL,
                   iteration.number = re[[9]], model = "complete genotyping model")
  }
  names(result[[2]]) <- effectname
  return(result)
}
