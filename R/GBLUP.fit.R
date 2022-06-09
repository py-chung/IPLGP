#' Muti-trait GBLUP Model
#'
#' Built the muti-trait GBLUP model using the phenotypic and genotypic data of a
#' training population by 'mmer' from R package 'sommer'. Then, output the fitted
#' values of the training population.
#'
#' @param t1 vector. The phenotype of trait1. The missing value must be coded as NA.
#' The length of all triat must be the same.
#' @param t2 vector. The phenotype of trait2. The missing value must be coded as NA.
#' The length of all triat must be the same.
#' @param t3 vector. The phenotype of trait3. The missing value must be coded as NA.
#' The length of all triat must be the same.
#' @param t4 vector. The phenotype of trait4. The missing value must be coded as NA.
#' The length of all triat must be the same.
#' @param t5 vector. The phenotype of trait5. The missing value must be coded as NA.
#' The length of all triat must be the same.
#' @param geno matrix. An n*p matrix with n individuals and p markers of the
#' training population. The markers must be coded as 1, 0, or -1 for alleles AA,
#' Aa, or aa. The missing value must have been already imputed.
#' @param K matrix. An n*n matrix denotes the genomic relationship matrix of the
#' training population if geno is set to be NULL.
#' @param outcross logical. A logical variable, if outcross is set to be TRUE,
#' the crop is regarded as an outcross crop. The kinship matrix of dominance
#' effects are also considered in the model. The geno data must be given when
#' outcross being TRUE.
#'
#' @return
#' \item{fitted.value}{The fitted values.}
#' \item{fitted.A}{The additive effect part of fitted values.}
#' \item{fitted.D}{The dominance effect part of fitted values.}
#' \item{mu}{The average value of fitted values.}
#'
#'
#' @note
#' Due to restrictions on the use of the funtion 'mmer', if an unknown error occurs
#' during use, please try to input the phenotype data as the format shown in the
#' example.
#'
#' @export
#'
#' @seealso
#' \code{\link[sommer]{mmer}}
#'
#' @references
#'
#' Habier D, Fernando RL, Dekkers JCM. 2007. The impact of genetic relationship
#' information on genome-assisted breeding values. Genetics 177:2389-2397.
#'
#' VanRaden PM. 2008. Efficient methods to compute genomic predictions.
#' J Dairy Sci. 91:4414-4423.
#'
#' @examples
#' # generate simulated data
#' t1 <- rnorm(50,30,10)
#' t2 <- rnorm(50,10,5)
#' t3 <- rnorm(50,20,20)
#' t4 <- NULL
#' t5 <- NULL
#'
#' # run with the marker score matrix
#' geno.test <- matrix(sample(c(1, -1), 5000, replace = TRUE), 50, 100)
#' result1 <- GBLUP.fit(t1, t2, t3, t4, t5, geno = geno.test)
#' result1$fitted.value
#'
#' # run with the genomic relationship matrix
#' K.test <- geno.test%*%t(geno.test)/ncol(geno.test)
#' result2 <- GBLUP.fit(t1, t2, t3, t4, t5, K = K.test)
#' result2$fitted.value
GBLUP.fit <- function(t1, t2, t3, t4, t5, geno = NULL, K = NULL, outcross = FALSE){

  phetest <- c(length(t1), length(t2), length(t3), length(t4), length(t5))
  if(length(table(phetest[phetest != 0])) != 1){
    stop("Phenotype data error, please input number vectors with the same length.", call. = FALSE)
  }

  phe <- cbind(t1, t2, t3, t4, t5)
  datatry <- try(phe*phe, silent=TRUE)
  if(class(datatry)[1] == "try-error"){
    stop("Phenotype data error, please input number vectors with the same length.", call. = FALSE)
  }
  nt <- ncol(phe)

  if(is.null(geno) & is.null(K)){
    stop("One of the arguments 'geno' and 'K' must be assigned.", call. = FALSE)
  }

  if(!outcross[1] %in% c(0,1)){
    outcross <- FALSE
  } else {outcross <- outcross[1]}

  if(is.null(geno)){
    if(outcross){
      stop("The geno data must be given when outcross being TRUE.", call. = FALSE)
    }

    KA <- K
    rownames(KA) <- 1:nrow(KA)
    colnames(KA) <- 1:nrow(KA)
    id <- rownames(KA)

    dat <- data.frame(id,phe)
    r.variable <- ~sommer::vs(id, Gu = KA, Gtc = sommer::unsm(nt))
  } else {
    tg <- table(as.matrix(geno))
    if(FALSE %in% (names(tg) %in% c(1, 0, -1))){
      stop("Genotype data error, please check and fix.", call. = FALSE)
    }
    if(outcross){
      if(length(geno[geno == 0]) == 0){
        warning("No heterozygous in genotype data of outcross crop, it may cause estimation errors.")
      }

      K0 <- geno.d(geno)
      KA <- K0$KA
      KD <- K0$KD

      rownames(KA) <- 1:nrow(KA)
      colnames(KA) <- 1:nrow(KA)
      rownames(KD) <- paste0("D",1:nrow(KD))
      colnames(KD) <- paste0("D",1:nrow(KD))

      id1 <- rownames(KA)
      id2 <- rownames(KD)

      dat <- data.frame(id1,id2,phe)
      dat$id1 <- factor(id1)
      dat$id2 <- factor(id2)

      r.variable <- ~sommer::vs(id1, Gu = KA, Gtc = sommer::unsm(nt)) + sommer::vs(id2, Gu = KD, Gtc = sommer::unsm(nt))
    } else {
      rownames(geno) <- 1:nrow(geno)
      id <- rownames(geno)
      KA <- geno%*%t(geno)/ncol(geno)

      dat <- data.frame(id,phe)

      r.variable <- ~sommer::vs(id, Gu = KA, Gtc = sommer::unsm(nt))
    }
  }

  datatry <- try(KA%*%as.matrix(phe), silent = TRUE)
  if(class(datatry)[1] == "try-error" | NA%in%KA){
    stop("Input data error, please cheak your input data.", call. = FALSE)
  }

  if(nt == 1){
    t1 <- c(phe)
    fit <- sommer::mmer(t1~1,
                        random = r.variable,
                        rcov = ~sommer::vs(units, Gtc = sommer::unsm(nt)),
                        data = dat)

    tol <- 10^-5
    while(length(fit) == 0){
      cat("Try bigger tolparinv. tolparinv:", tol, "\n")
      fit <- sommer::mmer(t1~1,
                          random = r.variable,
                          rcov = ~sommer::vs(units, Gtc = sommer::unsm(nt)),
                          data = dat,
                          tolparinv = tol)
      tol <- tol*10
    }
  } else {
    fit <- sommer::mmer(cbind(t1, t2, t3, t4, t5)~1,
                        random = r.variable,
                        rcov = ~sommer::vs(units, Gtc = sommer::unsm(nt)),
                        data = dat)

    tol <- 10^-5
    while(length(fit) == 0){
      cat("Try bigger tolparinv. tolparinv:", tol, "\n")
      fit <- sommer::mmer(cbind(t1, t2, t3, t4, t5)~1,
                          random = r.variable,
                          rcov = ~sommer::vs(units, Gtc = sommer::unsm(nt)),
                          data = dat,
                          tolparinv = tol)
      tol <- tol*10
    }
  }

  if(outcross){
    u0 <- list()
    for(i in 1:nt){
      u0[[i]] <- fit$U$`u:id1`[[i]]+fit$U$`u:id2`[[i]]
    }
    uA <- fit$U$`u:id1`
    uD <- fit$U$`u:id2`

    fitted.D <- matrix(unlist(uD), ncol = nt)
    colnames(fitted.D) <- names(u0)
    fitted.D <- fitted.D[order(as.numeric(names((u0[[1]])))),]

    if(!is.null(row.names(phe))){
      row.names(fitted.D) <- row.names(phe)
    }
    fitted.D <- as.matrix(fitted.D)
  } else {
    u0 <- fit$U$`u:id`
    uA <- fit$U$`u:id`
    uD <- 0

    fitted.D <- NULL
  }

  fitted.value <- matrix(unlist(u0), ncol = nt)
  fitted.A <- matrix(unlist(uA), ncol = nt)
  colnames(fitted.value) <- names(u0)
  colnames(fitted.A) <- names(u0)

  mu <- fit$fitted[1,]
  fitted.value <- fitted.value+matrix(mu, nrow(fitted.value), nt, byrow = TRUE)
  fitted.value <- fitted.value[order(as.numeric(names((u0[[1]])))),]
  fitted.A <- fitted.A[order(as.numeric(names((u0[[1]])))),]
  if(!is.null(row.names(phe))){
    row.names(fitted.value) <- row.names(phe)
    row.names(fitted.A) <- row.names(phe)
    names(mu) <- row.names(phe)
  }
  fitted.value <- as.matrix(fitted.value)
  fitted.A <- as.matrix(fitted.A)

  return(list(fitted.value = fitted.value, fitted.A = fitted.A, fitted.D = fitted.D, mu = mu))
}

