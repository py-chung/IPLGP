#' Standardize Phenotypic Values
#'
#' Standardize the phenotypic values of all the target traits from a training
#' population. Then, output the standardized phenotypic values, the mean vector,
#' and the standard deviation vector of the target traits.
#'
#' @param phe matrix. An n*t matrix with n individuals and t traits, denotes the
#' phenotypic values. The missing value must be coded as NA.
#'
#' @return
#' \item{standardize.phe}{An n*t matrix contains the standardized phenotypic values.}
#' \item{mu}{A vector with length t contains the averages of the phenotypic values
#' of the t target traits.}
#' \item{sd}{A vector with length t contains the standard deviations of the phenotypic
#' values of the t target traits.}
#'
#' @export
#'
#' @examples
#' # generate simulated data
#' phe.test <- data.frame(trait1 = rnorm(50,30,10), trait2 = rnorm(50,10,5), trait3 = rnorm(50,20,20))
#'
#' # run and output
#' result <- phe.sd(phe.test)
#' result
phe.sd <- function(phe){

  mean0 <- function(x){
    x1 <- mean(x, na.rm = TRUE)
    return(x1)
  }
  sd0 <- function(x){
    x1 <- stats::sd(x, na.rm = TRUE)
    return(x1)
  }

  phe <- as.matrix(phe)
  datatry <- try(phe%*%t(phe), silent = TRUE)
  if(class(datatry)[1] == "try-error"){
    stop("Phenotype data error, please cheak your phenotype data.", call. = F)
  }

  mu1 <- apply(phe, 2, mean0)
  sd1 <- apply(phe, 2, sd0)
  sd2 <- sd1
  sd1[sd1 == 0] <- 1

  phe1 <- t((t(phe)-mu1)*(1/sd1))
  return(list(standandize.phe = phe1, mu = mu1, sd = sd2))
}
