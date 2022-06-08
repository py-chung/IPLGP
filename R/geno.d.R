#' Generate the Genetic Design Matrix with dominance Effect
#'
#' Input the commonly used additive effect genetic design matrix to generate the
#' design matrix and kinship matrix of additive and dominance effects respectively.
#'
#' @param geno matrix. An n*p matrix denotes the commonly used additive effect
#' genetic design matrix of the training population.
#' @param AA number or character. The code denote alleles AA in the geno data.
#' @param Aa number or character. The code denote alleles Aa in the geno data.
#' @param aa number or character. The code denote alleles aa in the geno data.
#'
#' @return
#' \item{genoA}{An n*p matrix denote additive effects, and the markers are coded
#' as 1, 0, or -1 for alleles AA, Aa, or aa.}
#' \item{genoD}{An n*p matrix denote dominance effects, and the markers are coded
#' as 0.5, -0.5, or 0.5 for alleles AA, Aa, or aa.}
#' \item{KA}{An n*n matrix denote the kinship matrix of individuals with additive
#' effects. Whitch is caculated by genoA.}
#' \item{KD}{An n*n matrix denote the kinship matrix of individuals with dominance
#' effects. Whitch is caculated by genoD.}
#'
#' @export
#'
#' @references
#'
#' Cockerham, C. C., 1954. An extension of the concept of partitioning
#' hereditary variance for analysis of covariances among relatives When
#' epistasis is present. Genetics 39: 859–882.
#'
#' @examples
#'
#' geno <- rbind(rep(1,10),rep(0,10),rep(-1,10),c(rep(1,5),rep(-1,5)),c(rep(-1,5),rep(1,5)))
#' geno
#'
#' geno2 <- geno.d(geno)
#'
#' geno2$genoD
#' geno2$KD
geno.d <- function(geno, AA = 1, Aa = 0, aa = -1){
  geno <- as.matrix(geno)

  geno[geno == AA[1]] <- 1
  geno[geno == Aa[1]] <- 0
  geno[geno == aa[1]] <- -1
  geno[!geno %in% c(1, 0, -1)] <- 1

  genod <- geno
  genod <- 0.5-abs(genod)

  colnames(geno) <- paste(1:ncol(geno), "a", sep = "")
  colnames(genod) <- paste(1:ncol(genod), "d", sep = "")

  kd <- genod%*%t(genod)/ncol(genod)

  ka <- geno%*%t(geno)/ncol(geno)

  return(list(genoA = geno, genoD = genod, KA = ka, KD = kd))
}
