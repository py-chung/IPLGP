#' Search For A Subset With The Highest D-score
#'
#' Search for an optimal subset of the candidate individuals such that it
#' achieves the highest D-score by genetic algorithm (GA).
#'
#' @param K matrix. An n*n matrix denotes the genomic relationship matrix of
#' the n candidate individuals, where n > 4.
#' @param size integer. An integer denotes the size of the subset, note that
#' 3 < size < n.
#' @param keep vector. A vector indicates those candidate individuals which
#' will be retained in the subset before the search. The length of keep must
#' be less than size.
#' @param n0 integer. An integer indicates the number of chromosomes
#' (solutions) in the genetic algorithm, note that n0 > 3.
#' @param mut integer. An integer indicates the number of mutations in the
#' genetic algorithm, note that mut < size.
#' @param cri integer. An integer indicates the stopping criterion, note that
#' cri < 1e+06. The genetic algorithm will stop if the number of iterations
#' reaches cri.
#' @param console logical. A logical variable, if console is set to be TRUE,
#' the searching process will be shown in the R console.
#'
#' @return
#' \item{subset}{The optimal subset with the highest D-score.}
#' \item{D.score}{The D.score of the optimal subset.}
#' \item{time}{The number of iterations.}
#'
#' @export
#'
#' @references
#'
#' Chung PY, Liao CT. 2020. Identification of superior parental lines for
#' biparental crossing via genomic prediction. PLoS ONE 15(12):e0243159.
#'
#' Ou JH, Liao CT. 2019. Training set determination for genomic selection.
#' Theor Appl Genet. 132:2781-2792.
#'
#' @examples
#' # generate simulated data
#' geno.test <- matrix(sample(c(1, -1), 600, replace = TRUE), 20, 30)
#' K.test <- geno.test%*%t(geno.test)/ncol(geno.test)
#'
#' # run with no specified individual
#' result1 <- GA.Dscore(K.test, 6, cri = 1000, console = TRUE)
#' result1
#'
#' # run with some specified individuals
#' result2 <- GA.Dscore(K.test, 6, keep = c(1, 5, 10), cri = 1000, console = TRUE)
#' result2
GA.Dscore <- function(K, size, keep=c(), n0 = size, mut = 3, cri = 10000, console = FALSE){

  datatry <- try(K%*%K, silent=TRUE)
  if(class(datatry)[1] == "try-error" | length(K)<25 | NA%in%K){
    stop("Kinship matrix error, please cheak your kinship matrix.", call. = FALSE)
  }

  np <- nrow(K)
  p0 <- 1:np

  if(size > np | size < 4){
    stop("'size' argument error, must be 4 or more and less than number of the candidate individuals.", call. = FALSE)
  }

  if(length(keep) > 0){
    if(max(keep) > np){
      stop("'keep' argument error, element in keep vector not in the candidate set.", call. = FALSE)
    }
  }

  if(mut > size | mut < 1){
    stop("'mut' argument error, must be 1 or more and less than number of the candidate individuals.", call. = FALSE)
  }

  if(n0 < 4){
    stop("'n0' argument error, must be 4 or more.", call. = FALSE)
  }

  p <- p0
  if(length(keep) > 0){p <- p0[-keep]}
  if(!console[1] %in% c(0,1) | length(console) > 1){console <- FALSE}

  ans <- matrix(0, n0, size)
  for(i in 1:n0){
    ans[i,] <- c(keep, sort(sample(p, size-length(keep), replace = FALSE)))
  }
  time <- 1
  max <- -1e-10
  improve <- 1e+10
  improvetime <- 0

  if(console){
    cat("time", "max D-score\t", "improvement time", "\n", sep = "\t")
  }
  while(time < 1e+06 & improvetime < cri){
    d <- c()
    for(i in 1:n0){
      d[i] <- det(K[ans[i,], ans[i,]])
    }
    ans <- cbind(ans, d)
    ans <- ans[order(ans[, (size+1)]),]
    improve <- ans[n0, (size+1)]-max
    if(abs(improve) < 1e-05){improvetime <- improvetime+1
    } else {improvetime <- 0}
    max <- ans[n0, (size+1)]
    result <- ans[n0,]

    # selection step
    select1 <- floor(9*n0/10)
    select2 <- max-ans[1:select1, (size+1)]
    if(0%in%select2){select2 <- rep(1, select1)}
    select2 <- select2/sum(select2)
    select3 <- sample(1:select1, 2, prob=select2)
    ans <- ans[-select3, -(size+1)]
    ans.c <- ans
    if(class(ans.c)[1] != "matrix"){ans.c <- matrix(ans.c, 1, length(ans.c))}
    if(length(keep) > 0){ans.c <- ans.c[, -(1:length(keep))]}


    # crossover step
    cross1 <- sample(p, 1)
    cross2 <- sample(1:(n0-2), 2, replace = FALSE)
    cross3 <- ans.c[cross2[1], ans.c[cross2[1],] < cross1]
    cross4 <- ans.c[cross2[2], !ans.c[cross2[2],] < cross1]
    cross5 <- c(cross3, cross4)
    if(length(cross5) < (size-length(keep))){
      cross6 <- sample(p[!p%in%cross5], (size-length(cross5)-length(keep)))
      cross5 <- sort(c(cross5, cross6))
    } else if (length(cross5) > (size-length(keep))){
      cross5 <- sort(sample(cross5, (size-length(keep))))
    }
    ans.c <- rbind(cross5, ans.c)

    # mutation step
    for(i in 1:(n0-1)){
      mutation1 <- ans.c[i,]
      mutation2 <- sample(1:(size-length(keep)), mut, replace = FALSE)
      mutation1 <- mutation1[-mutation2]
      mutation3 <- sample(p[!p%in%mutation1], mut)
      mutation4 <- sort(c(mutation1, mutation3))
      if(i < (n0-1)){ans.c[i,] <- mutation4
      } else {ans.c <- rbind(ans.c, mutation4)}
    }
    ans <- ans.c
    if(length(keep) > 0){
      ans <- cbind(matrix(keep, nrow(ans.c), length(keep), byrow = TRUE), ans.c)}

    time <- time+1
    t0 <- 1000
    if(cri<2000){t0 <- 100}
    if(console){
      if(time%%t0 == 0 | improvetime == cri){
        cat(as.character(c(time, max, improvetime)), "\n", sep = "\t")}
    }
  }
  D.score <- result[(size+1)]
  subset <- result[1:size]
  return(list(subset = subset, D.score = D.score, time = time))
}

