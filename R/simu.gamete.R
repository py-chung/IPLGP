#' Simulate The Genotypic Values Of A Gamete
#'
#' Generate the genotypic values of a gamete from the genotypic data of its
#' parents by Monte Carlo simulation. The recombination rate is caculate by
#' Haldane’s mapping function.
#'
#' @param marker data frame. A p*4 data frame whose first column indicates the
#' chromosome number to which a marker belongs; second column indicates the
#' position of the marker in centi-Morgan (cM); and 3rd and 4th columns
#' indicates the genotype of the marker (numeric or character).
#'
#' @return
#' The SNP sequence of gamete.
#'
#' @export
#'
#' @references
#'
#' Haldane J.B.S. 1919. The combination of linkage values and the calculation
#' of distance between the loci for linked factors. Genetics 8: 299–309.
#'
#' @examples
#' # generate simulated data
#' marker.test <- data.frame(c(1,1,1,1,1,2,2,2,2,2),c(10,20,30,40,50,10,20,30,40,50),
#' c("A","T","C","G","A","A","G","A","T","A"),c("A","A","G","C","T","A","G","T","T","A"))
#'
#' # run
#' simu.gamete(marker.test)

simu.gamete <- function(marker){

  if(ncol(marker) != 4 | NA%in%marker){
    stop("Marker data error, please cheak your marker data.", call. = FALSE)
  }

  marker <- marker[order(marker[,1], marker[,2]),]
  distance <- marker[,2]

  datatry <- try(sum(distance), silent = TRUE)
  if(class(datatry)[1] == "try-error"){
    stop("Marker data error, please cheak your marker data.", call. = FALSE)
  }

  distance.c <- distance[-1]-distance[-length(distance)]
  distance.r <- (1-exp(-2*distance.c/100))/2
  distance.r[distance.r<0] <- 0.5

  index <- sample(3:4, 1)
  marker1 <- marker[,index]
  marker2 <- marker[,(7-index)]

  marker.new <- marker1[1]
  for(i in 2:nrow(marker)){
    x <- stats::runif(1)
    if(x < distance.r[i-1]){
      marker.new[i] <- marker2[i]
      D <- marker1
      marker1 <- marker2
      marker2 <- D
    } else { marker.new[i] <- marker1[i] }
  }

  mn <- table(marker[,1])
  mn1 <- c()
  for(i in 1:length(mn)){
    mn1 <- c(mn1, 1:mn[i])
  }

  names(marker.new) <- paste(marker[,1],mn1,sep = "_")
  return(marker.new)
}

