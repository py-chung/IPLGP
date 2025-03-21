#' QTL Hotspot
#'
#' This function generates both numerical and graphical summaries of the QTL
#' hotspot detection in the genomes, including information about the flanking
#' markers of QTLs.
#'
#' @param DataQTL data.frame. A data frame with 5 columns for QTL information.
#' The columns represent the serial number of QTLs, the trait names, the
#' chromosome numbers, the left flanking marker positions(in cM) of QTLs, and
#' the right flanking marker positions(in cM) of QTLs.
#' @param DataCrop data.frame. A data frame with 3 columns for chromosome
#' information consists of the chromosome names, the center positions(in cM)
#' and the lengths of chromosomes.
#' @param ScanStep numeric. A value for the length(cM) of every bin.
#' @param NH integer. A value for the number of spurious hotspots in the
#' proposed method.
#' @param NP integer. A value for permutation times to calculate the
#' threshold.
#' @param save.pdf logical. When set to TRUE, the PDF file of plots will be
#' saved in the working directory instead of being displayed in the console.
#'
#' @return
#' \item{EQF}{The expected QTL frequency(EQF) in every bin per chromosome.}
#' \item{P.threshold}{The EQF thresholds for proposed method.}
#' \item{Q.threshold}{The EQF thresholds for the Q method.}
#' \item{nHot}{The numbers of detected hotspots per chromosome for the proposed
#' method and Q method.}
#'
#' Graphical outputs for visualizing the summarized results includes the
#' expected QTL frequency of scan steps, and the composition of QTLs for
#' different traits in the detected hotspots.
#'
#' @note
#' This program may generate a large amount of graphical output. To manage this,
#' it's recommended to save the results in a PDF file using the "save.pdf"
#' argument.
#'
#' @export
#'
#' @references
#'
#' Wu, P.-Y., M.-.H. Yang, and C.-H. KAO 2021 A Statistical Framework for QTL
#' Hotspot Detection. G3: Genes, Genomes, Genetics: jkab056. <doi: 10.1093/g3journal/jkab056>
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "QHOTexample.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' result <- Qhot(QTL.example, crop.example, 5, 20, 100, save.pdf = FALSE)
Qhot <- function(DataQTL, DataCrop, ScanStep = 1, NH = 100, NP = 1000, save.pdf = TRUE){
  DataQTL <- data.frame(DataQTL)
  if(ncol(DataQTL) != 5 | TRUE %in% is.na(DataQTL)){
    stop("Input data DataQTL error, please check and fix.", call. = FALSE)
  }
  colnames(DataQTL) <- c("X", "Trait", "chr", "L", "R")

  if(!is.numeric(DataQTL[, 4]) | !is.numeric(DataQTL[, 5])){
    stop("Input data DataQTL error, please check and fix.", call. = FALSE)
  }

  if(length(table(DataQTL[, 1])) != nrow(DataQTL)){
    DataQTL <- DataQTL[order(DataQTL[, 3], DataQTL[, 4]),]
    DataQTL[, 1] <- 1:nrow(DataQTL)
  }

  DataCrop <- data.frame(DataCrop)
  if(ncol(DataCrop) != 3 | TRUE %in% is.na(DataCrop)){
    stop("Input data DataCrop error, please check and fix.", call. = FALSE)
  }
  colnames(DataCrop) <- c("CHR", "Center.cM.", "Length.cM.")

  if(!is.numeric(DataCrop[, 2]) | !is.numeric(DataCrop[, 3]) | length(table(DataCrop[, 1])) != nrow(DataCrop)){
    stop("Input data DataCrop error, please check and fix.", call. = FALSE)
  }

  croptest <- TRUE %in% (!names(table(DataQTL[, 3])) %in% DataCrop[, 1])
  for(i in 1:nrow(DataCrop)){
    datatest <- DataQTL[DataQTL[, 3] == DataCrop[i, 1],]
    if(nrow(datatest) > 0){
      croptest <- c(croptest, DataCrop[i, 3] < max(datatest[, 5]))
    }
  }
  if(TRUE %in% croptest){
    stop("DataCrop does not match DataQTL, please check and fix.", call. = FALSE)
  }

  if(!is.numeric(ScanStep) | length(ScanStep) > 1 | min(ScanStep) < 0){
    stop("Parameter ScanStep error, please input a positive number.", call. = FALSE)
  }

  if(!is.numeric(NH) | length(NH) > 1 | min(NH) < 0){
    stop("Parameter NH error, please input a positive integer.", call. = FALSE)
  }

  if(!is.numeric(NP) | length(NP) > 1 | min(NP) < 0){
    stop("Parameter NP error, please input a positive integer.", call. = FALSE)
  }

  if(!save.pdf[1] %in% c(0,1) | length(save.pdf) > 1){save.pdf <- TRUE}


  if(save.pdf){grDevices::pdf("QHOTplot.pdf")}
  Data.freq <- function(DataQTL, ScanStep){

    Len.DataQTL <- nrow(DataQTL)
    Data.h <- c()
    Na.h <- c()

    X.max <- ceiling(max(DataQTL)/ScanStep)*ScanStep

    Break <- seq(0, X.max, by = ScanStep)

    for(j in 1:Len.DataQTL){

      M1 <- DataQTL[j, 1]
      M2 <- DataQTL[j, 2]
      Break1 <- Break[max(which(Break <= M1)):min(which(Break >= M2))]


      if(length(Break1) <= 2){
        Data.h <- c(Data.h, 1)
        Na.h <- c(Na.h, Break1[1])
      } else {
        Temp.break <- Break[min(which(Break > M1)):max(which(Break < M2))]
        Len.QTL <- diff(c(M1, M2))
        Data.h <- c(Data.h, diff(sort(c(Temp.break, M1, M2)))/Len.QTL)
        Na.h <- c(Na.h, Break1[1:length(Break1)-1])
      }

    }
    Data.h <- c(Data.h,rep(0, times = length(Break)-1))
    Na.h <- c(Na.h, Break[1:length(Break)-1])
    Freq.h <- tapply(Data.h, Na.h, sum)
    Break.h1 <- levels(cut(1, breaks = Break, right = FALSE))
    data.frame(BreakPoint = Break.h1, Left = Break[1:length(Break)-1],
               Right = Break[2:length(Break)], ExpectedNumber = Freq.h)
  }

  PHot <- function(Data, NP, NH){
    L <- dim(Data)[1]
    Trait.no <- dim(Data)[2]
    t1 <- Sys.time()
    TV.p <- c()
    ps <- c()
    Data2 <- Data
    for(k in 1:NP){
      for(j in 1:Trait.no){
        Index <- sample(1:L, L)
        Data2[, j] <- Data[Index, j]
      }
      tmp <- sort(apply(Data2, 1, sum), decreasing = TRUE)
      ps<-c(ps, tmp)
    }
    ps2 <- matrix(ps, byrow = TRUE, nrow = NP)
    TV.p <- apply(ps2, 2, sort)[round(NP*0.95), 1:NH]

    t2 <- Sys.time()
    Output <- list(time = t2-t1, TV.p = TV.p)
    invisible(Output)

  }

  NoQTL <- table(DataQTL[, 3])
  Na.chr <- names(NoQTL)

  Na.chr.crop <- factor(DataCrop[, 1])
  Index2 <- Na.chr.crop%in%Na.chr
  DataCrop2 <- DataCrop[Index2,]
  Na.chr.crop2 <- Na.chr.crop[Index2]
  Index3 <- order(Na.chr.crop2)
  DataCrop3 <- DataCrop2[Index3,]
  ChrL <- DataCrop3[, 3]

  Na.trait <- levels(factor(as.vector(DataQTL[, 2])))
  N.trait <- length(Na.trait)

  DataF <- vector("list", length(NoQTL))
  DataF.P <- c()
  DataF.Q <- c()
  L.ALL <- ceiling(ChrL/ScanStep)
  L.ALL2 <- sum(L.ALL)


  for(k in 1:length(NoQTL)){

    DataQTL2 <- DataQTL[DataQTL[, 3] == Na.chr[k],4:5]
    DataF[[k]] <- Data.freq(DataQTL2, ScanStep)
    L.break <- dim(DataF[[k]])[1]
    for(t in 1:N.trait){

      DataQTL.t <- DataQTL[DataQTL[, 3] == Na.chr[k] & DataQTL[, 2] == Na.trait[t], 4:5]
      if(nrow(DataQTL.t) != 0){
        Temp <- Data.freq(DataQTL.t, ScanStep)$ExpectedNumber
        L.diff <- L.break-length(Temp)
        Temp2 <- c(Temp, rep(0, L.diff))
        } else {Temp2 <- 0}
      DataF[[k]] <- data.frame(DataF[[k]], Temp2)
      names(DataF[[k]])[ncol(DataF[[k]])] <- Na.trait[t]

    }

    for(tt in 1:dim(DataQTL2)[1]){
      Temp <- Data.freq(DataQTL2[tt,], ScanStep)$ExpectedNumber
      DataF.Q <- c(DataF.Q, Temp[Temp > 0])
    }


    DataF.P <- rbind(DataF.P, DataF[[k]][, 4:(4+N.trait)])
  }


  DataF.tmp <- data.frame(matrix(0, nrow = (L.ALL2-dim(DataF.P)[1]), ncol = (N.trait+1)))
  names(DataF.tmp) <- names(DataF.P)
  DataF.P2 <- rbind(DataF.P, DataF.tmp)

  Output <- PHot(DataF.P2[, 2:(N.trait+1)], NP, NH)
  TV.P2 <- Output$TV.p

  TV.Q <- c()
  for(i in 1:NP){
    Index.Q <- sample(1:L.ALL2, length(DataF.Q), replace = TRUE)
    Index.Q2 <- c(Index.Q, 1:L.ALL2)
    DataF.Q2 <- c(DataF.Q, rep(0,L.ALL2))
    TV.Q <- c(TV.Q, max(tapply(DataF.Q2, Index.Q2, sum)))
  }
  TV.Q2 <- sort(TV.Q)[ceiling(NP*0.95)]


  DataHotspot <- vector("list", length(NoQTL))
  TableH <- c()
  #i=1
  for(i in 1:length(NoQTL)){

    Data.Freq <- DataF[[i]][, 4]
    Data.Left <- DataF[[i]][, 2]
    DataHotspot[[i]] <- vector("list", (NH+1))
    for(nh in 1:NH){
      Index.P <- Data.Freq > TV.P2[nh]
      DataHotspot[[i]][[nh]] <- Data.Left[Index.P]
      TableH <- c(TableH,sum(Index.P))
    }
    Index.Q <- Data.Freq > TV.Q2
    DataHotspot[[i]][[(NH+1)]] <- Data.Left[Index.Q]
    TableH <- c(TableH, sum(Index.Q))
  }

  TableH2 <- matrix(TableH, byrow = FALSE, nrow = (NH+1))

  Ymax <- max(DataF.P2[, 1])

  x1 <- max(ChrL)*(-0.05)
  x2 <- max(ChrL)*1.1
  y1 <- max(Ymax)*(-0.4)
  y2 <- max(Ymax)*1.1


  for(i in 1:length(NoQTL)){

    plot(0, 1, type = "n", bty = 'n', xaxt = "n", yaxt = "n",
         xlim = c(x1, x2), ylim = c(y1, y2), main = " ", xlab = " ", ylab = " ")

    yw <- max(Ymax)*0.2

    xb <- 0
    yb <- -yw
    xt <- ChrL[i]
    yt <- 0
    graphics::rect(xb, yb, xt, yt, col = "gray40", border = 1)


    Hotspot <- DataHotspot[[i]]
    Y1 <- (seq(yb*0.4, yb*0.6, length = 2))
    color <- c("Navy", "red")
    for(k in 1:2){
      Hot.k <- NH:(NH+1)
      X1 <- Hotspot[[Hot.k[k]]]
      if(length(X1) != 0){
        X2 <- X1+0.5*ScanStep
        Y2 <- Y1[k]
        graphics::points(X2, rep(Y2, length(X2)), pch = 20, col = color[k])
      }
    }


    Na <- paste("Chr", Na.chr[i], sep = "")
    graphics::text(0, (yt+yb)*0.5, Na, pos = 2, col = 1, cex = 1)

    ChrC <- DataCrop3[, 2]
    if(is.na(ChrC[i]) != TRUE){
      graphics::points(ChrC[i], yb*1.1, col = 1, pch = 17)}

    by.x <- ceiling(diff(seq(0, xt, length = 20))[1]/5)*5
    graphics::axis(1, seq(0, xt, by = by.x), labels = seq(0, xt, by = by.x), cex.axis = 0.75,
         las = 0, tck = -0.005, pos = yb, col = 1, col.axis = "gray40", mgp = c(3, 0.005, 0))

    MarkerPos <- c(as.matrix(DataQTL[DataQTL[, 3] == Na.chr[i], 4:5]))
    graphics::axis(1, MarkerPos, labels = FALSE, cex.axis = 0.6, las = 0, tck = -0.02,
         pos = 0, col = "green", col.axis = 1, mgp = c(3, 0.05, 0))


    ylim <- ceiling(max(Ymax)/5)*5
    by.y <- ceiling(diff(seq(0, ylim,length = 20))[1]/5)*5

    graphics::axis(2, seq(0, ylim, by = by.y), labels = seq(0, ylim, by = by.y), cex.axis = 0.7,
         las = 1, tck = -0.005, pos = 0, col = "gray40", col.axis = "gray40", mgp = c(3, 0.5, 0))

    by.y.qtl <- seq(y2*0.05, y2*0.95, length=NoQTL[i])

    count.t <- 0
    for(j in 1:N.trait){
      DataQTL2 <- DataQTL[DataQTL[, 3] == Na.chr[i] & DataQTL[, 2] == Na.trait[j], 4:5]
      if(nrow(DataQTL2) > 0){
        Index.L <- order(DataQTL2[, 1])
        DataQTL3 <- DataQTL2[Index.L,]
        No.trait <- nrow(DataQTL3)
        count.t <- count.t+No.trait
        by.y2 <- by.y.qtl[(count.t-No.trait+1):count.t]
        L3 <- DataQTL3[, 1]
        R3 <- DataQTL3[, 2]
        graphics::segments(L3, by.y2, R3,by.y2, col = grDevices::rainbow(N.trait, start = 0.1, end = 0.9, v = 0.8)[j], lty = j)
        graphics::points(L3,by.y2, pch = j, col = grDevices::rainbow(N.trait, start = 0.1, end = 0.9, v = 0.8)[j], lty = j)
        graphics::points(R3,by.y2, pch = j, col = grDevices::rainbow(N.trait, start = 0.1, end = 0.9, v = 0.8)[j], lty = j)
      }
    }

    DataQTL.chr <- factor(as.vector(DataQTL[DataQTL[, 3] == Na.chr[i], 2]), levels = Na.trait)
    DataQTL.chr2 <- table(DataQTL.chr)
    DataQTL.chr3 <- round(DataQTL.chr2/length(DataQTL.chr), 3)*100

    by.y.chr <- seq(yt, y2, length = (N.trait+2))
    y.pt <- by.y.chr[2:(N.trait+1)]
    x.pt <- rep(xt+(0.02*x2), N.trait)

    graphics::points(x.pt, y.pt, pch = 1:N.trait,
           col = grDevices::rainbow(N.trait, start = 0.1, end = 0.9, v = 0.8)[1:N.trait])

    Na.x <- paste(DataQTL.chr2, " ( ", DataQTL.chr3, "% )", sep = "")
    graphics::axis(4, y.pt, labels = Na.x, cex.axis = 0.8, las = 1, tck = -0.005,
         pos = xt+(0.05*x2), col = graphics::par()$bg, col.axis = "navy", mgp = c(3, 0.005, 0))

    graphics::axis(4, yt,labels = "Total", cex.axis = 1, las= 1, tck = -0.005,
         pos = xt+(0.05*x2), col = graphics::par()$bg, col.axis = "navy", mgp = c(3, 0.005, 0))

    Na.x.total <- paste(NoQTL[i], " ( ",100,"% )", sep = "")
    graphics::axis(4, 0.5*(yb+yt), labels = Na.x.total, cex.axis = 0.75, las = 1, tck = -0.005,
         pos = xt+(0.05*x2), col = graphics::par()$bg, col.axis = "navy", mgp = c(3, 0.005, 0))

    Left <- DataF[[i]][, 2]
    Freq <- DataF[[i]][, 4]

    if((max(Left)+ScanStep) < xt){
      x.last <- (max(Left)+ScanStep)
    } else {x.last <- xt}
    graphics::lines(c(Left, x.last), c(Freq, 0), type = "s", lwd = 1, col = "gray20")

    Na.submain <- paste("All flanking markers of QTLs in the chromosome ", Na.chr[i], sep = "")
    graphics::axis(1,0.5*ChrL[i], labels = Na.submain, cex.axis = 1, las = 0, tck = -0.005,
         pos = 0.5*(yb+y1), col = graphics::par()$bg, col.axis = 1, mgp = c(3, 0.005, 0))

    Hot.max <- sort(TableH2[, i], decreasing = TRUE)[1]

    if(Hot.max > 0){
      Index.H <- order(TableH2[, i], decreasing = TRUE)[1]
      Pos.hot <- DataHotspot[[i]][[Index.H]]

      for(nHot in 1:Hot.max){

        plot(0,1, type = "n", bty = 'n', xaxt = "n", yaxt = "n",
             xlim = c(x1, x2), ylim = c(y1, y2), main = " ", xlab = " ", ylab = " ")

        yw <- max(Ymax)*0.2

        xb <- 0
        yb <- -yw
        xt <- ChrL[i]
        yt <- 0
        graphics::rect(xb, yb, xt, yt, col = "gray40", border = 1)

        graphics::rect(Pos.hot[nHot], yb, Pos.hot[nHot]+ScanStep, yt, col = "gray50", border = "gray80")
        graphics::rect(Pos.hot[nHot], yt, Pos.hot[nHot]+ScanStep, y2, col = "gray99", border = "gray95")

        Hotspot <- DataHotspot[[i]]
        Y1 <- (seq(yb*0.4, yb*0.6, length = 2))
        color <- c("navy", "red")
        for(k in 1:2){
          Hot.k <- NH:(NH+1)
          X1 <- Hotspot[[Hot.k[k]]]
          if(length(X1) != 0){
            X2 <- X1+0.5*ScanStep
            Y2 <- Y1[k]
            graphics::points(X2, rep(Y2, length(X2)), pch = 20, col = color[k])
          }
        }

        Na <- paste("Chr", Na.chr[i], sep = "")
        graphics::text(0, (yt+yb)*0.5, Na, pos = 2, col = 1, cex = 1)
        ChrC <- DataCrop3[, 2]
        if(is.na(ChrC[i]) != TRUE){
          graphics::points(ChrC[i], yb*1.1, col = 1, pch = 17)}
        by.x <- ceiling(diff(seq(0, xt, length = 20))[1]/5)*5
        graphics::axis(1, seq(0, xt, by = by.x), labels = seq(0, xt, by = by.x), cex.axis = 0.75,
             las = 0, tck = -0.005, pos = yb, col = 1, col.axis = "gray40", mgp = c(3, 0.005, 0))
        graphics::segments(0, 0, xt, 0, col = 1)
        ylim <- ceiling(max(Ymax)/5)*5
        by.y <- ceiling(diff(seq(0, ylim, length = 20))[1]/5)*5
        graphics::axis(2, seq(0, ylim, by = by.y), labels = seq(0, ylim, by = by.y), cex.axis = 0.7,
             las = 1, tck = -0.005, pos = 0, col = "gray40", col.axis = "gray40", mgp = c(3, 0.5, 0))

        DataQTL.nHot <- DataQTL[DataQTL[, 3] == Na.chr[i], 4:5]
        Index.LH <- DataQTL.nHot[, 1] <= (Pos.hot[nHot]+ScanStep)
        Index.RH <- DataQTL.nHot[, 2] >= Pos.hot[nHot]
        Index.LRH <-Index.LH&Index.RH

        graphics::axis(4, yt, labels = "Total", cex.axis = 1, las = 1, tck = -0.005,
             pos = xt+(0.01*x2), col = graphics::par()$bg, col.axis = "navy", mgp = c(3, 0.005, 0))

        Na.x.total <- paste(sum(Index.LRH), " ( ", 100, "% )", sep = "")
        graphics::axis(4, 0.5*(yb+yt), labels = Na.x.total, cex.axis = 0.8, las = 1, tck = -0.005,
             pos = xt+(0.05*x2),col = graphics::par()$bg, col.axis = "navy", mgp = c(3, 0.005, 0))
        DataQTL.nHot2 <- factor(as.vector(DataQTL[DataQTL[, 3] == Na.chr[i], 2][Index.LRH]), levels = Na.trait)

        DataQTL.nHot3 <- table(DataQTL.nHot2)
        DataQTL.nHot4 <- round(DataQTL.nHot3/sum(DataQTL.nHot3), 3)*100
        graphics::points(x.pt, y.pt, pch = 1:N.trait,
               col = grDevices::rainbow(N.trait, start = 0.1, end = 0.9, v = 0.8)[1:N.trait])
        Na.x <- paste(DataQTL.nHot3, " ( ", DataQTL.nHot4, "% )", sep = "")
        graphics::axis(4, y.pt, labels = Na.x, cex.axis = 0.8, las = 1, tck = -0.005,
             pos = xt+(0.05*x2), col = graphics::par()$bg, col.axis = "navy", mgp = c(3, 0.005, 0))


        by.y.qtl <- seq(y2*0.05, y2*0.95, length = sum(Index.LRH))

        count.t <- 0
        for(j in 1:N.trait){

          DataQTL2 <- DataQTL[DataQTL[, 3] == Na.chr[i] & DataQTL[, 2] == Na.trait[j], 4:5]
          Index.LH <- DataQTL2[, 1]<=(Pos.hot[nHot]+ScanStep)
          Index.RH <- DataQTL2[, 2]>=Pos.hot[nHot]
          Index.LRH <- Index.LH&Index.RH
          if(sum(Index.LRH) > 0){
            DataQTL3 <- DataQTL2[Index.LRH,]
            Index.L <- order(DataQTL3[, 1])
            DataQTL4 <- DataQTL3[Index.L,]
            No.trait <- nrow(DataQTL4)
            count.t <- count.t+No.trait
            by.y2 <- by.y.qtl[(count.t-No.trait+1):count.t]
            L3 <- DataQTL3[, 1]
            R3 <- DataQTL3[, 2]
            graphics::segments(L3, by.y2, R3, by.y2, col = grDevices::rainbow(N.trait, start = 0.1, end = 0.9, v = 0.8)[j], lty = j)
            graphics::points(L3, by.y2, pch = j, col = grDevices::rainbow(N.trait, start = 0.1, end = 0.9, v = 0.8)[j], lty = j)
            graphics::points(R3, by.y2, pch = j, col = grDevices::rainbow(N.trait, start = 0.1, end = 0.9, v = 0.8)[j], lty = j)
          }
        }


        Left <- DataF[[i]][, 2]
        Freq <- DataF[[i]][, 4]

        if((max(Left)+ScanStep) < xt){
          x.last <- (max(Left)+ScanStep)
        } else {x.last <- xt}
        graphics::lines(c(Left, x.last), c(Freq, 0), type = "s", lwd = 1, col = "gray20")

        Index.hh <- DataF[[i]]$Left == Pos.hot[nHot]

        Na.submain <- paste(DataF[[i]]$BreakPoint[Index.hh], "Expected QTL frequency =",
                          round(DataF[[i]]$ExpectedNumber[Index.hh], 2), sep = " ")
        graphics::axis(1, 0.5*ChrL[i],labels = Na.submain, cex.axis = 1, las = 0, tck = -0.005,
             pos = 0.5*(yb+y1), col = graphics::par()$bg, col.axis = 1, mgp = c(3, 0.005, 0))
      }
    }
  }

  graphics::layout(matrix(c(1, 2), 1, 2, byrow = TRUE), width = c(1, 1))

  plot(0, 1, type = "n", bty = 'n', xaxt = "n", yaxt = "n", main = " ", xlab = " ", ylab = " ")

  Na.legd <- c(paste("Proposed method n= ", NH, sep = ""), "Q method" )
  graphics::legend("left", Na.legd, pch = 20, col = color, bty = "n")

  plot(0, 1, type = "n", bty = 'n', xaxt = "n", yaxt = "n", main = " ", xlab = " ", ylab = " ")

  Na.legd <- graphics::legend("left", Na.trait, pch = 1:N.trait, lty=1:N.trait,
                  col = grDevices::rainbow(N.trait, start = 0.1, end = 0.9, v = 0.8)[1:N.trait], bty = "n")

  names(DataF) <- paste("Chr", Na.chr, sep = " ")
  names(TV.P2) <- paste("n= ", 1:NH, sep = "")
  names(TV.Q2) <- c("Q method")
  rownames(TableH2) <- c(paste("Proposed method n= ", 1:NH,sep = ""), "Q method")
  colnames(TableH2) <- paste("Chr.", Na.chr, sep = "")
  if(save.pdf){grDevices::dev.off()}

  return(list(EQF = DataF, P.threshold = TV.P2, Q.threshold = TV.Q2, nHot = TableH2))
}
