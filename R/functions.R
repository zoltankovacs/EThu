# descriptive statistic ###########################################################################################################################
doDescrStat <- function(Data, VarName, fileName = "", range = 1.5){
  Ind <- which(colnames(Data) == VarName)
  AnovaData <- data.frame(cbind(Data[,Ind], as.matrix(Data$sensors)))
  colnames(AnovaData)[1] <- "groups"
  model <- aov(groups ~ ., data = AnovaData)
  fileName <- paste0(VarName, fileName, "_AnovaTable.txt")
  PdfComment(c("Results of Anova saved as",
               fileName,
               "descriptive evaluation results:"))
  #write out the Anova table to txt
  AovModelSum <- summary(model)
  AovModelSumTable <- file(paste0("results/", fileName), "wt")
  sink(file=AovModelSumTable, append = TRUE, type = c("output"), split = FALSE)
  cat("\n\n\n")
  print(paste0("ANOVA results of ", VarName, ":"))
  print(levels(Data[,Ind]))
  cat("-----------------------------------------\n")
  print(AovModelSum)
  cat("-----------------------------------------\n")
  sink()
  close(AovModelSumTable)
  
  #a CA érzékeny a beoldási időre
  SignInd <- which(AovModelSum[[1]]$`Pr(>F)` <= 0.05)
  outAll<-NA
  for (i in 1:ncol(Data$sensors)){
    sensName <- colnames(Data$sensors)[i]
    main = paste0("Boxplot of ", VarName, " - sensor ", sensName)
    xlab = ""
    ylab = paste0(sensName, " intenzitas")
    if(any(SignInd == i)){
      col.sub = "red"
    }else{
      col.sub = "black" 
    }
    sub = paste0("p=", round(AovModelSum[[1]]$`Pr(>F)`[i], 4))
    out <- boxplot(Data$sensors[,i] ~ Data[,Ind], xlab=xlab, ylab=ylab, main=main, 
                   col=unique(Data$numRep[,Ind]), sub=sub, col.sub=col.sub, las=2, range = range)
    outAll <- c(outAll, list(sensor = out))
    names(outAll)[length(outAll)] <- colnames(Data$sensors)[i]
  }
  invisible(outAll[-1])
}#Eof func

doOutDetect <- function(Data, VarName, fileName = ""){
  Ind <- which(colnames(Data) == VarName)
  AnovaData <- data.frame(cbind(Data[,Ind], as.matrix(Data$sensors)))
  colnames(AnovaData)[1] <- "groups"
  model <- aov(groups ~ ., data = AnovaData)
  fileName <- paste0(VarName, fileName, "_AnovaTable.txt")
  
  AovModelSum <- summary(model)
  AovModelSumTable <- file(paste0("results/", fileName), "wt")
  close(AovModelSumTable)
  
  #a CA érzékeny a beoldási időre
  SignInd <- which(AovModelSum[[1]]$`Pr(>F)` <= 0.05)
  outAll<-NA
  for (i in 1:ncol(Data$sensors)){
    sensName <- colnames(Data$sensors)[i]
    main = paste0("Boxplot of ", VarName, " - sensor ", sensName)
    xlab = ""
    ylab = paste0(sensName, " intenzitas")
    if(any(SignInd == i)){
      col.sub = "red"
    }else{
      col.sub = "black" 
    }
    out <- boxplot(Data$sensors[,i] ~ Data[,Ind],plot=F)
    outAll <- c(outAll, list(sensor = out))
    names(outAll)[length(outAll)] <- colnames(Data$sensors)[i]
  }
  invisible(outAll[-1])
}#Eof func


doSensorBars <- function(Data, VarName, wherLeg = "topright"){
  Ind <- which(colnames(Data) == VarName)
  
  Ave <- ddply(data.frame(Data$sensors), .(Data[,Ind]), colwise(mean))
  SDs <- ddply(data.frame(Data$sensors), .(Data[,Ind]), colwise(sd))
  
  AveP <- Ave[,-1]+SDs[,-1]
  AveN <- Ave[,-1]-SDs[,-1]
  
  
  Nc = ncol(Ave)-1
  Nr = nrow(Ave)
  xlim=c(0,Nc*Nr)
  ylim=range(rbind(AveP, AveN))
  for (i in 1:Nr){
    if(!i==1){par(new=T)}
    col=unique(Data$numRep[,Ind])[i]
    plot(seq(i, xlim[2], Nr), Ave[i,-1], type="h", xlim=xlim, ylim=ylim, lwd=4, col=col, axes=F, xlab="", ylab="")
    par(new=T)
    for (j in 1:ncol(Ave)){
      points(rep(seq(i, xlim[2], Nr)[j],2), c(AveP[i,j], AveN[i,j]), pch='-', type="l", lwd=2, col="gray70")
    }
    points(seq(i, xlim[2], Nr), AveP[i,], pch='-', cex=1.5, lwd=3, col="gray70")
    points(seq(i, xlim[2], Nr), AveN[i,], pch='-', cex=1.5, lwd=3, col="gray70")
    par(new=F)
  }
  abline(h=0, lty=2, col="gray70")
  at = colMeans(rbind(seq(1, xlim[2], Nr), seq(i, xlim[2], Nr)))
  axis(1, at, labels = colnames(Data$sensors), las=2)
  axis(2)
  main=paste0("Average and +- SD of sensors' intensity by ", VarName)
  title(main=main, xlab="sensors", ylab = "intensity")
  legend(wherLeg, col=unique(Data$numRep[,Ind]), lty=1, legend=unique(Data[,Ind]), bty = "n")
}

DetOutsFromBoxplotList <- function(doDescrStat_List, Data){
  AllOuts <- NA
  for(i in 1: length(doDescrStat_List)){
    Outs <- doDescrStat_List[[i]]$out
    if (length(Outs)>0){
      AllOuts <- c(AllOuts, which(Data$sensors[,i] %in% Outs))
    }
  }
  invisible(unique(sort(AllOuts[-1])))
}

# euclidean distance ########################################
doEu_dsitance <- function(DataAll){
  res<-as.data.frame(matrix(0,nrow=length(unique(DataAll$SmplName)),
                            ncol=length(unique(DataAll$SmplName))))
  names(res)<-unique(DataAll$SmplName)
  rownames(res)<-unique(DataAll$SmplName)
  for (i in unique(DataAll$SmplName) ) # do by sample
  {
    for (j in unique(DataAll$SmplName) ) # do by sample
    {
      temp<-0
      for (z in names(DataAll$sensors)) # by sensor
      {
        tmp<-t.test(DataAll$sensors[DataAll$SmplName==i,z] , DataAll$sensors[DataAll$SmplName==j,z])
        if (tmp$p.value <= 0.05) temp<-temp + diff(tmp$estimate)^2
      } 
      res[i,j]<-sqrt(temp)
    }
  }
res
  }

plotEuDistances <- function(DataAll){
  layout(matrix(c(1,2),ncol=2,nrow=1),widths=c(8,2),heights=8)
  par(mar=c(8,8,1,1),bg='white')
  res<-11 # a skála felosztása
  DIFF<-doEu_dsitance(DataAll)
  
  image(as.matrix(DIFF),xaxt='n',yaxt='n',
        col=grey.colors(10*res,start=1,end=0),cex.main=0.7,
        main="")
  axis(1,at=seq(0,1,length.out=length(unique(DataAll$SmplName))),labels=unique(DataAll$SmplName),las=2)  
  axis(2,at=seq(0,1,length.out=length(unique(DataAll$SmplName))),labels=unique(DataAll$SmplName),las=2)  
  #távolságok
  for (i in 1:dim(DIFF)[1])
  {
    for (j in 1:dim(DIFF)[2])
    {
      if (i<=j) next
      text(seq(0,1,length.out=dim(DIFF)[1])[i],
           seq(0,1,length.out=dim(DIFF)[1])[j],round(DIFF[i,j]),font=2,col=round(1-DIFF[i,j]/max(DIFF)))
    }
  }
  par(mar=c(3,3,3,1))
  image(matrix(seq(min(DIFF,na.rm=T),max(DIFF,na.rm=T),length.out=res*10),ncol=res*10,nrow=1),
        col=grey.colors(10*res,start=1,end=0),xlab="",ylab="",xaxt="n",yaxt="n")
  axis(2,at=seq(0,1,length.out=res),
       labels=substr(as.character(seq(min(DIFF,na.rm=T),ceiling(max(DIFF,na.rm=T)/10)*10,length.out=res)),1,7),
       cex.axis=1,line=0,las=2)
  par(mar=c(5.1, 4.1, 4.1, 2.1),bg='transparent') # reset
layout(1)
  }


# descriptive statistic END #######################################################################################################################

# data corrections ################################################################################################################################
flying_correctorD <- function(DataAll){
  for (i in unique(DataAll$repeats)){
    Ind <- which(DataAll$repeats == i)
    tmp <- DataAll$sensors[Ind,]
    tmpCorr <- t(apply(tmp, 1, function(x) x - colMeans(tmp)))
    DataAll$sensors[Ind,] <- tmpCorr
  }# E for i
  out <- DataAll
}# E func
# data corrections END ############################################################################################################################




### PCA ###########################################################################################################################################
# calc PCA model
#' @title calc PCA model
#' @description calc PCA model
#' @details XXX Here the details of how the folder should be named, with 
#' separators etc.
#' @param rawData data.frame having $sensors
#' @param day optional argument useful if the experiment was performed in different days
#' @param nameChange optional argument if the name of the groups has to be renamed
#' @return a table with additional columns
#' @export
calcPCA <- function(rawData, center = TRUE, scale. = FALSE){
  matrix.pca <- prcomp(rawData$sensors, center = center, scale. = scale.)
  matrix.pc <- predict(matrix.pca)          
  loadings <- matrix.pca$rotation
  out <- list(PCAmodel = matrix.pca, scores = matrix.pc, loadings = loadings)
}

# get out the PCA variances and SDs
PCAvariance <- function(PCAmodel){
  Nrs <- 1:ncol(summary(PCAmodel$PCAmodel)$importance)
  Standard_deviation <- paste0("PC", Nrs,": SD=",format(summary(PCAmodel$PCAmodel)$importance[1,], digits = 2)," ","(", summary(PCAmodel$PCAmodel)$importance[2,]*100,"%", ")", sep="")
  Proportion_of_Variance <- paste0("PC", Nrs, " - ", summary(PCAmodel$PCAmodel)$importance[2,]*100,"%", sep="")
  summaryTable <- rbind(Standard_deviation, Proportion_of_Variance)
  colnames(summaryTable) <- colnames(summary(PCAmodel$PCAmodel)$importance)
  summaryTable
}

plotScreeplot <- function(PCAmodel, main = "", NrPCs = 3, PCAvariance = NA){
  screeplot(PCAmodel, type = "lines", main = main)
  if (!any(is.na(PCAvariance))){
    legend("topright", cex = 1, xjust = 1, yjust = 1, legend = PCAvariance[1, 1:NrPCs], box.col = "black")
  }
}

plotPcaScores <- function(PCAmodel, PCs = 1:2, ...){
  xlab <- PCAvariance(pcaMod)[2, PCs[1]]
  ylab <- PCAvariance(pcaMod)[2, PCs[2]]
  plot(PCAmodel$scores[,PCs], xlab = xlab, ylab = ylab, ...)
}

### PCA END #######################################################################################################################################

### LDA ###########################################################################################################################################
getLdaVars <- function(LDAmodel, round = 2){
  Vars <- round((LDAmodel$svd)^2/sum((LDAmodel$svd)^2)*100, round)
}
### LDA END #######################################################################################################################################

PdfComment <- function(comment = "HELLO", xjust=0.5, ...){
  plot(1:100, type="n", axes=F, xlab = "", ylab ="")
  legend("center", legend=comment, bty = "n", xjust=xjust, ...)
}

