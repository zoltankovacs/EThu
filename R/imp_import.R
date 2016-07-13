findMaxNrElements <- function(fileNames, NrAdjust, sep = const_sep){
  maxNr <- max(unlist(lapply(strsplit(fileNames, sep), length)))
  maxNrInd <- which.max(unlist(lapply(strsplit(fileNames, sep), length)))
  maxNr <- maxNr -1 
  maxNr <- maxNr + NrAdjust
  return(list(maxNr=maxNr, maxNrInd=maxNrInd))
  } #Eof

createSingleLineList <- function(singlePath){
  pathSplit <- unlist(strsplit(singlePath, "/"))[-1]
  last <- pathSplit[length(pathSplit)]
  
  hashSplit <- unlist(strsplit(last, const_userPrefSep))
   if (length(hashSplit) == 1) { # if there is no '#' present
     configID <- substr(last, 1, nchar(last) - const_nrRestCharCsv)
     pathSplit <- pathSplit[-length(pathSplit)]
    } else {
      configID <- substr(hashSplit[2], 1, nchar(hashSplit[2]) - const_nrRestCharCsv)
      userInput <- unlist(strsplit(hashSplit[1], const_elemSep))
      pathSplit <- c(pathSplit[-length(pathSplit)], userInput)
    }
  pathSplit <- c(pathSplit, paste("configID", configID, sep = "@"))
  
  sepSplit <- strsplit(pathSplit, "@")
  
  ID <- unlist(lapply(sepSplit, function(x) x[1]))
  value <- unlist(lapply(sepSplit, function(x) x[2]))
  return(list(ID=ID, value=value))
} #Eof

reFactor <- function(df) {
  for (i in 1: ncol(df)) {
    if (!is.factor(df[,i])) {
      df[,i] <- as.factor(df[,i])
    }
  }
  return(df)
}#eof

createInfoTable <- function(fileNames, NrAdjust = 1){
  maxNrList <- findMaxNrElements(fileNames, NrAdjust)
  InfoTable <- as.data.frame(matrix(NA, length(fileNames), maxNrList$maxNr))
  colNameInfTable <- createSingleLineList(fileNames[maxNrList$maxNrInd])$ID
  colnames(InfoTable) <- colNameInfTable
  
  for (f in 1:length(fileNames)){
    a <- createSingleLineList(fileNames[f])
    Ind <- match(colNameInfTable, a$ID)
    InfoTable[f, ] <- a$value[Ind]
  }
  InfoTable <- reFactor(InfoTable)
  return(InfoTable)
} #Eof

askSep <- function(cns){
  nrOk <- FALSE
  tmp <- data.frame(classVariables = cns)
  while (!nrOk){
    cat("Please provide the index of the class variable to be used for phyical separation of the data, provide 0 for no separation\n\n")
    print(tmp)
    cat("\n\n")
    a <- readLines(n = 1)
    a <- as.numeric(a)
    if (all(is.numeric(a)) & length(a) == 1 & a <= length(cns)){
         nrOk = TRUE
    }
  } # E while
  return(a)
} #Eof

getHeader_old <- function(rawCsvFile, tz = "EST", yearCorrect = 2016){
  ColumnNames <- c("ScannerID", "type", "PixelWidth","SystemTemp", "DetectorTemp", "Humidity", "NrMeasPoint", "NrRep", "Gain", "MeasTime", "absTime")
  # infoTable <- data.frame(matrix(NA, nrow = 1, ncol = length(ColumnNames)))
  # colnames(infoTable) <- ColumnNames
  ScannerID <- data.frame(as.character(rawCsvFile[11,2]), stringsAsFactors = F)
  type <- data.frame(as.numeric(as.character(rawCsvFile[13, 2])), stringsAsFactors = FALSE)
  if (type == "1"){type <- "Hadamard"} else {type <- "Column"}
  type <- data.frame(type, stringsAsFactors = FALSE)
  PixelWidth <- data.frame(as.numeric(as.character(rawCsvFile[16, 2])))
  SystemTemp <-data.frame(as.numeric(as.character(rawCsvFile[3,2]))/100)
  DetectorTemp <- data.frame(as.numeric(as.character(rawCsvFile[4,2]))/100)
  Humidity <- data.frame(as.numeric(as.character(rawCsvFile[5,2]))/100)
  NrMeasPoint <- data.frame(as.numeric(as.character(rawCsvFile[17,2])))
  NrRep <- data.frame(as.numeric(as.character(rawCsvFile[18,2])))
  Gain <- data.frame(as.numeric(as.character(rawCsvFile[19,2])))
  MeasTime <- data.frame(as.numeric(as.character(rawCsvFile[20,2])))
  absTime <- as.POSIXct(strptime(as.character(rawCsvFile[1,2]), "%d/%m/%Y @ %H:%M:%S"), tz)
  if(!difftime(absTime, "1996-06-15 10:59:51 EST", tz = tz) > 3650){
    absTime <- absTime-difftime(cut.POSIXt(absTime, "year"), paste0(yearCorrect, "-01-01"))
  }
  absTime <- as.character(absTime)
  absTime <- data.frame(absTime, stringsAsFactors = FALSE)
  infoTable <- cbind(ScannerID, type, PixelWidth, SystemTemp, DetectorTemp, Humidity, NrMeasPoint, NrRep, Gain, MeasTime, absTime)
  colnames(infoTable) <- ColumnNames
  invisible(infoTable)
}

getHeader <- function(singleFilename, tz = "EST", yearCorrect = 2016){
  con <- file(singleFilename, open="r")
  char <- readLines(con, n=19)
  close(con)
  char <- unlist(lapply(strsplit(char, ","), function(x) x[2]))
  ScannerID <- char[10]
  type <- if (char[12] == "1"){type <- "Hadamard"} else {type <- "Column"}
  PixelWidth <- as.numeric(char[15])
  SystemTemp <- as.numeric(char[4])/100
  DetectorTemp <- as.numeric(char[5])/100
  Humidity <- as.numeric(char[6])/100
  NrMeasPoint <- as.numeric(char[16])
  NrRep <- as.numeric(char[17])
  Gain <- as.numeric(char[18])
  MeasTime <- as.numeric(char[19])
  absTime <- as.POSIXct(strptime(as.character(char[2]), "%d/%m/%Y @ %H:%M:%S"), tz)
  if(!difftime(absTime, "1996-06-15 10:59:51 EST", tz = tz) > 3650){
    absTime <- absTime-difftime(cut.POSIXt(absTime, "year"), paste0(yearCorrect, "-01-01"))
  }
  absTime <- as.character(absTime)
  #
  out<- data.frame(ScannerID, type, PixelWidth, SystemTemp, DetectorTemp, Humidity, NrMeasPoint, NrRep, Gain, MeasTime, absTime, stringsAsFactors = FALSE)
  return(out)
} #eof

getAbsRefSmplSpect_old <- function(rawCsvFile, wlsRnd = 4){
  wavelength <- as.numeric(as.character(rawCsvFile[-(1:21),1]))
  absorbance <- t(data.frame(as.numeric(as.character(rawCsvFile[-(1:21),2]))))
  ref <- t(data.frame(as.numeric(as.character(rawCsvFile[-(1:21),3]))))
  smpl <- t(data.frame(as.numeric(as.character(rawCsvFile[-(1:21),4]))))
  AbsRefSmplSpect <- rbind(absorbance, ref, smpl)
  rownames(AbsRefSmplSpect) <- c("absorbance", "reference", "sample")
  colnames(AbsRefSmplSpect) <- paste0("X", round(wavelength, wlsRnd))
  invisible(AbsRefSmplSpect)
} #eof

getAbsRefSmplSpect <- function(singleFilename){
  a <- read.csv(singleFilename, skip=19)
  nir <- t(a[,2:4])
  wls <- paste("X", a[,1], sep="")
  colnames(nir) <- wls
  rownames(nir) <- c("absorbance", "reference", "sample")
  return(nir)
} #eof

getNrOfCol <- function(singlePath, nrColAdd){
  tmp <- read.csv(singlePath)
  headerTable <- getHeader(tmp)
  AbsRefSmplSpect <- getAbsRefSmplSpect(tmp)
  Nr <- ncol(headerTable) + ncol(AbsRefSmplSpect) + nrColAdd
  colNames <- c(colnames(headerTable), colnames(AbsRefSmplSpect))
  return(list(Nr=Nr, colNames=colNames))
}#Eof

makeListLayout <- function(infoTable, ColInd_s, colInd_c) {
  uniq_s <- unique(infoTable[, ColInd_s])
  leScanner <- length(uniq_s)
  leConf <- levConf <-  NULL
  for (i in 1:leScanner){
    ind <- which(infoTable[,ColInd_s] == uniq_s[i])
    tmp <- infoTable[ind,]
    a <- levels(tmp[, colInd_c])
    b <- length(a)
    levConf <- c(levConf, a)
    leConf <- c(leConf, b)
  } # end for i
  
  return(list(leScanner=leScanner, uniqS=as.character(uniq_s), levConf = levConf, leConf=leConf))
} #eof

mainF <- function(scanIdCol = const_scanIdCol, confIdCol = const_confIdCol, NrAdjust = 1){
  fileNames <- list.files(const_rawdataFolder, full.names = TRUE, recursive = T, pattern = const_fileExtension)
  infoTable <- createInfoTable(fileNames, NrAdjust)
  print(length(fileNames)); wait()
  infoTableRed <- infoTable[, which(!colnames(infoTable) %in% c(scanIdCol, confIdCol))]
  sepID <- askSep(colnames(infoTable))
  Ind_s <- which(colnames(infoTable)  == scanIdCol)
  Ind_c <- which(colnames(infoTable)  == confIdCol)
  
  if (sepID == 0){
    uniq_s <- unique(infoTable[, Ind_s])
    # for (us in uniq_s){
    for (us in uniq_s[1]){
      Ind <- which(infoTable[, Ind_s] == us)
      selinfTable <- infoTable[Ind, ]
      uniq_c <- unique(selinfTable[, Ind_c])
      # for (uc in uniq_c){
      for (uc in uniq_c[1]){
        Ind2 <- which(selinfTable[, Ind_c] == uc)
        print(length(Ind2)); wait()
        a <- getNrOfCol(fileNames[Ind[Ind2]][1], ncol(infoTable)-2)
        nrOfCol <- a$Nr
        colNames <- c(colnames(infoTableRed), a$colNames)
        colDfAbs <- colDfRef <- colDfSmpl <- as.data.frame(matrix(NA, length(Ind2), nrOfCol))
        print(nrow(colDfAbs)); wait()
        
        colnames(colDfAbs) <- colnames(colDfRef) <- colnames(colDfSmpl) <- colNames
        fileNamesSelect <- fileNames[Ind[Ind2]]
        print(1:length(fileNamesSelect)); wait()
        for (i in 1:length(fileNamesSelect)){
        #  tmp <- read.csv(fileNamesSelect[i])
          Header   <- getHeader(fileNamesSelect[i])
          spectra <- getAbsRefSmplSpect(fileNamesSelect[i])
          addHeader <- infoTableRed[Ind[Ind2][i], ]
          colDfAbs[i,] <- cbind(addHeader, Header, t(spectra[1,]))
          colDfRef[i,] <- cbind(addHeader, Header, spectra[2,, drop=FALSE])
          colDfSmpl[i,] <- cbind(addHeader, Header, spectra[3,, drop=FALSE])
        } #Efor i
        # here should save the giveb confog's data
        # before turn posix and character to factor
      } #Efor uc
    } #Efor us
    
    
  } else {
    
  } #Eif
  return(list(colDfAbs=colDfAbs, colDfRef=colDfRef, colDfSmpl=colDfSmpl))
} #Eof

