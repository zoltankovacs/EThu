findMaxNrElements <- function(fileNames, NrAdjust, sep = const_sep) {
  maxNr <- max(unlist(lapply(strsplit(fileNames, sep), length)))
  maxNrInd <- which.max(unlist(lapply(strsplit(fileNames, sep), length)))
  maxNr <- maxNr -1 
  maxNr <- maxNr + NrAdjust
  return(list(maxNr=maxNr, maxNrInd=maxNrInd))
} #Eof

createSingleLineList <- function(singlePath) {
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
  ID <- c(ID, const_clNameConsec)
  value <- unlist(lapply(sepSplit, function(x) x[2]))
  value <- c(value, NA)
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

reCharacter <- function(df) {
  for (i in 1: ncol(df)) {
    if (!is.character(df[,i]) & !is.numeric(df[,i])) {
      df[,i] <- as.character(df[,i])
    }
  }
  return(df)
}#eof

createInfoTable <- function(fileNames, scanIdCol = const_scanIdCol, confIdCol = const_confIdCol, NrAdjust = 1) {
  maxNrList <- findMaxNrElements(fileNames, NrAdjust)
  infoTable <- as.data.frame(matrix(NA, length(fileNames), maxNrList$maxNr + 1))# to add 1 for the consecutives
  colNameInfTable <- createSingleLineList(fileNames[maxNrList$maxNrInd])$ID
  colnames(infoTable) <- colNameInfTable
  for (f in 1:length(fileNames)) {
    a <- createSingleLineList(fileNames[f])
    Ind <- match(colNameInfTable, a$ID)
    infoTable[f, ] <- a$value[Ind]
  }
  infoTable <- reFactor(infoTable)
  allChar <- apply(infoTable[,-ncol(infoTable)], 1, function(x) paste0(x, collapse="-"))
  a <- rle(allChar)
  grouping <- a$lengths
  consecsAll <- NULL
  for (g in 1:length(grouping)) {
    consecs <- seq(1, grouping[g])
    consecsAll <- c(consecsAll, consecs)
  } #E for
  infoTable[,ncol(infoTable)] <- consecsAll
  return(infoTable)
} #Eof

askSep <- function(infoTableRed, doPar) {
  cns <- colnames(infoTableRed)[-ncol(infoTableRed)]
  naInd <- which(is.na(infoTableRed[-ncol(infoTableRed)]), arr.ind = TRUE)
  if (nrow(naInd) != 0){
    naInd <- unique(naInd[,2])
    cns <-  cns[-naInd]
  }
  nrOk <- FALSE
  tmp <- data.frame(classVariables = cns)
  while (!nrOk) {
    cat("Please provide the index of the class variable to be used for phyical separation of the data, provide 0 for no separation.\n\n")
    print(tmp)
    a <- readLines(n = 1)
    a <- as.numeric(a)
    if (all(is.numeric(a)) & length(a) == 1 & a <= length(cns)) {
      nrOk = TRUE
    }
  } # E while
  if (a == 0) {
  	msg <- "without separation... "
  } else {
  	msg <- paste("separated by ", cns[a], "... ", sep="")
  }
  if (doPar) {dpadd <- "(parallel) "} else {dpadd <- "(seriell) "}
  cat(paste("Importing ", nrow(infoTableRed), " datafiles ", dpadd, msg,  sep=""))
  return(list(sepID = a, sepChar = cns[a]))
} #Eof

getHeader <- function(singleFilename, tz = "EST", yearCorrect = 2016) {
  con <- file(singleFilename, open="r")
  char <- readLines(con, n=19)
  close(con)
  char <- unlist(lapply(strsplit(char, ","), function(x) x[2]))
  ScannerID <- char[10]
  type <- if (char[12] == "1") {type <- "Hadamard"} else {type <- "Column"}
  PixelWidth <- as.numeric(char[15])
  SystemTemp <- as.numeric(char[4])/100
  DetectorTemp <- as.numeric(char[5])/100
  Humidity <- as.numeric(char[6])/100
  NrMeasPoint <- as.numeric(char[16])
  NrRep <- as.numeric(char[17])
  Gain <- as.numeric(char[18])
  MeasTime <- as.numeric(char[19])
  absTime <- as.POSIXct(strptime(as.character(char[2]), "%d/%m/%Y @ %H:%M:%S"), tz)
  if(!difftime(absTime, "1996-06-15 10:59:51 EST", tz = tz) > 3650) {
    absTime <- absTime-difftime(cut.POSIXt(absTime, "year"), paste0(yearCorrect, "-01-01"))
  }
  absTime <- as.character(absTime)
  #
  out<- data.frame(ScannerID, type, PixelWidth, SystemTemp, DetectorTemp, Humidity, NrMeasPoint, NrRep, Gain, MeasTime, absTime, stringsAsFactors = FALSE)
  return(out)
} #eof

getAbsRefSmplSpect <- function(singleFilename) {
  a <- utils::read.csv(singleFilename, skip=19)
  nir <- t(a[,2:4])
  wls <- paste("X", a[,1], sep="")
  colnames(nir) <- wls
  rownames(nir) <- const_nameSpectList
  return(nir)
} #eof

getNrOfCol <- function(singlePath, nrColAdd) {
  headerTable <- getHeader(singlePath)
  AbsRefSmplSpect <- getAbsRefSmplSpect(singlePath)
  Nr <- ncol(headerTable) + ncol(AbsRefSmplSpect) + nrColAdd
  colNames <- c(colnames(headerTable), colnames(AbsRefSmplSpect))
  return(list(Nr=Nr, colNames=colNames))
}#Eof

makeListLayout <- function(infoTable, ColInd_s, colInd_c) {
  nrSpectDf <- const_nrSpectDf
  uniq_s <- unique(infoTable[, ColInd_s])
  leScanner <- length(uniq_s)
  leConf <- levConf <-  NULL
  spectList <- vector("list", nrSpectDf)
  names(spectList) <- const_nameSpectList
  # first get the dimensions of the list
  for (i in 1:leScanner) {
    ind <- which(infoTable[,ColInd_s] == uniq_s[i])
    tmp <- infoTable[ind,]
    a <- levels(tmp[, colInd_c])
    b <- length(a)
    levConf <- c(levConf, a)
    leConf <- c(leConf, b)
  } # end for i
  outList <- vector("list", leScanner)
  names(outList) <- as.character(uniq_s)
  ao <- NULL
  for (i in 1: leScanner) {
    a <- rep(i, leConf[i])
    ao <- c(ao, a)
  }
  # create the list
  for (i in 1:leScanner) {
    outList[[i]] <- lapply(vector("list", leConf[i]), function(x) c(x, spectList))
    # Names <- levConf[(sum(leConf[1:i])-leConf[i]+1):sum(leConf[1:i])]
    names(outList[[i]]) <- as.character(levConf[which(ao == i)])
  }#efor
  return(outList)
  #return(list(leScanner=leScanner, uniqS=as.character(uniq_s), levConf = levConf, leConf=leConf))
} #eof

makeListLayout_configIDs_parallel <- function(uniq_c) {
	nrSpectDf <- const_nrSpectDf
	nameSpectList <- const_nameSpectList
	#
	spectList <- vector("list", nrSpectDf) # the spectList is the Innermost list, containing the three types of spectra
	names(spectList) <- nameSpectList
	configList <- lapply(vector("list", length(uniq_c)), function(x) c(x, spectList)) # make a list with length of unique configIDs and have in each three spectra lists
	names(configList) <- uniq_c
	return(configList)
} # EOF

dataImport_inner_old <- function(infoTable, infoTableRed, fileNames, scanIdCol, confIdCol) {
  Ind_s <- which(colnames(infoTable)  == scanIdCol)
  Ind_c <- which(colnames(infoTable)  == confIdCol)
  outList <- makeListLayout(infoTable, Ind_s, Ind_c)
  uniq_s <- unique(infoTable[, Ind_s])
  for (sid in 1: length(uniq_s)) {
    Ind <- which(infoTable[, Ind_s] == uniq_s[sid])
    selinfTable <- infoTable[Ind, ]
    uniq_c <- unique(selinfTable[, Ind_c])
    for (cid in 1: length(uniq_c)) {
      Ind2 <- which(selinfTable[, Ind_c] == uniq_c[cid])
      a <- getNrOfCol(fileNames[Ind[Ind2]][1], ncol(infoTable)-2)
      nrOfCol <- a$Nr
      colNames <- c(colnames(infoTableRed), a$colNames)
      colDfAbs <- colDfRef <- colDfSmpl <- as.data.frame(matrix(NA, length(Ind2), nrOfCol))
      colnames(colDfAbs) <- colnames(colDfRef) <- colnames(colDfSmpl) <- colNames
      fileNamesSelect <- fileNames[Ind[Ind2]]
      for (i in 1:length(fileNamesSelect)) {
        Header   <- getHeader(fileNamesSelect[i])
        spectra <- getAbsRefSmplSpect(fileNamesSelect[i])
        addHeader <- infoTableRed[Ind[Ind2][i], ]
        colDfAbs[i,] <- cbind(addHeader, Header, t(spectra[1,]))
        colDfRef[i,] <- cbind(addHeader, Header, spectra[2,, drop=FALSE])
        colDfSmpl[i,] <- cbind(addHeader, Header, spectra[3,, drop=FALSE])
      } #Efor i
      outList[[sid]][[cid]][[1]] <- colDfAbs
      outList[[sid]][[cid]][[2]] <- colDfRef
      outList[[sid]][[cid]][[3]] <- colDfSmpl
      # before turn posix and character to factor
    } #Efor uc
  } #Efor us
  return(outList)
} # Eof

dataImport_inner <- function(infoTable, infoTableRed, fileNames, scanIdCol, confIdCol, doPar) {
  Ind_s <- which(colnames(infoTable)  == scanIdCol)
  Ind_c <- which(colnames(infoTable)  == confIdCol)
  uniq_s <- unique(infoTable[, Ind_s])
  # check if we have to register a parallel backend
  if (doPar) {
  	registerParallelBackend() # own custom function
  } else {
  	foreach::registerDoSEQ() # from package foreach, is registering a serial backend -- so NO parallel computations
  }
  if (doPar) {
		parStackList <- foreach(sid= 1: length(uniq_s)) %dopar% { # going through the single scanner IDs
  	  		Ind <- which(infoTable[, Ind_s] == uniq_s[sid])
			selinfTable <- infoTable[Ind, ] # is the infoTable that only contains a single scannerID
			uniq_c <- unique(selinfTable[, Ind_c])
			configList <- makeListLayout_configIDs_parallel(uniq_c) # here make the collection of lists for the configIDs
			##	
			for (cid in 1: length(uniq_c)) { # going through the single config IDs within a single scanner
			  Ind2 <- which(selinfTable[, Ind_c] == uniq_c[cid])
			  a <- getNrOfCol(fileNames[Ind[Ind2]][1], ncol(infoTable)-2)
			  nrOfCol <- a$Nr
			  colNames <- c(colnames(infoTableRed), a$colNames)
			  colDfAbs <- colDfRef <- colDfSmpl <- as.data.frame(matrix(NA, length(Ind2), nrOfCol))
			  colnames(colDfAbs) <- colnames(colDfRef) <- colnames(colDfSmpl) <- colNames
			  fileNamesSelect <- fileNames[Ind[Ind2]]
			  for (i in 1:length(fileNamesSelect)) {
				Header   <- getHeader(fileNamesSelect[i])
				spectra <- getAbsRefSmplSpect(fileNamesSelect[i])
				addHeader <- infoTableRed[Ind[Ind2][i], ]
				colDfAbs[i,] <- cbind(addHeader, Header, t(spectra[1,]))
				colDfRef[i,] <- cbind(addHeader, Header, spectra[2,, drop=FALSE])
				colDfSmpl[i,] <- cbind(addHeader, Header, spectra[3,, drop=FALSE])
			  } #Efor i
			  configList[[cid]][[1]] <- colDfAbs
			  configList[[cid]][[2]] <- colDfRef
			  configList[[cid]][[3]] <- colDfSmpl
			} #Efor uc
			out <- configList # the result from the parallel
  	  	} # end parallel for i
  	  	names(parStackList) <- as.character(uniq_s) # give it back the names
  	  	return(parStackList)
  } else { # so we do NOT do it in parallel
	outList <- makeListLayout(infoTable, Ind_s, Ind_c) # the list for the serial processing
  	for (sid in 1: length(uniq_s)) {
    	Ind <- which(infoTable[, Ind_s] == uniq_s[sid])
   		selinfTable <- infoTable[Ind, ]
    	uniq_c <- unique(selinfTable[, Ind_c])
    	for (cid in 1: length(uniq_c)) {
      		Ind2 <- which(selinfTable[, Ind_c] == uniq_c[cid])
		    a <- getNrOfCol(fileNames[Ind[Ind2]][1], ncol(infoTable)-2)
			nrOfCol <- a$Nr
		    colNames <- c(colnames(infoTableRed), a$colNames)
			colDfAbs <- colDfRef <- colDfSmpl <- as.data.frame(matrix(NA, length(Ind2), nrOfCol))
			colnames(colDfAbs) <- colnames(colDfRef) <- colnames(colDfSmpl) <- colNames
			fileNamesSelect <- fileNames[Ind[Ind2]]
			for (i in 1:length(fileNamesSelect)) {
		        Header   <- getHeader(fileNamesSelect[i])
        		spectra <- getAbsRefSmplSpect(fileNamesSelect[i])
		        addHeader <- infoTableRed[Ind[Ind2][i], ]
        		colDfAbs[i,] <- cbind(addHeader, Header, t(spectra[1,]))
		        colDfRef[i,] <- cbind(addHeader, Header, spectra[2,, drop=FALSE])
        		colDfSmpl[i,] <- cbind(addHeader, Header, spectra[3,, drop=FALSE])
      		} # end for i
		    outList[[sid]][[cid]][[1]] <- colDfAbs
		    outList[[sid]][[cid]][[2]] <- colDfRef
		    outList[[sid]][[cid]][[3]] <- colDfSmpl
      		# before turn posix and character to factor
      	} # end for cid
	} # end for sid
	return(outList)
 } # end else if dopar
} # Eof

checkMainInput <- function(doPar, scanIdCol, confIdCol) {
	if (!all(is.logical(doPar)) | length(doPar) != 1) {
		stop("Please provide either 'TRUE' or 'FALSE' to the argument 'doPar'. Thank you.", call.=FALSE)
	}
	if (scanIdCol == "def") {
		assign("scanIdCol", const_scanIdCol, pos=parent.frame(n=1))
	}
	if (confIdCol == "def") {
		assign("confIdCol", const_confIdCol, pos=parent.frame(n=1))
	}
} # EOF

#' @title Import Data
#' @description Imports data from the folder-structure in the 'rawdata' folder
#' @details XXX Here the details of how the folder should be named, with 
#' separators etc.
#' @param doPar Logical. If datafile import should be done in parallel.
#' @param scanIdCol The standard colum name for the scanner Id.
#' @param confIdCol  The standard colum name for the configuration Id.
#' @param NrAdjust Defaults to 1 XXX
#' @return A (big) list.
#' @export
dataImport <- function(doPar=TRUE, scanIdCol="def", confIdCol="def", NrAdjust = 1) {
  checkMainInput(doPar, scanIdCol, confIdCol) # is assigning variables here !!!
  fileNames <- list.files(const_rawdataFolder, full.names = TRUE, recursive = T, pattern = const_fileExtension)
  infoTable <- createInfoTable(fileNames, NrAdjust)
  infoTableRed <- infoTable[, which(!colnames(infoTable) %in% c(scanIdCol, confIdCol))]
  infoTableRed <- reCharacter(infoTableRed)
  sepID <- askSep(infoTableRed, doPar)
  Ind_s <- which(colnames(infoTable)  == scanIdCol)
  Ind_c <- which(colnames(infoTable)  == confIdCol)
  outList <- makeListLayout(infoTable, Ind_s, Ind_c)
  if (sepID$sepID == 0) {
    outList <- dataImport_inner(infoTable, infoTableRed, fileNames, scanIdCol, confIdCol, doPar)
  } else {
    cN <- sepID$sepChar # this is the splitting col name
    indcN <- which(colnames(infoTable) == cN)
    a <- levels(infoTable[,indcN])
    outList <- vector("list", length(a))
    names(outList) <- a
    for (i in 1:length(a)) {
      ind <- which(infoTable[, indcN] == a[i])
      infoTableSel <- infoTable[ind, ]
      infoTableSelRed <- infoTableSel[, which(!colnames(infoTableSel) %in% c(scanIdCol, confIdCol))]
      infoTableSelRed <- reCharacter(infoTableSelRed)
      selectionList <- dataImport_inner(infoTableSel, infoTableSelRed, fileNames, scanIdCol, confIdCol, doPar)
      outList[[i]] <- selectionList
    } # Efor i
  } #Eif
  cat("ok.\n")
  return(outList)
  #return(list(colDfAbs=colDfAbs, colDfRef=colDfRef, colDfSmpl=colDfSmpl))
} #Eof
