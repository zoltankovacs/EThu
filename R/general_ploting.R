makeHeatColors <- function(vector = vector, startCol = "blue", endCol = "red", s = 1, v = 1, decreasing = FALSE){
  vector0 <- vector
  if (any(is.na(vector0))){
    naIndex <- which(is.na(vector0))
    vector <- vector0[-naIndex]
  }
  n <- length(unique(vector))
  colList <- as.matrix(c(0, 1/6, 2/6, 3/6, 4/6, 5/6), ncol=1)
  Colors <- c("red", "yellow", "green", "cyan", "blue", "magenta")
  rownames(colList) <- Colors
  startColInd <- which(rownames(colList) == startCol)
  endColInd <- which(rownames(colList) == endCol)
  vector
  if (length(startColInd) == 0 || length(endColInd) == 0){
    print("Please select from colors: red, yellow, green, cyan, blue, magenta")
    
  } else {  
    Cols <- rainbow(n, s = s, v = v, start = as.numeric(colList[startColInd,1]), end = as.numeric(colList[endColInd,1]))
    
    Colors <- rep(NA, n)
    
    if (!decreasing){
      for (i in 1: n){
        actLevel <- sort(unique(vector))[i]
        actLevelindex <- which(vector == actLevel)
        Colors[actLevelindex] <- Cols[i]
      }
    } else {
      for (i in 1: n){
        actLevel <- sort(unique(vector), decreasing = TRUE)[i]
        actLevelindex <- which(vector == actLevel)
        Colors[actLevelindex] <- Cols[i]
      }
    }
    if (any(is.na(vector0))){
      vector0[-naIndex] <- Colors
      Colors <- vector0
    }
    return(Colors)
  }
} #Eof

chngGroupCol <- function(dataset, classVar) {
	classVarInd <- which(colnames(dataset) == classVar)
	numRepInd <- which(colnames(dataset$numRep) == classVar)
	classVarCont <- dataset[,classVarInd]
	newCols <- dataset$numRep[, numRepInd]
	cat(paste(1:length(levels(classVarCont)), levels(classVarCont), "\n"))
	colPalet <- colors()
	nrOk <- "st"
	while (!nrOk == "end") {
    cat("Please provide the index of the class variable to be used for phyical separation of the data, provide 0 for no separation.\n\n")
    a <- nrOk <- readLines(n = 1)
	if (nrOk != "end") {
		a <- strsplit(a, "-")[[1]]
		Nr <- as.numeric(a[1])
		color <- a[2]
		if (!any(colPalet == color)){
			cat("Please provide a valid color name from the following list.\n")
			print(colPalet)
		} else {
			grToChange <- levels(classVarCont)[Nr]
			newCols[which(classVarCont == grToChange)] <- color
		}# Eif
		dataset$numRep[, numRepInd] <- newCols
	} # if
  } # E while
  return(dataset)
} #Eof