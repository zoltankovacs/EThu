
impD <- function(NrFile = 1) { # import the raw txt file(s), NrFile: which file to import
  files <- list.files(paste0(getwd(), "/rawdata"), full.names = TRUE, pattern = "*.txt")
  filesShort <- list.files(paste0(getwd(), "/rawdata"), full.names = FALSE, pattern = "*.txt")
  print(paste0("The file: '", filesShort[NrFile], "' was imported"))
  return(read.table(files[NrFile]))
} #Eof

getFolderName <- function() {
  return(unlist(strsplit(getwd(),split=c("/")))[length(unlist(strsplit(getwd(),split=c("/"))))])
} #Eof

getNames <- function(rawData, nameChange = FALSE) {
  smplName <- unlist(lapply(strsplit(rownames(rawData), "_"), function(x) paste(x[1:length(x)-1], collapse ="_")))
  smplPos <- as.numeric(unlist(lapply(strsplit(rownames(rawData), "_"), function(x) x[length(x)])))
  repeats <- as.numeric(apply(data.frame(rle(smplName)$lengths), 1, function(x) seq(1:x)))
  print("The following samples are in the set:")
  smplName <- as.factor(smplName)
  smplNameLev <- levels(smplName)
  print(smplNameLev)
  if (nameChange) {
	for (i in 1:length(smplNameLev)) {
		cat(paste0("\nPlease provide the new name for ", smplNameLev[i], " (use '-' to not change): "))
		a <- readLines(n = 1)
		if (!a == "-") {
			levels(smplName)[i] <- a
		}
	} #Efor
  } #Eof
  return(cbind(data.frame(smplName), data.frame(smplPos, repeats)))
} #Eof

#' @title add variables to ET data
#' @description Add variables
#' @details XXX Here the details of how the folder should be named, with 
#' separators etc.
#' @param rawData ET raw data imported from the txt file
#' @param day optional argument useful if the experiment was performed in different days
#' @param nameChange optional argument if the name of the groups has to be renamed
#' @return a table with additional columns
#' @export
mDataStr <- function(rawData, day = 1, nameChange = FALSE){ # make data structure
  data <- data.frame(getNames(rawData, nameChange = nameChange))
  if (!is.na(day)) {
	data <- cbind(data, day)
  }
  data$sensors <- rawData
  return(data)
} #Eof

chngNames <- function(charVect, toWhat){
  NewNames <- as.factor(charVect)
  levels(NewNames) <- toWhat
  return(NewNames)
} #Eof

addNumRep <- function(data){
  sInd <- which(colnames(data) == "sensors")
  numRep <- data[,-sInd]
  for(i in 1:ncol(numRep)){
    if (length(unique(numRep[,i])) > 9) {
		numRep[,i] <- makeHeatColors(numRep[,i], startCol = "red", endCol = "blue")
	} else {
		numRep[,i] <- as.numeric(as.factor(numRep[,i]))
	} #Eif
  } #Efor
  dataa <- data[, -sInd]
  dataa$numRep <- data.frame(numRep)
  dataa$sensors <- data[, sInd]
  return(dataa)
} #Eof


# addVars <- function(rawData, day = 1) {
#   colnames(rawData) <- paste("X", colnames(rawData), sep = "_")
#   header <- getNames(rawData)
#   return(cbind(header, day, rawData))
# } #Eof

importStructureData <- function(NrFile = 1, day = 1, nameChange = FALSE) {
	rawData <- impD(NrFile)
	dataStr <- mDataStr(rawData, nameChange = nameChange)
	dataStrNr <- addNumRep(dataStr)
	return(dataStrNr)
} #Eof



















