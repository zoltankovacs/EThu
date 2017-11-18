ssc <- function(dataset, varName, smplName, includ = TRUE) {
	varInd <- which(colnames(dataset) == varName)
	if (length(varInd) < 1) stop(paste("There is now such variable:", varName))
	smplInd <- which(dataset[,varInd] == smplName)
	if (length(smplInd) < 1) stop(paste0("There is now such sample name: ", varName, ", in variable: ", varName))
	if (includ){
		dataset <- dataset[smplInd, ]
	} else {
		dataset <- dataset[-smplInd, ]
	}
	return(dataset)
} #Eof

print_datasetInfo <- function(dataset) {
	colNames <- colnames(dataset)
	colNames <- colNames[2:which(colNames == "numRep")-1]
		cat(paste0("\nVariables and there levels:\n"))
		for (i in 1: length(colNames)) {
			uniques <- unique(as.character(dataset[,i]))
			if (length(uniques) > 10) {
				uniques <- c(uniques[1:10], "...")
			} # Eif
			cat(c(paste0("\n", colNames[i], ": "), uniques))
		} #Efori
} #Eof


do_ddply <- function(dataset, VectorToAve, NonNumCols = NA){
  if (any(is.na(NonNumCols))){
	out <- ddply(data.frame(dataset), .(VectorToAve), colwise(mean))
  } else {
	NumPart <- ddply(data.frame(dataset[ , -NonNumCols]), .(VectorToAve), colwise(mean))
	#NonNumPart <- NA
	for (i in sort(unique(VectorToAve))){
		Ind <- min(which(VectorToAve == i))
		if (!exists("NonNumPart")){
			NonNumPart <- dataset[Ind, NonNumCols]
		} else {
			NonNumPart <- rbind(NonNumPart, dataset[Ind, NonNumCols])
		}
	}#E for i
	
	rownames(NonNumPart) <- rownames(NumPart)
  out <- cbind(NumPart[,1], NonNumPart[,], NumPart[,-1])
  }
  colnames(out)[1] <- "ToAve"
  out
}#E do_ddplyZ
