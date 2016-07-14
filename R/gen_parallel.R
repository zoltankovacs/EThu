checkHaveParallel <- function() {
	a <- foreach::getDoParName()
	if (is.null(a)) { return(FALSE) }
	if (a == "doSEQ") { return(FALSE) }
	return(TRUE)
} # EOF

# this to be used
registerParallelBackend <- function() {
	nrWorkers <- parallel::detectCores(logical=FALSE)
	haveParallel <- checkHaveParallel()
	if (!haveParallel) {
		if (is.na(nrWorkers)) {
			doParallel::registerDoParallel()
		} else {
			if (is.numeric(nrWorkers) & (length(nrWorkers) == 1)) {
				doParallel::registerDoParallel(nrWorkers)		
			} else {
				stop("Please provide a length one numeric as the number of worker processes in the settings file. Thank you very much.", call.=FALSE)
			}
		} # end else
	} # end if !haveParallel
} # EOF
