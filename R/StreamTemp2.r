###################################
#Functions for routing LPF-GUESS Runoff
#
#Created by Jerad Hoy
#Date Created: 7/7/2015
#Date Last Modified:
#
#Sourced into streamNetwork.r file to do actual routing
#

#' Stream Temperature
#'
#' Routes runoff through a stream network using various physical relationships.
#'
#' @param 
#' @param catchments
#' @param Rsurf
#' @param Rsub
#' @param spinUpCycles
#' @param spinUpYears
#' @param debugMode
#' @param by
#' @param widthCoeffs
#' @param manningN
#' @param slopeMin
#' @param aCoeffCoeff
#' @param outputExtraVars
#'
#' @return List of matrices of routing variable timeseries.
#'
#' @examples
#' flow <- RouteWater(edges=edgesInBounds, catchments=catchmentsInBounds, Rsurf=surfaceRunoff,  Rsub=subsurfRunoff, spinUpCycles=gwSpinUpCycles, spinUpYears=spinUpYears, debugMode=F, by=timeStep, widthCoeffs=streamWidthCoeffs, manningN=manningN, slopeMin=slopeMin, aCoeffCoeff=aCoeffCoeff)
#'
#' @name StreamTemp
#' @export

# Routes surface and subsurface water through river network
StreamTempCpp <- function(edges, catchments, RsurfSnow, RsurfNoSnow, Tair, simFlow, defaults=setupList, by="month", outputExtraVars=T, debugMode=T, runStart=1, runStop=NULL, K=.1, etaInt=10, outFile=NULL, prof=NULL){

	if(!is.null(outFile)){
		out <- file(outFile)
		sink(out, append=T, type=c("output", "message"), split=F)
	}

	idField <- defaults$edgeIdField
	nextDownField <- defaults$edgeNextDownField
    
	edges <- as.data.frame(edges@data)

    # Order edges by Shreve order so calculation is in right order
    edges <- edges[order(edges[, defaults$edgeOrderField]),]

	lengthKM <- edges[, defaults$edgeLengthField] * 120

	edges[,idField] <- as.character(edges[, idField])
	edges[,nextDownField] <- as.character(edges[, nextDownField])


    # Set the timeLength to run simulation
	runStop <- ifelse(is.null(runStop), nrow(simFlow[[1]]), runStop)

	RsurfSnow <- RsurfSnow[runStart:runStop,match(edges[, idField], colnames(RsurfSnow))]

	RsurfNoSnow <- RsurfNoSnow[runStart:runStop,match(edges[, idField], colnames(RsurfNoSnow))]

	annualTmean <- rep(getAnnualTmean(Tair), times=1, each=ifelse(by == "month", 12, 365))

	Tair <- Tair[runStart:runStop,match(edges[, idField], colnames(Tair))]

	simFlow <- lapply(simFlow, function(x) x[runStart:runStop,match(edges[, idField], colnames(x))])	

    timeLength <- length(runStart:runStop)
    # Create seed matrix to use for storing data in results

	Tsnow <- .1

    if(by == "day"){
		# Set the velocity conversion factor for a daily timestep
		vConvFactor <- 60*60*24/1000
	}

	# Set days in month to use for monthly-timestep velocity conversion factor
	daysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30, 31)
	if(by == "month"){
		##############Assumes start month is january
		vConvFactor <- rep(daysInMonth*60*60*24/1000, length.out=timeLength)
		velocities <- simFlow$v * vConvFactor
	}

	

	orders <- edges[, defaults$edgeOrderField] 
	ids <- edges[, idField]

	start <- as.numeric(format(Sys.time(), "%s"))

	parentList <- lapply(edges[,idField], function(x) {which(edges[,nextDownField] == x)})


    print("About to Run")

	if(!is.null(prof)){
		Rprof(prof, line.profiling=T, memory.profiling=T)
	}

	####Rearrage all matrices in correct stream order
	####Make parent list into numeric column indices

	#print(paste0("Type: ", typeof(), " Class: ", class()))
	lapply(ls(), function(x) {print(paste0(x," - Type: ", typeof(get(x)), " Class: ", class(get(x))))})
	lapply(ls(), function(x) {print(paste0(x," - is.null?: ", is.null(get(x))))})

	temp <- streamTempLoop(timeLength=timeLength, edgeIDs=ids, orders=orders, velocities=velocities, lengths=lengthKM, RsurfSnow=as.matrix(RsurfSnow), RsurfNoSnow=as.matrix(RsurfNoSnow), flowqSub=as.matrix(simFlow$qSub), flowqOut=as.matrix(simFlow$qOut), flowqIn=as.matrix(simFlow$qIn), flowsRiv=as.matrix(simFlow$sRiv), annualTmean=annualTmean, by=by, parentList=parentList, K=K, Tair=as.matrix(Tair))

	colnames(temp[[1]]) <- colnames(RsurfSnow)
	rownames(temp[[1]]) <- rownames(RsurfSnow)


	if(!is.null(prof)){
		Rprof(NULL)
	}

	if(!is.null(outFile)){
		close(out)
	}

	return(temp)
}



getAnnualTmean <- function(tMeanFrame){
	results <- c()
    for(i in 1:(nrow(tMeanFrame)/12)){
		#print(1:12+12*(i-1))
		results <- c(results, mean(colMeans(tMeanFrame[1:12+12*(i-1),])))
	}
	return(results)
}
