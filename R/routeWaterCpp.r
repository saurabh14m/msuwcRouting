###################################
#Functions for routing LPF-GUESS Runoff
#
#Created by Jerad Hoy
#Date Created: 7/7/2015
#Date Last Modified:
#
#Sourced into streamNetwork.r file to do actual routing
#


#' Routes surface and subsurface water through river network
#'
#'
#' @param dat dataframe of guage data to be cleaned and filled
#' @param simStartDate date of the beginning of the simulation
#' @param simEndDate date of the end of the simulation
#'
#' @return Dataframe of stream gauge data with missing values filled with NAs
#'
#' @examples
#' cleanDat(dat, simStartDate, simEndDate)
#'
#' @name RouteWater
#' @export

RouteWaterCpp <- function(edges, catchments, Rsurf, Rsub, defaults=setupList, spinUpCycles=0, spinUpYears=10, debugMode=F, by="day", widthCoeffs=c(.3, .6), manningN=.07, slopeMin=.01, aCoeffCoeff=3, outputExtraVars=T, etaInt=10){ 

	edges <- as.data.frame(edges@data)
    edges <- edges[order(edges[, edgeOrderField]),]

	idField <- defaults$edgeIdField
	nextDownField <- defaults$edgeNextDownField

    edges <- AssignContribArea(edges, catchments)
    edges <- AssignBfWidth(edges, widthCoeffs[1], widthCoeffs[2])
    edges <- CorrectEdgeSlopes(edges, slopeMin)
    edges <- AssignAcoeff(edges, catchments, aCoeffCoeff)

    LengthKM <- edges[, defaults$edgeLengthField] * 120
	slopes <- edges[, defaults$edgeSlopeField]


    # Order edges by Shreve order so calculation is in right order
	edges[,idField] <- as.character(edges[,idField])
	edges[,nextDownField] <- as.character(edges[, nextDownField])

	Rsurf <- Rsurf[,match(edges[, idField], colnames(Rsurf))]
	Rsub <- Rsub[,match(edges[, idField], colnames(Rsub))]


    # Set the timeLength to avoid re-executing nrow function
    timeLength <- nrow(Rsurf) 
    # Create seed matrix to use for storing data in results

    # Set days in month to use for monthly-timestep velocity conversion factor
    daysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30, 31)
    
    if(by == "day"){
		# Set the velocity conversion factor for a daily timestep
		vConvFactor <- 60*60*24/1000
    }
    
	if(by == "month"){
		# Set velocity conversion factor based on month of timestep
		#ASSSUMES that we are starting in january, need to fix to be flexible with start dates
		vMonthFactor <- rep(60*60*24*daysInMonth/1000, length.out=timeLength)
	}

	parentList <- lapply(edges[,idField], function(x) {which(edges[,nextDownField] == x)})
	orders <- edges[, defaults$edgeOrderField]
	ids <- edges[, idField]

	
    print("About to Run")

	start <- as.numeric(format(Sys.time(), "%s"))

	flow <- routeWaterLoop(timeLength=timeLength, edgeIDs=ids, orders=orders, streamLengths=lengthKM, streamWidths=, streamSlopes=slopes, RsurfSnow=as.matrix(RsurfSnow), RsurfNoSnow=as.matrix(RsurfNoSnow), flowqSub=as.matrix(simFlow$qSub), flowqOut=as.matrix(simFlow$qOut), flowqIn=as.matrix(simFlow$qIn), flowsRiv=as.matrix(simFlow$sRiv), annualTmean=annualTmean, by=by, parentList=parentList, K=K, Tair=as.matrix(Tair))

    return(flow)
}
