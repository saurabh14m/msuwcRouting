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

RouteWater <- function(edges, catchments, Rsurf, Rsub, spinUpCycles=0, spinUpYears=10, debugMode=F, by="day", widthCoeffs=c(.3, .6), manningN=.07, slopeMin=.01, aCoeffCoeff=3, outputExtraVars=T, etaInt=10){ 

	edges <- as.data.frame(edges@data)


    edges <- AssignContribArea(edges, catchments)
    edges <- AssignBfWidth(edges, widthCoeffs[1], widthCoeffs[2])
    edges <- CorrectEdgeSlopes(edges, slopeMin)
    edges <- AssignAcoeff(edges, catchments, aCoeffCoeff)
    edges$LengthKM <- edges[, edgeLengthField] * 120



    # Order edges by Shreve order so calculation is in right order
    edges <- edges[order(edges[, edgeOrderField]),]
	edges[,edgeIdField] <- as.character(edges[,edgeIdField])

	Rsurf <- Rsurf[,match(edges[, edgeIdField], colnames(Rsurf))]
	Rsub <- Rsub[,match(edges[, edgeIdField], colnames(Rsub))]


    # Set the timeLength to avoid re-executing nrow function
    timeLength <- nrow(Rsurf) 
    # Create seed matrix to use for storing data in results


    seedMatrix <- matrix(0, nrow=timeLength, ncol=ncol(Rsurf),
        dimnames=list(rownames(Rsurf), edges[, edgeIdField]))



    results <- list(
        qIn = seedMatrix,
        qOut = seedMatrix,
        sRiv = seedMatrix,
        sSub = seedMatrix, 
        qSub = seedMatrix, 
        v = seedMatrix,
        h = matrix(1, nrow=timeLength, ncol=ncol(Rsurf),
        dimnames=list(rownames(Rsurf), edges[, edgeIdField]))
    )
    rm(seedMatrix)

	print(timeLength)
	print(dim(results$sRiv))
    
    if(by == "day"){
	# Set the velocity conversion factor for a daily timestep
	vConvFactor <- 60*60*24/1000
    }
    
	if(by == "month"){
	# Set velocity conversion factor based on month of timestep
		vConvFactor <- 60*60*24*daysInMonth[as.numeric(format(zoo::as.yearmon(rownames(Rsurf[timeStep,])), "%m"))]/1000
	}

    # Set days in month to use for monthly-timestep velocity conversion factor
    daysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30, 31)

    cycleCount <- 0
	stepsLooped <- 0
	
    print("About to Run")

	start <- as.numeric(format(Sys.time(), "%s"))

    while(cycleCount <= spinUpCycles){

		###### BEGIN ROUTING
		# Loop through each day/timestep
		#print(paste("cycleCount:", cycleCount))

		spinUpSteps <- spinUpYears*ifelse(by == "day", 365, 12)
		print(timeLength)
		print(cycleCount)
		print(spinUpCycles)

		if(cycleCount < spinUpCycles){
			stepsToLoop <- spinUpSteps
		} else {
			stepsToLoop <- timeLength
			spinUpSteps <- stepsToLoop
		}
			
		for(timeStep in 1:stepsToLoop){
			#print(paste("Cycle:", cycleCount, "Timestep:",  timeStep))

			# Start time to use for fucntion time-estimates



			# Loop through edges in river network
			for(i in 1:nrow(edges)){
			#print(i)

				# Set hydroID of edges so don have to keep subsetting
				hydroID <- edges[i,edgeIdField]
				
				if(edges[i, edgeOrderField] == 1){
					# Set qIn to 0 for edges of Order 1
					qIn <- 0
				} else {
					# Sum qOut of parent edges for edges of Order > 1
					qIn <- sum(results$qOut[timeStep,
						   edges[edges[, edgeNextDownField] == edges[i, edgeIdField], edgeIdField]])
				}


				# Set rS(surface runoff), rSub(sub runoff)
				rS <- Rsurf[timeStep, hydroID]
				rSub <- Rsub[timeStep, hydroID]

				# Set fixed channel dimensions
				len <- edges$LengthKM[i]
				width <- edges$bfWidth[i]

				if(timeStep == 1){
					#print(by)
					# On the first timestep, sRiv and sSub (storage) is 0, height needs to be initialized
					print(spinUpSteps)
					print(hydroID)
					print(dim(results$sRiv))

					sRiv <- results$sRiv[spinUpSteps, hydroID]
					height <- results$h[spinUpSteps, hydroID]
					sSub <- results$sSub[spinUpSteps, hydroID]
					#print(sSub)

				} else {

					# Set sRiv, sSub from previous day
					sRiv <- results$sRiv[timeStep-1, hydroID]
					sSub <- results$sSub[timeStep-1, hydroID]

					# Calculate height form sRiv and channel dimensions
					height <- (sRiv * vConvFactor*1000)/(len * 1000 * width)
				}

				# Manning's Eqation used to calculate velocity
				v <- (((height*width)/(2*height+width))^(2/3) * edges[i, "Slope2"]^.5)/manningN 
				# Set velocity caps, upper cap may be unncessary, need to calibrate values to get realisic velocity
				v <- ifelse(v <= .01, .01, v)	

				# Convert velocity form m/s to km/timestep (km/day)
				v <- v*vConvFactor

				# Caluclate groundwater discharge
				# Ignores current timestep subsurface runoff, assuming that groundwater movement is too slow
				qSub <- ((sSub+rSub)/edges$aCoeff[i])^(1/.5)
				# Could base it off stream dimensions
				qLoss <- 0
				
				if(qIn > 1e5){
					stop("Ridiculous Values!!!")
				}
				
				# Assumes that l/v <= 1, if not, need different routing scheme 
				# Delta t is always one, so l/v can work
				qOut <- if(len/v < 1){
					sRiv + (1 - len/v)*qIn + (1 - len/(2*v))*(rS + qSub)
				} else {
					(v/len)*sRiv + 0*qIn + (v/(2*len))*(rS + qSub) 
				}


				# Store values in results list
				results$sRiv[timeStep, hydroID] <- sRiv + rS + qIn + qSub - qOut - qLoss
				results$sSub[timeStep, hydroID] <- sSub + rSub - qSub
				results$qOut[timeStep, hydroID] <- qOut

				if(outputExtraVars){
					results$qIn[timeStep, hydroID] <- qIn
					results$qSub[timeStep, hydroID] <- qSub
					results$v[timeStep, hydroID] <- v/vConvFactor
					results$h[timeStep, hydroID] <- height
				}
				
				if(debugMode){
					# DebugMode in order to debug problems
					print(" ")
					print(paste("Cycle:", cycleCount, "Timestep:",timeStep, "Edge:", hydroID, "Order:", edges[i, edgeOrderField]))
					print(paste("qIn =", qIn))
					print(paste("qOut=",qOut))
					print(paste("qSub=",qSub))
					print(paste("sRiv =",sRiv))
					print(paste("rS =", rS))
					print(paste("rSub =",rSub))
					print(paste("sSub =",sSub))
					print(paste("len =", len))
					print(paste("width =", width))
					print(paste("height =", height))
				}


			}
			print(stepsLooped)

			stepsLooped <- stepsLooped + 1



			if((stepsLooped %% etaInt) == 0){
				print(spinUpCycles)
				print(spinUpYears)
				print(timeLength)
				print(stepsLooped)
				timeElapsed <- as.numeric(format(Sys.time(), "%s"))-start
				print(paste("Seconds elapsed:", timeElapsed))

				secondsToFinish <- (timeElapsed/stepsLooped * (timeLength + spinUpCycles*spinUpYears*ifelse(by == "day",365, 12) - stepsLooped))

				print(paste("ETA:", 
					#round((Sys.time()-start) * (timeLength*(spinUpCycles1+1)length((timeStep+1):timeLength)/10*(spinUpCycles+1), 2), 
					round(secondsToFinish), 
					"seconds or", 
					#round((Sys.time()-start)*length((timeStep+1):timeLength)/10/60*(spinUpCycles+1), 2), 
					round(secondsToFinish/60, digits=2), 
					"minutes. Will finish at", 
					secondsToFinish + Sys.time()
				))
			}
		}
		cycleCount <- cycleCount + 1
    }	
    return(results)
}
