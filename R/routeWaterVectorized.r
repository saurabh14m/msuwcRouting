RouteWaterVectorized <- function(edges, catchments, Rsurf, Rsub, spinUpCycles=0, spinUpYears=10, debugMode=F, by="day", widthCoeffs=c(.3, .6), manningN=.07, slopeMin=.01, aCoeffCoeff=3, outputExtraVars=T){
    
    edges <- AssignContribArea(edges, catchments)
    edges <- AssignBfWidth(edges, widthCoeffs[1], widthCoeffs[2])
    edges <- CorrectEdgeSlopes(edges, slopeMin)
    edges <- AssignAcoeff(edges, catchments, aCoeffCoeff)
    edges$LengthKM <- edges@data[, edgeLengthField] * 120



    # Order edges by Shreve order so calculation is in right order
    edges <- edges@data[order(edges@data[, edgeOrderField]),]
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

    
    if(by == "day"){
	# Set the velocity conversion factor for a daily timestep
	vConvFactor <- 60*60*24/1000
    }
    
    # Set days in month to use for monthly-timestep velocity conversion factor
    daysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30, 31)

    cycleCount <- 0

	#set unique edge orders
	orders  <- sort(unique(edges[, edgeOrderField]))



    print("About to Run")
    while(cycleCount <= spinUpCycles){

		###### BEGIN ROUTING
		# Loop through each day/timestep
		#print(paste("cycleCount:", cycleCount))

		spinUpSteps <- spinUpYears*ifelse(by == "day", 365, 12)

		if(cycleCount < spinUpCycles){
			stepsToLoop <- spinUpSteps
		} else {
			stepsToLoop <- timeLength
		}
			







			startTimeCounter <- Sys.time()

		for(timeStep in 1:stepsToLoop){
			#print(paste("Cycle:", cycleCount, "Timestep:",  timeStep))

			# Start time to use for fucntion time-estimates
			if(cycleCount == 0 & timeStep == 1){
			start <- Sys.time()
			}

			if(by == "month"){
			# Set velocity conversion factor based on month of timestep
			#vConvFactor <- 60*60*24*daysInMonth[as.integer(format(zoo::as.yearmon(rownames(Rsurf[timeStep,])), "%m"))%%12 + 1]/1000
			vConvFactor <- 60*60*24*daysInMonth[as.numeric(format(zoo::as.yearmon(rownames(Rsurf[timeStep,])), "%m"))]/1000
			
			}

			for(edgeToProcess in 1:nrow(edges)){

				#print("got here")
				if(edges[edgeToProcess,edgeOrderField] != 1){

					#qIn <- sum(results$qOut[timeStep,
					#	   edges[edges[, edgeNextDownField] == edges[edgeToProcess, edgeIdField], edgeIdField]])
				}
			}
			

			next
			for(ord in 1:length(orders)){



				# Set hydroID of edges so don have to keep subsetting
				hydroID <- edges[i,edgeIdField]

				# Set rS(surface runoff), rSub(sub runoff)
				rS <- Rsurf[timeStep, hydroID]
				rSub <- Rsub[timeStep, hydroID]

				# Set fixed channel dimensions
				len <- edges$LengthKM[i]
				width <- edges$bfWidth[i]

				if(timeStep == 1){
					#print(by)
					# On the first timestep, sRiv and sSub (storage) is 0, height needs to be initialized
					#sRiv <- tail(results$sRiv[1:(stepsToLoop*ifelse(by == "day",365, 12)), hydroID], n=1)
					sRiv <- results$sRiv[spinUpSteps, hydroID]
					height <- results$h[spinUpSteps, hydroID]
					sSub <- results$sSub[spinUpSteps, hydroID]
					print(sSub)

					#sRiv <- 0
					#height <- 1
					#sSub <- 0
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
				#v <- ifelse(v >= 3.0, 3.0, v)	
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
				#qOut <- sRiv + (1-len/v)*qIn + (1 - len/(2*v))*(rS + qSub)
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
					#results$qIn[timeStep, hydroID] <- qIn
					#results$qSub[timeStep, hydroID] <- qSub
					#results$v[timeStep, hydroID] <- v/vConvFactor
					#results$h[timeStep, hydroID] <- height
				}
				
				if(debugMode){
					# DebugMode in order to debug problems
					print(" ")
					print(paste("Day:",timeStep, "Edge:", hydroID, "Order:", edges[i, edgeOrderField]))
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



			if(cycleCount == 0 & timeStep == 10){

			print(paste("ETA:", 
					#round((Sys.time()-start) * (timeLength*(spinUpCycles1+1)length((timeStep+1):timeLength)/10*(spinUpCycles+1), 2), 
					round((Sys.time()-start) * (timeLength + spinUpCycles*spinUpYears*ifelse(by == "day",365, 12))/10), 
					"seconds or", 
					#round((Sys.time()-start)*length((timeStep+1):timeLength)/10/60*(spinUpCycles+1), 2), 
					round((Sys.time()-start) * (timeLength + spinUpCycles*spinUpYears*ifelse(by == "day",365, 12))/10/60), 
					"minutes. Will finish at", 
					(Sys.time()-start)*length((timeStep+1):timeLength)/10*(spinUpCycles+1)+Sys.time()
				))
			}
		}
		cycleCount <- cycleCount + 1
    }	
	print(Sys.time() - startTimeCounter)
    return(results)
}
