
# Routes surface and subsurface water through river network
StreamTemp <- function(edges, catchments, RsurfSnow, RsurfNoSnow, Tair, simFlow, defaults=setupList, by="month", outputExtraVars=T, debugMode=T, runStart=1, runStop=NULL, K=.1, etaInt=10, outFile=NULL, prof=NULL){

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

    # Set the timeLength to run simulation
	runStop <- ifelse(is.null(runStop), nrow(simFlow[[1]]), runStop)

	RsurfSnow <- RsurfSnow[runStart:runStop,match(edges[, idField], colnames(RsurfSnow))]

	RsurfNoSnow <- RsurfNoSnow[runStart:runStop,match(edges[, idField], colnames(RsurfNoSnow))]

	annualTmean <- rep(getAnnualTmean(Tair), times=1, each=ifelse(by == "month", 12, 365))


	Tair <- Tair[runStart:runStop,match(edges[, idField], colnames(Tair))]

	simFlow <- lapply(simFlow, function(x) x[runStart:runStop,match(edges[, idField], colnames(x))])	

	#########
	## Need to subset all of simflow for faster computation
	#####


    timeLength <- length(runStart:runStop)
    # Create seed matrix to use for storing data in results


    seedMatrix <- matrix(0, nrow=timeLength, ncol=ncol(RsurfSnow),
        dimnames=list(rownames(RsurfSnow), edges[, idField]))



    results <- list(
        Twater = seedMatrix,
        TwaterLocal = seedMatrix,
        TwaterLocalWarmed = seedMatrix,
        TwaterQin = seedMatrix,
        TwaterSriv = seedMatrix,
        kRsQsub= seedMatrix,
        kQIn = seedMatrix,
        kQsRiv = seedMatrix
    )

	tracemem(results)

    rm(seedMatrix)

	Tsnow <- .1

    if(by == "day"){
		# Set the velocity conversion factor for a daily timestep
		vConvFactor <- 60*60*24/1000
	}

	# Set days in month to use for monthly-timestep velocity conversion factor
	if(by == "month"){
		##############Assumes start month is january
		daysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30, 31)
		vConvFactor <- rep(daysInMonth*60*60*24/1000, length.out=timeLength)
		velocities <- simFlow$v * vConvFactor
	}
	

    

	stepsLooped <- 0
	start <- as.numeric(format(Sys.time(), "%s"))

	edges <- data.table(edges)
	setkeyv(edges, idField)

	orders <- edges[, defaults$edgeOrderField] 
	ids <- edges[, idField]

	print(class(edges))
	a <- edges[defaults$edgeOrderField > 1, idField]
	parentList <- lapply(a, function(x) {as.character(edges[NextDownField <- 2 == x, idField])})
	names(parentList) <- a
	rm(a)

    print("About to Run")

	if(!is.null(prof)){
		Rprof(prof, line.profiling=T, memory.profiling=T)
	}

	for(timeStep in 1:timeLength){


		# Loop through edges in river network
		for(i in 1:length(ids)){
			#print(i)

			# Set hydroID of edges so don have to keep subsetting
			hydroID <- ids[i]
			ord = orders[i]


			
			v = velocities[timeStep, hydroID] ##m/s ##NEED to convert to m/month if month

			len = lengthKM[i]
			
			qSnow = RsurfSnow[timeStep, hydroID]
			qNoSnow = RsurfNoSnow[timeStep, hydroID]

			qGw = simFlow$qSub[timeStep, hydroID]
			qOut = simFlow$qOut[timeStep, hydroID]
			qIn = simFlow$qIn[timeStep, hydroID]

			Tgw = annualTmean[timeStep]+1.5 ###often 1-2*C higher than mean annual air temperature of a region, annual time series

			if(timeStep > 1){
				TwaterOld = results$Twater[timeStep - 1, hydroID] 
				sRiv = simFlow$sRiv[timeStep - 1, hydroID]

			} else {
				TwaterOld = 0
				sRiv = 0
			}

			TairLocal = Tair[timeStep, hydroID]

			if(by == "month"){
				TairLag = TairLocal
				lamda = 1
			} else {
				print(" by == day ??")
				TairLag #### Need to make some sort of lag for daily air temp
				lamda ####
			}

			##
			## NEED TO CHECK that qOut - qLocal equals qIn
			##

			qLocal = (qSnow+qGw+qNoSnow)
			#print(qLocal)




			# temperature of the water before effects of air and upstream temperature
			if(qLocal > 0){
				TwaterLocal = ((Tsnow*qSnow)+(Tgw*qGw)+(lamda*TairLag)*qNoSnow)/qLocal
			} else {
				TwaterLocal = 0
			}
			

			#print("Caculated TwaterLocal")




			#print(ord)
			if(ord == 1){

				#TwaterInitial = TwaterLocal
				TwaterUpstream = 0

			} else {
				#May need to calculated weighted temperatures
				qInUpstream = simFlow$qOut[timeStep, parentList[[hydroID]]]

				TwaterUpstream = results$Twater[timeStep, parentList[[hydroID]]]

				#print(TwaterIn)
				#print(qInUpstream)
				#print(TwaterUpstream)
				#print(length(TwaterUpstream))

				if(any(sapply(c(qInUpstream, TwaterUpstream), function(x) is.nan(x)))){
					stop()
				}
				
				if(length(TwaterUpstream) == 0 | length(qInUpstream) == 0){
					TwaterUpstream = 0

				} else if(sum(qInUpstream) == 0){

					TwaterUpstream = 0

				} else {
					
					TwaterUpstream = sum(qInUpstream*TwaterUpstream)/sum(qInUpstream)

				}
				#print(TwaterUpstream)

				#TwaterInitial = (TwaterUpstream(qOut - qLocal) + TwaterLocal*qLocal)/qOut
			}




			if(len/v < 1){

				rSqSubTT = len/(2*v) ##
				sRivTT = len/(2*v)

				if(ord > 1){
					###Get dimensionless factor, multiply times time to get TT?
					qInTT  = len/v ##gives fraction of timestep
				} else {
					qInTT = 0
				}

			} else {
				print("len/v > 1????")
				qInTT = 0
				rSqSubTT = 1-v/(2*len)
				sRivTT = 1 - v/len
			}

			#reachTravelTime = simFlow$TT[timeStep, hydroID] ## TT (hour) travel time of water through the subbasain 

			## Need to modify k ##

			#TairLocal = 5

			##A function of travel time and shading, given a large enough travel time, it should be close to one

			#if(TairLocal > 0){
			if(T){
				TwaterLocalWarmed = TwaterLocal + (TairLocal - TwaterLocal)*ifelse(K*rSqSubTT > 1, 1, K*rSqSubTT)
				TwaterQin = TwaterUpstream + (TairLocal - TwaterUpstream)*ifelse(K*qInTT > 1, 1, K*qInTT)
				TwaterSriv = TwaterOld + (TairLocal - TwaterOld)*ifelse(K*sRivTT > 1, 1, K*sRivTT)


				####OLD Twater code
				#Twater = TwaterInitial + (TairLocal - TwaterInitial) * K * reachTravelTime
			} else {
				epsilon = 2.5 ####Set to 2.5*C before sensitivity testing

				TwaterLocalWarmed = TwaterLocal + ((TairLocal + epsilon) - TwaterLocal)*K*rSqSubTT

				TwaterQin = TwaterUpstream + ((TairLocal + epsilon) - TwaterUpstream)*K*qInTT
		
				TwaterSriv = TwaterOld + ((TairLocal + epsilon) - TwaterOld)*K*sRivTT
				###Twater = TwaterInitial + ((TairLocal + epsilon) - TwaterInitial) * K * (reachTravelTime)
			}

			if(qIn == 0 & qLocal == 0 & sRiv == 0){
				Twater = 0
			} else {
				Twater = (TwaterQin*qIn + TwaterLocalWarmed*qLocal + TwaterSriv*sRiv)/(qIn+qLocal+sRiv)
			}

			Twater = ifelse(Twater < 0, .1, Twater)

			# Store values in results list
			results$Twater[timeStep, hydroID] = Twater


			if(debugMode == T | !is.null(outFile)){
				try(c(
					print(paste("Day:",timeStep, "Edge:", hydroID, "Order:", ord)),

					print(paste0("v = ", v)),
					print(paste0("len = ", len)),
					print(paste0("len/v = ", len/v)),
					print(paste0("lamda = ", lamda)),

					print(paste0("qSnow = ", qSnow)),
					print(paste0("qGw = ", qGw)),
					print(paste0("qNoSnow = ", qNoSnow)),
					print(paste0("qIn = ", qIn)),
					print(paste0("qLocal = ", qLocal)),
					print(paste0("sRiv = ", sRiv)),

					print(paste0("rSqSubTT= ", rSqSubTT*vConvFactor*1000/60/60)),
					print(paste0("sRivTT = ", sRivTT*vConvFactor*1000/60/60)),
					print(paste0("qInTT = ", qInTT*vConvFactor*1000/60/60)),
					print(paste0("Tsnow = ", Tsnow)),
					print(paste0("TairLocal = ", TairLocal)),
					print(paste0("Tgw = ", Tgw)),
					print(paste0("TairLag = ", TairLag)),
					print(paste0("TwaterLocal = ", TwaterLocal)),
					print(paste0("TwaterOld = ", TwaterOld)),
					print(paste0("TwaterUpstream = ", TwaterUpstream)),

					print(paste0("TwaterQin = ", TwaterQin)),
					print(paste0("TwaterLocalWarmed = ", TwaterLocalWarmed)),
					print(paste0("TwaterSriv = ", TwaterSriv)),

					print(paste("Twater = ", Twater))
				))
			}

			if(is.nan(Twater) | is.na(Twater)){
				warning("NaN or NA found!!")
			}

			if(outputExtraVars){
				results$TwaterLocal[timeStep, hydroID] = TwaterLocal
				results$TwaterLocalWarmed[timeStep, hydroID] = TwaterLocalWarmed
				results$TwaterQin[timeStep, hydroID] = TwaterQin
				results$TwaterSriv[timeStep, hydroID] = TwaterSriv
				results$kRsQsub[timeStep, hydroID] = K*rSqSubTT
				results$kQIn[timeStep, hydroID] = K*qInTT
				results$kQsRiv[timeStep, hydroID] = K*sRivTT
			}
			

			if(Twater > 50) stop()

		}

		stepsLooped = stepsLooped + 1


        if((stepsLooped %% etaInt) == 0){

			timeElapsed = as.numeric(format(Sys.time(), "%s")) - start
            print(paste("Seconds elapsed:", timeElapsed))

		secondsToFinish = (timeElapsed/stepsLooped * (timeLength - stepsLooped))

		message(paste("ETA:",
			round(secondsToFinish),
			"seconds or",
			round(secondsToFinish/60, digits=2),
			"minutes. Will finish at",
			secondsToFinish + Sys.time()
			))
		}
    }

	if(!is.null(prof)){
		Rprof(NULL)
	}

	if(!is.null(outFile)){
		close(out)
	}

	return(results)
}




evalCpp('
		#include <Rcpp.h>
		using namespace Rcpp;
		NumericMatrix xx(5, 4)
')

cppFunction('int add(int x, int y, int z) {
			  int sum = x + y + z;
			    return sum;
}')

evalCpp('NumericMatrix(5,5)(1,2)')
evalCpp('NumericMatrix(5,5)(1,(1,3))')

cppFunction('NumericVector meow(List parentList, int n){
			NumericVector a = parentList[n];
			int len = a.size();
			return len;
}')


meow(list(c(1,59,88), c(22,55), numeric(0)), 1)


#####create parent list of numers actually
