###################################
#Functions for routing LPF-GUESS Runoff
#
#Created by Jerad Hoy
#Date Created: 7/7/2015
#Date Last Modified:
#
#Sourced into streamNetwork.r file to do actual routing
#
makeHydrographs <- function(flowData, gaugeData, var=1, precip=NULL, spack=NULL, saveGraphs=F, plotTogether=F, interact=T, plotStats=T, plotSnowpack=F, plotPrecip=F, plotSeason=F, plotAnnual=F, dataMin=50, yLabel="Flow (m3/s)"){

    
    gauges <- gaugeData[names(gaugeData) %in% colnames(flowData[[var]])]

    if(length(gauges) >= 1){

		if(plotTogether){
			par(mfrow=c(length(gauges), 1))
		} else {
			#par(.pardefault)
		}
		for(i in 1:length(gauges)){

			if(saveGraphs){
				png(paste(plotDir, "/",gsub(" ", "_", gsub("[[:punct:]]", "", gauges@data[1, "SITENAME"])), ".png", sep=""))
			}

			if(timeStep == "month"){
				dates <- zoo::as.Date(zoo::as.yearmon(rownames(flowData[[var]])))
			} else {
				dates <- zoo::as.Date(rownames(flowData[[var]]))
			}
			if(plotSeason){

				par(mfrow=c(4,1), mar=c(3,4.9,4,2))
				for(season in 1:4)

				print("Plotting Seasonal")
				flow1 <- sumSeasonal(data.frame(flowData[[var]][, names(gauges)[i]]))
				flow1 <- flow1[seq(season, to=length(flow1), by=4)]
				gaugeFlow <- sumSeasonal(data.frame(gauges[[i]][,2], row.names=gauges[[i]][,1]))
				gaugeFlow <<- gaugeFlow[seq(season, to=length(gaugeFlow), by=4)]

				if(length(which(!is.na(gaugeFlow[,2]))) < dataMin){
				   print(paste("Only", length(which(!is.na(gaugeFlow[,2]))), "seasons of data, skipping gauge", i))
				   next
				}

				print("Plotting Seasonal")
				plot(seq(zoo::as.Date(simStartDate), length.out=length(flow1), by="year"), flow1, type="l", lwd=2, cex=1.3, cex.lab=1.4,  cex.axis=1.3, col="red", xlab="",  ylab=paste0(yLabel), ylim=c(0, max(c(max(flow1, na.rm=T), max(gaugeFlow, na.rm=T)))))
				print("Guage data")

				lines(seq(zoo::as.Date(simStartDate), length.out=length(gaugeFlow), by="year"), gaugeFlow, lwd=2)
				lines(gaugeFlow)
				print("Plotting Seasonal")
				#gofs <- CalcGOFs(flow1, gaugeFlow)
				#legend("topright", legend=paste(rownames(gofs), gofs), ncol=2, cex=.7)
				title(paste(c("MAM", "JJA", "SON", "DJF")[season], "for",  gaugesInBounds@data[gaugesInBounds@data[,8] == as.numeric(names(gauges[i])),"SITENAME"][1]), cex.main=1.3)

			} else if(plotAnnual){

				print("Plotting Annual")
				flow1 <- data.frame(flowData[[var]][, names(gauges)[i]])
				flow1 <- aggregate(flow1, list(substring(rownames(flow1), 4)), FUN=sum)
				gaugeFlow <- data.frame(gauges[[i]][,2], row.names=gauges[[i]][,1])
				gaugeFlow <- aggregate(gaugeFlow, list(substring(rownames(gaugeFlow), 4)), FUN=sum)

				if(length(which(!is.na(gaugeFlow[,2]))) < dataMin){
					print(paste("Only", length(which(!is.na(gaugeFlow[,2]))), "years of data, skipping gauge", i))
					next
				}
				print("got here")
				print(dim(flow1))
						plot(flow1, type="l", col="red", lwd=2, cex=1.3, cex.lab=1.3, cex.axis=1.3, xlab="",  ylab=expression(paste("Flow (m"^"3", "/s)")), ylim=c(0, max(flow1[,2])*1.05))
				print("got here")
						lines(gaugeFlow, lwd=2)
				print("got here")

			} else {    

				if(length(which(!is.na(gauges[[i]][,2]))) < dataMin){
				   print(paste("Only", length(which(!is.na(gauges[[i]][,2]))), "months of data, skipping gauge", i))
				   next
				}

				print("Plotting Monthly")
				plot(dates, flowData[[var]][, names(gauges)[i]], type="l", col="red", lwd=2, cex=1.3, cex.lab=1.3, cex.axis=1.3, xlab="",  ylab=paste0(yLabel))
				#plot(dates, flowData[[var]][, names(gauges)[i]], type="l", col="red", lwd=2, cex=1.3, cex.lab=1.3, cex.axis=1.3, xlab="",  ylab="Flow (m3/s)")
				lines(zoo::as.Date(gauges[[i]][,1]), gauges[[i]][,2], lwd=2)

			}

			if(plotSnowpack){
				par(new=T)
				plot(dates, snowpack, col="cyan", type="l", lty=3, lwd=2,axes=F, xlab=NA, ylab=NA)
			}
			if(plotPrecip){
				par(new=T)
				plot(dates, precip, col="cyan", type="l", lty=3, lwd=2,axes=F, xlab=NA, ylab=NA)
			}

			abline(0, 0)
			if(!plotSeason){
				title(gaugesInBounds@data[gaugesInBounds@data[,8] == as.numeric(names(gauges[i])),"SITENAME"][1], cex.main=1.3)
				legend("topleft", col=c("red", "black"), legend=c("Modelled LPJ-Guess Data", "Observations"), lty=1)
			}

			if(plotStats){
				gofs <<- CalcGOFs(flowData, gauges[i])
				gofs <<- gofs[c(1,4,6,7,9,10,11,17),]
				
				legend("topright", legend=paste(names(gofs), gofs), ncol=2, cex=1)
			}

			if(saveGraphs){
				dev.off()
			}

			print(paste("Gauge ", which(names(gaugeData) == names(gauges[i]))))

			if(interact){
				if(i < length(gauges)){
					answer <- readline("Press Enter To Continue")
					if(answer == ""){
					next
					} else {
					break
					}
				}
			}
			
		}
    } else {
		print("No gauges to plot!")
		print("Enter edgeIDs to plot! If plotting multiple, seperate by spaces.")
		edgeIdList <- readline()

		if(edgeIdList != ""){
			edgeIdList <- strsplit(edgeIdList, " ")[[1]]
			if(plotTogether){
				par(mfrow=c(length(edgeIdList), 1))
			} else {
				par(.pardefault)
			}
			for(i in 1:length(edgeIdList)){
				if(saveGraphs){
					png(paste(plotDir, "/",edgeIdList[i], ".png", sep=""))
				}

				if(timeStep == "month"){
					dates <- zoo::as.Date(zoo::as.yearmon(rownames(flowData[[var]])))
				} else {
					dates <- zoo::as.Date(rownames(flowData[[var]]))
				}

				plot(dates, flowData[[var]][, edgeIdList[i]], type="l", col="red", xlab="",  ylab="Flow (m/s)")
				abline(0, 0)

				title(paste("Edge", edgeIdList[i]))

				if(saveGraphs){
					dev.off()
				}
				if(interact){
					if(i < length(edgeIdList)){
						answer <- readline("Press Enter To Continue")
						if(answer == ""){
							next
						} else {
							break
						}
					}
				}
			}
	    }
    }
}

plotHydroVar <- function(flowData, gauges, hydroVar, edgeIdList=NULL, saveGraphs=F, plotTogether=F, interact=T){

    gauges <- gauges[names(gauges) %in% colnames(flowData[[hydroVar]])]

    if(length(gauges) >= 1){
		if(plotTogether){
			par(mfrow=c(length(gauges), 1))
		} else {
			par(.pardefault)
		}
		for(i in 1:length(gauges)){
			if(saveGraphs){
				png(paste(plotDir, "/",gsub(" ", "_", gsub("[[:punct:]]", "", gaugesInBounds@data[1, "SITENAME"])), ".png", sep=""))
			}

			if(timeStep == "month"){
				dates <- zoo::as.Date(zoo::as.yearmon(rownames(flowData[[hydroVar]])))
			} else {
				dates <- zoo::as.Date(rownames(flowData[[hydroVar]]))
			}
			plot(dates, flowData[[hydroVar]][, names(gauges)[i]], type="l", col="red", xlab="",  ylab=hydroVar)

			abline(0, 0)
			title(gaugesInBounds@data[gaugesInBounds@data[,8] == as.numeric(names(gauges[i])),"SITENAME"])
			legend("topleft", col=c("red", "black"), legend=c("Routed LPJ-Guess Runoff", "Gauge Data"), lty=1)

			if(saveGraphs){
				dev.off()
			}
			if(interact){
				if(i < length(gauges)){
					answer <- readline("Press Enter To Continue")
					if(answer == ""){
						next
					} else {
						break
					}
				}
			}
			
		}
    } else {
		if(is.null(edgeIdList)){
			print("No gauges to plot!")
			print("Enter edgeIDs to plot! If plotting multiple, seperate by spaces.")
			edgeIdList <- readline()
			if(edgeIdList != ""){
				edgeIdList <- strsplit(edgeIdList, " ")[[1]]
			}
		}
		if(plotTogether){
			par(mfrow=c(length(edgeIdList), 1))
		} else {
			par(.pardefault)
		}
		for(i in 1:length(edgeIdList)){
			if(saveGraphs){
				png(paste(plotDir, "/",edgeIdList[i], ".png", sep=""))
			}

			if(timeStep == "month"){
				dates <- zoo::as.Date(zoo::as.yearmon(rownames(flowData[[hydroVar]])))
			} else {
				dates <- zoo::as.Date(rownames(flowData[[hydroVar]]))
			}

			plot(dates, flowData[[hydroVar]][, edgeIdList[i]], type="l", col="red", xlab="",  ylab="Flow (m/s)")
			abline(0, 0)

			title(paste("Edge", edgeIdList[i]))

			if(saveGraphs){
				dev.off()
			}
			if(interact){
				if(i < length(edgeIdList)){
					answer <- readline("Press Enter To Continue")
					if(answer == ""){
					next
					} else {
					break
					}
				}
			}
		}
    }
}





CalcGOFs <- function(flowData, gauges, var=1){
    
    gauges <- gauges[names(gauges) %in% colnames(flowData[[var]])]

    flowData[[var]] <- flowData[[var]][1:nrow(gauges[[1]]), ]
     
    if(nrow(flowData[[var]]) != nrow(gauges[[1]])){
		print("rows don't match!")
		for(g in 1:length(gauges)){
			daysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
			out <- c()
			for(row in 1:(nrow(flow[[var]])-1)){
				print(row)
				print(col)
				out <- c(out, seq(flow[[var]][row, col], flow[[var]][row+1, col], length.out=daysInMonth[as.numeric(format(zoo::as.yearmon(rownames(flow[[var]])[row]), "%m"))]))

			}
			out <- c(out, seq(flow[[var]][row, col], 0, length.out=daysInMonth[as.numeric(format(zoo::as.yearmon(rownames(flow[[var]])[row]), "%m"))]+1))
		}
    }

    d <- c()

    for(i in 1:length(gauges)){
		gaugeGof <- tryCatch(gof(flowData[[var]][, names(gauges)[i]], gauges[[i]][,2], na.rm=T), error=function(e){NULL})
		if(is.null(gaugeGof)){
			next
		}
		d <- cbind(d,gof(flowData[[var]][, names(gauges)[i]], gauges[[i]][,2], na.rm=T))
    }

    #colnames(d) <- names(gauges)
    return(d)
}


makeTaylorDiagrams <- function(flowData, gauges, saveGraphs=F, interact=F){
	if(saveGraphs){
	    png(paste(plotDir, "/",gsub(" ", "_", gsub("[[:punct:]]", "", gaugesInBounds@data[1, "SITENAME"])), ".png", sep=""))
	}
colors <- rainbow(length(gauges))
	taylor.diagram(flowData[[var]][, names(gauges)[1]], gauges[[1]][,2], col=colors[1], pch=1)
	for(j in 2:length(gauges)){
	    print("plotting point")
	    taylor.diagram(flowData[[var]][, names(gauges)[j]], gauges[[j]][,2], add=T, col=colors[j], pch=length(gauges)%%25)
	}

	#title(@data[gaugesInBounds@data[,8] == as.numeric(names(gauges[i])),"SITENAME"])

	if(saveGraphs){
	    dev.off()
	}

}

sumSeasonal <- function(dat){
    i <- 1
    sumdData <- c()
        while(i <= (nrow(dat)-2)){
	    if(as.numeric(format(zoo::as.yearmon(rownames(dat)[i]), "%m")) %in% c(3,6,9,12)){
		sumdData <- c(sumdData, sum(dat[i:(i+2),1]))
		i <- i + 3
	    } else {
		i <- i + 1
	    }
    }
    return(sumdData)
}

