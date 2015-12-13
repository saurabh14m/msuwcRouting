# Convert to m/m and correct negative and zero slopes in edges and set to a slopeMin value
CorrectEdgeSlopes <- function(edges, slopeMin=.01){
   slopes <- edges[, edgeSlopeFieldDeg]/120000
   slopes[slopes <= 0] <- slopeMin 
   edges$Slope2 <- slopes
   return(edges)
}



# Goes down edges and sums upstream contributing area, used in later calculations
AssignContribArea <- function(edges, catchments){
    edges <- edges[order(edges[, edgeOrderField]),]

    edges$ContribArea <- NA
    
    for(i in 1:nrow(edges)){
        
        if(edges[i,edgeOrderField] == 1){
            edges[i,"ContribArea"]  <- catchments@data[which(catchments@data[,catchIdField] == edges[i, edgeIdField]), catchAreaField]

        } else {
	    

            edges[i,"ContribArea"] <- sum(edges[edges[, edgeNextDownField] == edges[i,edgeIdField],"ContribArea"]) + catchments@data[which(catchments@data[,catchIdField] == edges[i, edgeIdField]), catchAreaField] 
        }
        
    }

    edges$ContribArea <- edges$ContribArea*14400

    return(edges)
}

# Assigning bankfull depth, could be used for flood situation, but not currently used
#AssignBfDepth <- function(edges, a, b){

#    edges$bfDepth <- NA
    
#    for(i in 1:nrow(edges)){
#        edges[i,"bfDepth"] <- a*(edges[i,]$ContribArea)^b
#    }
#    return(edges)
#}


# Assigns bankfull width, but used as width for channel dimensions
AssignBfWidth <- function(edges, a, b){

    edges$bfWidth <- NA
    
    for(i in 1:nrow(edges)){
        edges[i,"bfWidth"] <- a*(edges[i,]$ContribArea)^b
    }
    return(edges)
}

# Assigns an 'a' coefficient for use in non-linear groundwater storage-discharge
AssignAcoeff <- function(edges, catchments, coeff){

    edges$aCoeff <- NA
    
    for(i in 1:nrow(edges)){
	edges[i,"aCoeff"] <- coeff * catchments@data[which(catchments@data[,catchIdField] == edges[i, edgeIdField]), catchAreaField] * 14400
    }

    return(edges)

}

