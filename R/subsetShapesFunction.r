##################################
#Functions for subsetting catchments and edges
#
#Created by Jerad Hoy
#Date Created: 7/7/2015
#Date Last Modified:
#
#Sourced into streamNetwork.r file to do actual routing
#

#' Get Parents
#'
#' Returns the parents of an edge or catchment given it's id
#'
#' @param shapeFrame a spatialFrame of edges or catchments
#' @param id ID of the feature that you are trying to find the parents of
#' @param idField the field that shape IDs are stored in. Defaulted to edgeIdField, a parameter set in the beginning of routeWater.r
#' @param nextDownField the field that nextDownId of the shapes are stored in. This is a field created automatically in ArcGIS processing. Defaulted to edgeNextDownField, a parameter set in the beginning of routeWater.r
#'
#' @return the parents of the edge being used
#'
#' @examples
#' getParents(edgesInBounds, "34342")
#'
#' @name getParents
#' @export
getParents <- function(shapeFrame, id, idField=edgeIdField, nextDownField=edgeNextDownField){
            return(shapeFrame[shapeFrame[, nextDownField] == id, idField])
}

#' Get Order
#'
#' Returns the order of an edge or catchment given it's id
#'
#' @param shapeFrame a spatialFrame of edges or catchments
#' @param id ID of the feature that you are trying to find the order of
#' @param idField the field that shape IDs are stored in. Defaulted to edgeIdField, a parameter set in the beginning of routeWater.r
#' @param orderField the field that stream order of the shapes are stored in. This is a field created automatically in ArcGIS processing. Defaulted to edgeNextDownField, a parameter set in the beginning of routeWater.
#'
#' @return the order of the edge being used
#'
#' @examples
#' getOrder(edgesInBounds, "34342")
#'
#' @name getOrder
#' @export
getOrder <- function(shapeFrame, id, idField=edgeIdField, orderField=edgeOrderField){
	if(typeof(shapeFrame) == "S4"){
		shapeFrame <- shapeFrame@data
	}

	return(shapeFrame[shapeFrame[, idField] == id, orderField])
}

#' Find All Parents
#'
#' Recursive fucntion that returns all the parents of an edge or catchment given it's ID
#'
#' @param shapeFrame a spatialFrame of edges or catchments
#' @param id ID of the feature that you are trying to return all the parents of
#' @param idField the field that shape IDs are stored in. Defaulted to edgeIdField, a parameter set in the beginning of routeWater.r
#' @param orderField the field that stream order of the shapes are stored in. This is a field created automatically in ArcGIS processing. Defaulted to edgeNextDownField, a parameter set in the beginning of routeWater.
#'
#' @return the order of the edge being used
#'
#' @examples
#' getOrder(edgesInBounds, "34342")
#'
#' @name findAllParents
#' @export
findAllParents <- function(shapeFrame, ID, createArray=T, idField=edgeIdField, nextDownField=edgeNextDownField, orderField=edgeOrderField){
	#print(typeof(shapeFrame))
	if(typeof(shapeFrame) == "S4"){
		shapeFrame <- shapeFrame@data
	}

	if(createArray == T){
		IDarray <<- c()
	}

	for(i in 1:length(ID)){

		if(shapeFrame[shapeFrame[, idField] == ID[i], orderField] > 1){

			#parents <- getParents(shapeFrame, ID[i], idField, nextDownField)
			parents <- shapeFrame[shapeFrame[, nextDownField] == ID[i], idField]
			#print(parents)
			#print(IDarray)
			IDarray <<- c(IDarray, parents)
			#idIndex <- 1
			#print(IDarray)

			#readline()
			for(j in 1:length(parents)){

				if(shapeFrame[shapeFrame[, idField] == parents[j], orderField] > 1){

					findAllParents(shapeFrame, parents[j], createArray=F, idField, nextDownField)
				}
			}
		}
	}
	return(IDarray)
}


# Subesets catchments by HUC10 codes, ignores NA HUC10 codes
#' Fill Data Gaps
#'
#' Fills gaps in daily stream gaugedata with NA so that statistics can later be calculated with data. Function only used by getGaugeData
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
#' @name cleanDat
#' @export
GetShapesInBounds <- function(shapeFrame, hucCodes, by="HUC"){

    if(shapeFrame@class[1] == "SpatialLinesDataFrame"){
	hucField <- edgeHucField
	orderField <- edgeOrderField
	nextDownField <- edgeNextDownField
	idField <- edgeIdField
    } else {
	hucField <- catchHucField
	orderField <- catchOrderField
	nextDownField <- catchNextDownField
	idField <- catchIdField
    }
    print(idField)


    shapesInHuc <- shapeFrame[shapeFrame@data[, hucField] %in% hucCodes,]
    print(nrow(shapesInHuc))
    
    maxOrderId <- shapesInHuc@data[which(shapesInHuc@data[, orderField] == max(shapesInHuc@data[, orderField])), idField]

    print(maxOrderId)

    maxOrderId <- shapesInHuc@data[shapesInHuc@data[, nextDownField] == shapesInHuc@data[shapesInHuc@data[, idField] == maxOrderId , nextDownField], idField]

    print(maxOrderId)



    shapeIDsInBounds <- c(maxOrderId, findAllParents(shapeFrame, maxOrderId, createArray=T, idField, nextDownField, orderField))



    if(nrow(shapesInHuc) != length(shapeIDsInBounds)){
	warning(paste("Some headwater shapes not in selected HUC. Adding", length(shapeIDsInBounds)-nrow(shapesInHuc), "shapes."))
    }

	shapesInBounds <- shapeFrame[shapeFrame@data[, idField] %in% shapeIDsInBounds,]

    return(shapesInBounds)
}


#' Fill Data Gaps
#'
#' Fills gaps in daily stream gaugedata with NA so that statistics can later be calculated with data. Function only used by getGaugeData
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
#' @name cleanDat
#' @export
GetShapesById <- function(shapeFrame, IDs){

    if(shapeFrame@class[1] == "SpatialLinesDataFrame"){
	hucField <- edgeHucField
	orderField <- edgeOrderField
	nextDownField <- edgeNextDownField
	idField <- edgeIdField
    } else {
	hucField <- catchHucField
	orderField <- catchOrderField
	nextDownField <- catchNextDownField
	idField <- catchIdField
    }

    shapeIDsInBounds <- c(IDs, findAllParents(shapeFrame, IDs, createArray=T, idField, nextDownField, orderField))

    shapesInBounds <- shapeFrame[shapeFrame@data[, idField] %in% shapeIDsInBounds,]


    return(shapesInBounds)
}
