#' @title Construct species distribution maps models
#' 
#' @description Constructs mean, standard deviation, and quantile (0.025, 0.5 and 0.975) maps for models calculated using \code{\link{uniSpace}} and \code{\link{ppSpace}} 
#' 
#' @param modelSpace An object of class \code{ppSpace} or of class \code{uniSpace}.
#' @param dims A vector of length 2 defining the number of pixels to use as rows and columns to define the map.
#' @param type Either "mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode" or "space". Defines the map to be drawn. If "space", maps the mean of the spatial field.
#' @param sPoly A spatial polygon to isolate the region of interest. If none is given, a map is drawn for the entire region covered by the mesh. 
#' 
#' @importFrom INLA inla.mesh.projector
#' @importFrom INLA inla.mesh.project
#' @importFrom INLA inla.stack.index
#' @importFrom raster raster
#' @importFrom raster mask
#' @importFrom raster xmin
#' @importFrom raster ymin
#' @importFrom raster xmax
#' @importFrom raster ymax
#'
#' @keywords hplot
#' 
#' @export
#' 
mapSpace <- function(modelSpace, dims, 
                       type = c("mean", "sd", "0.025quant", 
                                "0.5quant", "0.975quant",
                                "mode","space"), sPoly = NULL){
  ### General check
  if(length(type) > 1){
    stop("Only one type should be defined")
  }
  
  ### Define map basis
  if(is.null(sPoly)){
    mapBasis <- inla.mesh.projector(attributes(modelSpace)$mesh,
                                    dims = dims,
                                    crs = attributes(modelSpace)$mesh$crs)
  }else{
    mapBasis <- inla.mesh.projector(attributes(modelSpace)$mesh,
                                    dims = dims,
                                    xlim = c(xmin(sPoly), xmax(sPoly)),
                                    ylim = c(ymin(sPoly), ymax(sPoly)),
                                    crs = attributes(modelSpace)$mesh$crs)
  }
  
  
  if(type=="space"){
    
    mapPred <- inla.mesh.project(mapBasis, 
                                 modelSpace$summary.random[['i']][['mean']])
  }else{
  
    ### Find the mesh edges on which predictions should be made
    ID <- inla.stack.index(attributes(modelSpace)$Stack, tag="pred")$data
  
    ### Calculate prediction
    mapPred <- inla.mesh.project(mapBasis, 
                               modelSpace$summary.fitted.values[[type]][ID])
  
  }
  ### Transform map into a raster
  mapRaster <- raster(t(mapPred[,ncol(mapPred):1]),
                      xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                      ymn = min(mapBasis$y), ymx = max(mapBasis$y))
  
  ### Return
  return(mapRaster)
}