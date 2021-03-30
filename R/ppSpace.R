#' @title Spatial point process model
#' @name ppSpace
#'
#' @description Spatial point process model using INLA. This function is essentially a specialized wrapper over \code{inla}
#'
#' @param formula A formula that only relates the response \code{y} and some (or all) of the explanatory variables \code{X}. A paricularity of the is formula is that the response has to be defined as \code{y}.
#' @param sPoints A \code{SpatialPoint*} object that includes the sample location of the modelled species.
#' @param ppWeight An object of class \code{ppWeight}
#' @param explanaMesh An object of class \code{explanaMesh}
#' @param smooth A single value ranging between 0 and 2 passed to \code{inla.spde2.pcmatern}. It defines the smoothness of the Matern SPDE model. Default is set at 2.
#' @param prior.range A vector of length 2, with (range0, Prange) specifying that P(ρ < ρ_0) = p_ρ, where ρ is the spatial range of the random field. If Prange is NA, then range0 is used as a fixed range value. Default is c(0.05, 0.01).
#' @param prior.sigma A vector of length 2, with (sigma0, Psigma) specifying that P(σ > σ_0) = p_σ, where σ is the marginal standard deviation of the field. If Psigma is NA, then sigma0 is used as a fixed range value.  Default is c(1, 0.01).
#' @param many Logical. Whether the data in \code{sPoints} is large or not. See details. Default is \code{FALSE}.
#' @param fix A vector with the name of variables in the model that should be fixed to a given value when doing predictions. These values are used to map the intensities across the study area for a given value. Currently, the maximum of each variable is used as the fixed value, but it should be made more flexible in the future for example for playing more easily with climate change scenarios. Default is \code{NULL}, meaning no variables are fixed.
#' @param sboffset A character string with the name of the variable in the raster stack that should be used as an offset to scaled down the integration weights according to the level of effort across the study region. See details for further explanations. Default is \code{NULL}.
#' @param orthoCons Set to \code{TRUE} to force all the variance to go into the fixed effects. Sets constraints to have spatial field orthogonal to predictors. Experimental and currently not working...
#' @param \dots Arguments passed to \code{inla}
#'
#' @details 
#' 
#' If the argument \code{many = TRUE}, the estimation and the prediction will be carried out solely at the mesh edges, whereas when \code{many = FALSE} the estimation will be carried out at the mesh edges and at the sampled location. When the number of samples is very large (e.g. tens of thousands of samples or more) using \code{many = TRUE} can be much more computationally efficient. However, there is a precision trade-off. When \code{many = TRUE}, each sample is associated to an edge and the model is constructed using the number of samples associated to an edge as an importance value. In doing so, some spatial precision is lost at the expense of speed.
#' 
#'  The sampling bias offset argument \code{sboffset} is used to scaled down the weights (w) obtained from the dual mesh using a variable representing effort. This variable has to be a layer in the raster stack given for the predictors Specifically, values in the raster layer given will be summed for each polygon in the dual mesh to summarize the effort for each polygon. The extraction is made exact by using the \href{https://CRAN.R-project.org/package=exactextractr}{exactextractr} package. Once summed, values for each polygon (e) are 1) scaled with the weights, 2) rescaled between 0 and 1 and 3) multiplied with the original weights ((e/w) / max(e/w)) * w to adjust the weights in the integration mesh. This is an adaptation from Simpson et al. (2016). Note that polygons from the dual mesh that are partially overlapping the region of interest will get the weight associated with their area overlapping the study region and the effort considered is the effort associated with this overlapping area.
#'
#' @return
#' 
#' An object of class \code{ppSpace} that includes a model output, which is the model output of INLA.
#' 
#' In addition, it includes a series of attributes:
#' 
#'	  \item{\code{formula}}{The formula used to construct the model}
#'	  \item{\code{sPoints}}{A \code{SpatialPointDataFrame} object that includes the sample location and associated data of the modelled species.}
#'	  \item{\code{XEst}}{A matrix with all the explanatory variables used to construct the model. If there were factors in the original set of explanatory variables \code{X}, in \code{XEst}, they were decomposed into dummy variables. The values in \code{XEst} are the one from the sampled location.}
#'	  \item{\code{XPred}}{A matrix with all the explanatory variables used to construct the model. If there were factors in the original set of explanatory variables \code{X}, in \code{XPred}, they were decomposed into dummy variables. The values in \code{XPred} were gathered at the mesh edges. When \code{many = TRUE}, the values in \code{XPred} are exactly the same as the values in \code{XEst}}
#'	  \item{\code{mesh}}{An object of class \code{inla.mesh}. It is the mesh used to construct the model.}
#'	  \item{\code{Stack}}{An object of class \code{inla.data.stack}. It is a stack object constructed internally.}
#'
#' @references 
#' 
#' Simpson, D. Illian, J. B., Lindgren, F. Sørbye, S. H. and Rue, H. 2016. Going off grid: computationally efficient inference for log-Gaussian Cox processes. Biometrika, 103(1): 49-70 \url{https://doi.org/10.1093/biomet/asv064}
#'
#' @importFrom sp coordinates
#' @importFrom raster rasterFromXYZ
#' @importFrom sp SpatialPoints
#' @importFrom raster extract
#' @importFrom INLA inla
#' @importFrom INLA inla.spde2.pcmatern
#' @importFrom INLA inla.spde.make.A
#' @importFrom INLA inla.spde.make.index
#' @importFrom INLA inla.stack
#' @importFrom INLA inla.stack.data
#' @importFrom INLA inla.stack.A
#' @importFrom Matrix Diagonal
#' @importFrom stats model.matrix
#' @importFrom stats model.frame
#' @importFrom sf st_as_sf
#' @importFrom exactextractr exact_extract
#' 
#' @export
#' 
#' @keywords models
ppSpace <- function(formula,
                    sPoints, 
                    ppWeight, 
                    explanaMesh, 
                    smooth = 2,
                    prior.range = c(0.05, 0.01),
                    prior.sigma = c(1, 0.01), 
                    many = FALSE,
                    fix = NULL,
                    sboffset = NULL,
                    orthoCons = FALSE,
                    ...){

  #==============
  ### Basic check
  #==============
  if(as.character(formula[[2]]) != "y"){
    stop("'y' should be used to define the response variable")
  }
  
  ### Check if the mesh in ppWeight and dataPred are the same
  if(!identical(attributes(ppWeight)$mesh$graph, explanaMesh$mesh$graph)){
    stop("'ppWeight' and 'explanaMesh' were constructed using different mesh")
  }
  
  #================
  ### Basic objects
  #================
  nsPoints <- length(sPoints)
  nEdges <- explanaMesh$mesh$n
  xy <- coordinates(sPoints)
  colnames(xy) <- c("x","y") 

  #==============
  ### Define SPDE
  #==============
  SPDE <- inla.spde2.pcmatern(mesh=explanaMesh$mesh, alpha=smooth,
                              prior.range=prior.range,
                              prior.sigma=prior.sigma)

  #========================================================
  ### Rescale weights if sampling bias offset is included
  #========================================================
  
  if(!is.null(sboffset)){
    polys <- gIntersection(explanaMesh$sPoly,attributes(ppWeight)$dmesh,byid=TRUE)
    e <- exact_extract(explanaMesh$X[[sboffset]], 
                       st_as_sf(polys), 
                       fun = function(values, coverage){
                         sum(values * coverage, na.rm = TRUE)
                       },progress = FALSE)
    k <- ppWeight > 0
    ppWeight[k] <- ppWeight[k] * ((e/ppWeight[k])/max(e/ppWeight[k]))
  }
  
  
  #==========================
  ### Define response objects
  #==========================
  if(many){
    ### Aggregate spatial data
    xyDF <- as.data.frame(xy)
    spaceAgg <- aggData(xyDF, explanaMesh$mesh)
    
    ### Pseudo-absences are the number of edges on the mesh
    ### Occurences are the number of points
    yPP <- spaceAgg$Freq
    
    ### weight associated to pseudo-absences (w) and occurrences (0)
    ePP <- ppWeight[spaceAgg$space]
  }else{
    ### Pseudo-absences are the number of edges on the mesh
    ### Occurences are the number of points
    yPP <- rep(0:1, c(nEdges, nsPoints))
    
    ### weight associated to pseudo-absences (w) and occurrences (0)
    ePP <- c(ppWeight, rep(0, nsPoints))
  }
  
  #===========================
  ### Define covariate objects
  #===========================
  ### Organize data into a data.frame
  refData <- as.data.frame(values(explanaMesh$X))
  
  Xfactor <- unlist(lapply(explanaMesh$Xmesh, is.factor))
  if(any(Xfactor)){
    for(i in which(Xfactor)){
      refData[,i] <- as.factor(refData[,i])
      levels(refData[,i]) <- levels(explanaMesh$Xmesh[,i])
    }
  }
  
  refData <- data.frame(y = 1, refData)
  
  ### Organize X so that it follows the formula
  Xorg <- model.matrix(formula,model.frame(formula, 
                                           data = refData, 
                                           na.action = NULL))[,-1,drop=FALSE]
  
  ### Construct a brick out of Xorg
  xyXorg <- cbind(coordinates(explanaMesh$X),Xorg)
  Xbrick <- rasterFromXYZ(xyXorg)
  if(nlayers(Xbrick) == 1){
    Xbrick <- brick(Xbrick)  
  }
  names(Xbrick) <- colnames(Xorg)
  
  if(many){
    ### Extract covariate values for model estimation
    locEst <- SpatialPoints(coords = explanaMesh$mesh$loc[,1:2])
    XEst <- extract(Xbrick, locEst)
    XPred <- XEst
  }else{
    ### Extract covariate values for model estimation
    meshxy <- rbind(explanaMesh$mesh$loc[,1:2],xy)
    rownames(meshxy) <- 1:nrow(meshxy)
    
    locEst <- SpatialPoints(coords = meshxy)
    XEst <- extract(Xbrick, locEst)
    
    ### Extract covariate values for model prediction
    locPred <- SpatialPoints(coords = explanaMesh$mesh$loc[,1:2])
    XPred <- explanaMesh$Xmesh
  }
  
  #==================================================
  ### Fix given predictors
  #==================================================
  if(!is.null(fix)){
    m <- match(fix, colnames(XPred))
    v <- apply(XPred[, m, drop=FALSE], 2, max, na.rm = TRUE)
    XPred[, m] <- rep(v, each = nrow(XPred))
  }
  
  #==================================================
  ### Define projection matrix and build stack object
  #==================================================
  if(many){
    ### Projection matrix
    ProjInfer <- inla.spde.make.A(explanaMesh$mesh,
                                  loc = explanaMesh$mesh$loc[spaceAgg$space,])
    
    IDSpace <- inla.spde.make.index('i', nEdges)
  
    ### Build stack objects
    StackEst <- inla.stack(data = list(y = yPP, e = ePP), 
                           A = list(1, ProjInfer), 
                           effects = list(c(list(Intercept = 1), 
                                               asplit(XEst, 2)),
                                          IDSpace), 
                           tag = "est")
    
    StackPred <- inla.stack(data = list(y = NA, e = NA),
                            A = list(1, ProjInfer), 
                            effects = list(c(list(Intercept = 1), 
                                                asplit(XPred, 2)), 
                                           IDSpace),
                            tag = "pred")
  }else{
    #--------------------
    ### Projection matrix
    #--------------------
    ### For inference
    ProjInfer <- inla.spde.make.A(explanaMesh$mesh, xy)
    
    ### For integration
    ProjInter <- Diagonal(nEdges, rep(1, nEdges))
    
    ### Combine both projection matrix
    A <- rbind(ProjInter, ProjInfer)

    #----------------------
    ### Build stack objects
    #----------------------
    StackEst <- inla.stack(data = list(y = yPP, e = ePP), 
                           A = list(1, A), 
                           effects = list(c(list(Intercept = 1), 
                                               asplit(XEst, 2)), 
                                          list(i = 1:nEdges)), 
                           tag = "est")
    
    StackPred <- inla.stack(data = list(y = NA, e = NA),
                                A = list(1, ProjInter), 
                                effects = list(c(list(Intercept = 1), 
                                                    asplit(XPred, 2)), 
                                               list(i = 1:nEdges)),
                                tag = "pred")
  }
  
  Stack <- inla.stack(StackEst, StackPred)
  
  #===============
  ### Build models
  #===============
  
  ### Make variables names explicit in formula
  X <- paste(colnames(XEst),collapse=" + ")
  fixed <- paste("y ~ 0 + Intercept +",X)
  
  if(orthoCons){
    # build constraints
    XX = cbind(rep(1, nrow(XEst)), XEst)
    Q = qr.Q(qr(XX))
    AA <- as.matrix(t(Q)%*%ProjInfer)
    ee <- rep(0, ncol(XX))
    formule <- formula(paste(fixed,"f(i, model=SPDE, extraconstr = list(A = AA, e = ee))", sep=" + "))
  }else{
    formule <- formula(paste(fixed,"f(i, model=SPDE)", sep=" + "))
  }
  
  model <- inla(formule, family = "poisson", 
                data = inla.stack.data(Stack),
                control.predictor = list(A = inla.stack.A(Stack), 
                                         link = 1),
                E = inla.stack.data(Stack)$e, ...)
  
  nameRes <- names(model)
  
  #===============
  ### Return model
  #===============
  ### add a series of attributes to the result object
  if(many){
    attributes(model) <- list(formula = formula,
                              sPoints = sPoints,
                              XEst = XEst,
                              XPred = XPred,
                              mesh = explanaMesh$mesh,
                              Stack = Stack)
  }else{
    attributes(model) <- list(formula = formula,
                              sPoints = sPoints,
                              XEst = XEst,
                              XPred = XPred,
                              mesh = explanaMesh$mesh,
                              Stack = Stack)
  }
  
  names(model) <- nameRes
  
  class(model) <- "ppSpace"
  
  return(model)
}
