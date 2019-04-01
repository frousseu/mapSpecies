#' @title Univariate spatial model
#' @name uniSpace
#'
#' @description Spatial model for presence absence data using INLA. This function is essentially a sophisticated wrapper over \code{inla}
#'
#' @param formula A formula that only relates the response \code{y} and some (or all) of the explanatory variables \code{X}. A paricularity of the is formula is that the response has to be defined as \code{y}.
#' @param spdf A SpatialPointsDataFrame with a vector of 0 and 1 in the data.frame part \code{stack} (or \code{brick}) combining all of the explanatory variables to consider.
#' @param explanaMesh An object of class \code{explanaMesh}
#' @param family A character string describing the error distribution to be used when constructing the model. 
#' @param link A character string describing the link function to be used when constructing the model.
#' @param smooth A single value ranging between 0 and 2 passed to \code{inla.spde2.pcmatern}. It defines the smoothness of the Matern SPDE model. Default is set at 2.
#' @param prior.range A vector of length 2, with (range0, Prange) specifying that P(ρ < ρ_0) = p_ρ, where ρ is the spatial range of the random field. If Prange is NA, then range0 is used as a fixed range value. Default is c(0.05, 0.01).
#' @param prior.sigma A vector of length 2, with (sigma0, Psigma) specifying that P(σ > σ_0) = p_σ, where σ is the marginal standard deviation of the field. If Psigma is NA, then sigma0 is used as a fixed range value.  Default is c(1, 0.01).
#' @param \dots Arguments passed to \code{inla}
#'
#' @details 
#'
#' The underlying model used by this function is a generalized spatial linear model using the cloglog link function. The idea to use the cloglog link function instead of another link function (e.g. logit or probit) is that it is more flexible.
#'
#' @importFrom INLA inla.spde2.pcmatern
#' @importFrom sp coordinates
#' @importFrom raster rasterFromXYZ
#' @importFrom sp SpatialPoints
#' @importFrom raster extract
#' @importFrom INLA inla.spde.make.A
#' @importFrom INLA inla.stack
#' @importFrom INLA inla
#'
#' @export
#'
#' @keywords models
uniSpace <- function(formula, spdf,  explanaMesh,
                    family = NULL,
                    link = NULL,
                    smooth = 2,
                    prior.range = c(0.05, 0.01),
                    prior.sigma = c(1, 0.01), ...){

  
  #==============
  ### Basic check
  #==============
  if(is.null(family) | is.null(link)){
    stop("Either 'family' or 'link' (or both) need to be specified")
  }
  
  #================
  ### Basic objects
  #================
  nsmpl <- length(spdf)
  nEdges <- explanaMesh$mesh$n

  #==============
  ### Define SPDE
  #==============
  SPDE <- inla.spde2.pcmatern(mesh=explanaMesh$mesh,
                              alpha=smooth,
                              prior.range=prior.range,
                              prior.sigma=prior.sigma)

  #==========================
  ### Define response objects
  #==========================
  ### Coordinates use for the estimation
  xy <- coordinates(spdf)

  #===========================
  ### Define covariate objects
  #===========================
  ### Organize data into a data.frame
  ### Note that y here is bogus, it will be removed later
  refData <- data.frame(y = 1, values(explanaMesh$X))
  colnames(refData)[1] <- names(spdf)

  ### Organize X so that it follows the formula
  Xorg <- model.matrix(formula,model.frame(formula,
                                           data = refData,
                                           na.action = NULL))[,-1]

  ### Construct a brick out of Xorg
  xyXorg <- cbind(coordinates(explanaMesh$X),Xorg)
  Xbrick <- rasterFromXYZ(xyXorg)

  ### Extract covariate values for model estimation
  locEst <- SpatialPoints(coords = xy)
  XEst <- extract(Xbrick, locEst)

  ### Extract covariate values for model prediction
  locPred <- SpatialPoints(coords = explanaMesh$mesh$loc[,1:2])
  XPred <- explanaMesh$Xmesh

  #=====================
  ### Construct A matrix
  #=====================
  ### For estimation
  AEst <- inla.spde.make.A(mesh = explanaMesh$mesh, loc = xy)

  ### For prediction
  APred<-inla.spde.make.A(mesh = explanaMesh$mesh)

  #====================================
  ### Build stack object for estimation
  #====================================
  ### For estimation
  resp <- unlist(spdf@data)
  names(resp) <- NULL
  resp <- list(y = resp)
  names(resp) <- names(spdf)

  StackEst <- inla.stack(data = resp, A = list(AEst, 1),
                         effects = list(list(i = 1:nEdges),
                                             list(Intercept = 1,
                                             X = XEst)),
                         tag = "est")

  ### For prediction
  respNA <- list(y = NA)
  names(respNA) <- names(spdf)

  StackPred <- inla.stack(data = respNA, A = list(APred, 1),
                          effects = list(list(i = 1:nEdges),
                                         list(Intercept = 1,
                                              X = XPred)),
                          tag = "pred")

  ### Combine both stack objects
  Stack <- inla.stack(StackEst, StackPred)

  #===============
  ### Build models
  #===============
  formule <- formula(paste(names(spdf) ,"~ 0 + Intercept + X + f(i, model=SPDE)"))

  model <- inla(formule, family = family,
                data = inla.stack.data(Stack),
                control.link = list(link = link),
                control.predictor = list(A = inla.stack.A(Stack),
                                         link = 1),
                E = inla.stack.data(Stack)$e, ...)

  ### Return model
  res <- list(model = model, Stack = Stack, mesh = explanaMesh$mesh)
  class(res) <- "spatialPA"

  return(res)
}