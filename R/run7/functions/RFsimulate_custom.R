RFsimulate_custom <- function (
    model, x, y = NULL, z = NULL, T = NULL, grid=NULL,
    distances, dim, data, given = NULL, err.model,
    params, err.params, n = 1, RFopt = NULL, ...
) {
  if (!missing(model) && is(model, "RFfit"))
    stop("To continue with the output of 'RFfit' use 'predict' or give the components separately")
  mc <- as.character(deparse(match.call()))
  
  ### create local RFopt if not provided ####################################
  if (is.null(RFopt)) {
    if (!missing(distances) && length(distances) > 0) {
      RFoptOld <- RandomFields:::internal.rfoptions(
        xyz_notation = length(y) != 0,
        expected_number_simu = n, ...,
        general.spConform = FALSE,
        RELAX = !missing(model) && is(model,"formula")
      )
    } else {
      RFoptOld <- RandomFields:::internal.rfoptions(
        xyz_notation = length(y) != 0,
        expected_number_simu = n,
        RELAX = !missing(model) && is(model,"formula")
      )
    }
    RFopt <- RFoptOld[[2]]
  }
  
  if (n > 2 && RFopt$internal$examples_reduced) {
    message("number of simulations reduced")
    n <- 2
  }
  
  cond.simu <- !missing(data) && !is.null(data) 
  reg <- RFopt$registers$register
  
  ### simulate from stored register ########################################
  mcall <- as.list(match.call(expand.dots=FALSE))
  if (length(mcall) == 1 ||
      (length(mcall) == 2 && !is.null(mcall$n)) ||
      (length(mcall) == 3 && !is.null(mcall$n) && "..." %in% names(mcall))) {
    if (cond.simu) {
      stop("repeated performance of conditional simulation not programmed yet")
    } else {
      res <- rfDoSimulate(
        n = n, reg = reg, spConform = RFopt$general$spConform
      )
      if (RFopt$general$returncall) attr(res, "call") <- mc
      attr(res, "coord_system") <- .Call(C_GetCoordSystem, reg,
                                         RFopt$coords$coord_system,
                                         RFopt$coords$new_coord_system)
      return(res)
    }
  }
  
  ### preparations #########################################################
  stopifnot(!missing(model) && !is.null(model))
  model.orig <- model
  model <- PrepareModel2(model, params=params, ...)
  err.model <- if (missing(err.model)) NULL else PrepareModel2(err.model, params=err.params, ...)
  
  ### conditional simulation ###############################################
  if (cond.simu) {
    if (isSpObj(data)) data <- sp2RF(data)
    stopifnot(missing(distances) || is.null(distances))
    res <- switch(GetProcessType(model),
                  RPgauss = 
                    rfCondGauss(model=model.orig, x=x, y=y, z=z, T=T,
                                grid=grid, n=n, data=data, given=given,
                                err.model=err.model,
                                predict_register = MODEL_PREDICT,
                                ...),
                  stop(GetProcessType(model),
                       ": conditional simulation of the process not programmed yet")
    )
  } else { ## unconditional simulation ####
    if(!is.null(err.model))
      warning("error model is unused in unconditional simulation")
    if (exists(".Random.seed") && !is.na(RFopt$basic$seed)) {
      .old.seed <- .Random.seed
      on.exit(set.seed(.old.seed), add = TRUE)
    }
    rfInit(model=list("Simulate",
                      setseed=eval(parse(text="quote(set.seed(seed=seed))")),
                      env=.GlobalEnv, model), 
           x=x, y=y, z=z, T=T,
           grid=grid, distances=distances, dim=dim, reg=reg, RFopt=RFopt,
           y.ok = FALSE)
    if (n < 1) return(NULL)
    res <- rfDoSimulate(n=n, reg=reg, spConform=FALSE)
  }
  
  ## output: RFsp   #################################
  if ((!cond.simu || (!missing(x) && length(x) != 0)) &&
      RFopt$general$spConform) {
    info <- RFgetModelInfo(if (cond.simu) MODEL_PREDICT else reg, level=3)
    if (length(res) > 1e7) {
      message("Too big data set (>1e7 entries). Returned with spConform=FALSE")
      return(res)
    }
    prep <- prepare4RFspDataFrame(info, RFopt, x, y, z, T, grid,
                                  coordnames = attributes(res)$coordnames)
    attributes(res)$varnames <- attributes(res)$varnames
    res <- conventional2RFspDataFrame(
      data=res, coords=prep$coords,
      gridTopology=prep$gridTopology, n=n,
      vdim=info$vdim, T=info$loc$T,
      vdim_close_together=RFopt$general$vdim_close_together
    )
    if (is.raster(x)) {
      res <- raster::raster(res)
      raster::projection(res) <- raster::projection(x)
    }
  }
  
  if (RFopt$general$returncall) attr(res, "call") <- mc
  attr(res, "coord_system") <- .Call(C_GetCoordSystem,
                                     as.integer(reg),
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)
  return(res)
}
