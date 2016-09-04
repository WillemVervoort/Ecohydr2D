## copied from Hydromad: 
## https://github.com/josephguillaume/hydromad/blob/master/R/options.R
##
## Willem Vervoort

.defaultecohydro2dOptions <- function() {
  list(
    # ----------------
    # Grid
    # ---------------
    NX = 12,
    NY = 1,
    NRBL = 2,    #number of river blocks
    IXR = c(1,12), # location river blocks in X and Y dir
    IYR = c(1,1),
    DELX = c(20,seq(20,200,by=20),20),
    DELY = 100,   
    ITIM = 1,
    # --------------------
    # Groundwater
    # --------------------
    # Depth of water table:
    bottom = -25,   #meters
    GWthreshold = 0.0, #m 
    DELTcrit = 21,   #days
    dslope_x = 0, # m difference in x-direction
    dslope_y = 0, # m difference in y-direction
    # -------------------
    ########################
    #river: 
    ########################
    RES = 100
  )
}
    
#3 hydromad says:    
## code below copied from lattice

ecohydro2d.getOption <- function(name)
{
  .ecohydro2dEnv$options[[name]]
}

ecohydro2d.options <- function(...)
{
  ## this would have been really simple if only form allowed were
  ## lattice.options("foo", "bar") and
  ## lattice.options(foo=1, bar=2). But it could also be
  ## lattice.options(foo=1, "bar"), which makes some juggling necessary
  
  new <- list(...)
    if (is.null(names(new)) && length(new) == 1 && is.list(new[[1]])) new <- new[[1]]
  old <- .ecohydro2dEnv$options
  
  ## if no args supplied, returns full options list
  if (length(new) == 0) return(old)
  
  nm <- names(new)
  if (is.null(nm)) return(old[unlist(new)]) ## typically getting options, not setting
  isNamed <- nm != "" ## typically all named when setting, but could have mix
  if (any(!isNamed)) nm[!isNamed] <- unlist(new[!isNamed])
  
  ## so now everything has non-"" names, but only the isNamed ones should be set
  ## everything should be returned, however
  
  retVal <- old[nm]
  names(retVal) <- nm
  nm <- nm[isNamed]
  
  ## this used to be
  
  ## modified <- updateList(retVal[nm], new[nm])
  ## .LatticeEnv$lattice.options[names(modified)] <- modified
  
  ## but then calling lattice.options(foo = NULL) had no effect
  ## because foo would be missing from modified.  So, we now do:
  
  updateList <- function (x, val) {
    if (is.null(x)) x <- list()
    utils::modifyList(x, val)
  }
  .ecohydro2dEnv$options <- updateList(old, new[nm])
  
  ## return changed entries invisibly
  invisible(retVal)
}