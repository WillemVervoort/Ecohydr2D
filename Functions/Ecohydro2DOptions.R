## copied from Hydromad: 
## https://github.com/josephguillaume/hydromad/blob/master/R/options.R
##

.defaultEcoHydro2DOptions <- function() {
  list(
    # ----------------
    # Grid
    # ---------------
    NX = 12,
    NY = 1,
    NRBL = 0,    #number of river blocks
    DELX = rep(10,NX),
    DELY = 100,   
    ITIM = 1,
    # --------------------
    # Groundwater
    # --------------------
    # Depth of water table:
    init_heads = rep(2,NX), #meters
    bottom = -25,   #meters
    bottommatrix = matrix(bottom,NX,NY),
    GWthreshold = 0.0, #m 
    DELTcrit = 21,   #days
    # -------------------
    ########################
    #river: 
    ########################
    RES = NULL,
    criver = c(RES,0.1), #resistance (m) # WV says Increase this to limit leakage from river
    Ariver <- rep(100,NRBL), # river area (m^2) # this is a guess
#    if (NRBL > 0) {
      riverheads <- rep(2,NRBL) # water head river (m). 
 #   } else {
  #    hriver=0 # should have some value, otherwise: formula breaks off
   # }
  )
}
    
#3 hydromad says:    
## code below copied from lattice

EcoHydro2D.getOption <- function(name)
{
  .EcoHydro2DEnv$options[[name]]
}

EcoHydro2D.options <- function(...)
{
  ## this would have been really simple if only form allowed were
  ## lattice.options("foo", "bar") and
  ## lattice.options(foo=1, bar=2). But it could also be
  ## lattice.options(foo=1, "bar"), which makes some juggling necessary
  
  new <- list(...)
    if (is.null(names(new)) && length(new) == 1 && is.list(new[[1]])) new <- new[[1]]
  old <- .EcoHydro2DEnv$options
  
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
  .EcoHydro2DEnv$options <- updateList(old, new[nm])
  
  ## return changed entries invisibly
  invisible(retVal)
}