# function for water balance
wb <- function(t,State,pars,input) {
  return(with(as.list(c(State,pars)), {
    Q <- k1*S
    ET <- (1-exp(-k2*S))*5
    dS <- input(t) - ET - Q
    list(c(dS,ET,Q))}))
}

# pars are: vegpar.m, soilpar.m, model, delta.m, Z.prev, Z.m, fullday
# In is rainfall and maybe gwdepths?
# last_t_soil.sat is yini, maybe GWdepths is also yini


#WB_fun <- function(vegpar.m, In, fullday, last_t_soil.sat, 
#                   soilpar.m=soilpar, 
#                   model = model_in, 
#                   Zmean.m,
#                   deltat.m = deltat,
#                   Z.prev = (GWdepths[max(1,(fullday-1)),j]*-100),
#                   Z.m = GWdepths[fullday,j]*-100) {

WB_fun <- function(t, State, pars, input) {
  return(with(as.list(c(State,pars)), {

    # create internal parameters
    n <- soilpar.m$n
    Zr <- vegpar.m$Zr
    
    # Rainfall into infiltration
    Phi<-In
    surfoff <- 0
    # check if soil bucket full, create runoff
    if(Phi > (1 - last_t_soil.sat)){
      Phi <- 1 - last_t_soil.sat
      surfoff <- (In - Phi)*n*Zr
    } 
    # Call the loss model
    if (vegpar.m$DR==T) {
      loss <- do.call(model,list(s = last_t_soil.sat, 
                                 soilpar = soilpar.m, vegpar = vegpar.m,
                                 Z = Z.m, Zmean = Zmean.m, Z.prev = Z.prev))
      #print(loss[3])
    } else {
      loss <- do.call(model,list(s = last_t_soil.sat, Z = Z.m,
                                 soilpar = soilpar.m, vegpar = vegpar.m,
                                 Z.prev = Z.prev))
    }
    
  }
}
                   
                   
  n <- soilpar.m$n
  Zr <- vegpar.m$Zr
  
  # Rainfall into infiltration
  Phi<-In
  surfoff <- 0
  # check if soil bucket full, create runoff
  if(Phi > (1 - last_t_soil.sat)){
    Phi <- 1 - last_t_soil.sat
    surfoff <- (In - Phi)*n*Zr
  } 
  # Call the loss model
  if (vegpar.m$DR==T) {
    loss <- do.call(model,list(s = last_t_soil.sat, 
                               soilpar = soilpar.m, vegpar = vegpar.m,
                               Z = Z.m, Zmean = Zmean.m, Z.prev = Z.prev))
    #print(loss[3])
  } else {
    loss <- do.call(model,list(s = last_t_soil.sat, Z = Z.m,
                               soilpar = soilpar.m, vegpar = vegpar.m,
                               Z.prev = Z.prev))
  }
  # update the water balance
  # ss.t is temporary ss
  ss.t <- last_t_soil.sat + Phi - loss[1]/(n*Zr)*deltat.m
  # build in safety if capillary flux fills soil with dynamic water table
  if (ss.t > 1) {  # very wet conditions
    soil.sat <- 1
    qcap <- (loss[4] - (ss.t - 1)*n*Zr/deltat.m)*deltat.m
  } else { # normal conditions
    soil.sat<-ss.t
    qcap <- loss[4]*deltat.m
  }
  static_stress<- zeta(s=soil.sat,s_star=vegpar.m$s_star,s_w=vegpar.m$s_w,q=vegpar.m$q)
  Tg <- loss[3]*deltat.m
  #print(loss[3])
  Ts <- loss[2]*deltat.m
  Tt <- Tg + Ts
  Leakage <- ifelse(loss[5]>0,loss[5],0)*deltat.m
  
  smloss <- (loss[1]*deltat.m  - (qcap - loss[4]*deltat.m)) #increase with difference between originally calculated and new (pot.reduced) capillary flux
  GWrech <- Leakage-loss[3]*deltat.m-qcap
  #print(qcap)
  out <- list(s = soil.sat, qcap = qcap, Tgw = Tg, Tsoil = Ts, Ttotal = Tt, 
              Leakage = Leakage, GWrech = GWrech, smloss = smloss, 
              static_stress = static_stress, surfoff = surfoff)
  
}
