# Adjusted version by Willem to include VvdZ2008 function. 


# Evaporation function following Teuling and Troch 2005
E_Teuling <- function(s,vegpar) {
    if (s <= vegpar$s_w) {
      beta_T <- 0
    } else {
      if (s <= vegpar$s_star & s > vegpar$s_w) {
        beta_T <- (s-vegpar$s_w)/(vegpar$s_star-vegpar$s_w)
      } else {
        beta_T <- 1
      }
    }
    E <- vegpar$fr*beta_T*(1-exp(-vegpar$c_T*vegpar$LAI))*vegpar$Ep
    return(E)  
    }

# Model VvdZ2008 reworked to fit framework
#This version has 3 limits
rho_new_1<-function(s,Z,soilpar,vegpar) {
      Zr <- vegpar$Zr
      hb <- soilpar$psi_s_bar*-10^4
      # introduced fix to makes sure s.lim bounded at 1
      s.lim <- ifelse(Z > hb,(Z/hb)^(-1/soilpar$b),1)
      G1 <- G(soilpar$b,hb,Z) # using the G function
      m <- soilpar$K_s/(soilpar$n*Zr*(exp(soilpar$beta*(1-s.lim))-1))
      m2 <- soilpar$K_s*G1/(soilpar$n*Zr)
      m1 <- soilpar$K_s*G1/(soilpar$n*Zr*
          (1-exp(soilpar$beta*(vegpar$s_star-s.lim))))
          #(exp(soilpar$beta*(soilpar$s_fc-vegpar$s_star))-1))
          
      eta<-vegpar$E_max/(soilpar$n*Zr)

      r <- 0
 
      if (s > s.lim) {
        Q <- m*(exp(soilpar$beta*(s-s.lim))-1)
        E <- eta
      }
      if (m1 < eta) {
      # define s_cr: s-critical
        s_cr <- m2/eta*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
          if (s > s_cr && s<=vegpar$s_star)
          {
            E <- eta*((s-s_cr)/(vegpar$s_star-s_cr))
            Q <- -m2*((s-s_cr)/(vegpar$s_star-s_cr))
          } else {
            if (s > vegpar$s_star && s <= s.lim) {
               Q <- -m1*(1-exp(soilpar$beta*(s-s.lim))) #(exp(soilpar$beta*(soilpar$s_fc-s))-1)
               E <- eta
            }
          }
      } else {
        # fixed problem with shallow water tables 26/10/2011 WV
          if (s > vegpar$s_star && s <= s.lim) {
            E <- eta
        # fixed problem with shallow water tables 26/10/2011 WV
        # used to be Q <- - m1*exp(soilpar$beta*(s-s.lim))
        Q <- eta - m1*exp(soilpar$beta*(s-s.lim))
          }
      }
      L <- Q + E
      qcap <- ifelse(Q < 0, abs(Q),0)
      return(c(c(L,E,0,qcap,Q))
}




# Model Eagleson capillary flux equation (Eagleson model)
Fun_E <- function(s,Z,soilpar,vegpar) {
  hb <- -10^4*soilpar$psi_s_bar
  b1 <- 2+3/soilpar$b
  a1 <- 1+(3/2)/(b1-1)
  H1 <- a1*(hb/(Z-vegpar$Zr))^b1
  Q <- soilpar$K_s*(s^(2*soilpar$b+3)- H1) # bottom flux
    E <- sapply(s,E_Teuling,vegpar=vegpar)
    qcap <- ifelse(Q < 0,abs(Q),0) # capillary flow component
    L <- Q+E   # total losses
#   if (L < 0) L <- 0 # remove as this means bucket never fills from capillary flow
    return(c(L,E,0,qcap,Q))
}
  # Feedback model needs function Q_E (bottom flux at the root zone) #equation 6 VvdZ2009
  Q_E <- function(s,Z,soilpar,vegpar) {
    hb <- -10^4*soilpar$psi_s_bar
    b1 <- 2+3/soilpar$b
    a1 <- 1+(3/2)/(b1-1)
    H1 <- a1*(hb/(Z-vegpar$Zr))^b1
    Q <- soilpar$K_s*(s^(2*soilpar$b+3)- H1)
    return(Q)
  }

  # Deep root functions              #VvdZ 2009 eq.7
  f <- function(x,c1) exp(-c1*x/100)
  U <- function(z1,z2,Z,c1) integrate(f,z2,z1,c1)$value/integrate(f,0,Inf,c1)$value   #eq.8
  Rc_B <- function(z1,c1,z2) U(z1,z2=z2,Z=z1,c1=c1)*z2/(z1-z2)                 #eq.7

# Feedback DR Model Vervoort and van der Zee (2009)
Fun_FB <- function(s,soilpar,vegpar,Z) {
    
# Calculate the Transpiration components
    vegpar$E_max <-(1-exp(-vegpar$c_T*vegpar$LAI))*vegpar$Ep #equation 2 VvdZ2009
    Rc <- Rc_B(Z,vegpar$c1,vegpar$Zr)
    fr <- U(z1=vegpar$Zr,z2=0,Z=Z,c1=vegpar$c1)
    T_Zr <- U(z1=vegpar$Zr,z2=0,Z=Z,c1=vegpar$c1)*E_Teuling(s,vegpar=vegpar) #=fr*B(s)*Emax, eq3 V.etal2009
    T_DR <- min(Rc*vegpar$E_max,vegpar$E_max-T_Zr) #eq3 Verv.2009
#    T_DR <- min(Rc,1-T_Zr)*vegpar$E_max
 # Qbottom part
    k <- do.call(Q_E,list(s=s,Z=Z,soilpar=soilpar,vegpar=vegpar))    #k is the bottom flux at root zone. Positive means leakage. Negative k means capillary rise
    qcap <- ifelse(k < 0,abs(k),0)
    L <- k+T_Zr # Loss of soil moisture. k is leakage, T_Zr is uptake by plants.
#    if (L < 0) L <- 0 # remove as this means bucket never fills from capillary flow
    return(c(L, T_Zr, T_DR, qcap,k)) # returns the total losses, the T loss from soil (T_Zr) and the T from groundwater (T_DR)
}                        

zeta <- function(s,s_star,s_w,q){
    min(max(((s_star-s)/(s_star-s_w))^q,0),1)}
