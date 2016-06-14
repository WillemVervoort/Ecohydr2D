# New vegetation input file with all vegs
# Vegmodelling Willem Vervoort April 2009
# ##########################################
Veg <- function(vtype,soilpar) {
   E_w<-0.01 # cm/day  (Evaporation rate at wilting point)
   k<-0.5  # index of plant resistance to water stress (page 739 second column, Porporato et al., 2001)
                      
  if (vtype == "Grass") {
     #................................................
     # Vegetation 1 (Grass)
     # paspalum secateum F-I and R-I, 2004
     # Gives only parameters
     # Willem Vervoort 1/07/05
     #.................................................

     Zr<-50 # soil depth (cm)        !! = 40cm in F-I&R-I 2004
     delta <- 0.1  # cm interception depth (Fernandez-Illascez & Rodriguez-Iturbe, Table 1)

     # calculate s_w and s_star
   
      psi_ss<--0.09          # soil matric potential at stomatal closure
     s_star<-(psi_ss/soilpar$psi_s_bar)^(-1/soilpar$b)
     psi_sw<--3 # MPa

      #psi_sw<--1.5
      s_w<-(psi_sw/soilpar$psi_s_bar)^(-1/soilpar$b)
      E_max<-0.33 # cm/day

     #
     # Colonisation parameters
     # Table 4, F-I and R-I, 2004
     R1<-0.6 # named this R1, since rainfall is already called R
     S<-7.9
     m1<-0.1  # mortality rate
     colpar <- list(a=R1,b=S,m=m1)

     LAI <- 1.5    # Roberts et al. 2000 assuming Kc = 0.8
     c_T <- 0.45 # Tarrawarra dataset parameter in Teuling and Troch 2005
     fr <- 1
     Ep <- 0.65 # cm/day BOM data average Moree
     
     #new groundwater uptake function
     # defines the fraction of roots close to the groundwater
     c1 <- NULL
	 fs <- NULL
     
     # vegmodelling param
     AboveGround_PF <- 0.6 # above ground partioning factor how much carbon goes to roots
     Leaf_PF <- 1.0  # leaf partioning (how much of above ground goes to leaves)
     SLA <- 10  # Specific leaf area m^2/kg?
     Mort <- 0.015  # 68 days Tjoelker et al. 2005
     fc.coef <- 0.9 #Lu et al. 2001
     veg.psi <- 0.03  # Haxeltine et al
     dorm<-TRUE
     dorm.cutoff <- 0.4
     DR <- FALSE     
     
  }

  if (vtype == "TreesDR") {
       #......................
     # Vegetation 2 (trees)
     # Eucalyptus tree generic
     # Gives only parameters
     # Willem Vervoort 1/07/05
     #............................

     Zr<-100 # root depth
     delta <- 0.2  # cm interception rate (Laio et al., 2001b page 748 for trees)

     psi_ss<--0.12  # MPa Table 1, F-I & R-I 2004
     s_star<-(psi_ss/soilpar$psi_s_bar)^(-1/soilpar$b)
     psi_sw<--7 # MPa (Eucalyptus Camaldulensis Whitehead and Beadle 2004)
     s_w<-(psi_sw/soilpar$psi_s_bar)^(-1/soilpar$b)

     E_max<-0.5 # 0.2 cm/day  Cooke et al 1998 0.3 cm/day Beadle and Whitehead 2004

     #
     # Colonisation parameters
     # table 4 F-I and R-I, 2004
     T1<-0.5
     U<-7.6
     m2<-0.1  # mortality rate
     colpar <- list(a=T1,b=U,m=m2)
     
     # Parameters for root uptake Teuling and Troch
     LAI <- 0.9 # 3.2 maximum Peel et al. 2005 Eastern Australia 0.9 Huntley et al 2001 (table 1) 2.5 Whitehead and Beadle Table 4, native forest
     c_T <- 0.45 # Mid value from WAVES
     fr <- 1.0 # Stick to 1 with a box model
     Ep <- 0.65 # cm/day BOM data average Goondiwindi
#     Ep <- 0.9 # cm/day BOM data summer Moree

     #new groundwater uptake function
     # fraction of roots close to the groundwater
     c1 <- 1.0
     fs <- 0.25
   
     # vegmodelling param   
     AboveGround_PF <- 0.35  # above ground partioning factor how much carbon goes to roots
     Leaf_PF <- 0.5  # leaf partioning (how much of above ground goes to leaves)
     SLA <- 3.5  # Specific leaf area 3.6 - 10 Whitehead and Beadle 3.5 Schulze et al. 2006
     Mort <- 0.0015 # leaves live 2 years
     fc.coef <- 0.9 #Lu et al. 2001
     veg.psi <- 0.038  # Haxeltine et al   # probably can leave that constant
     dorm<-FALSE
     dorm.cutoff <- 0.1    # minimum leaf area index
     DR<- TRUE
  }
  
    if (vtype == "TreesNoDR") {
       #......................
     # Vegetation 2 (trees)
     # Eucalyptus tree generic
     # Gives only parameters
     # Willem Vervoort 1/07/05
     #............................

     Zr<-100 # root depth
     delta <- 0.2  # cm interception rate (Laio et al., 2001b page 748 for trees)

     psi_ss<--0.12  # MPa Table 1, F-I & R-I 2004
     s_star<-(psi_ss/soilpar$psi_s_bar)^(-1/soilpar$b)
     psi_sw<--5 # MPa (Eucalyptus Camaldulensis Whitehead and Beadle 2004)
     s_w<-(psi_sw/soilpar$psi_s_bar)^(-1/soilpar$b)

     E_max<-0.5 # 0.2 cm/day  Cooke et al 1998 0.3 cm/day Beadle and Whitehead 2004

     #
     # Colonisation parameters
     # table 4 F-I and R-I, 2004
     T1<-0.5
     U<-7.6
     m2<-0.1  # mortality rate
     colpar <- list(a=T1,b=U,m=m2)
     
     # Parameters for root uptake Teuling and Troch
     LAI <- 0.9 # 3.2 maximum Peel et al. 2005 Eastern Australia 0.9 Huntley et al 2001 (table 1) 2.5 Whitehead and Beadle Table 4, native forest
     c_T <- 0.45 # Mid value from WAVES
     fr <- 1.0 # Stick to 1 with a box model
     Ep <- 0.65 # cm/day BOM data average Goondiwindi
#     Ep <- 0.9 # cm/day BOM data summer Moree

     #new groundwater uptake function
     # fraction of roots close to the groundwater
     c1 <- NULL
     fs <- NULL
     # vegmodelling param   
     AboveGround_PF <- 0.35  # above ground partioning factor how much carbon goes to roots
     Leaf_PF <- 0.5  # leaf partioning (how much of above ground goes to leaves)
     SLA <- 3.5  # Specific leaf area 3.6 - 10 Whitehead and Beadle 3.5 Schulze et al. 2006
     Mort <- 0.0015 # leaves live 2 years
     fc.coef <- 0.9 #Lu et al. 2001
     veg.psi <- 0.038  # Haxeltine et al   # probably can leave that constant
     dorm<-FALSE
     dorm.cutoff <- 0.1    # minimum leaf area index
     DR <- FALSE     
  }

  if (vtype == "Lignum") {
     #................................................
     # Vegetation 3 (Lignum)
     # Suggested values F. van Ogtrop & W. Vervoort
     # Gives only parameters
     # Willem Vervoort 20/05/2009
     #.................................................

     Zr<-100 # soil depth (cm)
     delta <- 0.2  # cm interception depth (Fernandez-Illascez & Rodriguez-Iturbe, Table 1)

     # calculate s_w and s_star

      psi_ss<--0.12
     s_star<-(psi_ss/soilpar$psi_s_bar)^(-1/soilpar$b)
     psi_sw<--5 # MPa

      #psi_sw<--1.5
      s_w<-(psi_sw/soilpar$psi_s_bar)^(-1/soilpar$b)
      E_max<-0.33 # cm/day

     #
     # Colonisation parameters
     # not important to be changed later, place holder
     # Table 4, F-I and R-I, 2004
     R1<-0.6 # named this R1, since rainfall is already called R
     S<-7.9
     m1<-0.1  # mortality rate
     colpar <- list(a=R1,b=S,m=m1)

     # This is still important
     LAI <- 2.03 # starting LAI (Roberts et al. 2000)
     c_T <- 0.45 # Tarrawarra dataset parameter in Teuling and Troch 2005
     fr <- 1    # no deep roots
     Ep <- 0.65 # cm/day BOM data average Moree

     #new groundwater uptake function
     # defines the fraction of roots close to the groundwater
     c1 <- NULL
	fs <- NULL

     # vegmodelling param
     AboveGround_PF <- 0.65 # quite high  Capon et al
     Leaf_PF <- 0.3 # very low Capon et al
     SLA <- 25 # very high Capon et al. in press: 25 is that correct?
     Mort <- 0.015 # Leaves less than 1 year
     fc.coef <- 0.9 #Lu et al. 2001
     veg.psi <- 0.03  # Haxeltine et al
     dorm <- TRUE
     dorm.cutoff <- 0.5
     DR <- FALSE     
  }
  
if (vtype == "Bare") {
     #................................................
     # Bare Soil
     # added for 2D project
     # Willem Vervoort 12/06/13
     #.................................................

     Zr<-25 # soil depth (cm)
     delta <- 0  # cm interception depth 

     # calculate s_w and s_star
   
      psi_ss<--0.01          # soil matric potential at stomatal closure, very small
     s_star<-(psi_ss/soilpar$psi_s_bar)^(-1/soilpar$b)
     psi_sw<--1.5 # MPa

      #psi_sw<--1.5
      s_w<-(psi_sw/soilpar$psi_s_bar)^(-1/soilpar$b)
      E_max<-0.33 # cm/day

     #  This is irrelevant for bare soil
     # Colonisation parameters
     # Table 4, F-I and R-I, 2004
     R1<-NULL # named this R1, since rainfall is already called R
     S<-NULL
     m1<-NULL  # mortality rate
     colpar <- list(a=R1,b=S,m=m1)

     LAI <- 1    # using LAI 1 for bare soil
     c_T <- 0.45 # Tarrawarra dataset parameter in Teuling and Troch 2005
     fr <- 1
     Ep <- 0.65 # cm/day BOM data average Moree
     
     #new groundwater uptake function
     # defines the fraction of roots close to the groundwater
     c1 <- NULL
	fs <- NULL
     
     # vegmodelling param not needed
     AboveGround_PF <- NULL # above ground partioning factor how much carbon goes to roots
     Leaf_PF <- NULL  # leaf partioning (how much of above ground goes to leaves)
     SLA <- NULL  # Specific leaf area m^2/kg?
     Mort <- NULL  # 68 days Tjoelker et al. 2005
     fc.coef <- NULL #Lu et al. 2001
     veg.psi <- NULL  # Haxeltine et al
     dorm<-FALSE
     dorm.cutoff <- 0.4
     DR <- FALSE     
     
  }


return(list(E_max=E_max,E_w=E_w,k=k,Zr=Zr,delta=delta,s_w=s_w,s_star=s_star,
          LAI=LAI,c_T=c_T,fr=fr,Ep=Ep,c1=c1,fs=fs,AboveGround_PF=AboveGround_PF,
          Leaf_PF=Leaf_PF,SLA=SLA,Mort=Mort,fc.coef=fc.coef,veg.psi=veg.psi,dorm=dorm,
          dorm.cutoff=dorm.cutoff,colpar=colpar,DR=DR))
}
