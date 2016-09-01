\name{Soil_cpp}
\alias{Soil_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
This function generates a list of soil parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Based on a choice of generic soil types, this function generates a list of the necessary soil parameters for ecohydrological modelling using the extensions on the Rodriguez-Iturbe model and Vervoort and van der Zee ecohydrological models.
}
\usage{
Soil_cpp(stype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stype}{
%%     ~~Describe \code{stype} here~~
character string describing the generic soil type. Currently implemented are: "Medium Light Clay", "S Clay Loam", "M Clay", "H Clay", "C Sand", "Loamy Sand".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function only generates the default list for each vegetation, individual parameters in the list can be redefined directly by assigning new values to elements of the list.
}
\value{
%%  ~Describe the value returned
A set of parameters specific for all soil types
\item{n }{porosity (-)}  
\item{K_s }{Hydraulic conductivity in cm/day}
\item{b }{Parameter "b" in Campbell's water retention function}
\item{hb }{bubbling pressure (cm)}
\item{psi_s_bar }{bubbling pressure (kPa)}
\item{s_fc }{soil saturation at field capacity (0 - 1)}  
\item{s_h }{soil saturation at hygroscopic point (0 - 1)}
\item{beta }{parameter in soil retention function Laio et al (2001)}
\item{beta1 }{beta parameter in Eagleson function (Eagleson (1978))}
\item{a1 }{parameter in Eagleson function (Eagleson (1978))}
\item{spec_y}{specific yield calculated as n - s_fc}
}
\references{
%% ~put references to the literature/web site here ~
Vervoort and van der Zee (2008)
Vervoort and van der Zee (2009)
Vervoort and van der Zee (2012)
Laio et al. (2001)
Eagleson (1978)
}
\author{
%%  ~~who you are~~
R.Willem Vervoort willemvervoort@gmail.com
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{Veg_cpp},\link{E_Teuling_cpp}} 

}
\examples{
soilpar <- Soil_cpp("S Clay Loam")
str(soilpar)
names(soilpar)
# assign a value
vegpar$K_s <- 0.5
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
