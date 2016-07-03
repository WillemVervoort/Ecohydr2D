\name{Veg_cpp}
\alias{Veg_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
This function generates a list of vegetation parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Based on a choice of vegetations, this function generates a list of the necessary vegetation parameters for ecohydrological modelling using the extensions on the Rodriguez-Iturbe model.
}
\usage{
Veg_cpp(vtype, soilpar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vtype}{
%%     ~~Describe \code{vtype} here~~
character string describing the vegetation type. Currently implemented are: "Grass", "TreesDR" (with deep roots), "Bare", "Lignum", "TreesNoDR" (without deep roots).
}
  \item{soilpar}{
%%     ~~Describe \code{soilpar} here~~
A list of the necessary soil parameters generated using Soil_cpp()
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function only generates the default list for each vegetation, individual parameters in the list can be redefined directly by assigning new values to elements of the list.
}
\value{
%%  ~Describe the value returned
A set of parameters that are the same across all vegetation
\item{E_w }{Evaporation at wilting point, minimum evaporation}  \item{k }{a parameter describing ....}
And a set of vegetation specific parameters
\item{Zr }{The root depth of the vegetation in cm}
\item{delta }{The rainfall interception depth cm}
\item{s_star }{The soil saturation value of incipient stomatal closure (0 - 1)}
\item{s_w }{The soil saturation at wilting point (0 - 1)}
\item{Emax }{The maximum evaporation/transpiration rate in cm, for the Rodriguez-Iturbe mode}
\item{LAI }{The leaf area index (-)}
\item{c_T }{The value of radiation conversion in the Teuling & Troch model}
\item{fr }{deprecated, The fraction of roots, generally set to 1, Vervoort and van der Zee, 2009}
\item{Ep }{Potential Evaporation rate in cm}
\item{c1 }{parameter in the new root distribution function (0.5 - 3) Vervoort and van der Zee, 2012}
\item{fs }{parameters in the new root distribution function Vervoort and van der Zee, 2012}
\item{DR }{a logical value whether or not to model deep roots}
\item{q }{vegetation stress value to calculate static stress}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Rodriguez-Iturbe et al 
Rodriguez-Iturbe and Porporato 
Rodriguez-Iturbe and Fernandez-Illascez 
Vervoort and van der Zee, 2008 
Vervoort and van der Zee, 2009  
Vervoort and van der Zee, 2012 
}
\author{
%%  ~~who you are~~
R. Willem Vervoort willemvervoort@gmail.com
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{Soil_cpp},\link{E_Teuling_cpp}} 
}
\examples{
soilpar <- Soil_cpp("S Clay Loam")
vegpar <- Veg_cpp("TreesDR",soilpar)
str(vegpar)
names(vegpar)
# assign a value
vegpar$Zr <- 50
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
