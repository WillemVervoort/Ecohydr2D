\name{E_Teuling_cpp}
\alias{E_Teuling_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
utitity function that calculates actual evaporation based on LAI and soil saturation.
%%  ~~function to do ... ~~
}
\description{
This function take LAI  and soil saturation and calculates the actual ET. This based on Teuling and Troch (2005). It is mostly an internal function but exported for testing reasons
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
E_Teuling_cpp(s, vegpar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
%%     ~~Describe \code{s} here~~
soil saturation value (range 0 - 1)
}
  \item{vegpar}{
%%     ~~Describe \code{vegpar} here~~
a list of vegetation parameters defined by Vegfun() or Vegfun_cpp(). The important ones here are: LAI, Ep, C_T and fr
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a single value, which is the value of the actual evaporation at that saturation level
}
\references{
%% ~put references to the literature/web site here ~
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
\code{\link{?Vegpar_cpp}, \link{?Soilpar_cpp}}
}
\examples{
soilpar <- Soil_cpp("S Clay Loam")
vegpar <- Veg_cpp("Grass",soilpar)
E_Teuling_cpp(0.5,vegpar)

# show for a range of values
s <- seq(0, 1, length = 100)
plot(s, E_Teuling_cpp(s, vegpar), type = "l", col = "red", lwd = 2,
	xlab = "soil saturation", ylab = "actual transpiration")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }% use one of  RShowDoc("KEYWORDS")
\keyword{ function }% __ONLY ONE__ keyword per line
