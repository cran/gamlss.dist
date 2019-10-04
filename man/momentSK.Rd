\name{momentSK}
\alias{momentSK}
\alias{centileSK}
\alias{centileSkew}
\alias{centileKurt}
\alias{theoCentileSK}
\alias{plotCentileSK}
\alias{boundary}
\alias{fEGB2_1}
\alias{fEGB2_2}
\alias{fJSU}
\alias{fSEP3}
\alias{fST3_1}
\alias{fST3_2}
\alias{SKmoment_col}
\alias{SKmoment_gray}
\alias{SKmomentBoth}
\alias{checkMomentSK}
\alias{cEGB2_1}
\alias{tEGB2_1}
\alias{cEGB2_2}
\alias{tEGB2_2}
\alias{cEGB2_1Data}
\alias{cJSU}
\alias{tJSU}
\alias{cSB}
\alias{tSB}
\alias{cSEP3}
\alias{tSEP3}
\alias{cSHASH}
\alias{tSHASH}
\alias{cST3_1}
\alias{tST3_1}
\alias{cST3_2}
\alias{tST3_2}
\alias{SKcentile_col}
\alias{SKcentile_gray}
\alias{SKcentileBoth}
\alias{checkCentileSK}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample and theoretical Moment and Centile Skewness and Kurtosis Functions}
\description{
The  functions \code{momentSK()}, \code{centileSK()},  \code{centileSkew()} and \code{centileKurt()},  calculate sample statistics related to skewness and kurtosis.  The function \code{theoCentileSK()} calculates the theoretical centile statistics from a given \code{gamlss.family} distribution. The \code{plotCentileSK()} plots the theoretical centile skewness and kurtosis against \code{p} (see below). 

The function \code{checkMomentSK()} can be use to check (a) whether the moment skewness and kurtosis of a fitted model are modelled adequantly (the residuals of the model are used). (b)  whether a given sample display skewness or kurtosis. 
}
\usage{
momentSK(x)
centileSK(x, cent = c(1, 25))
centileSkew(x, cent = 1)
centileKurt(x, cent = 1)

theoCentileSK(fam = "NO", p = 0.01, ...)
plotCentileSK(fam = "NO", plotting = c("skew", "kurt", "standKurt"), 
             add = FALSE, col = 1, lty = 1, lwd = 1, ylim = NULL, ...)
             
checkMomentSK(x, add = FALSE, bootstrap = TRUE, no.bootstrap = 99, 
               col.bootstrap = "lightblue", pch.bootstrap = 21, 
               asCharacter = TRUE, col.point = "black", pch.point = 4, 
               lwd.point = 2, cex.text = 1.5, col.text = "black", 
               show.legend = TRUE) 
               
checkCentileSK(x, type = c("central", "tail"), add = FALSE, 
              bootstrap = TRUE, no.bootstrap = 99, 
              col.bootstrap = "lightblue", pch.bootstrap = 21, 
              asCharacter = TRUE, col.point = "black", pch.point = 4, 
              lwd.point = 2, cex.text = 1.5, col.text = "black", 
              show.legend = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data vector or gamlss model}
  \item{cent}{the centile required}
  \item{type}{For centile skewness and kurtosis only whether "central" (default) or "tail")}
  \item{fam}{A gamlss distribution family}
  \item{plotting}{what to plot}
  \item{add}{whether to add the line to the existing plot}
  \item{col}{the colour of the line}
  \item{lty}{the type of the line}
  \item{lwd}{the width of the line}
  \item{ylim}{the y limit of the graph}
  \item{p}{the value determiming the centile skewness or kurtosis}
  \item{\dots}{additional  arguments pass to \code{theoCentileSK()} function i.e. the values of the distribution parameters}
  \item{bootstrap}{whether a plot of the bootstrap skewness and kurtosis measures should be added in the plot}
  \item{no.bootstrap}{the number of boostrap skewness and kurtosis measures}
  \item{col.bootstrap}{the coloue for boostraps}
  \item{pch.bootstrap}{the point type of  boostraps}
  \item{asCharacter}{whether to plot the estimated  skewness and kurtosis measure as character or as point}
  \item{col.point}{the colour of the skewness and kurtosis measure}
  \item{pch.point}{the  point type  of the skewness and kurtosis measure}
  \item{lwd.point}{the width of the plotted line}
  \item{cex.text}{the size of the text}
  \item{col.text}{the colour of the text}
  \item{show.legend}{whether to show the legent}
}
\details{
Those function calculate sample moment and centile skewness and kurtosis statistics and theoretical  centile values for a specific distribution.
}
\value{Different functions produce different output:
The function \code{momentSK()} produce:
\item{mom.skew:}{sample moment skewness}
\item{trans.mom.skew:}{sample transformed moment skewness}
\item{mom.kurt:}{sample moment kurtosis}
\item{excess.mom.kurt:}{sample excess moment kurtosis}
\item{trans.mom.kurt:}{sample ransformed moment excess kurtosis}
\item{jarque.bera.test:}{the value of the Jarque-bera test for testing whether skewness and excess kurtosis are zero or not}

 The function \code{centileSK()} produces:
 \item{S0.25:}{sample centile central skewness}
 \item{S0.01:}{sample centile tail skewness}
 \item{trans.S0.25:}{sample centile transformed central skewness}
 \item{trans.S0.01:}{sample centile transformed tail skewness}
 \item{K0.01:}{sample centile kurtosis}
 \item{standK0.01:}{standardise centile kurtosis, (\code{K0.01/3.449})}
 \item{exc.K0.01:}{excess centile kurtosis, (\code{K0.01-3.449})}
 \item{trans.K0.01:}{transfored excess centile kurtosis, (exc.K0.01/(1+abs(exc.K0.01))}

The function \code{centileSkew()} for a given argument \code{p} produces:
 \item{p:}{the value determiming the centile skewness}
 \item{Sp:}{sample centile skewness at \code{p}}
 \item{tSp:}{sample transformed centile skewness at \code{p}}
 

The function \code{centileKurt()} for a given argument \code{p} produces:
 \item{p}{the value determiming the centile kurtosis}
 \item{Kp}{sample  centile kurtosis at \code{p}}
 \item{sKp}{sample  standardise centile kurtosis at \code{p}}
 \item{ex.Kp:}{sample excess centile kurtosis at \code{p}}
 \item{teKp:}{sample transformed excess centile kurtosis at \code{p}}


The function \code{theoCentileSK} for a given \code{gamlss.family}  produces:
 \item{IR}{the interquartile range of the distribution}
 \item{SIR}{the semi interquartile range of the distribution}
 \item{S_0.25}{the central skewness of the distribution}
 \item{S_0.01:}{the tail skewness of the distribution}
 \item{K_0.01:}{the centile kurtosis of the distribution}
 \item{sK_0.01:}{the standardise centile kurtosis of the distribution}

}
\references{
Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC. An older version can be found in \url{http://www.gamlss.com/}.
}
\author{
Mikis Stasinopoulos, Bobert Rigby, Gillain Heller and Fernanda De Bastiani. 
}



\seealso{\code{\link{gamlss.family}}}
\examples{
Y <- rSEP3(1000)
momentSK(Y)
centileSK(Y)
centileSkew(Y, cent=20)
centileKurt(Y, cent=30)

theoCentileSK("BCCG", mu=2, sigma=.2, nu=2)
plotCentileSK(fam="BCCG",  mu=2, sigma=.2, nu=2)

checkMomentSK(Y)
checkCentileSK(Y)
#checkCentileSK(Y, type="tail")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}% use one of  RShowDoc("KEYWORDS")
