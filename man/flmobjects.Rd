\name{flmobjects}

\alias{print.flm}
\alias{summary.flm}

\title{Print and Summary of \code{flm} Objects}

\usage{
\method{print}{flm}(x, \dots)
\method{summary}{flm}(object, \dots)
}

\arguments{
\item{x}{An object of class \code{"flm"}.}
\item{object}{An object of class \code{"flm"}.}
\item{...}{Additional arguments to be passed to other methods.}
}

\value{
The functions return the number of profile factors in the functional linear model, 
the number of runs, the criterion, the objective value of the final design,
the number of iterations of the coordinate exchange algorithm to get to the final
design, and the computational elapsed time in hours:minutes:seconds. 
}

\description{
Print and Summary of objects of class \code{"flm"}.
}

\note{
For examples see \link[fdesigns]{pflm}.
}

\author{
Damianos Michaelides <\email{dm3g15@soton.ac.uk}>, Antony Overstall, Dave Woods
}










