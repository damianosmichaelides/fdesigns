\name{fglmobjects}

\alias{print.fglm}
\alias{summary.fglm}

\title{Print and Summary of \code{fglm} Objects}

\usage{
\method{print}{fglm}(x, \dots)
\method{summary}{fglm}(object, \dots)
}

\arguments{
\item{x}{An object of class \code{"fglm"}.}
\item{object}{An object of class \code{"fglm"}.}
\item{...}{Additional arguments to be passed to other methods.}
}

\value{
The functions return the number of profile factors in the functional generalised linear 
model, the number of runs, the criterion, the objective value of the final design,
the number of iterations of the coordinate exchange algorithm to get to the final
design, the method of approximation of the expectation of the objective function, the family 
distribution and link function, and the computational elapsed time in hours:minutes:seconds. 
}

\description{
Print and Summary of objects of class \code{"fglm"}.
}

\note{
For examples see \link[fdesigns]{pfglm}.
}

\author{
Damianos Michaelides <\email{dm3g15@soton.ac.uk}>, Antony Overstall, Dave Woods
}








