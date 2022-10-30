\name{summary.flm}

\alias{summary.flm}

\title{Print and Summary of flm Objects}

\usage{
## S3 method for class "flm"
print(obj, ...)
## S3 method for class "flm"
summary(obj, ...)
}

\arguments{
\item{obj}{An object of class "flm".}
\item{...}{Additional arguments to be passed to other methods.}
}

\value{
The functions return the number of profile factors in the functional linear model, 
the number of runs, the criterion, the objective value of the final design,
the number of iterations of the coordinate exchange algorithm to get to the final
design, and the computational elapsed time in hours:minutes:seconds. 
}

\description{
Print and Summary of the flm Object.
}

\examples{
For examples see \link[fdesigns]{pflm}.
}










