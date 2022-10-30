\name{summary.fglm}

\alias{summary.fglm}

\title{Print and Summary of fglm Objects}

\usage{
## S3 method for class "fglm"
print(obj, ...)
## S3 method for class "fglm"
summary(obj, ...)
}

\arguments{
\item{obj}{An object of class "fglm".}
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
Print and Summary of the fglm Object.
}

\examples{
For examples see \link[fdesigns]{pfglm}.
}










