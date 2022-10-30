\name{P}

\alias{P}

\title{Compute Profile Factor Polynomials}

\usage{
P(x, deg)
}

\arguments{
\item{x}{A coefficient matrix from the basis expansion of a profile factor. The name passed needs
to match the name of the profile factor in startd in \link[fdesigns]{pflm} and \link[fdesigns]{pfglm}.}
\item{deg}{The degree of the polynomial effect for the profile factor.}
}

\value{
The function returns an attributes list. The list contains the polynomial coefficient matrix,
the argument x, and the argument deg.
}

\description{
A function to be used in the formula argument in \link[fdesigns]{pflm} and \link[fdesigns]{pfglm}. 
It computes polynomials of the basis functions of profile factors to be passed 
in functional models and find optimal designs. 
}

\details{
The function \link[fdesigns]{P} is intented to be used in the formula argument in the 
function \link[fdesigns]{pflm} and \link[fdesigns]{pfglm}. 

In the traditional linear models, polynomial effects of factors are handled using the 
function \link[base]{I}. However, profile factors are expanded using basis functions
and the coefficients are a matrix instead of a vector. 

In other words, the function \link[fdesigns]{P} is an extension to the function \link[fdesigns]{I},
but for functional models with profile factors. 
}

\examples{
## Example 1:

d <- list(x1=matrix(runif(48), nrow=12))
attach(d)
example1 <- P(x = x1, deg = 2)
detach(d)

## Example 2: Use the function in a formula
## This example involves finding a D-optimal design for a functional linear model of n=20 runs
## depending on two profile factors. In addition to the main effects, the model includes the
## interaction of the profile factors and the quadratic effect of the second profile factor. 
## The quadratic effect is identified in the formula argument using the \link[fdesigns]{P} function. 

example2 <- pflm(formula = ~ x1 + x2 + x1:x2 + P(x2, 2), nsd = 1, mc.cores = 1,
  npf = 2, tbounds = c(0, 1), nruns = 20, startd=NULL, dx = c(2, 2),
  knotsx = list(c(0.33, 0.66), c(0.25, 0.50, 0.75)), 
  pars = c("power", "power", "bspline", "bspline"), db = c(1, 1, 2, 2), 
  knotsb = list(c(), c(), c(0.25, 0.50, 0.75), c(0.25, 0.50, 0.75)), 
  criterion = "D", lambda = 1, tol = 0.0001, dlbound = -1, dubound = 1, 
  progress = FALSE)
}










