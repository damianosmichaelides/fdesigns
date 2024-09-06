\name{fdesigns-package}
\alias{fdesigns-package}
\alias{fdesigns}
\docType{package}
\title{
\packageTitle{fdesigns}
}
\description{
\packageDescription{fdesigns}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{fdesigns}
The most important functions are \link[fdesigns]{pflm} and \link[fdesigns]{pfglm} which can be used to identify optimal designs for functional linear and generalised linear models, respectively, using the coordinate exchange algorithm.}
\author{
\packageAuthor{fdesigns}

Maintainer: \packageMaintainer{fdesigns}
}
\references{
Michaelides, D. (2023). Design of experiments for models involving profile factors (Doctoral dissertation, University of Southampton).
}

\examples{
## Example 1:
## This example involves finding an A-optimal design for a functional linear model of 4 runs
## depending on one profile factor. The settings of the profile factor are represented by a 
## B-spline basis of degree zero and a single knot at (0.5). The single functional parameter 
## is represented with a linear power series basis. Five random starts are chosen.

example1 <- pflm(formula = ~ x1, nsd = 5, mc.cores = 1, npf = 1,
  tbounds = c(0, 1), nruns = 4, startd = NULL,  dx = c(0), 
  knotsx = list(c(0.5)), pars = c("power"), db = c(1), 
  knotsb = list(c()), criterion = "A", lambda = 0, 
  dlbound = -1, dubound = 1, tol = 0.0001, progress = FALSE)

print(example1) ##  prints the output of example1.
##
## The number of profile factors is: 1
##
## The number of runs is: 4
##
## The objective criterion is: A-optimality
##
## The objective value is: 8.75
##
## The number of iterations is: 6
##
## The computing elapsed time is: 00:00:00


## Example 2:
## This example involves finding an A-optimal design for a functional logistic
## model of 12 runs depending on one profile factor. The settings of the profile 
## factor are represented by a B-spline basis of degree zero and a three interior knots 
## at (0.25, 0.50, 0.75). The single functional parameter is represented with a linear 
## power series basis. The method of approximation is Monte Carlo with the prior 
## specified by the function prmc. Three random starts are chosen.

set.seed(100) ## Set seed to achieve reproducibility.

prmc <- function(B,Q) {
  matrix(rnorm(B*Q, mean=0, sd=sqrt(2)), nrow=B, ncol=Q)
}
## A function which specifies the prior. This function returns a 
## B by Q matrix of randomly generated values from the prior 
## distribution for the model parameters.

\donttest{
example2 <- pfglm(formula = ~ 1 + x1, nsd = 3, mc.cores = 1, npf = 1, 
             tbounds = c(0,1), nruns = 12, startd = NULL, 
             dx = c(0), knotsx = list(c(0.25,0.50,0.75)), 
             pars = c("power"), db = c(1), knotsb = list(c()), 
             lambda = 0, criterion = "A", family = binomial,
             method=c("MC"), level = 6, B = 10000, prior = prmc,
             dlbound = -1, dubound = 1, tol = 0.0001, progress = TRUE)

print(example2) ##  prints the output of example1.
##
## The number of profile factors is: 1
##
## The number of runs is: 12
##
## The objective criterion is: A-optimality
##
## The objective value is: 20.23283
##
## The number of iterations is: 5
##
## The method of approximation is: MC
##
## The family distribution and the link function are: binomial and logit
##
## The computing elapsed time is: 00:00:12
}
}
