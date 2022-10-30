\name{pflm}

\alias{pflm}

\title{Optimal designs for functional linear models 
using the coordinate exchange algorithm}

\usage{
pflm(formula, nsd = 1, mc.cores = 1, npf, tbounds, 
  nruns, startd = NULL, dx, knotsx, pars, db, 
  knotsb = NULL, lambda = 0, criterion = c("A", "D"), 
  tol = 1e-04, dlbound = -1, dubound = 1, progress = FALSE)
}

\arguments{
\item{formula}{Object of type formula, to create the model equation. Elements need to match the list names for startd. Main effects are called using the names of the profile factors in startd, interactions are called using the names of the profile factors in startd seperated with :, and polynomial effects are called using
the function \link[fdesigns]{P}. Scalar factors are called using the same way and degree and knots through the arguments dx and knotsx are used to specify the scalar factors. A scalar factor is equivalent to a profile factor with degree 0 and no interior knots.}
\item{nsd}{The number of starting designs. The default entry is 1.}
\item{mc.cores}{The number of cores to use. The option is initialized from environment variable MC_CORES if set. Must be at least one, and for parallel computing at least two cores are required. The default entry is 1.}
\item{npf}{The total number of (profile) factors in the model.}
\item{tbounds}{A time vector of length 2, representing the boundaries of time, i.e., 0 and T.}
\item{nruns}{The number of runs of the experiment.}
\item{startd}{Representing the starting design but if NULL then random designs are automatically generated. 
It should be a list of length nsd, and each component should be a list of length npf.}
\item{dx}{A vector of length npf, representing the degree of B-spline basis functions for the functions of the
functional factors. A scalar factor must have a zero degree entry.}
\item{knotsx}{A list of length npf, with every object in the list representing the knot vectors of 
each functional factor. A Scalar factor must have no interior knots, i.e., an empty knot vector.}
\item{pars}{A vector of length equal to the total terms in formula, representing the basis choice for 
the (functional) parameters. Entries should be "power" or "bspline". A scalar parameter is represented through a "power" basis.}
\item{db}{A vector of length equal to the total terms in formula, representing the degree of the basis for
the (functional) parameters. For power series basis the degree is: 1 for linear, 2 for quadratic, etc. A scalar parameter must have degree 0.}
\item{knotsb}{A list of length equal to the total terms in formula, representing the knot vector of each 
(functional) parameter. For parameters represented by a power series basis, the knot vector should be empty or NULL.}
\item{lambda}{Smoothing parameter to penalise the complexity of the functions of the profile factors. The default value is 0, i.e., no penalty.}
\item{criterion}{The choice of objective function. Currently there are two available choices:
1. A-optimality (criterion = "A");
2. D-optimality (criterion = "D").}
\item{tol}{The tolerance value in the optimisation algorithm. Default value is 0.0001.}
\item{dlbound}{The design's lower bound. The default lower bound is -1.}
\item{dubound}{The design's upper bound. The default upper bound is 1.}
\item{progress}{If TRUE, it returns the progress of iterations from the optimisation process. The default entry is FALSE.}
}

\value{
The function returns an object of class "flm" which is a list with the following components:
\item{objval}{The objective value of the final design found from pflm.} 
\item{design}{The final design found from pflm. The final design is a list of length equal to the number
of profile factors, exactly as the starting design startd.}
\item{nits}{The total number of iterations needed to identify the final design.} 
\item{time}{The computational elapsed time in finding the final design.} 
\item{startd}{If starting designs were passed as an argument in pflm, then this is the starting design from the argument startd that led to the final design. If no starting designs were passed to pflm, this is the starting design generated randomly by pflm that led to the final design.}
\item{tbounds}{The argument tbounds.}
\item{npf}{The argument npf.}
\item{criterion}{The argument criterion.}
\item{nruns}{The argument nruns.}
\item{formula}{The argument formula.}
\item{dx}{The argument dx.}
\item{knotsx}{The argument knotsx.}
\item{lambda}{The argument lambda.}
\item{dbounds}{A vector of length 2, containing the arguments dlbound and dubound.}
\item{bestrep}{A scalar value indicating the repetition that led to the final design.}
\item{allobjvals}{A vector of length equal to nsd, representing the objective value from all of the repetitions.}
\item{alldesigns}{A list of length equal to nsd of all the final designs. Each component of the list is a list of length equal to npf representing the final design in each repetition of the coordinate exchange algorithm.}
\item{allstartd}{If starting designs were passed as an argument in pflm, then this is the argument. If no starting designs were passed to pflm, this is the starting designs generated randomly by pflm.}
}

\description{
Optimal designs for functional linear models for which the functional factors are represented
 as B-spline basis functions and the functional parameters are represented as power series basis functions
 or as B-spline basis functions.
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

plot(example1) ##  plots the functions of the profile factors.

## Example 2:
## This example involves finding a D-optimal design for a functional linear model of n=20 runs
## depending on two profile factors. In addition to the main effects, the model includes the
## interaction of the profile factors and the quadratic effect of the second profile factor. 
## The settings of the profile factors are represented by B-spline basis of quadratic degrees
## and knots at (0.33, 0.66) and (0.25, 0.50, 0.75). The functional parameters are represented 
## with linear power basis and quadratic B-spline basis with knots at (0.25, 0.50, 0.75). 
## The complexity of the designs is penalised with the smoothing value equal to 1.

example2 <- pflm(formula = ~ x1 + x2 + x1:x2 + P(x2, 2), nsd = 1, mc.cores = 1,
  npf = 2, tbounds = c(0, 1), nruns = 20, startd = NULL, dx = c(2, 2),
  knotsx = list(c(0.33, 0.66), c(0.25, 0.50, 0.75)), 
  pars = c("power", "power", "bspline", "bspline"), db = c(1, 1, 2, 2), 
  knotsb = list(c(), c(), c(0.25, 0.50, 0.75), c(0.25, 0.50, 0.75)), 
  criterion = "D", lambda = 1, tol = 0.0001, dlbound = -1, dubound = 1, 
  progress = FALSE)

print(example2) ##  prints the output of example2.
##
## The number of profile factors is: 2
##
## The number of runs is: 20
##
## The objective criterion is: D-optimality
##
## The objective value is: 0.05706758
##
## The number of iterations is: 6
##
## The computing elapsed time is: 00:00:17

plot(example2) ## plots the functions of profile factors.
}










