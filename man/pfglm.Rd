\name{pfglm}

\alias{pfglm}

\title{Optimal designs for functional generalised linear 
models using the coordinate exchange algorithm}

\usage{
pfglm(formula, nsd = 1, mc.cores = 1, npf, tbounds, nruns, startd = NULL, 
        dx, knotsx, pars, db, knotsb = NULL, lambda = 0, criterion = c("A","D"),
        family, method = c("quadrature", "MC"), level = NULL, B = NULL, prior, 
        dlbound = -1, dubound = 1, tol = 0.0001, progress = FALSE)
}

\arguments{
\item{formula}{Object of type formula, to create the model equation. Elements need to match the list names for startd. Main effects are called using the names of the profile factors in startd, interactions are called using the names of the profile factors in startd seperated with :, and polynomial effects are called using
the function \link[fdesigns]{P}. Scalar factors are called using the same way and degree and knots through the arguments dx and knotsx are used to specify the scalar factors. A scalar factor is equivalent to a profile factor with degree 0 and no interior knots.}
\item{nsd}{The number of starting designs. The default entry is 1.}
\item{mc.cores}{The number of cores to use. The option is initialized from environment variable MC_CORES if set. Must be at least one, and for parallel computing at least two cores are required. The default entry is 1.}
\item{npf}{The total number of (functional) factors in the model.}
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
\item{family}{Specifies the error distribution and the link function of the functional generalised linear model. It can be the name of a family in the form of a character string or an R family function; see \link[stats]{family} for details. Currently, the methodology is implemented only for the binomial family with the logit link, i.e., family = binomial(link = "logit"), and the Poisson family with the log link, i.e., family = poisson(link = "log").}
\item{method}{A character argument specifying the method of approximation of the expectation of the objective function with respect to a prior distribution of the parameters. Currently there are two available choices:
1. Deterministic quadrature approximation (method = "quadrature");
2. Stochastic Monte Carlo approximation (method = "MC").}
\item{level}{An optional argument that specifies the accuracy level in the quadrature approximation. It is the number of points in each dimension. If NULL and method = "quadrature", then it defaults to 5. A high value of level may increase the computation time; especially for complicated models. If the model is complicated, i.e., several profile factors or interactions and polynomials, prefer to use method = "MC".}
\item{B}{An optional argument that specifies the size of the Monte Carlo samples. If NULL and method = "MC", then it defaults to 10000. For method = "quadrature", B is computer automatically according to the dimensionality of the functional model and the level argument.}
\item{prior}{An argument to specify the prior distribution.
For method = "MC", it should be a function of two arguments B and Q. Both arguments are integers. The value of B corresponds to the argument B, and the value of Q represents the total number of basis functions of the functional parameters. The function must generate a matrix of dimensions B by Q, that contains a random sample from the prior distribution of the parameters.
For method = "quadrature", normal and uniform prior distribution for the parameters are allowed. For a normal prior distribution, the argument prior needs to be a list of length 2, with the entries named "mu" for the prior mean and "sigma2" for the prior variance-covariance matrix. The prior mean can be a scalar value that means all parameters have the same prior mean, or a vector of prior means with length equal to the number of parameters in the functional model. The prior variance-covariance can be a scalar value that means all parameters have a common variance, or a vector of prior variances with length equal to the number of parameters in the functional model, or a square matrix with the number of rows and columns equal to he number of parameters in the functional model. For a uniform prior distribution, the argument prior needs to be a list of a single entry named "unifbound" for the lower and upper bounds of the prior distribution. The bounds can be a vector of length 2 that means all parameters have the same bounds, or a matrix with the number of rows equal to 2 and the number of columns equal to the number of parameters in the functional model.}
\item{dlbound}{The design's lower bound. The default lower bound is -1.}
\item{dubound}{The design's upper bound. The default upper bound is 1.}
\item{tol}{The tolerance value in the optimisation algorithm. Default value is 0.0001.}
\item{progress}{If TRUE, it returns the progress of iterations from the optimisation process. The default entry is FALSE.}
}

\value{
The function returns an object of class "flm" which is a list with the following components:
\item{objective.value}{The objective value of the final design found from pfglm.} 
\item{design}{The final design found from pfglm. The final design is a list of length equal to the number
of profile factors, exactly as the starting design startd.}
\item{n.iterations}{The total number of iterations needed to identify the final design.} 
\item{time}{The computational elapsed time in finding the final design.} 
\item{startd}{If starting designs were passed as an argument in pfglm, then this is the argument startd. 
        If no starting designs were passed to pfglm, this is the starting design generated randomly by pfglm.}
\item{tbounds}{The argument tbounds.}
\item{npf}{The argument npf.}
\item{criterion}{The argument criterion.}
\item{nruns}{The argument nruns.}
\item{formula}{The argument formula.}
\item{family}{A vector of length equal to 2, containing the family and the link function.}
\item{method}{The argument method.}
\item{B}{The argument B.}
\item{prior}{The argument prior.}
\item{dx}{The argument dx.}
\item{knotsx}{The argument knotsx.}
\item{lambda}{The argument lambda.}
\item{dbounds}{A vector of length 2, containing the arguments dlbound and dubound.}
\item{bestrep}{A scalar value indicating the repetition that led to the final design.}
\item{allobjvals}{A vector of length equal to nsd, representing the objective value from all of the repetitions.}
\item{alldesigns}{A list of length equal to nsd of all the final designs. Each component of the list is a list of length equal to npf representing the final design in each repetition of the coordinate exchange algorithm.}
\item{allstartd}{If starting designs were passed as an argument in pfglm, then this is the argument. If no starting designs were passed to pfglm, this is the starting designs generated randomly by pfglm.}
}

\description{
Optimal designs for functional generalised linear models for which the functional factors are represented
 as B-spline basis functions and the functional parameters are represented as power series basis functions
 or as B-spline basis functions.
}

\examples{
## Example 1:
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

example1 <- pfglm(formula = ~ 1 + x1, nsd = 3, mc.cores = 1, npf = 1, 
             tbounds = c(0,1), nruns = 12, startd = NULL, 
             dx = c(0), knotsx = list(c(0.25,0.50,0.75)), 
             pars = c("power"), db = c(1), knotsb = list(c()), 
             lambda = 0, criterion = "A", family = binomial,
             method=c("MC"), level = 6, B = 10000, prior = prmc,
             dlbound = -1, dubound = 1, tol = 0.0001, progress = TRUE)

print(example1) ##  prints the output of example1.
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

plot(example1) ##  plots the functions of the profile factors.

## Example 2:
## This example involves finding a A-optimal design for a functional logistic
## model of 12 runs depending on one profile factor. The settings of the profile 
## factor are represented by a B-spline basis of degree zero and a three interior knots 
## at (0.25, 0.50, 0.75). The single functional parameter is represented with a linear 
## power series basis. The method of approximation is Quadrature with Normal prior 
## distribution, with all parameters having mean 0 and variance 2. 
## Five random starts are chosen.

## The number of profile factors is: 1
## 
## The number of runs is: 12
## 
## The objective criterion is: A-optimality
## 
## The objective value is: 1.8554e+24
## 
## The number of iterations is: 28
## 
## The method of approximation is: quadrature
## 
## The family distribution and the link function are: binomial and logit
## 
## The computing elapsed time is: 00:00:05

plot(example2) ##  plots the functions of the profile factors.


}










