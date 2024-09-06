\name{plot.flm}

\alias{plot.flm}

\title{Plot of \code{flm} Objects}

\usage{
\method{plot}{flm}(x, \dots)
}

\arguments{
\item{x}{An object of class \code{"flm"}.}
\item{...}{Additional arguments to be passed to other methods.}
}

\value{
The function returns the question:
"Which profile factor to plot?".
The answer needs to be an integer to specify the profile factor for which 
to plot the optimal functions in every run. The value needs to be between 
1 and the argument npf from the function \link[fdesigns]{pflm}. 

After that, the function returns n.runs (see \link[fdesigns]{pflm}) plots of the optimal functions 
of the profile factor indicated. The x-axis represents the time, the y-axis represents the
values of the function of the profile factor, and the title indicated the number of run of 
each plot. 
}

\description{
Plot of \code{"flm"} objects. For the choice of a profile factor, the optimal functions are plotted.
}

\note{
For examples see \link[fdesigns]{pflm}.
}

\author{
Damianos Michaelides <\email{dm3g15@soton.ac.uk}>, Antony Overstall, Dave Woods
}








