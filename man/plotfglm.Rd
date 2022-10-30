\name{plot.fglm}

\alias{plot.fglm}

\title{Plot fglm Objects}

\usage{
## S3 method for class "fglm"
plot(obj, ...)
}

\arguments{
\item{obj}{An object of class "fglm".}
\item{...}{Additional arguments to be passed to other methods.}
}

\value{
The function returns the question:
"Which profile factor to plot?".
The answer needs to be an integer to specify the profile factor for which 
to plot the optimal functions in every run. The value needs to be between 
1 and the argument npf from the function \link[fdesigns]{pfglm}. 

After that, the function returns n.runs (see \link[fdesigns]{pfglm}) plots of the optimal functions, 
of the profile factor indicated. The x-axis represents the time, the y-axis represents the
values of the function of the profile factor, and the title indicated the number of run of 
each plot. 
}

\description{
Plot fglm Objects. For the choice of a profile factor, the optimal functions are plotted.
}

\examples{
For examples see \link[fdesigns]{pflm}.
}










