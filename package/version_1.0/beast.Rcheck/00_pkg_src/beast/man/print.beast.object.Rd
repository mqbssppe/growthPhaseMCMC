

\name{print.beast.object}
\alias{print.beast.object}
\title{
Print function
}
\description{
This function prints a summary of objects returned by the \code{beast} function.
}
\usage{
\method{print}{beast.object}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{beast.object}, which is returned by the \code{beast} function.
}
  \item{...}{
        ignored.
}
}
\details{
The function prints a summary of the most probable number (MAP) of change-points per time-series in the form of a table, as well as a list containing the MAP number of change-points and the corresponding locations (posterior medians) per time-series.
}
\author{
Panagiotis Papastamoulis
}

