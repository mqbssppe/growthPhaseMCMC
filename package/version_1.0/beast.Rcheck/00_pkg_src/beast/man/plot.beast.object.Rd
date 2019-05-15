

\name{plot.beast.object}
\alias{plot.beast.object}
\title{
Plot function
}
\description{
This function plots objects returned by the \code{beast} function. All output is diverted to a pdf file provided in the \code{fileName} argument.
}
\usage{
\method{plot}{beast.object}(x, fileName, width, height, pointsize, ylab, xlab, timeScale, myPal, ...)
}
\arguments{
  \item{x}{
An object of class \code{beast.object}, which is returned by the \code{beast} function.
}
  \item{fileName}{
Name of the output pdf file.
}
  \item{width}{
Width of pdf file. Default: 9
}
  \item{height}{
Height pdf file. Default: 6
}
  \item{pointsize}{
Pointsize. Default: 12
}
  \item{ylab}{
\eqn{y}-axis label. Default: x.
}
  \item{xlab}{
\eqn{x}-axis label. Default: t.
}
  \item{timeScale}{
A multiplicative-factor which will be used to scale the \eqn{x}-axis labels. For example, if time-points correspond to 10-minute periods, then setting \code{timeScale = 1/6} will make the \eqn{x}-axis labels correspond to hours. Default: \code{timeScale = 1}.
}
  \item{myPal}{
Vector of colors that will be used to produce the plot with all time-series overlayed. If the distinct values of the inferred numbers of change-points is less than 10,  the \code{Set1} pallete of the \code{RColorBrewer} library is used. Otherwise, the user has to manually define the colors.
}
  \item{...}{
        ignored.
}
}
\details{
The function will produce a plot with all time-series coloured according to the corresponding number of change-points. In addition, it will generate individual plots per time-series displaying the observed data with boxplots which summarize the posterior distribution of change-points locations, conditionally on the most probable number of change-points.
}
\author{
Panagiotis Papastamoulis
}

