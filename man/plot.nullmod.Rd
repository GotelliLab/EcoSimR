% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/null_model_engine.R
\name{plot.nullmod}
\alias{plot.nullmod}
\title{plot a histogram null model}
\usage{
\method{plot}{nullmod}(x, ...)
}
\arguments{
\item{x}{the null model object to plot.}

\item{...}{Other variables to be passed on to base plotting.}
}
\description{
Plot a null model object.
}
\details{
The "hist" plot type is common to all EcoSimR modules. The blue histogram represents the NRep values of the metric for the simulated assemblages. The red vertical line represents the metric value for the real assemblage. The two pairs of vertical dashed black lines represent the one-tailed (long dash) and two-tailed (short dash) 95% confidence exact confidence intervals of the simulated data.
}

