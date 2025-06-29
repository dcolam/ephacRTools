% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imagingTools.R
\name{mergeSEandImg}
\alias{mergeSEandImg}
\title{Merge together the imaging-results into the Column Data of the SE}
\usage{
mergeSEandImg(se, df_img, tableType = "pa")
}
\arguments{
\item{se}{SummarizedExperiment Object with the Ephys-Data}

\item{df_img}{DataFrame with imaging results returned by prepareImgDF()}
}
\value{
A dataframe
}
\description{
Merge together the imaging-results into the Column Data of the SE
}
