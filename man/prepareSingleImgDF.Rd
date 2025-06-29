% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imagingTools.R
\name{prepareSingleImgDF}
\alias{prepareSingleImgDF}
\title{Prepare Imaging-results tables from Cluster-Analysis SQLite databases}
\usage{
prepareSingleImgDF(
  pathDB,
  analysis = c("pa", "coloc"),
  id_cols = c("Date", "Plate_ID", "Well", "Image_ID", "Channel_Name", "Selection",
    "Selection_Area"),
  num_cols = c("Area", "Mean", "IntDen"),
  scale_num = FALSE,
  scale_cols = NULL,
  scale_fun = function(x) as.numeric(scale(x, TRUE, TRUE))
)
}
\arguments{
\item{pathDB}{Path to SQlite-DB}

\item{analysis}{pa or coloc, which table to extract}

\item{id_cols}{Columns which metadata about the image and measurement}

\item{num_cols}{Numeric Columns about the particle metrices}

\item{scale_num}{Boolean parameter to include scaled numeric metrices (Default FALSE)}

\item{scale_fun}{Scale function passed through the numeric columns}
}
\value{
A dataframe
}
\description{
Prepare Imaging-results tables from Cluster-Analysis SQLite databases
}
