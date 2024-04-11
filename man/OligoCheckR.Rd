% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OligoCheckR.R
\name{OligoCheckR}
\alias{OligoCheckR}
\title{Check and keep similar oligonucleotides}
\usage{
OligoCheckR(OligoFound, filtGaps = TRUE, nOcc, plot = TRUE, outpath)
}
\arguments{
\item{OligoFound}{(mandatory) it should be the output from OligoFindR() function, the dataframe with the oligonucleotides and their characteristics.}

\item{filtGaps}{(mandatory, default: TRUE) if TRUE, the oligonucleotides with gaps are removed from the analysis.}

\item{nOcc}{(mandatory) the number of times the oligonucleotide is found through the diffenrent clusters. The theorical maximum is the number of cluster
(meaning the oligonucleotide is found everywhere and well conserved). The higher is nOcc, the stronger should be the oligonucleotide.}

\item{plot}{(optional, default: TRUE) if TRUE, creates an overview of the most abundant oligonucleotides and their location on a graph, using ggplot2.
The picture is a new variable created named "oligOverview" (print(oligOverview) to display it).}

\item{outpath}{(optional) the path to the FOLDER in which the function outputs will be stored.}
}
\value{
Another dataframe of the remaining and most abundant oligonucleotides.
}
\description{
This function allows to look for the most abundant similar oligonucleotides from each consensus to keep
the most redundants which are above a threshold defined by the user.
}
\note{
Example made using default settings
To save the dataframe (as .Rdata) in a given folder (here "JellyOutputs")
add the following argument: outpath = "path_to/JellyOutputs"
}
\examples{
\dontrun{
JellyOligoChecked <- OligoCheckR(JellyOligoFound, filtGaps = TRUE, nOcc = 56, plot = TRUE)
}
}