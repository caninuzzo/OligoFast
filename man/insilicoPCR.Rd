% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insilicoPCR.R
\name{insilicoPCR}
\alias{insilicoPCR}
\title{in silico PCR for primers coverage and investigation}
\usage{
insilicoPCR(
  DNAList,
  forward,
  reverse,
  maxdiffF = 2,
  maxdiffR = 2,
  targetReg,
  amplifData = TRUE,
  outpath
)
}
\arguments{
\item{DNAList}{(mandatory) a DNAStringSet object, such as the output from 'loadFASTAasDNA' function.}

\item{forward}{(mandatory) character, which corresponds to the forward primers. Note: the primers are given in the conventional 5'->3' orientation.}

\item{reverse}{(mandatory) character, which corresponds to the reverse primers. Note: the primers are given in the conventional 5'->3' orientation.}

\item{maxdiffF}{(mandatory set to 2 by default) number, which corresponds to the maximum mismatches allowed between sequences and forward primer.}

\item{maxdiffR}{(mandatory set to 2 by default) number, which corresponds to the maximum mismatches allowed between sequences and reverse primer.}

\item{targetReg}{(optional) the approximate target region to target the amplification within a specific region. Note: choose a wider region to include potential indels and gaps!}

\item{amplifData}{(optional, set to TRUE by default) If TRUE, informations about the amplification of each sequences will be gathered in a dataframe (PCROut).
It generates as well a visual plot of the size of the amplicon: AmplifOverview.}

\item{outpath}{(optional) the path to the FOLDER in which the function outputs will be stored.}
}
\value{
Different objects including: amplicons, non-resolutive amplicons, sequences not-amplified and a dataframe that sums up the amplification process.
}
\description{
This function allows to perform an in silico PCR on a given dataset of sequences using a given couple of primers. The amplification is based on characters (bases) only, none biochemical parameters are considered.
The number of maximum mismatches allowed can be set, ambiguous nucleotides are taken into account.
The function returns the amplicons produced, but as well a list of the non amplified sequences and the non resolutive amplicon.s.
The resolutive power of the amplicon for the given primers couple is also showed, and a size range of the amplicons can be displayed on a graph.
}
\note{
Example made using default settings
To save the outputs (as .Rdata) in a given folder (here "JellyOutputs")
add the following argument: outpath = "path_to/JellyOutputs"
}
\examples{
\dontrun{
insilicoPCR(JellyDNA, JellyOligoMatched$Forward, JellyOligoMatched$Reverse, maxdiffF=2, maxdiffR=2)
}
}
