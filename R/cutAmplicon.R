#' Truncate the amplicons resulting from insilicoPCR
#'
#' The goal of this function is simply to gather one part of the amplicon (or sequences) from a location to another (in bp).
#' @param amplicon (mandatory) DNAStringSet object: the amplicons or sequences to cut.
#' @param cut (mandatory) a numeric vector containing the start and the end (in bp) to define where to cut. Example: c(200,800) to keep only the region from 200bp to 800bp.
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @return a DNAStringSet object with the remaining parts of the amplicons or sequences.
#' @importFrom Biostrings subseq
#' @export
#' @examples
#' \dontrun{
#' JellyFrag_100_to_500 <- OligoTestR(JellyDNA, cut = c(100,500))
#' }
#' @note
#' To save the outputs (as .Rdata) in a given folder (here "JellyOutputs")
#' add the following argument: outpath = "path_to/JellyOutputs"

cutAmplicon <- function(amplicon, cut, outpath) {
  #
  if (cut[1] < cut[2] && cut[2] < min(width(amplicon))) {
    cutsubseq <- Biostrings::subseq(amplicon, start = cut[1], end = cut[2])
  }
  return(cutsubseq)
}
