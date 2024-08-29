#' in silico PCR for primers coverage and investigation
#'
#' This function allows to perform an in silico PCR on a given dataset of sequences using a given couple of primers. The amplification is based on characters (bases) only, none biochemical parameters are considered.
#' The number of maximum mismatches allowed can be set, ambiguous nucleotides are taken into account.
#' The function returns the amplicons produced, but as well a list of the non amplified sequences and the non resolutive amplicon.s.
#' The resolutive power of the amplicon for the given primers couple is also showed, and a size range of the amplicons can be displayed on a graph.
#' @param DNAList (mandatory) a DNAStringSet object, such as the output from 'loadFASTAasDNA' function.
#' @param forward (mandatory) character, which corresponds to the forward primers. Note: the primers are given in the conventional 5'->3' orientation.
#' @param reverse (mandatory) character, which corresponds to the reverse primers. Note: the primers are given in the conventional 5'->3' orientation.
#' @param maxdiffF (mandatory set to 2 by default) number, which corresponds to the maximum mismatches allowed between sequences and forward primer.
#' @param maxdiffR (mandatory set to 2 by default) number, which corresponds to the maximum mismatches allowed between sequences and reverse primer.
#' @param targetReg (optional) the approximate target region to target the amplification within a specific region. Note: choose a wider region to include potential indels and gaps!
#' @param amplifData (optional, set to TRUE by default) If TRUE, informations about the amplification of each sequences will be gathered in a dataframe (PCROut).
#' It generates as well a visual plot of the size of the amplicon: AmplifOverview.
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @return Different objects including: amplicons, non-resolutive amplicons, sequences not-amplified and a dataframe that sums up the amplification process.
#' @importFrom Biostrings vcountPattern vmatchPattern subseq DNAStringSet
#' @importFrom dplyr intersect
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' insilicoPCR(JellyDNA, JellyOligoMatched$Forward, JellyOligoMatched$Reverse, maxdiffF=2, maxdiffR=2)
#' }
#' @note
#' Example made using default settings
#' To save the outputs (as .Rdata) in a given folder (here "JellyOutputs")
#' add the following argument: outpath = "path_to/JellyOutputs"

insilicoPCR <- function(DNAList, forward, reverse, maxdiffF=2, maxdiffR=2, targetReg, amplifData=TRUE, outpath) {

  if (class(DNAList)!="DNAStringSet"){
    stop("Error: DNAList should be a list of DNA sequences belonging to the class DNAStringSet")
  }

  # check the presence of "NNNNNNNNNN" in sequences (10 successive N)
  if (length(grep("NNNNNNNNNN",as.character(DNAList)))>0){
    cat("\r", length(grep("NNNNNNNNNN",as.character(DNAList))),
        "sequences with 10 consecutive 'N' or more were removed from the analysis.")
    DNAList <- DNAList[-grep("NNNNNNNNNN",as.character(DNAList))]
  }

  cat(paste0("Performing in silico PCR with: 5'",
             forward,"3' and 5'",reverse,"3' (",
             maxdiffF, " and ", maxdiffR, " maximum mismatches allowed/primers respectively)\n"))

  # convert forward et reverse for matching with DNAList
  forward <- toupper(forward)
  reverse <- toupper(OligoFast::revCompCase(reverse))

  if(!missing(targetReg)){
    if (length(targetReg)==2 && targetReg[1]<targetReg[2]){
      DNAList <- Biostrings::subseq(DNAList, start = targetReg[1], end = targetReg[2])
    } else {
      stop("Error in the target region specified (targetReg argument...)")
    }
  }


  # match forward with sequences
  couF <- Biostrings::vcountPattern(forward, DNAList, max.mismatch = maxdiffF, fixed = F)
  couF0 <- Biostrings::vcountPattern(forward, DNAList, max.mismatch = 0, fixed = F)
  keeF <- which(couF>0)
  matF <- Biostrings::vmatchPattern(forward, DNAList, max.mismatch = maxdiffF, fixed = F)

  # match reverse with sequences
  couR <- Biostrings::vcountPattern(reverse, DNAList, max.mismatch = maxdiffR, fixed = F)
  couR0 <- Biostrings::vcountPattern(reverse, DNAList, max.mismatch = 0, fixed = F)
  keeR <- which(couR>0)
  matR <- Biostrings::vmatchPattern(reverse, DNAList, max.mismatch = maxdiffR, fixed = F)


  ### Prepare output dataframes:
  # full data :
  PCROut <- data.frame(matrix(nrow=length(DNAList), ncol=6))
  colnames(PCROut) <- c("seq_identifier",
                         "amplicon_size", # if amplicon
                         "exact_F_matching", # if mismatchesF=0
                         "exact_R_matching", # if mismatchesR=0
                         "multipleMatchF",
                         "multipleMatchR")

  PCROut$seq_identifier <- names(DNAList)

  if (length(dplyr::intersect(keeF,keeR))>0) {
    for (e in dplyr::intersect(keeF,keeR)) {
      if (couF[e]==1 && couR[e]==1) {
        PCROut$amplicon_size[e] <- unlist(matR@ends[e]) - (unlist(matF@ends[e])-matF@width0[e]+1) + 1
      }
      if (couF[e]>1) {
        PCROut$multipleMatchF[e] <- couF[e]-1
      } else {
        PCROut$multipleMatchF[e] <- 0
      }
      if (couR[e]>1) {
        PCROut$multipleMatchR[e] <- couR[e]-1
      } else {
        PCROut$multipleMatchR[e] <- 0
      }
    }
  }
  # indicates if exact matching
  PCROut$exact_F_matching[which(couF0==1)] <- TRUE
  PCROut$exact_R_matching[which(couR0==1)] <- TRUE

  PCROut$exact_F_matching[which(is.na(PCROut$exact_F_matching))] <- FALSE
  PCROut$exact_R_matching[which(is.na(PCROut$exact_R_matching))] <- FALSE

  # gather amplicon produced (if there are)
  if (length(which(PCROut$amplicon_size>0))>0) {
    amplified <- NULL
    for (e in which(PCROut$amplicon_size>0)) {
      # when x mismatches allowed, some last nucleotides of reverse oligo can be missing at the end of the sequence
      # deal with it without stopping the process with an error:
      if (length(DNAList[[e]])<unlist(matR@ends[e])) {
        amplified <- append(amplified,
                            as.character(Biostrings::subseq(DNAList[[e]],
                                                            start = (unlist(matF@ends[e])-matF@width0[e]+1)+1,
                                                            end = length(DNAList[[e]]))))
      } else {
        amplified <- append(amplified,
                            as.character(Biostrings::subseq(DNAList[[e]],
                                                            start = (unlist(matF@ends[e])-matF@width0[e]+1)+1,
                                                            end = unlist(matR@ends[e]))))
      }
    }
    names(amplified) <- PCROut$seq_identifier[which(PCROut$amplicon_size>0)]
    amplified <- Biostrings::DNAStringSet(amplified)
  } else {
    cat("No sequences were amplified with the given primers")
  }

  # gather duplicated amplicons:
  if (length(which(duplicated(amplified)==TRUE))>0) {
    non_resolutive_amplicon <- amplified[which(duplicated(amplified)==TRUE)]
  }

  # gather sequences not amplified
  not_amplified <- as.character(DNAList[which(is.na(PCROut$amplicon_size))])
  names(not_amplified) <- PCROut$seq_identifier[which(is.na(PCROut$amplicon_size))]
  not_amplified <- Biostrings::DNAStringSet(not_amplified)

  cat(paste0(length(amplified),
      " amplicons obtained from ",
      length(DNAList),
      " sequences. [ Mean length of ",
      round(mean(amplified@ranges@width)),
      "pb, min. size: ",
      min(amplified@ranges@width),
      "pb max. size: ",
      max(amplified@ranges@width),
      "pb ]"))

  cat("\n Resolution of the primers couple: ",
      length(amplified)-length(non_resolutive_amplicon),
      "/",length(amplified)," amplicons.")

  cat("\nNote: If amplicons size are significatively different, consider curation process with CutAmplicon() function.")

  # return outputs
  if (amplifData==TRUE){
    AmplifOverview <<- ggplot(PCROut, aes(x = amplicon_size)) +
                       geom_histogram(aes(y=after_stat(count)),
                                      binwidth = 1, color = "black", fill = "lightgreen", alpha = 0.7)

    PCROut <<- PCROut
    amplified <<- amplified
    non_resolutive_amplicon <<- non_resolutive_amplicon
    not_amplified <<- not_amplified
  } else {
    PCROut <<- PCROut
    amplified <<- amplified
    non_resolutive_amplicon <<- non_resolutive_amplicon
    not_amplified <<- not_amplified
  }

  if (!missing(outpath)){
    PCR_results <- list(PCROut = PCROut,
                        amplified = amplified,
                        non_resolutive_amplicon = non_resolutive_amplicon,
                        not_amplified = not_amplified)
    save(PCR_results, file = paste0(outpath,"/insilicoPCR_output.Rdata"))
  }

}
