#' Build consensus sequences for each cluster.
#'
#' This function aims to build a consensus sequence for each cluster of aligned sequences according to parameters
#' set by the user.
#' @param AlignedClust (mandatory) the output from ShuffleAndAlign() which is a list of clusters of aligned sequences.
#' @param ntThresHi (mandatory, set to 0.9 by default) a nucleotide at the given position will be kept in the consensus sequence
#' if its observation frequency (0->1) is equal or higher to this threshold.
#' @param ntThresLo (mandatory, set to 0.3 by default) if the frequency of the most observed nucleotide at a given position is
#' below this threshold, then it will result as an ambiguous base (i.e. 'N') for this position.
#' Moreover, when ignGap=FALSE, if the observation frequency of gaps (i.e. '-') for this position is ABOVE this threshold, then the
#' corresponding nucleotide will be put in lowercase in the consensus sequence to highlight its weakness.
#' @param ignGap (mandatory, set to FALSE by default) consider or not the gaps (i.e. '-') in the sequence alignments. To get more
#' accurate results, ignGap should be FALSE (as default) but in some cases if the user wants to remove gaps it is possible by setting TRUE.
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @param suppData (mandatory, set to FALSE by default) if TRUE, an object is created: a list of nucleotides frequency matrixes for each cluster
#' @return A list of consensus sequences (character) resulting from each clusters.
#' @import Biostrings
#' @export
#' @examples
#' \dontrun{
#' JellyConsensus <- BuildConsensus(JellyAlignments, ntThresHi=0.9, ntThresLo=0.3, ignGap=FALSE, suppData=FALSE)
#' }
#' @note
#' Example made using default settings
#' To save the different consensus (as .Rdata) in a given folder (here "JellyOutputs")
#' add the following argument: outpath = "path_to/JellyOutputs"

BuildConsensus <- function(AlignedClust, ntThresHi=0.9, ntThresLo=0.3, ignGap=FALSE, outpath, suppData=FALSE){

  # rise exception and stop function if AlignedClust is not the list of msa DNA Alignment from ShuffleAndAlign function
  if (class(AlignedClust[[1]])!="MsaDNAMultipleAlignment" && class(AlignedClust)!="list") {
    stop("The object provided for 'AlignedClust' argument must be the complete list output from 'ShuffleAndAlign' function")
  }

  # start time counter to estimate function progression
  t1<- Sys.time()

  # build list to store supplementary data (frequency matrixes) if required
  if (isTRUE(suppData)){
    freqData <- list()
  }

  # -> Building consensus sequences:
  ConsensusClust <- list()
  time_storage <- NULL
  for (e in 1:length(AlignedClust)) {
    # time estimation start
    ti <- Sys.time()

    # setting variables to use and frequency matrix
    mat <- Biostrings::consensusMatrix(AlignedClust[[e]])
    mat <- OligoFast::DegToNucl(mat)
    freqMat <- mat/colSums(mat)
    cons <- NULL

    # -> The user wants to ignore gaps:
    if (isTRUE(ignGap)){
      freqMat["-",] <- 0
      freqMat <- sweep(freqMat, 2, colSums(freqMat), "/")
    }

    # -> Start iteration position by position within the frequency matrix:
    # let's consider X the nucleotide with the highest observation frequency
    for (i in 1:ncol(freqMat)) {
      ntEval <- sort(freqMat[,i], decreasing = T)

      # (1) X is higher than ntThresHi:
      if (max(ntEval) >= ntThresHi) {
        cons[i] <- names(which.max(freqMat[,i]))

        # (2) X is lower than ntThresLo:
      } else if (max(ntEval) < ntThresLo) {
        cons[i] <- "N"

        # (3) X is between ntThresHi and ntThresLo:
      } else {
        ntObs <- names(which(round(cumsum(ntEval),2)<1))[-which(names(ntEval)=="-")] # here round() avoid an issue

        # (3.A.) X is associated to 1 nucleotide and a gap:
        if (length(ntObs)==1) {
          cons[i] <- ntObs

          # (3.B.) X is associated to 2 nucletoides and (potentially) a gap:
        } else if (length(ntObs)==2) {
          cons[i] <- as.character(IUPAC_COMPLETE_MAP[which(names(IUPAC_COMPLETE_MAP)==
                                                             paste(ntObs,collapse = ""))])

          # (3.C.) X is associated to 3 nucleotides and (potentially) a gap:
        } else if (length(ntObs)==3) {
          cons[i] <- as.character(IUPAC_COMPLETE_MAP[which(names(IUPAC_COMPLETE_MAP)==
                                                             paste(ntObs,collapse = ""))])

          # (3.D.) the nucleotide diversity is considered too high to establish a consensus, set N:
        } else {
          cons[i] <- "N"
        }
      }

      # -> Gap treatment using lowercase:
      if (as.numeric(ntEval["-"]) >= ntThresLo && ignGap==FALSE) {
        cons[i] <- tolower(cons[i])
      }
    } # end i loop

    # -> SuppData: fill the list with freqMat of each cluster
    if (isTRUE(suppData)) {
      freqData[[e]] <- freqMat
    }

    ConsensusClust[[e]] <- cons

    # end processing time for storage
    tf <- Sys.time()

    #### Displaying time estimation with progressive loading bar ####
    while (e<length(AlignedClust)) {
      time_storage[e] <- as.numeric(round(difftime(tf, ti, units = "mins"),2))
      # bar
      fillBar <- round(e/length(AlignedClust)*60,0)
      # display
      cat("\r [",
          paste0(rep("#",fillBar),collapse = ""),
          paste0(rep("_",(60-fillBar)),collapse = ""),
          "] Building consensus for cluster ", e, "/", length(AlignedClust),". Remaining time estimation: ~",
          round(mean(time_storage)*(length(AlignedClust)-e),2),
          "min.     ")
      flush.console()
      break
    }
    if (e==length(AlignedClust)) {
      cat("\r [",
          paste0(rep("#",61),collapse=""),
          "] Job finished.",
          rep(" ",35),
          "\n")
      flush.console()
    }
    #### end loading bar code ####

  } #end e loop

  # -> SuppData: save final list in R environment
  if (isTRUE(suppData)) {
    freqData <<- freqData
  }

  # -> oupath: save final objects in a folder ?
  if (!missing(outpath)){
    save(ConsensusClust, file = paste0(outpath,"/Consensus_out_Hi",
                                 gsub("\\.", "", format(ntThresHi, nsmall = 1)),
                                 "Lo", gsub("\\.", "", format(ntThresLo, nsmall = 1)),
                                 "Gap", substr(as.character(ignGap),1,1),".Rdata"))
  }

  # stop time counter and display exec. time:
  t2 <- Sys.time()
  cat(paste0("Total function execution time: ",
             as.numeric(round(difftime(t2, t1, units = "mins"),2)),
             " minutes"))

  return(ConsensusClust)

}
