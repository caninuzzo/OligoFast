#' Look for suitable oligonucleotides to create primers pairs
#'
#' This function allows to check if, among the oligonucleotides previously kept, some of them could constitute pairs to
#' amplify gene fragments suitable with the size expected for the final amplicons.
#' @param OligoChecked (mandatory) it must be the output from OligoCheckR() function, the dataframe with the oligonucleotides kept and their full characteristics.
#' @param ampSize (mandatory) a numeric vector is expected, with the minimum and maximum size range of the amplicon. Example: ampSize = c(400,600) to select amplicon from 400pb to 600pb (primers included in the size).
#' @param target (optional) a numeric vector is expected, if you want to look for primers couples only in a particular region within the sequences.
#' You should indicate the values as a numeric vector: (for example c(400,800) to target ONLY a variable region located around 400 and 800pb).
#' @param GCtail (mandatory, default: TRUE) if TRUE then the oligonucleotides presenting 3 or more G, C (and other bases corresponding) on their last nucleotides at the 3' tail will be removed from the analysis.
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @return Another dataframe of the matching oligonucleotides couples with the given parameters. In this one, all the nucleotides are given in the 5'->3' orientation. There are sorted in the 'probably most promising' order, according to 'occ_score' field value.
#' @import dplyr
#' @importFrom stringr str_count
#' @export
#' @examples
#' \dontrun{
#' JellyOligoMatched <- OligoMatchR(JellyOligoChecked, ampSize = c(300,500), GCtail = TRUE)
#' }
#' @note
#' Example made using default settings
#' To save the dataframe (as .Rdata) in a given folder (here "JellyOutputs")
#' add the following argument: outpath = "path_to/JellyOutputs"

OligoMatchR <- function(OligoChecked, ampSize, target, GCtail = TRUE, outpath) {

  # rise exception and stop function if OligoChecked is not the expected object
  if (class(OligoChecked$oligos)!="character" && class(OligoChecked)!="data.frame") {
    stop("The object provided for 'OligoChecked' argument must be the complete list output from 'OligoCheckR' function")
  }
  # rise exception and stop function if ampSize is not defined
  if (missing(ampSize) || is.null(ampSize) || length(ampSize) != 2) {
    stop("Argument 'ampSize' is missing or NULL. Please enter a range to define expected amplicon size (ex.: c(200,400)).")
  }

  # start time counter
  t1 <- Sys.time()

  # filter only target area if defined:
  if (!missing(target) && class(target) == "numeric" && target[1] < target[2]) {
    OligoChecked <- OligoChecked[!OligoChecked$meanStartPos < (target[1]+20) & !OligoChecked$meanEndPos > (target[2]-20),]
  }

  OligoMatched <- data.frame()

  ### DEV FASTER ALTERNATIVE
  # when i is set on one oligo, first filter all the other oligo which can match for the right size. then do the calculations (time saver)
  for (i in 1:nrow(OligoChecked)) {
    fmrOligoChecked <- OligoChecked[(i+1):nrow(OligoChecked),]
    fmrOligoChecked <- fmrOligoChecked %>%
                    filter(meanStartPos >= OligoChecked$meanEndPos[i] - ampSize[2] &
                           meanStartPos <= OligoChecked$meanEndPos[i] - ampSize[1] |
                           meanEndPos >= OligoChecked$meanStartPos[i] + ampSize[1] &
                           meanEndPos <= OligoChecked$meanStartPos[i] + ampSize[2])

    # progressive bar WITHOUT TIME ESTIMATION, but processing infos
    while (i>1 && i < nrow(OligoChecked)) {
      fillBar <- round(i/nrow(OligoChecked)*60,0)
      cat("\r [",
          paste0(rep("#",fillBar),collapse = ""),
          paste0(rep("_",(60-fillBar)),collapse = ""),
          "]  Confronting oligo #", i, "( /", nrow(OligoChecked), ") with", nrow(fmrOligoChecked), " compatible ones.")
      flush.console()
      break
    }
    if (nrow(fmrOligoChecked)>0){
      for (j in 1:nrow(fmrOligoChecked)) {
        # determine reading sense
        ReadSens <- OligoChecked$meanStartPos[i] - fmrOligoChecked$meanStartPos[j]

        if (ReadSens < 0) {
          ampLength <- fmrOligoChecked$meanEndPos[j] - OligoChecked$meanStartPos[i] + 1

          iline <- nrow(OligoMatched)
          OligoMatched[iline+1,1] <- as.character(OligoChecked$oligos[i])
          OligoMatched[iline+1,2] <- as.numeric(OligoChecked$meanStartPos[i])
          OligoMatched[iline+1,3] <- as.character(fmrOligoChecked$oligos[j])
          OligoMatched[iline+1,4] <- as.numeric(fmrOligoChecked$meanEndPos[j])
          OligoMatched[iline+1,5] <- ampLength
          OligoMatched[iline+1,6] <- as.numeric(mean(OligoChecked$wbs[i], fmrOligoChecked$wbs[j]))
        }

        if (ReadSens > 0) {
          ampLength <- OligoChecked$meanEndPos[i] - fmrOligoChecked$meanStartPos[j] + 1

          iline <- nrow(OligoMatched)
          OligoMatched[iline+1,1] <- as.character(fmrOligoChecked$oligos[j])
          OligoMatched[iline+1,2] <- as.numeric(fmrOligoChecked$meanStartPos[j])
          OligoMatched[iline+1,3] <- as.character(OligoChecked$oligos[i])
          OligoMatched[iline+1,4] <- as.numeric(OligoChecked$meanEndPos[i])
          OligoMatched[iline+1,5] <- ampLength
          OligoMatched[iline+1,6] <- as.numeric(mean(fmrOligoChecked$wbs[j], OligoChecked$wbs[i]))
        }
      }
    }# end j

    # display stop
    if (i == nrow(OligoChecked)) {
      cat("\r [",
          paste0(rep("#",61),collapse=""),
          "] Job finished.",
          rep(" ",30),
          "\n")
      flush.console()
    }
    #### end loading bar code ####


  }# end i

  # tidy up:
  names(OligoMatched)=c("Forward","StartPos","Reverse","EndPos","AmpliconSize","wbs")
  # reverse complement keeping up/low cases
  OligoMatched$Reverse <- OligoFast::revCompCase(OligoMatched$Reverse)

  # remove GCtail at the 3' extremities:
  if (GCtail == TRUE){
    OligoMatched <- OligoMatched[!stringr::str_count(substr(OligoMatched$Forward, nchar(OligoMatched$Forward)-4, nchar(OligoMatched$Forward)), "(G|C|Y|R|K|S|g|c|y|r|k|s)") >= 3
                                   & !stringr::str_count(substr(OligoMatched$Reverse, nchar(OligoMatched$Reverse)-4, nchar(OligoMatched$Reverse)), "(G|C|Y|R|K|S|g|c|y|r|k|s)") >= 3,]
  }


  # check after GCtail:
  if (nrow(OligoMatched)>0) {

    # sort the primers output 'OligoMatched' by their probability of occurrences in sequences
    # add 2 new fields for index calculation
    OligoMatched$cumFreqF <- rep(NA,nrow(OligoMatched))
    OligoMatched$cumFreqR <- rep(NA,nrow(OligoMatched))
    # reattribute associated cumFreq values
    for (i in 1:nrow(OligoMatched)) {
      OligoMatched$cumFreqF[i] <- OligoChecked[grep(paste0("^",OligoMatched$Forward[i],"$"),
                                                    OligoChecked$oligos),"cumFreq"]
      OligoMatched$cumFreqR[i] <- OligoChecked[grep(paste0("^",revCompCase(OligoMatched$Reverse[i]),"$"),
                                                    OligoChecked$oligos),"cumFreq"]
    }
    # reorder according to (cumFreqF+cumFreR)*(1-wbs) <-> probability of being found with strong nucleotides
    OligoMatched$wbs[is.na(OligoMatched$wbs)] <- 0
    OligoMatched$occ_score <- (OligoMatched$cumFreqF + OligoMatched$cumFreqR)*(1-OligoMatched$wbs)
    OligoMatched <- OligoMatched[order(-OligoMatched$occ_score),]

    # -> oupath: save final objects in a folder ?
    if (!missing(outpath)) {
      if (!missing(target)) {
        save(OligoMatched, file = paste0(outpath,"/OligoMatched_out_size",
                                         ampSize[1],"-",ampSize[2],
                                         "targ",target[1],"-",target[2],
                                         "GCtail",substr(as.character(GCtail),1,1),".Rdata"))
      }
      save(OligoMatched, file = paste0(outpath,"/OligoMatched_out_size",
                                       ampSize[1],"-",ampSize[2],
                                       "GCtail",substr(as.character(GCtail),1,1),".Rdata"))
    }


    t2 <- Sys.time()
    cat(paste0("Total function execution time: ",
               as.numeric(round(difftime(t2, t1, units = "mins"),2)),
               " minutes"))

    return(OligoMatched)

  } else {
    cat("No primers couples found... Try to modify some parameters in previous steps...")
  }

}
