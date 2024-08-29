#' Test a set of different primers on a dataset
#'
#' This function allows to test a set of different primers on a dataset of DNA sequences in order to sort them by their amplification ability.
#' @param DNAList (mandatory) a DNAStringSet object, such as the output from 'loadFASTAasDNA' function.
#' @param PrimersCouples (mandatory) a set of different primers (output from the 'OligoMatchR' function) OR a vector composed by the 2 primers (in 5'->3' conventional orientation) such as c("FORWARD","REVERSE").
#' @param mmCalc (mandatory, set to 2 by default) the number of mismatches used in the primers scoring calculation and the amplicon resolution showed.
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @return A dataframe of the matching oligonucleotides couples with the given parameters. In this one, all the nucleotides are given in the 5'->3' orientation.
#' @importFrom Biostrings vcountPattern vmatchPattern subseq
#' @importFrom dplyr intersect
#' @export
#' @examples
#' \dontrun{
#' OligoTestR(JellyDNA, JellyOligoMatched, mmCalc = 2)
#' }
#' @note
#' Example made using default settings
#' To save the dataframe (as .Rdata) in a given folder (here "JellyOutputs")
#' add the following argument: outpath = "path_to/JellyOutputs"


OligoTestR <- function(DNAList, PrimersCouples, mmCalc = 2, outpath) {

  # start time counter
  t1 <- Sys.time()
  cat("Processing...")

  # check the presence of "NNNNNNNNNN" in sequences (10 successive N)
  if (length(grep("NNNNNNNNNN",as.character(DNAList)))>0){
    cat("\r", length(grep("NNNNNNNNNN",as.character(DNAList))),
        "sequences with 10 consecutive 'N' or more were removed from the analysis.")
    DNAList <- DNAList[-grep("NNNNNNNNNN",as.character(DNAList))]
  }


  # start process if PrimersCouples is a dataframe (classic output from previous step)
  if (class(PrimersCouples)=="data.frame") {
    # if ignGap=TRUE & filtGap=FALSE: residual "-" may be present.
    # find them and replace by "N" to avoid error in process
    if (length(grep("-",PrimersCouples$Forward))>0) {
      PrimersCouples$Forward <- gsub("-","N",PrimersCouples$Forward)
      cat("primers with '-' gaps detected, replaced by 'N'\n")
    }
    if (length(grep("-",PrimersCouples$Reverse))>0) {
      PrimersCouples$Reverse <- gsub("-","N",PrimersCouples$Reverse)
      cat("primers with '-' gaps detected, replaced by 'N'\n")
    }

    # prepare OligoResults to be filled
    OligoResults <- data.frame(matrix(nrow=nrow(PrimersCouples),ncol=18))
    colnames(OligoResults) <- c("Forward",
                                "Reverse",
                                "F_0_mm",
                                "F_1_mm",
                                "F_2_mm",
                                "F_3_mm",
                                "F_4_mm",
                                "F_5_mm",
                                "R_0_mm",
                                "R_1_mm",
                                "R_2_mm",
                                "R_3_mm",
                                "R_4_mm",
                                "R_5_mm",
                                "multipleMatchF",
                                "multipleMatchR",
                                "amplicon_resolution",
                                "warning")
    # start loop for each primers couples to test (i)
    for (i in 1:nrow(PrimersCouples)){
      # fill OligoResults and set variables
      rownames(OligoResults)[i] <- paste0("couple_",i)
      rownames(PrimersCouples)[i] <- paste0("couple_",i)
      forward <- toupper(PrimersCouples$Forward[i])
      reverse <- toupper(OligoFast::revCompCase(PrimersCouples$Reverse[i]))
      OligoResults$Forward[i] <- PrimersCouples$Forward[i]
      OligoResults$Reverse[i] <- PrimersCouples$Reverse[i]
      # process pattern recognition with different mismatches (j)
        for (j in 5:0) {
         # fill with sequences number for each mismatches thresholds set for F & R
         OligoResults[i,3+j] <- length(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)==1))
         OligoResults[i,9+j] <- length(which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)==1))
         # fill if multiples matches occur
         if (length(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)>1))) {
           OligoResults[i,15] <- paste0("TRUE from ", j, " mm (",
                                        length(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F) > 1))
                                        ,")")
         }
         if (length(which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)>1))) {
           OligoResults[i,16] <- paste0("TRUE from ", j, " mm (",
                                        length(which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F) > 1))
                                        ,")")
         }
         # Calculate the number of different amplicons when mismatches number is mmCalc
         if (j==mmCalc){
           # keep the sequences where both primers are detect and give a positive size amplicon
           fmr_amp <- dplyr::intersect(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)==1),which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)==1))
           fmrAmpF <- Biostrings::vmatchPattern(forward, DNAList[fmr_amp], max.mismatch = j, fixed = F)
           fmrAmpR <- Biostrings::vmatchPattern(reverse, DNAList[fmr_amp], max.mismatch = j, fixed = F)

           # to avoid error if reverse match before forward :
           startF <- unlist(fmrAmpF@ends)-nchar(forward)+1
           endR <- unlist(fmrAmpR@ends)
           compFR <- which(startF > endR)
           overlapFR <- which(((startF+nchar(forward))-(endR-(nchar(reverse))))>0)
           # when the result is negative, remove it:
           if (length(compFR) > 0) {
             fmr_amp <- fmr_amp[-compFR]
             startF <- startF[-compFR]
             endR <- endR[-compFR]
           }
           # when the primers overlap, remove it:
           if (length(overlapFR) > 0) {
             fmr_amp <- fmr_amp[-overlapFR]
             startF <- startF[-overlapFR]
             endR <- endR[-overlapFR]
           }
           # when, due to mmCalc>0, endR is over real sequence end, adjust it to avoid error:
           if (length(which(width(DNAList[fmr_amp])<endR))>0){
            endR[which(width(DNAList[fmr_amp])<endR)] <- endR[which(width(DNAList[fmr_amp])<endR)] - mmCalc
           }
           # keep the amplicon only if there are for the mmCalc value given
           if (length(fmr_amp)>0) {
             ampliCalc <- Biostrings::subseq(DNAList[fmr_amp],
                              start = startF,
                              end = endR)
             # detect potential 'outliers' / i.e. introns/exons
             if (length(which(ampliCalc@ranges@width < (mean(ampliCalc@ranges@width)-20) |
                              ampliCalc@ranges@width > (mean(ampliCalc@ranges@width)+20)))>0) {
               OligoResults[i,18] <- paste0(length(which(ampliCalc@ranges@width < (mean(ampliCalc@ranges@width)-20) |
                                                           ampliCalc@ranges@width > (mean(ampliCalc@ranges@width)+20)))
                                            ," potential introns/exons")
             }
             # give an amplification ability / resolution:
             OligoResults[i, 17] <- paste0(length(DNAList[fmr_amp]) -
                                             length(which(duplicated(ampliCalc) == TRUE)),
                                           "/", length(DNAList[fmr_amp]))
           } else {
             # fmr_amp = integer(0) for the given j==mmCalc
             # give an amplification ability / resolution:
             OligoResults[i, 17] <- NA
           }
        } # end j==mmCalc
      } # end j loop
      # progressive bar WITHOUT TIME ESTIMATION, but processing infos
      while (i>1 && i < nrow(PrimersCouples)) {
        fillBar <- round(i/nrow(PrimersCouples)*60,0)
        cat("\r [",
            paste0(rep("#",fillBar),collapse = ""),
            paste0(rep("_",(60-fillBar)),collapse = ""),
            "]  Testing primers couple #", i, " /", nrow(PrimersCouples))
        flush.console()
        break
      }
      if (i==nrow(PrimersCouples)){
        cat("\r [",
            paste0(rep("#", 60),collapse = ""),
            "]",
            paste0(rep(" ", 30)))
        flush.console()
      }
    } # end i time display

  # order dataframe according to amplification ability and resolution
  if (all(is.na(OligoResults$amplicon_resolution))) {
    stop("Problem in the primers couples. None are able to produce amplicons")
  }

  OligoResults <- OligoResults[order(-as.numeric(unlist(lapply(strsplit(OligoResults$amplicon_resolution, split = "/"), function(x) x[[1]][1]))),
                                      -OligoResults$F_0_mm,
                                      -OligoResults$R_0_mm),]
  # prepare outputs
  OligoResults <<- OligoResults
  OligoTested <<- PrimersCouples
  return_list <- list(OligoTested = OligoTested, OligoResults = OligoResults)

  } else {

    # if PrimersCouples is not used but instead it's a vector c("FORWARD","REVERSE")
    if (length(PrimersCouples)==2 && class(PrimersCouples)=="character") {
      # if ignGap=TRUE & filtGap=FALSE: residual "-" may be present.
      # find them and replace by "N" to avoid error in process
      if (length(grep("-",PrimersCouples[1]))>0){
        PrimersCouples[1] <- gsub("-","N",PrimersCouples[1])
        cat("Gap detected in the forward primer, replaced by 'N'")
      }
      if(length(grep("-",PrimersCouples[2]))>0) {
        PrimersCouples[2]<- gsub("-","N",PrimersCouples[2])
        cat("Gap detected in the reverse primer, replaced by 'N'")
      }

      # prepare OligoResults
      OligoResults <- data.frame(matrix(nrow=1,ncol=16))
      colnames(OligoResults) <- c("F_0_mm",
                                  "F_1_mm",
                                  "F_2_mm",
                                  "F_3_mm",
                                  "F_4_mm",
                                  "F_5_mm",
                                  "R_0_mm",
                                  "R_1_mm",
                                  "R_2_mm",
                                  "R_3_mm",
                                  "R_4_mm",
                                  "R_5_mm",
                                  "multipleMatchF",
                                  "multipleMatchR",
                                  "amplicon_resolution",
                                  "warning")
      # set variables
      forward <- toupper(PrimersCouples[1])
      reverse <- toupper(OligoFast::revCompCase(PrimersCouples[2]))
      # process pattern recognition with different mismatches (j)
      for (j in 5:0) {
        # fill with sequences number for each mismatches thresholds set for F & R
        OligoResults[1,1+j] <- length(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)==1))
        OligoResults[1,7+j] <- length(which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)==1))
        # fill if multiples matches occur
        if (length(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)>1))) {
          OligoResults[1,13] <- paste0("TRUE from ", j, " mm (",
                                       length(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F) > 1))
                                       ,")")
        }
        if (length(which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)>1))) {
          OligoResults[1,14] <- paste0("TRUE from ", j, " mm (",
                                       length(which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F) > 1))
                                       ,")")
        }
        # Calculate the number of different amplicons when mismatches number is mmCalc
        if (j==mmCalc){
          # keep the sequences where both primers are detect and give a positive size amplicon
          fmr_amp <- dplyr::intersect(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)==1),
                                      which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)==1))

          fmrAmpF <- Biostrings::vmatchPattern(forward, DNAList[fmr_amp], max.mismatch = j, fixed = F)
          fmrAmpR <- Biostrings::vmatchPattern(reverse, DNAList[fmr_amp], max.mismatch = j, fixed = F)

          # to avoid error if reverse match before forward
          startF <- unlist(fmrAmpF@ends)-nchar(forward)+1
          endR <- unlist(fmrAmpR@ends)
          compFR <- which(startF > endR)
          # when the result is negative remove it
          if (length(compFR) > 0) {
            fmr_amp <- fmr_amp[-which(startF > endR)]
            startF <- startF[-compFR]
            endR <- endR[-compFR]
          }
          # keep only the amplicons
          ampliCalc <- Biostrings::subseq(DNAList[fmr_amp],
                              start = startF,
                              end = endR)
          # detect potential 'outliers' / i.e. introns/exons
          if (length(which(ampliCalc@ranges@width < (mean(ampliCalc@ranges@width)-20) |
                           ampliCalc@ranges@width > (mean(ampliCalc@ranges@width)+20)))>0){
            OligoResults[1,16] <- paste0(length(which(ampliCalc@ranges@width < (mean(ampliCalc@ranges@width)-20) |
                                                        ampliCalc@ranges@width > (mean(ampliCalc@ranges@width)+20)))
                                         ," potential introns/exons")
          }
          # give an amplification ability / resolution:
          OligoResults[1,15] <- paste0(length(DNAList[fmr_amp])-length(which(duplicated(ampliCalc)==TRUE)),
                                       "/",length(DNAList[fmr_amp]))
        } # end j==mmCalc loop
      } # end j loop
      OligoResults <<- OligoResults
    } # end if
  } # end else

  # if outpath: save final objects in a folder ?
  if (!missing(outpath)) {
    save(return_list, file = paste0(outpath,"/OligoTestR_output.Rdata"))
  }

  # end time counter
  t2 <- Sys.time()
  cat(paste0("\n \r Total function execution time: ",
             as.numeric(round(difftime(t2, t1, units = "mins"),2)),
             " minutes"))
  flush.console()
}
