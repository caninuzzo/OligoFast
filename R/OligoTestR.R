#' Test a set of different primers on a dataset
#'
#' This function allows to test a set of different primers on a dataset of DNA sequences in order to sort them by their amplification ability.
#' @param DNAList (mandatory) a DNAStringSet object, such as the output from 'loadFASTAasDNA' function.
#' @param PrimersCouples (mandatory) a set of different primers (output from the 'OligoMatchR' function) OR a vector composed by the 2 primers (in 5'->3' conventional orientation) such as c("FORWARD","REVERSE").
#' @param mmCalc (mandatory, set to 2 by default) the number of mismatches used in the primers scoring calculation and the amplicon resolution showed.
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @return A dataframe of the matching oligonucleotides couples with the given parameters. In this one, all the nucleotides are given in the 5'->3' orientation.
#' @import Biostrings
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

  if (class(PrimersCouples)=="data.frame") {

    OligoResults <- data.frame(matrix(nrow=nrow(PrimersCouples),ncol=17))
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
                                "amplicon_resolution")

    for (i in 1:nrow(PrimersCouples)){
      rownames(OligoResults)[i] <- paste0("couple_",i)
      rownames(PrimersCouples)[i] <- paste0("couple_",i)
      forward <- toupper(PrimersCouples$Forward[i])
      reverse <- toupper(OligoFast::revCompCase(PrimersCouples$Reverse[i]))
      OligoResults$Forward[i] <- PrimersCouples$Forward[i]
      OligoResults$Reverse[i] <- PrimersCouples$Reverse[i]
        for (j in 5:0) {
         OligoResults[i,3+j] <- length(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)==1))
         OligoResults[i,9+j] <- length(which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)==1))
         if (length(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)>1))) {
           OligoResults[i,15] <- paste0("TRUE from ", j, " mm")
         }
         if (length(which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)>1))) {
           OligoResults[i,16] <- paste0("TRUE from ", j, " mm")
         }
         if (j==mmCalc){

           fmr_amp <- intersect(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)==1),which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)==1))

           fmrAmpF <- Biostrings::vmatchPattern(forward, DNAList[fmr_amp], max.mismatch = j, fixed = F)
           fmrAmpR <- Biostrings::vmatchPattern(reverse, DNAList[fmr_amp], max.mismatch = j, fixed = F)

           # to avoid error if reverse match before forward :
           ### EDIT 28.05.2024 [ startF <- (unlist(fmrAmpF@ends)-fmrAmpF@width0+1)+1 ] by :
           startF <- unlist(fmrAmpF@ends)-nchar(forward)+1
           endR <- unlist(fmrAmpR@ends)
           if (length(which(startF>endR))>0){
             fmr_amp <- fmr_amp[-which(startF>endR)]
           }


           ampliCalc <- Biostrings::subseq(DNAList[fmr_amp],
                            start = startF,
                            end = endR)
           OligoResults[i,17] <- paste0(length(DNAList[fmr_amp])-length(which(duplicated(ampliCalc)==TRUE)),
                                "/",length(DNAList[fmr_amp]))


         }
        }

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
    } # end i

  OligoResults <- OligoResults[order(-as.numeric(unlist(lapply(strsplit(OligoResults$amplicon_resolution, split = "/"), function(x) x[[1]][1]))),
                                      -OligoResults$F_0_mm,
                                      -OligoResults$R_0_mm),]
  OligoResults <<- OligoResults
  OligoTested <<- PrimersCouples

  return_list <- list(OligoTested = OligoTested, OligoResults = OligoResults)

  } else {

    # if PrimersCouples is not used but instead it's a vector c("FORWARD","REVERSE")
    if (length(PrimersCouples)==2 && class(PrimersCouples)=="character") {

      OligoResults <- data.frame(matrix(nrow=1,ncol=15))
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
                                  "amplicon_resolution")

      forward <- toupper(PrimersCouples[1])
      reverse <- toupper(OligoFast::revCompCase(PrimersCouples[2]))
      for (j in 5:0) {
        OligoResults[1,1+j] <- length(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)==1))
        OligoResults[1,7+j] <- length(which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)==1))
        if (length(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)>1))) {
          OligoResults[1,13] <- paste0("TRUE from ", j, " mm")
        }
        if (length(which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)>1))) {
          OligoResults[1,14] <- paste0("TRUE from ", j, " mm")
        }
        if (j==mmCalc){

          fmr_amp <- intersect(which(Biostrings::vcountPattern(forward, DNAList, max.mismatch = j, fixed = F)==1),which(Biostrings::vcountPattern(reverse, DNAList, max.mismatch = j, fixed = F)==1))

          fmrAmpF <- Biostrings::vmatchPattern(forward, DNAList[fmr_amp], max.mismatch = j, fixed = F)
          fmrAmpR <- Biostrings::vmatchPattern(reverse, DNAList[fmr_amp], max.mismatch = j, fixed = F)

          # to avoid error if reverse match before forward :
### EDIT 28.05.2024 [ startF <- (unlist(fmrAmpF@ends)-fmrAmpF@width0+1)+1 ] by :
          startF <- unlist(fmrAmpF@ends)-nchar(forward)+1
          endR <- unlist(fmrAmpR@ends)
          if (length(which(startF>endR))>0){
            fmr_amp <- fmr_amp[-which(startF>endR)]
          }

          ampliCalc <- Biostrings::subseq(DNAList[fmr_amp],
                              start = startF,
                              end = endR)
          OligoResults[1,15] <- paste0(length(DNAList[fmr_amp])-length(which(duplicated(ampliCalc)==TRUE)),
                                       "/",length(DNAList[fmr_amp]))


        }
      } # end j

      OligoResults[1,1] <- paste0(PrimersCouples[1], " / ", PrimersCouples[2])
      OligoResults <<- OligoResults

    } # end if
  } # end else


  # -> oupath: save final objects in a folder ?
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
