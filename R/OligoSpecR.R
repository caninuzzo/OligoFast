#' Check the primers specificity to the target sequences
#'
#' This function allows to check the specificity of a set of primers couples to a selection of sequences.
#' The specificity will be defined according to the ability to amplify -in silico- a target dataset of sequences
#' while failing to amplify another dataset of sequences to exclude.
#' @param PrimersCouples (mandatory) a set of different primers (output from the 'OligoMatchR' or 'OligoTestR' function) OR a vector composed by the 2 primers (in 5'->3' conventional orientation) such as c("FORWARD","REVERSE").
#' @param TargetSeq (mandatory) a DNAStringSet object (such as the output from 'loadFASTAasDNA' function) containing the sequences to target by the amplification with the primers couple(s).
#' @param AvoidSeq (mandatory) a DNAStringSet object (such as the output from 'loadFASTAasDNA' function) containing the sequences for which the amplification by the primers couple(s) should be avoided.
#' @param mmCalc (mandatory, set to 2 by default) the number of mismatches used in the primers used in the in silico PCR process.
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @return A dataframe showing the number of sequences amplified from each dataset and a specificity score of the primers couples.
#' @importFrom Biostrings vcountPattern
#' @importFrom dplyr intersect
#' @export
#' @examples
#' \dontrun{
#' OligoSpecR(OligoResults[1:10,], JellyDNA, NOJellyDNA, mmCalc = 2)
#' }
#' @note
#' Example made using default settings and using only the first 10 primers obtained from OligoTestR function.
#' To save the different consensus (as .Rdata) in a given folder (here "JellyOutputs")
#' add the following argument: outpath = "path_to/JellyOutputs"

OligoSpecR <- function(PrimersCouples, TargetSeq, AvoidSeq, mmCalc = 2, outpath){

  # start time counter
  t1 <- Sys.time()
  cat("Processing...")

  # check the presence of "NNNNNNNNNN" in 'target' sequences dataset (10 successive N)
  if (length(grep("NNNNNNNNNN",as.character(TargetSeq)))>0){
    cat("\r", length(grep("NNNNNNNNNN",as.character(TargetSeq))),
        "sequences with 10 consecutive 'N' or more were removed from dataset of target sequences.")
    TargetSeq <- TargetSeq[-grep("NNNNNNNNNN",as.character(TargetSeq))]
  }

  # check the presence of "NNNNNNNNNN" in 'avoid' sequences dataset(10 successive N)
  if (length(grep("NNNNNNNNNN",as.character(AvoidSeq)))>0){
    cat("\r", length(grep("NNNNNNNNNN",as.character(AvoidSeq))),
        "sequences with 10 consecutive 'N' or more were removed from the dataset of sequences to avoid.")
    AvoidSeq <- AvoidSeq[-grep("NNNNNNNNNN",as.character(AvoidSeq))]
  }

  if (class(PrimersCouples)=="data.frame") {

    OligoSpecResults <- data.frame(matrix(nrow=nrow(PrimersCouples),ncol=5))
    colnames(OligoSpecResults) <- c("Forward",
                                    "Reverse",
                                    "Target",
                                    "Avoid",
                                    "Specificity_score")
    rownames(OligoSpecResults) <- rownames(PrimersCouples)

    for (i in 1:nrow(PrimersCouples)) {
      OligoSpecResults$Forward[i] <- PrimersCouples$Forward[i]
      OligoSpecResults$Reverse[i] <- PrimersCouples$Reverse[i]
      forward <- toupper(PrimersCouples$Forward[i])
      reverse <- toupper(OligoFast::revCompCase(PrimersCouples$Reverse[i]))

      ### Test primers couples on TARGET sequences:
      fmr_amp <- dplyr::intersect(which(Biostrings::vcountPattern(forward, TargetSeq, max.mismatch = mmCalc, fixed = FALSE)==1),
                                  which(Biostrings::vcountPattern(reverse, TargetSeq, max.mismatch = mmCalc, fixed = FALSE)==1))

      # target amplification score
      OligoSpecResults$Target[i] <- paste0(length(TargetSeq[fmr_amp]),"/",length(TargetSeq))
      targSpec <- length(TargetSeq[fmr_amp])/length(TargetSeq)

      ### Test primers couples on sequences TO AVOID:
      fmr_amp <- dplyr::intersect(which(Biostrings::vcountPattern(forward, AvoidSeq, max.mismatch = mmCalc, fixed = FALSE)==1),
                                  which(Biostrings::vcountPattern(reverse, AvoidSeq, max.mismatch = mmCalc, fixed = FALSE)==1))

      # target amplification score
      OligoSpecResults$Avoid[i] <- paste0(length(AvoidSeq[fmr_amp]),"/",length(AvoidSeq))
      avoSpec <- length(AvoidSeq[fmr_amp])/length(AvoidSeq)

      ### Specificity score:
      OligoSpecResults$Specificity_score[i] <- targSpec/avoSpec

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
            paste0(rep("#",60),collapse = ""),
            paste0(rep(" ",30)),
            "]")
        flush.console()
      }

    } # end i
  } else {
    # if PrimersCouples is not used but instead it's a vector c("FORWARD","REVERSE")
    if (length(PrimersCouples)==2 && class(PrimersCouples)=="character") {

      OligoSpecResults <- data.frame(matrix(nrow=1,ncol=4))
      colnames(OligoSpecResults) <- c("primers_couple",
                                      "Target",
                                      "Avoid",
                                      "Specificity_score")

      cat("\r Testing the primers couple ", PrimersCouples[1], " / ", PrimersCouples[2])
      forward <- toupper(PrimersCouples[1])
      reverse <- toupper(OligoFast::revCompCase(PrimersCouples[2]))

      OligoSpecResults$primers_couples[1] <- paste0(PrimersCouples[1], " / ", PrimersCouples[2])

      ### Test primers couples on TARGET sequences:
      fmr_amp <- dplyr::intersect(which(Biostrings::vcountPattern(forward, TargetSeq, max.mismatch = mmCalc, fixed = FALSE)==1),
                                  which(Biostrings::vcountPattern(reverse, TargetSeq, max.mismatch = mmCalc, fixed = FALSE)==1))


      # target amplification score
      OligoSpecResults$Target[1] <- paste0(length(TargetSeq[fmr_amp]),"/",length(TargetSeq))
      targSpec <- length(TargetSeq[fmr_amp])/length(TargetSeq)

      ### Test primers couples on sequences TO AVOID:
      fmr_amp <- dplyr::intersect(which(Biostrings::vcountPattern(forward, AvoidSeq, max.mismatch = mmCalc, fixed = FALSE)==1),
                                  which(Biostrings::vcountPattern(reverse, AvoidSeq, max.mismatch = mmCalc, fixed = FALSE)==1))

      # target amplification score
      OligoSpecResults$Avoid[1] <- paste0(length(AvoidSeq[fmr_amp]),"/",length(AvoidSeq))
      avoSpec <- length(AvoidSeq[fmr_amp])/length(AvoidSeq)

      ### Specificity score:
      OligoSpecResults$Specificity_score[1] <- targSpec/avoSpec
    }
  }

  # end time counter
  t2 <- Sys.time()
  cat(paste0("\n \r Total function execution time: ",
             as.numeric(round(difftime(t2, t1, units = "mins"),2)),
             " minutes"))
  flush.console()

  if (!missing(outpath)){
    save(OligoSpecResults, file = paste0(outpath,"/OligoSpecR_output.Rdata"))
  }

  OligoSpecResults <<- OligoSpecResults
}
