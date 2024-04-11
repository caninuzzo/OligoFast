#' IUPAC_COMPLETE_MAP
#'
#' Creates a complete mapping of all nucleotides, as IUPAC_CODE_MAP from Biostrings, but including gaps as "-"
#' and specifications to attribute degenerated code to nucleotides
#' @usage Especially called in BuildConsensus function
#' @export
IUPAC_COMPLETE_MAP <- NULL
IUPAC_CODE_MAP <- Biostrings::IUPAC_CODE_MAP

for (i in 1:length(IUPAC_CODE_MAP)) {
  ncharIUPAC <- as.numeric(nchar(IUPAC_CODE_MAP[i]))
  if (ncharIUPAC==1){
    IUPAC_COMPLETE_MAP[names(IUPAC_CODE_MAP[i])] <- as.character(IUPAC_CODE_MAP[i])
  }

  if (ncharIUPAC>1){
    d <- apply(gtools::permutations(ncharIUPAC,ncharIUPAC,
                            unlist(strsplit(as.character(IUPAC_CODE_MAP[i]),split="")), repeats.allowed = F),
               1,paste,collapse="")
    IUPAC_COMPLETE_MAP[d] <- names(IUPAC_CODE_MAP[i])
  }
}

# add Gaps:
IUPAC_COMPLETE_MAP <- c(IUPAC_COMPLETE_MAP,"-"="-")
