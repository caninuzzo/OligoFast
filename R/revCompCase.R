#' Reverse complement (case sensitive)
#'
#' This function allows to create the reverse complement of an oligonucleotides keeping the characters in upper and lowercases (as well compatible with gaps and degenerated bases).
#' @param seq (mandatory) the sequence or oligonucleotides to set as reverse complement (character or DNAString object).
#' @param class (default: vector) the class of the reverse complement returned by the function (could be "vector" or "DNAString").
#' @return the reverse complement as a vector (chain of characters) or a DNAString object.
#' @importFrom Biostrings reverse DNAString
#' @export

revCompCase <- function(seq, class = "vector"){
  comp <- chartr("AGCTYRWSKMDVHBNagctyrwskmdvhbn-","TCGARYWSMKHBDVNtcgarywsmkhbdvn-",seq)
  revcomp <- Biostrings::reverse(comp)
  if (class=="DNAString"){
    revcomp <- Biostrings::DNAString(revcomp)
  }
  return(revcomp)
}
