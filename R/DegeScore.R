#' Calculate Degenerescence Score
#'
#' This function allows to calculate the degenerescence score for a given oligonucleotide or sequence.
#' (W, S, R, Y, K, M have a 2-deg. score ; 
#' B, D, H, V have a 3-deg. score ; 
#' N has a 4-deg. score).
#' @param x (mandatory) an oligonucleotide or sequence (character)
#' @return The degenerescence score (numeric) associated to the oligonucleotide/sequences provided as x.
#' @export

# simple function to calculate the degenerescence score of oligos
DegeScore <- function(x) {
  Dscore <- sum(
    sum(grepl("W",x,ignore.case = T) | 
          grepl("S",x,ignore.case = T) | 
          grepl("R",x,ignore.case = T) | 
          grepl("Y",x,ignore.case = T) | 
          grepl("K",x,ignore.case = T) | 
          grepl("M",x,ignore.case = T))*2,
    sum(grepl("B",x,ignore.case = T) | 
          grepl("D",x,ignore.case = T) | 
          grepl("H",x,ignore.case = T) | 
          grepl("V",x,ignore.case = T))*3,
    sum(grepl("N",x,ignore.case = T))*4)
  
  return(Dscore)
}
