#' Degenerated Nucleotide to Nucleotide 
#'
#' This function should not be used by an user, but it is required for BuildConsensus() function.
#' Indeed, when a degenerated nucleotide is found in a sequence, it is associated to a corresponding score
#' for the observation frequency matrix.
#' Calculate Degenerescence Score

#' @param x (mandatory) an oligonucleotide or sequence (character)
#' @return The 'fraction of real nucleotides' related to the given degenerated nucleotide.
#' @export
# When a degenerate base is found, apply its occurrence to the corresponding bases (A, T, C, G)
DegToNucl <- function(x){
  for (n in c("W","S","R","Y","K","M","B","D","H","V","N")){
    if (sum(x[n,])>0){
      # deg 2:
      if (n=="W"){
        x["A",which(x["W",]>0)] <- x["A",which(x["W",]>0)] + x["W",which(x["W",]>0)]/2
        x["T",which(x["W",]>0)] <- x["T",which(x["W",]>0)] + x["W",which(x["W",]>0)]/2
        x["W",which(x["W",]>0)] <- 0
      }
      if (n=="S"){
        x["C",which(x["S",]>0)] <- x["C",which(x["S",]>0)] + x["S",which(x["S",]>0)]/2
        x["G",which(x["S",]>0)] <- x["G",which(x["S",]>0)] + x["S",which(x["S",]>0)]/2
        x["S",which(x["S",]>0)] <- 0
      }
      if (n=="R"){
        x["A",which(x["R",]>0)] <- x["A",which(x["R",]>0)] + x["R",which(x["R",]>0)]/2
        x["G",which(x["R",]>0)] <- x["G",which(x["R",]>0)] + x["R",which(x["R",]>0)]/2
        x["R",which(x["R",]>0)] <- 0
      }
      if (n=="Y"){
        x["C",which(x["Y",]>0)] <- x["C",which(x["Y",]>0)] + x["Y",which(x["Y",]>0)]/2
        x["T",which(x["Y",]>0)] <- x["T",which(x["Y",]>0)] + x["Y",which(x["Y",]>0)]/2
        x["Y",which(x["Y",]>0)] <- 0
      }
      if (n=="K"){
        x["G",which(x["K",]>0)] <- x["G",which(x["K",]>0)] + x["K",which(x["K",]>0)]/2
        x["T",which(x["K",]>0)] <- x["T",which(x["K",]>0)] + x["K",which(x["K",]>0)]/2
        x["K",which(x["K",]>0)] <- 0
      }
      if (n=="M"){
        x["A",which(x["M",]>0)] <- x["A",which(x["M",]>0)] + x["M",which(x["M",]>0)]/2
        x["C",which(x["M",]>0)] <- x["C",which(x["M",]>0)] + x["M",which(x["M",]>0)]/2
        x["M",which(x["M",]>0)] <- 0
      }
      # deg 3:
      if (n=="B"){
        x["C",which(x["B",]>0)] <- x["C",which(x["B",]>0)] + x["B",which(x["B",]>0)]/3
        x["T",which(x["B",]>0)] <- x["T",which(x["B",]>0)] + x["B",which(x["B",]>0)]/3
        x["G",which(x["B",]>0)] <- x["G",which(x["B",]>0)] + x["B",which(x["B",]>0)]/3
        x["B",which(x["B",]>0)] <- 0
      }
      if (n=="D"){
        x["A",which(x["D",]>0)] <- x["A",which(x["D",]>0)] + x["D",which(x["D",]>0)]/3
        x["T",which(x["D",]>0)] <- x["T",which(x["D",]>0)] + x["D",which(x["D",]>0)]/3
        x["G",which(x["D",]>0)] <- x["G",which(x["D",]>0)] + x["D",which(x["D",]>0)]/3
        x["D",which(x["D",]>0)] <- 0
      }
      if (n=="H"){
        x["C",which(x["H",]>0)] <- x["C",which(x["H",]>0)] + x["H",which(x["H",]>0)]/3
        x["T",which(x["H",]>0)] <- x["T",which(x["H",]>0)] + x["H",which(x["H",]>0)]/3
        x["A",which(x["H",]>0)] <- x["A",which(x["H",]>0)] + x["H",which(x["H",]>0)]/3
        x["H",which(x["H",]>0)] <- 0
      }
      if (n=="V"){
        x["C",which(x["V",]>0)] <- x["C",which(x["V",]>0)] + x["V",which(x["V",]>0)]/3
        x["A",which(x["V",]>0)] <- x["A",which(x["V",]>0)] + x["V",which(x["V",]>0)]/3
        x["G",which(x["V",]>0)] <- x["G",which(x["V",]>0)] + x["V",which(x["V",]>0)]/3
        x["V",which(x["V",]>0)] <- 0
      }
      # deg 4:
      if (n=="N"){
        x["T",which(x["N",]>0)] <- x["T",which(x["N",]>0)] + x["N",which(x["N",]>0)]/4
        x["C",which(x["N",]>0)] <- x["C",which(x["N",]>0)] + x["N",which(x["N",]>0)]/4
        x["A",which(x["N",]>0)] <- x["A",which(x["N",]>0)] + x["N",which(x["N",]>0)]/4
        x["G",which(x["N",]>0)] <- x["G",which(x["N",]>0)] + x["N",which(x["N",]>0)]/4
        x["N",which(x["N",]>0)] <- 0
      }
    }
  }
  return(x)
}
