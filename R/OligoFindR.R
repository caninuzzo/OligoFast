#' Find and extract oligonucleotides
#'
#' This function allows to look for conserved oligonucleotides in each consensus sequences according to
#' the parameters given by the user (see args).
#' @param ConsensusClust (mandatory) it should be the output from BuildConsensus() function, a list of consensus sequences.
#' @param PriMin (mandatory, set to 18 by default) the minimum size (in bp) expected for the primers.
#' @param PriMax (mandatory, set to 22 by default) the maximum size (in bp) expected for the primers.
#' @param maxDeg (mandatory, set to 2 by default) the maximum of degenerescence allowed for the primers.
#' (The degenerescence score of W -which is A or T- is 2 ; the one of B -which can be T or C or G- is 3 ; N is 4 etc.).
#' @param GCrange (mandatory, set to c(40,60) by default) the range (in %) of %GC content within the oligo.
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @return A dataframe with the oligonucleotides fulfilling the given parameters and kept from each clusters.
#' @import data.table
#' @export
#' @examples
#' \dontrun{
#' JellyOligoFound <- OligoFindR(JellyConsensus, PriMin=18, PriMax=22, maxDeg=0, GCrange=c(40,60))
#' }
#' @note
#' Example made using default settings
#' To save the dataframe (as .Rdata) in a given folder (here "JellyOutputs")
#' add the following argument: outpath = "path_to/JellyOutputs"

OligoFindR <- function(ConsensusClust, PriMin=18, PriMax=22, maxDeg=2, GCrange=c(40,60), outpath) {

  # rise exception and stop function if subAli is not the list of msa DNA Alignment from ShuffleAndAlign function
  if (class(ConsensusClust[[1]])!="character" && class(ConsensusClust)!="list") {
    stop("The object provided for 'ConsensusCust' argument must be the complete list output from 'BuildConsensus' function")
  }

  # start time counter
  t1 <- Sys.time()

  # set variables to process function
  PriMin <- PriMin-1
  PriMax <- PriMax-1

  OligoFound <- data.table::data.table(subalignment=numeric(),
                                       startPos=numeric(),
                                       endPos=numeric(),
                                       oligos=character(),
                                       GC=numeric(),
                                       deg=numeric(),
                                       stringsAsFactors = F,
                                       row.names = NULL)

  time_storage <- NULL

  # Increase function speed:
  # if number of consensus to analyse <100: Nothing ; if 100-1000 : /10 ; if >1000 : /100. Then merge 'fragments'
  if (length(ConsensusClust)>100){
    frag <- 10
    print("enter in loop frag=10")
    if (length(ConsensusClust)>1000){
      frag <- 100
      print("enter in loop frag=100")
    }
    int_e <- length(ConsensusClust) %/% frag
    u1 <- seq(2, length(ConsensusClust), by = int_e)
    u1[1] <- 1
    u2 <- seq(from = int_e + 1, length(ConsensusClust), by = int_e)
    if (length(u1) != length(u2)){
      u2 <- append(u2,length(ConsensusClust))
    }
    print(paste0("u1:",u1))
    print(paste0("u2:",u2))
    OligoFoundList <- list()

    for (F in 1:length(u1)){
      # start time estimation
      tini <- Sys.time()
      # start loop to look for patterns in each consensus sequences
      OligoFound <- data.table::data.table(subalignment=numeric(),
                                           startPos=numeric(),
                                           endPos=numeric(),
                                           oligos=character(),
                                           GC=numeric(),
                                           deg=numeric(),
                                           stringsAsFactors = F,
                                           row.names = NULL)

      for (e in u1[F]:u2[F]) {
        print(paste0("enter in loop interval:",u1[F],":",u2[F]))

        # first display
        if (e==1){
          cat("Processing...")
        }

        # define start and end position for oligos detection that moves along the consensus from p1 to p1+p2
        startPos <- 1
        endPos <- length(ConsensusClust[[e]])-PriMin
        for (p1 in startPos:endPos) {
          for (p2 in PriMin:PriMax) {

            # evaluation according to criteria set previously
            xdeg <- OligoFast::DegeScore(ConsensusClust[[e]][p1:(p1+p2)])
            xGC <- round(sum(grepl("C",ConsensusClust[[e]][p1:(p1+p2)],ignore.case = T) |
                               grepl("G",ConsensusClust[[e]][p1:(p1+p2)],ignore.case = T))/length(ConsensusClust[[e]][p1:(p1+p2)])*100,2)

            # if criteria are respected, let's save the candidate oligos
            if (xdeg<=maxDeg && xGC > GCrange[1] && xGC < GCrange[2]) {
              # creates new line to add to OligoFound
              newL <- list("subalignment"=as.numeric(e),
                           "startPos"=as.numeric(p1),
                           "endPos"=as.numeric(p1+p2),
                           "oligos"=as.character(paste0(ConsensusClust[[e]][p1:(p1+p2)],collapse="")),
                           "GC"=as.numeric(xGC),
                           "deg"=as.numeric(xdeg))
              OligoFound <- data.table::rbindlist(list(OligoFound,newL))
            }

            # startPosition arrive at the end of the consensus sequence, avoid to produce error
            if ((p1+p2)>=length(ConsensusClust[[e]])) {
              break
            }
          }# end p1
        }# end p2

        # time estimation:
        tfin <- Sys.time()
        while (e<length(ConsensusClust)) {
          # bar
          fillBar <- round(e/length(ConsensusClust)*60,0)
          # display
          cat("\r [",
              paste0(rep("#",fillBar),collapse = ""),
              paste0(rep("_",(60-fillBar)),collapse = ""),
              "] Searching in consensus sequence #", e, "/", length(ConsensusClust))
          flush.console()
          break
        }
        if (e==length(ConsensusClust)) {
          cat("\r [",
              paste0(rep("#",61),collapse=""),
              "] Job finished.",
              rep(" ",20),
              "\n")
          flush.console()
        }
      } # end e
      OligoFoundList[[F]] <- OligoFound
    }# end F
    OligoFound <- do.call(rbind, OligoFoundList)

  } else { # less than 100 consensus sequences to analyse

    # start loop to look for patterns in each consensus sequences
    for (e in 1:length(ConsensusClust)) {
      # start time estimation
      tini <- Sys.time()

      # first display
      if (e==1){
        cat("Processing...")
      }

      # define start and end position for oligos detection that moves along the consensus from p1 to p1+p2
      startPos <- 1
      endPos <- length(ConsensusClust[[e]])-PriMin
      for (p1 in startPos:endPos) {
        for (p2 in PriMin:PriMax) {

          # evaluation according to criteria set previously
          xdeg <- OligoFast::DegeScore(ConsensusClust[[e]][p1:(p1+p2)])
          xGC <- round(sum(grepl("C",ConsensusClust[[e]][p1:(p1+p2)],ignore.case = T) |
                             grepl("G",ConsensusClust[[e]][p1:(p1+p2)],ignore.case = T))/length(ConsensusClust[[e]][p1:(p1+p2)])*100,2)

          # if criteria are respected, let's save the candidate oligos
          if (xdeg<=maxDeg && xGC > GCrange[1] && xGC < GCrange[2]) {
            # creates new line to add to OligoFound
            newL <- list("subalignment"=as.numeric(e),
                         "startPos"=as.numeric(p1),
                         "endPos"=as.numeric(p1+p2),
                         "oligos"=as.character(paste0(ConsensusClust[[e]][p1:(p1+p2)],collapse="")),
                         "GC"=as.numeric(xGC),
                         "deg"=as.numeric(xdeg))
            OligoFound <- data.table::rbindlist(list(OligoFound,newL))
          }

          # startPosition arrive at the end of the consensus sequence, avoid to produce error
          if ((p1+p2)>=length(ConsensusClust[[e]])) {
            break
          }
        }# end p1
      }# end p2

      # time estimation:
      tfin <- Sys.time()
      while (e<length(ConsensusClust)) {
        # bar
        fillBar <- round(e/length(ConsensusClust)*60,0)
        # display
        cat("\r [",
            paste0(rep("#",fillBar),collapse = ""),
            paste0(rep("_",(60-fillBar)),collapse = ""),
            "] Searching in consensus sequence #", e, "/", length(ConsensusClust))
        flush.console()
        break
      }
      if (e==length(ConsensusClust)) {
        cat("\r [",
            paste0(rep("#",61),collapse=""),
            "] Job finished.",
            rep(" ",20),
            "\n")
        flush.console()
      }
    } # end e
  } # end else

  t2 <- Sys.time()
  cat(paste0("Total function execution time: ",
             as.numeric(round(difftime(t2, t1, units = "mins"),2)),
             " minutes"))

  # -> oupath: save final objects in a folder ?
  if (!missing(outpath)){
    save(OligoFound, file = paste0(outpath,"/OligoFound_out_min",
                                   PriMin+1,"max",PriMax+1,"deg",maxDeg,
                                   "GC",GCrange[1],GCrange[2],".Rdata"))
  }

  return(OligoFound)

}
