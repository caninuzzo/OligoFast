#' Check and keep similar oligonucleotides
#'
#' This function allows to look for the most abundant similar oligonucleotides from each consensus to keep
#' the most redundants which are above a threshold defined by the user.
#' @param OligoFound (mandatory) it should be the output from OligoFindR() function, the dataframe with the oligonucleotides and their characteristics.
#' @param filtGaps (mandatory, default: TRUE) if TRUE, the oligonucleotides with gaps are removed from the analysis.
#' @param nOcc (mandatory) the number of times the oligonucleotide is found through the diffenrent clusters. The theorical maximum is the number of cluster
#' (meaning the oligonucleotide is found everywhere and well conserved). The higher is nOcc, the stronger should be the oligonucleotide.
#' @param plot (optional, default: TRUE) if TRUE, creates an overview of the most abundant oligonucleotides and their location on a graph, using ggplot2.
#' The picture is a new variable created named "oligOverview" (print(oligOverview) to display it).
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @return Another dataframe of the remaining and most abundant oligonucleotides.
#' @import dplyr
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' JellyOligoChecked <- OligoCheckR(JellyOligoFound, filtGaps = TRUE, nOcc = 56, plot = TRUE)
#' }
#' @note
#' Example made using default settings
#' To save the dataframe (as .Rdata) in a given folder (here "JellyOutputs")
#' add the following argument: outpath = "path_to/JellyOutputs"
OligoCheckR <- function(OligoFound, filtGaps = TRUE, nOcc, plot = TRUE, outpath) {

  # rise exception and stop function if OligoFound is not the dataframe from OligoFindR function
  if (class(OligoFound$oligos)!="character" && class(OligoFound)!="dataframe") {
    stop("The object provided for 'OligoFound' argument must be the dataframe output from 'OligoFindR' function")
  }

  # rise exception and stop function if nOcc is not defined
  if (missing(nOcc) || is.null(nOcc)) {
    stop("Argument 'nOcc' is missing or NULL. Please enter a numeric value.")
  }

  cat("Processing...")
  # start time counter
  t1 <- Sys.time()

  # remove the ones with gaps ('-') [optional but recommended]
  if (isTRUE(filtGaps)) {
    OligoFound <- OligoFound[-which(grepl("-",OligoFound$oligos)),]
  }

  # look for similar patterns from different subalignments (keep the ones over a threshold to define)
  occOligo <- data.frame(table(OligoFound$oligos))
  names(occOligo)[1]='oligos'
  occOligo$warnings <- vector(mode = "character", length = nrow(occOligo))

  ### DEV lowercase treatment
  oligoLow <- OligoFound[grep("[a-z]",OligoFound$oligos),]
  forafter <- NULL
  for (i in 1:nrow(oligoLow)){
    if (length(which(toupper(oligoLow$oligos[i]) == occOligo$oligos))>0){
      occOligo$Freq[which(toupper(oligoLow$oligos[i]) == occOligo$oligos)] <- occOligo$Freq[which(toupper(oligoLow$oligos[i]) == occOligo$oligos)] + 1
      occOligo$warnings[which(toupper(oligoLow$oligos[i]) == occOligo$oligos)] <- paste(occOligo$warnings[which(toupper(oligoLow$oligos[i]) == occOligo$oligos)],
                                                                                         paste0("Potential weakness: ",oligoLow$oligos[i]), sep = " ")
    forafter <- append(forafter,which(oligoLow$oligos[i] == OligoFound$oligos))
    }
  }
  ### DEV

  # primers that appears in more (or equal) to nOcc times
  occOligo <- occOligo[occOligo$Freq>=nOcc,]

  # convert to uppercase to match with oligos (warnings have already been set previously)
  if (length(forafter)>0){
    OligoFound$oligos[forafter] <- toupper(OligoFound$oligos[forafter])
  }

  OligoFound <- occOligo %>%
    merge(OligoFound, by='oligos') %>%
    distinct()

  OligoChecked <- OligoFound %>%
    group_by(oligos) %>%
    summarise(meanStartPos=mean(startPos),
              meanEndPos=mean(endPos),
              sdStartPos=sd(startPos),
              sdEndPos=sd(endPos))

  # to remove useless fields and keep minimum oligos for next fct
  OligoFound <- OligoFound %>% select(-c("subalignment","startPos","endPos"))

  OligoChecked <- OligoChecked %>%
    merge(OligoFound, by = 'oligos') %>%
    mutate(length=nchar(as.character(oligos))) %>%
    distinct()

  # returns image if the user wants it
   if (isTRUE(plot)) {
     oligOverview <<- ggplot(OligoChecked) +
       geom_point(mapping = aes(x=meanStartPos, y=Freq, color=Freq)) +
       scale_x_continuous(breaks = seq(0, max(OligoChecked$meanEndPos), by = 50)) +
       scale_y_continuous(breaks = OligoChecked$Freq, labels = OligoChecked$Freq) +
       theme(legend.position = "none")
  }

  # -> oupath: save final objects in a folder ?
  if (!missing(outpath)) {
    save(OligoChecked, file = paste0(outpath,"/OligoChecked_out_nocc",nOcc,".Rdata"))
  }

  # end time counter
  t2 <- Sys.time()
  cat(paste0("\r Total function execution time: ",
             as.numeric(round(difftime(t2, t1, units = "mins"),2)),
             " minutes"))

  return(OligoChecked)
}
