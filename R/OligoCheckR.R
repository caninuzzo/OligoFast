#' Check and keep similar oligonucleotides
#'
#' This function allows to look for the most abundant similar oligonucleotides from each consensus to keep
#' the most redundants which are above a threshold defined by the user.
#' @param OligoFound (mandatory) it should be the output from OligoFindR() function, the dataframe with the oligonucleotides and their characteristics.
#' @param filtGaps (mandatory, default: TRUE) if TRUE, the oligonucleotides with gaps are removed from the analysis.
#' @param pOcc (mandatory, default: 0.9) a threshold from 0 to 1. It is a percentage related to the number of clusters. The theorical maximum 1 means 100% so the oligonucleotides should appears in all the clusters to be kept.
#' with default setting, the oligonucleotides that appear in more than 90% of the clusters are kept. Naturally, the higher is pOcc, the stronger should be the oligonucleotide.
#' @param maxDeg (mandatory, default: 1) the maximum number of degenerated nucleotides allowed within the oligonucleotides.
#' As a reminder, the 'degeneration degree' of the nucleotides is dependent to the threshold fixed in the previous function for the argument of the same name.
#' @param plot (optional, default: TRUE) if TRUE, creates an overview of the most abundant oligonucleotides and their location on a graph, using ggplot2.
#' The picture is a new variable created named "oligOverview" (print(oligOverview) to display it).
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @return Another dataframe of the remaining and most abundant oligonucleotides.
#' @import dplyr
#' @import ggplot2
#' @importFrom pwalign stringDist
#' @importFrom Biostrings BStringSet
#' @export
#' @examples
#' \dontrun{
#' JellyOligoChecked <- OligoCheckR(JellyOligoFound, filtGaps = TRUE, nOcc = 56, plot = TRUE)
#' }
#' @note
#' Example made using default settings
#' To save the dataframe (as .Rdata) in a given folder (here "JellyOutputs")
#' add the following argument: outpath = "path_to/JellyOutputs"
OligoCheckR <- function(OligoFound, filtGaps = TRUE, pOcc = 0.9, maxDeg = 1, plot = TRUE, outpath) {

  ### Arguments check and formating:
  # rise exception and stop function if OligoFound is not the dataframe from OligoFindR function
  if (class(OligoFound$oligos)!="character" && class(OligoFound)!="dataframe") {
    stop("The object provided for 'OligoFound' argument must be the dataframe output from 'OligoFindR' function")
  }
  # from user maxDeg entry, adapt for maxDeg cluster threshold (nb bases different)
  if (maxDeg>4){
    stop("The value of 'maxDeg' is too high for optimal primers research. Please enter a value lower than 5.")
  }

  # store nb of clusters from previous fct:
  nClusts <- length(levels(as.factor(OligoFound$subalignment)))
  # start time counter
  t1 <- Sys.time()

  # Display function starting to process:
  cat("Processing... \n Merging occurrences of oligonucleotides and calculating weak bases score...")
  # remove the ones with gaps ('-') [optional but recommended]
  if (isTRUE(filtGaps)) {
    OligoFound <- OligoFound[!grep("-",OligoFound$oligos),]
  }

  ### lowercase characters treatment > merging > warning score > merging > adding Freq counts
  # look for similar patterns from different subalignments (keep the ones over a threshold to define)
  occOligo <- data.frame(table(OligoFound$oligos))
  names(occOligo)[1]='oligos'
  occOligo$wbs <- vector(mode = "character", length = nrow(occOligo))

  # separate oligos which include lowercase char. & oligos with full uppercase char.
  oligoLow <- OligoFound[grep("[a-z]",OligoFound$oligos),]
  oligoUP <- OligoFound[!grep("[a-z]",OligoFound$oligos),]


  # if lowercases characters (FALSE when ignGap=T in consensus step)
  # calculate occurrences of each oligos and introduce wbs score
  if (nrow(oligoLow)>0) {
    occOligoLow <- data.frame(table(oligoLow$oligos))
    names(occOligoLow)[1] <- "oligos"
    occOligoUP <- data.frame(table(oligoUP$oligos))
    names(occOligoUP)[1] <- "oligos"
    occOligoUP$wbs <- vector(mode = "character", length = nrow(occOligoUP))

    # keep only 1% for lowercases dataframe (less data to process, remove 'small events')
    occOligoUP <- occOligoUP %>% filter(Freq>round(0.01*nClusts))

  #if (nrow(oligoLow)>0){
    occOligoLow <- occOligoLow %>% filter(Freq>round(0.01*nClusts))
    # cumul Freq when lowercases and uppercases oligos match:
    for (i in occOligoLow$oligos){
      # do it only when occurrence is detected (less data to process)
      if (length(which(toupper(i) == occOligoUP$oligos))>0){
        cumFreq <- occOligoUP$Freq[which(toupper(i) == occOligoUP$oligos)] + occOligoLow[occOligoLow$oligos==i,"Freq"]
        occOligoUP$Freq[which(toupper(i) == occOligoUP$oligos)] <- cumFreq
        occOligoUP$wbs[which(toupper(i) == occOligoUP$oligos)] <- (nchar(gsub("[^a-z]", "", i))*occOligoLow[occOligoLow$oligos==i,"Freq"])/(nchar(as.character(occOligoUP$oligos[which(toupper(i)==occOligoUP$oligos)]))*cumFreq)
      }
    }
  } else {
    occOligoUP <- data.frame(table(oligoUP$oligos))
    names(occOligoUP)[1] <- "oligos"
    occOligoUP$wbs <- vector(mode = "character", length = nrow(occOligoUP))
    # keep only 1% for lowercases dataframe (less data to process, remove 'small events')
    occOligoUP <- occOligoUP %>% filter(Freq>round(0.01*nClusts))
    rm(oligoLow)
  }

  # preparing field
  occOligoUP$length <- apply(occOligoUP["oligos"],1,nchar)
  colnames(occOligoUP)[colnames(occOligoUP) == "Freq"] <- "cumFreq"

  ### if no degenerenated nucleotides are allowed then skip this process
  if (maxDeg>0) {
    # uppercases oligos only at this point - Display following of the process:
    cat("\n Summing the occurrences of similar oligonucleotides according to maxDeg threshold...\n")

    # clustering by oligos size to apply the right method with stringDist()
    # obsFreq <- old Freq
    occOligoUP$obsFreq <- occOligoUP$cumFreq

    OligoStorage <- list()
    OligoStored <- list()
    e <- 0

    # loop by size (from min to max oligo size <-> s) to get dissimilarity matrix with hamming method with maxDeg threshold:
    for (s in min(occOligoUP$length):max(occOligoUP$length)) {
      e <- e+1
      OligoStorage[[e]] <- occOligoUP %>% filter(length==s)
      disttest_fmr <- pwalign::stringDist(Biostrings::BStringSet(OligoStorage[[e]]$oligos), method = "hamming")
      gc()
      couples_fmr <- which(as.matrix(disttest_fmr)<=maxDeg & as.matrix(disttest_fmr)>0, arr.ind = TRUE)
      gc()
      rm(disttest_fmr)
      uniq_pairs_fmr <- couples_fmr[couples_fmr[,1] < couples_fmr[,2],]
      rm(couples_fmr)

      # cumul the freq of the oligo clustered together
      stor_list <- list()
      for (i in as.numeric(levels(as.factor(uniq_pairs_fmr)))){
        OligoStorage[[e]]$cumFreq[i] <- OligoStorage[[e]]$cumFreq[i] + sum(OligoStorage[[e]]$cumFreq[as.numeric(uniq_pairs_fmr[which(uniq_pairs_fmr[,1]==i),2])])
        if (length(as.numeric(uniq_pairs_fmr[which(uniq_pairs_fmr[,1]==i),2]))>0) {
          stor_list[[as.character(OligoStorage[[e]]$oligos[i])]] <- c(as.character(OligoStorage[[e]]$oligos[i]),as.character(OligoStorage[[e]]$oligos[as.numeric(uniq_pairs_fmr[which(uniq_pairs_fmr[,1]==i),2])]))
        }
      }
      # keep the most abundant oligos
      OligoStorage_fmr <- OligoStorage[[e]] %>% filter(cumFreq>=round(pOcc*nClusts))
      if (nrow(OligoStorage_fmr)>0) {
        stor_list <- stor_list[names(stor_list) %in% as.character(OligoStorage_fmr$oligos)]
        # in case of similar oligos without a consensus already established, creates it respecting maxDeg
        for (i in names(stor_list)){
          aj <- OligoStorage[[e]][OligoStorage[[e]]$oligos %in% stor_list[[i]],]
          aj$cumFreq <- max(OligoStorage[[e]]$cumFreq[which(OligoStorage[[e]]$oligos %in% stor_list[[i]])])
          suppressMessages(OligoStorage_fmr <- full_join(aj,OligoStorage_fmr))
        }
        OligoStored[[e]] <- OligoStorage_fmr %>% distinct(by=oligos, .keep_all = TRUE)
      } else { # end if nrow
        OligoStored[[e]] <- OligoStorage_fmr
      }
      gc()
    }

    # merge each df of the list to a single and big df:
    OligoStored <- Reduce(function(x, y) merge(x, y, all = TRUE), OligoStored)

    # if no primers found:
    if (nrow(OligoStored)==0){
      stop("Sorry no primers found with the defined parameters.
           Consider modifying them and retry. You can even modify some parameters from previous steps of the pipeline")
    }

    # remove "by" column
    OligoStored <- OligoStored[,1:5]

    ## get all the fields from OligoFound:
    OligoFound <- OligoStored %>%
      merge(OligoFound, by='oligos') %>%
      distinct()
    ## calculate mean start & end positions
    OligoChecked <- OligoFound %>%
      group_by(oligos) %>%
      summarise(meanStartPos=round(mean(startPos)),
                meanEndPos=round(mean(endPos)))

    # to remove useless fields and keep minimum oligos for next fct
    OligoFound <- OligoFound %>% select(-c("subalignment","startPos","endPos"))

    # to return the final objects: oligos with all useful data
    OligoChecked <- OligoChecked %>%
      merge(OligoFound, by = 'oligos') %>%
      distinct()

    ### if plot == TRUE, generate a ggplot:
    if (isTRUE(plot)) {

      # simple plot : bars reporting occurrences in fonction of the position along the sequences:
      oligOverview <<- ggplot(OligoChecked) +
        geom_rect(aes(xmin = meanStartPos, xmax = meanStartPos+length, ymin = 0, ymax = cumFreq), color = "black", alpha = 0.7) +
        scale_x_continuous(breaks = seq(from = 0, to = ceiling(max(OligoChecked$meanEndPos)/100)*100, by = 50),
                           limits = c(floor(min(OligoChecked$meanStartPos)/100)*100,ceiling(max(OligoChecked$meanEndPos)/100)*100+50)) +
        ylim(0,nClusts) +
        labs(title = "Potential oligonucleotides positions & occurrences along the sequences", x = "Estimated position (bp)", y = "Occurrences") +
        theme(legend.position = "none") +
        geom_hline(yintercept = round(pOcc*nClusts), linetype = "dashed", color = "red") +
        geom_hline(yintercept = nClusts, linetype = "solid", color = "green")
    } # end if plot==TRUE

  } else { # end if maxDeg>0
    # just do as it was done in v.0.0.1
    OligoStored <- occOligoUP %>% filter(cumFreq>=round(pOcc*nClusts))
    print(nrow(OligoStored))
    # if no primers found:
    if (nrow(OligoStored)==0) {
      stop("Sorry no primers found with the defined parameters.
           Consider modifying them and retry. You can even modify some parameters from previous steps of the pipeline")
    }

    ## get all the fields from OligoFound:
    OligoFound <- OligoStored %>%
      merge(OligoFound, by='oligos') %>%
      distinct()
    ## calculate mean start & end positions
    OligoChecked <- OligoFound %>%
      group_by(oligos) %>%
      summarise(meanStartPos=round(mean(startPos)),
                meanEndPos=round(mean(endPos)))

    # to remove useless fields and keep minimum oligos for next fct
    OligoFound <- OligoFound %>% select(-c("subalignment","startPos","endPos"))

    # to return the final objects: oligos with all useful data
    OligoChecked <- OligoChecked %>%
      merge(OligoFound, by = 'oligos') %>%
      distinct()


    ### plot :
    if (isTRUE(plot)) {
      oligOverview <<- ggplot(OligoChecked) +
        geom_rect(aes(xmin = meanStartPos, xmax = meanStartPos+length, ymin = 0, ymax = cumFreq), color = "black", alpha = 0.7) +
        scale_x_continuous(breaks = seq(from = 0, to = ceiling(max(OligoChecked$meanEndPos)/100)*100, by = 50),
                           limits = c(floor(min(OligoChecked$meanStartPos)/100)*100,ceiling(max(OligoChecked$meanEndPos)/100)*100+50)) +
        ylim(0,nClusts) +
        labs(title = "Potential oligonucleotides positions & occurrences along the sequences", x = "Estimated position (bp)", y = "Occurrences") +
        theme(legend.position = "none") +
        geom_hline(yintercept = round(pOcc*nClusts), linetype = "dashed", color = "red") +
        geom_hline(yintercept = nClusts, linetype = "solid", color = "green")
    } # end if plot==TRUE

  } # end else maxDeg==0

  # -> oupath: save final objects in a folder ?
  if (!missing(outpath)) {
    save(OligoChecked, file = paste0(outpath,"/OligoChecked_out_pocc",pOcc,"_maxDeg",maxDeg,".Rdata"))
  }

  # end time counter
  t2 <- Sys.time()
  cat(paste0("\r Total function execution time: ",
             as.numeric(round(difftime(t2, t1, units = "mins"),2)),
             " minutes"))

  return(OligoChecked)
}
