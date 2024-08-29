#' Shuffle the sequences, divide them into clusters and align them
#'
#' This function enables to create random clusters of sequences from the dataset and
#' align the sequences within each cluster according to the parameters given by the user.
#' To increase the robustness of the future results, it is also possible to iterate the creation of random cluster of sequences.
#' @param DNAList (mandatory) should be a DNAStringSet object (from loadFASTAasDNA function).
#' @param nseq (mandatory, set to 10 by default) the number of sequences in each cluster.
#' @param ite (mandatory, set to 10 by default) the number of iterations of the shuffling process.
#' Set 1 for no iteration, higher is the iteration number, stronger will be the results but higher will be the processing time.
#' @param outpath (optional) the path to the FOLDER in which the function outputs will be stored.
#' @return A list of clusters of aligned sequences (each cluster are MsaDNAMulipleAlignment object)
#' @importFrom Biostrings writeXStringSet DNAStringSet
#' @importFrom msa msa
#' @export
#' @note This function may take some time to process, especially for long sequences, big dataset or high number of iterations.
#' @examples
#' \dontrun{
#' JellyAlignments <- ShuffleAndAlign(JellyDNA, nseq = 10, ite = 10)
#' }
#' @note
#' Example made using default settings
#' To save the alignments list (as .Rdata) in a given folder (here "JellyOutputs")
#' add the following argument: outpath = "path_to/JellyOutputs"

ShuffleAndAlign <- function(DNAList, nseq = 10, ite = 10, outpath) {

  # rise exception and stop function if DNAList is not a DNAStringSet object
  if (class(DNAList)!="DNAStringSet") {
    stop("The object provided as DNA sequences must be a 'DNAStringSet' object!")
  }

  # start time counter to estimate function progression
  t1<- Sys.time()

  # WARNING ! If some sequences are similar:
  if (length(which(duplicated(DNAList)==TRUE)>0)){
    cat(paste0("WARNING! It seems that ", length(which(duplicated(DNAList)==TRUE)),
               " of the sequences provided are similar. Press Enter to continue or any other key to quit: "))

    # Pause the function and wait for user input
    user_input <- readline(prompt = "")

    # Check user input
    if (user_input == "") {
      cat("Continuing with the function...\n")
    } else {
      stop("Exiting the function...")
    }
  }

  # -> Shuffling step and iterations
  cat(paste0("Shuffling the sequences ", ite, " time.s. \n"))
  seq_shuffled <- replicate(ite,sample(DNAList, size = length(DNAList),replace = F))
  seq_shuffled <- do.call(c, seq_shuffled)

  # define the number of sub-alignments to process (K) from seq
  K <- floor(length(seq_shuffled)/nseq)
  time_storage <- NULL

  # -> Clustering and alignments of each clusters
  cat(paste0("Clustering: creating ",
             K,
             " clusters of ",
             nseq,
             " sequences from the ",
             length(seq_shuffled),
             " total shuffled sequences. \n"))

  # create storage list
  subAli <- list()

  # start loop to align the sequences from each clusters
  cat("Starting multiple sequences alignment for each clusters... \n")
  for (i in 1:K) {

    ti <- Sys.time()
    sink("fmre.txt")
    subAli[[i]] <- msa::msa(seq_shuffled[((nseq*i)-(nseq-1)):(nseq*i)], type = "dna",verbose = F)
    sink()
    file.remove("fmre.txt")
    tf <- Sys.time()

    #### Displaying time estimation with progressive loading bar ####
    while (i<K) {
      time_storage[i] <- as.numeric(round(difftime(tf, ti, units = "mins"),2))
      # bar
      fillBar <- round(i/K*60,0)
      # display
      cat("\r [",
          paste0(rep("#",fillBar),collapse = ""),
          paste0(rep("_",(60-fillBar)),collapse = ""),
          "] Aligning sequences from cluster ", i, "/", K,". Remaining time estimation: ~",
          round(mean(time_storage)*(K-i),2),
          "min.     ")
      flush.console()
      break
    }
    if (i==K) {
      cat("\r [",
          paste0(rep("#",61),collapse=""),
          "] Job finished.",
          rep(" ",30),
          "\n")
      flush.console()
    }
    #### end loading bar code ####

    # save sub-alignment as .fasta if desired
    if (!missing(outpath)) {
      Biostrings::writeXStringSet(Biostrings::DNAStringSet(subAli[[i]]),
                                  paste0(outpath,"/subalignment_cluster_",i,".fasta"),
                                  format = "fasta")
    }
  } # end i loop

  # save output as .RData file if asked
  if (!missing(outpath)) {
    save(subAli, file = paste0(outpath,"/ShuffleAndAlign_out_nseq",nseq,"ite",ite,".Rdata"))

    cat(paste0("The complete list of the aligned sequences clusters and each of them have been saved at: ",
               deparse(substitute(outpath)),
               " \n"))
  }

  # time counter (final)
  t2 <- Sys.time()
  cat(paste0("Total function execution time: ",
             as.numeric(round(difftime(t2, t1, units = "mins"),2)),
             " minutes"))

  return(subAli)
}
