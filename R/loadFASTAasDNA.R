#' FASTA loader
#'
#' Load a FASTA file and convert it to a DNA list (DNAStringSet object from Biostrings package).
#' @param fasta The path of the FASTA file to load in R.
#' @return A DNAStringSet object (list of DNA sequences).
#' @importFrom Biostrings readBStringSet DNAStringSet
#' @export
#' @examples
#' \dontrun{
#' JellyDNA <- loadFASTAasDNA("training_files/Cubozoa_and_Scyphozoa.fasta")
#' }

loadFASTAasDNA <- function(fasta){
  # Load a FASTA file and convert it to a DNAStringSet object.
  # (error will be produced here if the file is not a FASTA as excepted)
  xDNA <- Biostrings::readBStringSet(fasta)

  # convert Uracile to Thymine
  xDNA <- Biostrings::DNAStringSet(gsub("U","T",as.character(xDNA)))

  # return object and its description
  print(paste0("FASTA file correctly read and convert as a DNA sequences list"))
  print(summary(xDNA))
  return(xDNA)
}
