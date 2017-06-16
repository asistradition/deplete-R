library(ShortRead)

pathifyList <- function(listOfFiles, basePath) {
  returnList = c()
  for (name in listOfFiles) {
    returnList = c(returnList, file.path(basePath, name))
  }
  return(returnList)
}

loadAndTrim <- function(fileName, minQual = 25) {
  # Load Fastq file into a ShortRead object and Trim ends by phred score
  #
  # Args:
  #  fileName: A string containing the file path to open
  #  minQual: An integer with the minimum PHRED score to trim with
  #
  # Returns:
  #  A ShortRead object containing sequences that passed the filtering
  
  fastqSR <- readFastq(fileName)
  minPhred <- names(encoding(quality(fastqSR[1]))[minQual+1])
  fastqSR <- trimTailw(fastqSR, 2, minPhred, 2)
  
  return(fastqSR)
}
