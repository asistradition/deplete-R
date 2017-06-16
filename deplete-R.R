library(ShortRead)
library(seqLogo)

#Set these

WORKING.DIR <- "~"
FASTQ.PATH <- "."

MIN.QUALITY <- 25
MIN.READS <- 50

INVARIANT.SEED <- "gcaccaccgtcgccggcatt"
LENGTH.OF.DEGENERATE <- 5

NEG.FILE.NAMES <- c("FILE1",
                    "FILE2",
                    "FILE3")

EXP.FILE.NAMES <- c("FILE1",
                    "FILE2",
                    "FILE3")

setwd(WORKING.DIR)

pwmFromFastq <- function(negFileNames = NEG.FILE.NAMES, 
                         expFileNames = EXP.FILE.NAMES, 
                         filePath = FASTQ.PATH,
                         minQuality = MIN.QUALITY,
                         minReads = MIN.READS,
                         seedSequence = INVARIANT.SEED,
                         degenLen = LENGTH.OF.DEGENERATE) {
  
  # This is intented to take a set of FastQ files in negative
  # and in experimental data sets, find a degnerate loci defined by
  # an invariant seed sequence, and then build a position weight
  # matrix for that loci based on the differences between the control
  # and experimental data.
  # 
  # Basically it's a quick and dirty way to get Sequence Logos for
  # SELEX or plasmid depletion sequencing library experiments
  #
  # Args:
  #  negFileNames: A vector of strings with file names
  #  expFileNames: A vector of strings with file names
  #  filePath: A string with the path to the files in the first two arguments
  #  minQuality: An integer with the minimum PHRED score to filter reads
  #  minReads: An integer with the minimum number of reads required for any
  #            one specific degnerate sequence discovered
  #  seedSequence: A string with the invariant sequence to look at the 3'
  #                end of to find the degnerate loci
  #  degenLen: An integer with the length of the degenerate loci in bases
  #
  # Returns:
  #  A list of two  position-weight matrix (pwm) objects representing the
  #  PWMs of sequences enriched in the experimental over control, and depleted
  #  in the experimental over control
  
  source(file.path("lib","fileLoad.R"))
  source(file.path("lib","getDegenerateLoci.R"))
  
  #Build the workspace variables needed
  negFilePaths <- pathifyList(negFileNames, filePath)
  expFilePaths <- pathifyList(expFileNames, filePath)
  
  sequenceRegex <- paste(
    toupper(seedSequence),
    "[A-Z]{", toString(degenLen),
    "}", 
    sep="")
  
  #Basic QC
  
  qualityCheck <- qa(c(negFilePaths, expFilePaths))
  
  print(qualityCheck[["readCounts"]])
  print(qualityCheck[["baseCalls"]])
  
  #Load data and use grep to filter reads without the invariant requirement sequence
  
  negData <- sapply(X = negFilePaths, FUN = loadAndTrim, 
                    minQual = minQuality)
  expData <- sapply(X = expFilePaths, FUN = loadAndTrim, 
                    minQual = minQuality)
  
  negCounts <- sapply(X = negData, FUN = length)
  expCounts <- sapply(X = expData, FUN = length)
  
  negData <- sapply(X = negData, FUN = parseSRByRegex, 
                    regexExp = sequenceRegex)
  
  expData <- sapply(X = expData, FUN = parseSRByRegex, 
                    regexExp = sequenceRegex)
  
  #Narrow to the degenerate sequence of interest
  
  negData <- sapply(X = negData, FUN = getDegenerateLoci, 
                    searchSeq = seedSequence, degenLen = degenLen)
  
  expData <- sapply(X = expData, FUN = getDegenerateLoci, 
                    searchSeq = seedSequence, degenLen = degenLen)
  
  #Build frequency table
  
  
  negTables <- normalizeTable(
    sapply(X = sapply(X = negData, FUN = sread), FUN = table),
    minCount = minReads)
  expTables <- normalizeTable(
    sapply(X = sapply(X = expData, FUN = sread), FUN = table), 
    minCount = minReads)
  
  uniques <- unique(names(unlist(unname(c(negTables, expTables)))))
  
  negDF <- mergeDF(negTables, uniques)
  expDF <- mergeDF(expTables, uniques)
  
  diffs <- rowMeans(expDF, na.rm=TRUE)-rowMeans(negDF, na.rm=TRUE)
  
  return(list(makePWM(buildPWM(diffs[diffs > 0])),
              makePWM(buildPWM(diffs[diffs < 0]))))
}

pwms <- pwmFromFastq()

seqLogo(pwms[[1]])
seqLogo(pwms[[2]])
