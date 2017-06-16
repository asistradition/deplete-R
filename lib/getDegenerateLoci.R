library(ShortRead)

parseSRByRegex <- function(sro, regexExp) {
  # GREP a ShortRead object to find a specific regular expression
  # Both forward and reverse
  #
  # Args:
  #  sro: ShortRead object
  #  regexExp: A string containing the regular expression to grep
  #
  # Returns:
  #  A ShortRead object containing sequences which match the grep hit
  #  Adjusted to have the grep hit sense
  #
  # Notes:
  #  This is not safe against sequence result duplication in the
  #  event the grep hit is present in both F and R directions
  
  sro <- append(sro[grep(regexExp, sread(sro))], reverseComplement(sro)[grep(regexExp, sread(reverseComplement(sro)))])
  return(sro)
}

getDegenerateLoci <- function (sro, searchSeq, degenLen) {
  # Narrow a ShortRead object of sequencs to a specific degnerate loci
  # based on the presence of a seed sequence
  #
  # Args:
  #  sro: A ShortRead object containing sequences
  #  searchSeq: A string containing the invariant seed sequence to grep for
  #  degenLen: An integer with the length of the degnerate loci
  #
  # Returns:
  #  A ShortRead object narrowed to the specific degnerate loci of interest
  
  searchSeq <- toupper(searchSeq)
  regexF <- regexpr(searchSeq, sread(sro))
  
  startPad <- nchar(searchSeq)
  endPad <- startPad + degenLen - 1
  
  return(narrow(sro, start = (regexF+startPad), end=(regexF+endPad)))
}

normalizeTable <- function (table, minCount, totalCounts = NULL) {
  for (i in 1:length(table)) {
    table[[i]] <- table[[i]][table[[i]] >= minCount]
    if (is.null(totalCounts)) {
      table[[i]] <- table[[i]] / sum(table[[i]]) * 10^6
    }
    else {
      table[[i]] <- table[[i]] / totalCounts * 10^6
    }
  }
  return(table)
}

mergeDF <- function (list, rowNames, by = 0, all = TRUE, discardNA = TRUE) {
  
  # Merges data.frames row-wise and cleans up the resulting new data.frame
  #
  # Args:
  #  list: list of data.frame objects to merge together
  #  rowNames: A vector of the row names to use to start the initial data.frame
  #  by: 0 for row-wise and 1 for column-wise (passed to merge function)
  #  all: boolean for complete join (passed to merge function)
  #  discardNA: boolean for replacing all NA values with 0s
  #
  # Returns:
  #  A data.frame merged from the list of data.frames passed in
  
  merged <- data.frame(row.names = rowNames)
  for (i in 1:length(list)) {
    merged <- merge(merged, list[[i]], by = by, all = all)
    if (by == 0) {
      merged <- transform(merged, row.names = Row.names, Row.names = NULL)
    }
  }
  if (discardNA) {
    merged[is.na(merged)] <- 0
  }
  return(merged)
}

buildPWM <- function (seqs, chars = c("A", "T", "G", "C")) {
  
  # Calculates the position weight matrix of a named vector of sequences
  # using the read (normalized or raw) count.
  #
  # Args:
  #  seqs: Named vector with sequences as names and read count as numeric
  #  chars: Vector of sequence characters
  #
  # Returns:
  #  A data.frame of rows that are chars and columns that are positions,
  #  with numeric values for frequency of the char at that position
  #
  # Notes:
  #  chars need to cover all possible sequence characters that will be
  #  encountered or some downstream applications will fail because 
  #  frequencies won't total to 1
  #
  #  Sequences must all be the same length for the same reason as above
  
  pwm <- data.frame(row.names = chars)
  seqstrs <- names(seqs)
  totalValue <- sum(seqs)
  for (i in 1:nchar(seqstrs[1])) {
    activechars <- substring(seqstrs, i, i)
    activePWM <- data.frame(rep(0,length(chars)), row.names = chars)
    colnames(activePWM) <- paste("S", i, sep="")
    for (ach in chars) {
      activePWM[ach, 1] = sum(seqs[which(ach == activechars)])/totalValue
    }
    pwm <- mergeDF(list(pwm, activePWM), chars)
  }
  return (pwm)
}
