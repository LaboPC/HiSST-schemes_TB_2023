setwd("File_path_to/sdaA")
# read the fasta file
fasta_file <- readLines("sdaA_pseudo-db.fasta")

# initialize variables
current_seq <- ""
current_length <- 0
output <- ""

# iterate over each line in the file
for(line in fasta_file) {
  if(startsWith(line, ">")) { # a new sequence has started
    if(current_length >= 300) { # if the previous sequence had at least 300 nucleotides
      output <- paste(output, current_seq, sep = "\n")
    }
    current_seq <- line # store the new sequence header
    current_length <- 0 # reset the length
  } else {
    current_length <- current_length + nchar(line) # increment the length
    current_seq <- paste(current_seq, line, sep = "\n") # add the line to the sequence
  }
}

# check the last sequence
if(current_length >= 300) {
  output <- paste(output, current_seq, sep = "\n")
}

# write the output to a new file
writeLines(output, "sdaA_pseudo-db.fasta")

