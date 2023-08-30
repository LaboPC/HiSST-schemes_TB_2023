setwd("~/DATABASE/S_maltophilia/LOCUS/") # Change to your working directory file
LOCUS = "LOCUS" # Change "LOCUS" to the corresponding HiSST locus

# Create a data frame with the original sequence names and the corresponding new names
library(Biostrings)
Seq.txt <- readDNAStringSet(paste0(LOCUS,"_blast.txt"))

# Edit strain names
Seq.Name <- names(Seq.txt)
new_names <- gsub("strain ","",Seq.Name) # Delete 'strain' characters
new_names <- gsub(".*Steno","Steno",new_names) # Delete characters before 'Steno'

pattern.rm <- c(",.*", " chromosome.*"," genome.*"," complete.*"," DNA.*") # Delete characters after some characters
new_names. <- gsub(paste(pattern.rm,collapse="|"),"",new_names)

# Create a data frame with the sequence names and DNA sequences
Sequences <- data.frame(Names = new_names., Seq = paste(as.character(Seq.txt)))

# Write the data to the output file in FASTA format
writeLines(paste0(">", Sequences$Names, "\n", Sequences$Seq), con = paste0(LOCUS,".fasta"))
