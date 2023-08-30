
setwd("~/Example_Create_Paeruginosa_database/Step1_db_Pseudo/btuB/") # Change to your working directory file
LOCUS = "butB" # Change "LOCUS" to the corresponding HiSST locus

# Create a data frame with the original sequence names and the corresponding new names
Seq.xls <- read.delim(paste0("blastn_",LOCUS,"_pseudo-db.txt")) 
Seq.Name <- Seq.xls$Replicon.Name..Subject.
NewName <- sapply(Seq.Name, function(x) gsub(" \\(.+\\)", "",x))
NewName <- sapply(NewName, function(x) gsub("isolate ", "",x))
Sequences <- data.frame(Names = NewName, Seq = Seq.xls$Aligned.part.of.subject.sequence)

# Write the data to the output file in FASTA format
writeLines(paste0(">", Sequences$Names, "\n", Sequences$Seq), con = paste0(LOCUS,".fasta"))


# change the sequence names
#New_seq_name <- sapply(New_seq_name, function(x) sub(" - .*", "", x))

