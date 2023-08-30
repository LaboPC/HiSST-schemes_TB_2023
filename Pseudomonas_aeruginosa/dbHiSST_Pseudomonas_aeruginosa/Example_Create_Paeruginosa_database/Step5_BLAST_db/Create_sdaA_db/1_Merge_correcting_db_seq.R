
correct_db_nt <- function(LOCUS_HiSST, Limit_length){
  # Create a data frame with the original sequence names and the corresponding new names
  library(Biostrings)
  sequences <- readDNAStringSet(paste0("1_NCBI_db_", LOCUS_HiSST,".fasta"))
  
  #### 1_ Remove DNA sequences under 'n' bases ####
  
  # Calculate the lengths of each sequence
  sequence_lengths <- width(sequences)
  
  # Subset the DNAStringSet to remove sequences under 'Limit_length' bases
  filtered_sequences <- sequences[sequence_lengths >= Limit_length]
  
  # Check the filtered sequences
  print(filtered_sequences)
  
  
  #### 2_ Check whether the sequences are in same orientation. If they are not, perform a reverse complement operation ####
  
  library(Biostrings)
  library(foreach)
  library(doParallel)
  
  # Get the first sequence from the filtered sequences
  first_sequence <- filtered_sequences[[1]]
  
  # Set the number of CPU cores to use for parallel computation
  num_cores <- detectCores()
  
  # Initialize parallel backend
  registerDoParallel(num_cores)
  
  # Calculate the alignment scores of the original orientation of each sequence to the first sequence in parallel
  original_scores <- foreach(seq = filtered_sequences, .combine = c) %dopar% {
    Biostrings::pairwiseAlignment(first_sequence, seq, scoreOnly = TRUE)
  }
  
  # Calculate the alignment scores of the reverse complement orientation of each sequence to the first sequence in parallel
  reverse_complement_scores <- foreach(seq = filtered_sequences, .combine = c) %dopar% {
    rc_seq <- Biostrings::reverseComplement(seq)
    Biostrings::pairwiseAlignment(first_sequence, rc_seq, scoreOnly = TRUE)
  }
  
  # Perform reverse complement only if the reverse complement has a better alignment score
  corrected_sequences <- mapply(function(seq, original_score, rc_score) {
    if (rc_score > original_score) {
      reverseComplement(seq)
    } else {
      seq
    }
  }, filtered_sequences, original_scores, reverse_complement_scores)
  
  # Convert the corrected_sequences back to a DNAStringSet object
  corrected_sequences <- DNAStringSet(corrected_sequences)
  
  # Check the corrected sequences
  print(corrected_sequences)
  
  # Specify the file path for the output FASTA file
  #output_file <- "2_corrected_db_seq.fasta"
  
  # Write the corrected_sequences object to a FASTA file
  #writeXStringSet(corrected_sequences, file = output_file, format = "fasta")
  
  # Stop parallel backend
  stopImplicitCluster()
  
  
  #### 3_ Format Sequence names and combine 'corrected_sequences' with SST data base for Stenotrophomonas maltophilia ####
  
  # Create a data frame with the original sequence names and the corresponding new names
  library(Biostrings)
  #corrected_sequences <- readDNAStringSet("2_corrected_db_seq.fasta")
  
  # Edit strain names
  Seq.Name <- names(corrected_sequences)
  new_names <- c(sub("^[^ ]+ ", "", as.character(Seq.Name[2:length(Seq.Name)])))
  new_names <- c(sub(" ", "_", new_names, fixed = TRUE))
  
  #pattern.rm <- c(",.*", " chromosome.*"," genome.*"," complete.*"," DNA.*") # Delete characters after some characters
  #new_names. <- gsub(paste(pattern.rm,collapse="|"),"",new_names)
  
  # Load SST data base for Pseudomonas aeruginosa
  db_SST <- readDNAStringSet(paste0("SST_",LOCUS_HiSST,".fasta"))
  # Strain names
  db_SST.Name <- names(db_SST)
  db_SST.Name <- gsub("ID ","",db_SST.Name)
  db_SST.Name <- gsub(" ","_",db_SST.Name)
  
  # Create a data frame with the sequence names and DNA sequences
  Sequences <- data.frame(Names = c(db_SST.Name,new_names), 
                          Seq = paste(as.character(c(db_SST,corrected_sequences[2:length(corrected_sequences)]))))
  
  # Write the data to the output file in FASTA format
  writeLines(paste0(">", Sequences$Names, "\n", Sequences$Seq), con = paste0("3_db_nt_",LOCUS_HiSST,".fasta"))
}

#### Run function 'correct_db_nt' ####

setwd("~/Example_Create_Paeruginosa_database/Step5_BLAST_db/Create_sdaA_db/")

LOCUS_HiSST <- "sdaA"
# Subset the DNAStringSet to remove sequences under 310 bases
Limit_length <- 310

correct_db_nt(LOCUS_HiSST,Limit_length)
