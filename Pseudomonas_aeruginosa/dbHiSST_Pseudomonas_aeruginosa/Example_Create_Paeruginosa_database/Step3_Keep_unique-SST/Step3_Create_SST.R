

#### Function to create SST data bases ####

Create_database_per_SST <- function(input_files, HiSST_locus) {
  library(Biostrings)
  duplicate_sequences <- character()
  
  for (i in seq_along(input_files)) {
    locus <- readDNAStringSet(input_files[i])
    
    unique_sequences <- unique(locus)
    sequence_names <- names(locus)
    
    corresponding_names <- vector("list", length(unique_sequences))
    
    for (j in seq_along(unique_sequences)) {
      corresponding_names[[j]] <- paste(sequence_names[which(locus == unique_sequences[j])], collapse = ";")
    }
    
    df.locus <- data.frame(
      locus_SST = paste0("SST-", seq(1:length(unique_sequences))),
      Strains = paste(corresponding_names),
      Sequence = as.character(unique_sequences),
      stringsAsFactors = FALSE,
      row.names = seq(1:length(unique_sequences))
    )
    
    duplicates <- duplicated(df.locus$Sequence)
    if (any(duplicates)) {
      duplicate_sequences <- c(duplicate_sequences, paste("Duplicate sequences:", which(duplicates)))
    }
    
    write.table(df.locus, paste("SST_", HiSST_locus[i], "_data.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
    writeLines(paste0(">ID ",Genus_name,"_",Species_name," ", df.locus$locus_SST, "\n", df.locus$Sequence), con = paste0("SST_", HiSST_locus[i], ".fasta"))
  }
  
  if (length(duplicate_sequences) > 0) {
    print("Duplicate sequences:")
    print(duplicate_sequences)
  } else {
    print("No duplicate sequences found.")
  }
}


#### Change configurations below ####
setwd("~/Example_Create_Paeruginosa_database/Step3_Keep_unique-SST/")

HiSST_locus <- c("btuB","bvgS", "pheT", "sdaA")
input_files <- paste("~/Example_Create_Paeruginosa_database/Step2_keep_redond_loci/HiSST_",
                     HiSST_locus, ".fasta", sep = "")
# Pathogen names
Genus_name <- "Pseudomonas"
Species_name <- "aeruginosa"


#### Run the function 'Create_database_per_SST' ####

Create_database_per_SST(input_files,HiSST_locus)
