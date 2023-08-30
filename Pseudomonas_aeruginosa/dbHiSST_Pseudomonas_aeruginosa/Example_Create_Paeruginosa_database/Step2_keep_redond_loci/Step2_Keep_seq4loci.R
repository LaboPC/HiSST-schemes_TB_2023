library(Biostrings)

setwd("~/Example_Create_Paeruginosa_database/Step2_keep_redond_loci/")

# Define file paths and locus names
PATH.LOCUS <- "~/Example_Create_Paeruginosa_database/Step1_db_Pseudo/"
locus_names <- c("btuB", "bvgS", "pheT", "sdaA")

# Read DNA sequences for each locus
locus_seqs <- lapply(locus_names, function(locus) {
  readDNAStringSet(paste0(PATH.LOCUS,locus,"/", locus, ".fasta"))
})

# Get sequence names for each locus
locus_seq_names <- lapply(locus_seqs, names)

# Find common sequence names across all loci
common_seqs <- Reduce(intersect, locus_seq_names)

# Create a function to get common sequences for a given locus
get_common_seqs <- function(locus_seq) {
  idx <- match(common_seqs, names(locus_seq))
  locus_seq_common <- locus_seq[idx]
  return(locus_seq_common)
}

# Get common sequences for each locus
locus_common_seqs <- lapply(locus_seqs, get_common_seqs)

# Write common sequences to new files
for (i in seq_along(locus_names)) {
  writeXStringSet(locus_common_seqs[[i]], paste0("HiSST_", locus_names[i], ".fasta"))
}
