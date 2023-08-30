
library(Biostrings)
library(dplyr)

#### Change configurations below ####

# Change loci for each locus name of the HSST scheme
loci <- c("glnG","ribA", "tycC", "yvoA")
# Change the path file
file_path <- "~/path_to_your_file/"


#### Step 1: Combine fasta files in a data frame for each locus ####

#' The function 'combine_fasta' combines two fasta files with identical sequences and stores the sequence names in a new data frame.
#' The function uses the Biostrings package to read the fasta files and 
#' create data frames with columns for sequence names and sequence strings. 
#' The two data frames are then merged based on the sequence strings, 
#' and the resulting data frame is selected and renamed as needed. 

combine_fasta <- function(locus, SST_locus) {
  locus_seqs <- readDNAStringSet(locus)
  sst_seqs <- readDNAStringSet(SST_locus)
  
  locus_df <- data.frame(Strain = names(locus_seqs), Seq = as.character(locus_seqs))
  sst_df <- data.frame(SST = names(sst_seqs), Seq = as.character(sst_seqs))
  
  merged_df <- merge(locus_df, sst_df, by = "Seq")
  
  result_df <- merged_df %>% arrange(Strain) %>% select(Strain, SST, Seq)
  
  return(result_df)
}


#### ___ Run the script below ####
combined_data <- lapply(loci, function(locus) {
  combine_fasta(
    paste0(file_path, "Step2_keep_redond_loci/HiSST_", locus, ".fasta"),
    paste0(file_path, "Step3_Keep_unique-SST/SST_", locus, ".fasta")
  )
})

#### Step 2: Concatenate loci ####

setwd(paste0(file_path, "Step4_Assign_ST"))

locus_data <- lapply(seq_along(loci), function(i) {
  df <- data.frame(
    ID = combined_data[[i]]$Strain,
    SST = gsub(".*SST-", "", combined_data[[i]]$SST)
  )
  colnames(df)[2] <- loci[i]
  df
})

# Join all loci in a data frame
query <- Reduce(function(x, y) full_join(x, y, by = "ID"), locus_data)

query$combine_ST <- mapply(paste, sep = "-",query[2],query[3],query[4],query[5])


#### Step 3: Create HiSST data base for P. aeruginosa ####

# obtain unique HiSST profiles
unique_ST <- unique(query$combine_ST)
new.HiSST.ID <- data.frame(combine_ST = unique_ST, HiSST = seq(1:length(unique_ST)))
names(new.HiSST.ID)[1] <- "combine_ST"

# Combine all strain with HiSST profiles
HiSST_data <- full_join(query,new.HiSST.ID, by = "combine_ST")

#### Export HiSST profiles ####
# Keep only the unique STs
unique_df <- HiSST_data[!duplicated(HiSST_data$HiSST),]
# Create HiSST data frame
HiSST.profiles <- data.frame(HiSST = unique_df$HiSST,
                             unique_df[2],
                             unique_df[3],
                             unique_df[4],
                             unique_df[5])

write.table(HiSST_data, file = "HiSST_data_allstrains.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(HiSST.profiles, file = "new_version_HiSSTprofiles-St.txt", sep = "\t", row.names = FALSE, quote = FALSE)