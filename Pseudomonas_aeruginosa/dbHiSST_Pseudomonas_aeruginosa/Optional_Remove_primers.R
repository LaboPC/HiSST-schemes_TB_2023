library(Biostrings)

remove_primers <- function(input_files, forward_primers, reverse_primers) {
  if (length(input_files) != length(forward_primers) || length(input_files) != length(reverse_primers)) {
    stop("The number of input files does not match the number of primers.")
  }
  
  for (i in seq_along(input_files)) {
    input_file <- input_files[i]
    output_file <- paste0(tools::file_path_sans_ext(input_file), "_no_primers.fasta")
    
    REV.RC <- DNAStringSet(dada2:::rc(reverse_primers[i]))
    
    R1.flags <- paste0("-g ", forward_primers[i], ";max_error_rate=0.1", " -a ", REV.RC, ";max_error_rate=0.1") # 'max_error_rate=0.1' to set error tolerance at 10%
    
    system2(cut_path, args = c(R1.flags, "-n", 2,
                               "-o", output_file,
                               input_file))
  }
}

#### Change configurations below ####
setwd("path_to_your_working_directory/remove_primers/") # Change to your working directory path
input_files <- c("db_nt_bvgS.fasta", "db_nt_pheT.fasta", "db_nt_btuB.fasta", "db_nt_sdaA.fasta")
forward_primers <- c("ACGGCGACGARCTGTTGC", "GCGTGGACTTCTTCGACGC", "GCCAAGCCGTTCTTCTCCG", "ATCGTCGAGGACCGCACG")
reverse_primers <- c("GGCATGGTCGGCGTAACC", "GACAGCTCGCGGAACTTCG", "CAGGTTCTGCTCGCCGTC", "GTAGAGRTTGACCCAGTCGAGC")
#Cutadapt executable path
cut_path <- "path_to_cutadapt/cutadapt.exe" 

#### Run function to remove primers ####
remove_primers(input_files, forward_primers, reverse_primers)

