library(Biostrings)

remove_primers <- function(input_files, forward_primers, reverse_primers) {
  if (length(input_files) != length(forward_primers) || length(input_files) != length(reverse_primers)) {
    stop("The number of input files does not match the number of primers.")
  }
  
  for (i in seq_along(input_files)) {
    input_file <- input_files[i]
    output_file <- paste0(tools::file_path_sans_ext(input_file), "_no_primers.fasta")
    
    REV.RC <- DNAStringSet(dada2:::rc(reverse_primers[i]))
    
    R1.flags <- paste0("-g ", forward_primers[i], ";max_error_rate=0.2", " -a ", REV.RC, ";max_error_rate=0.2") # 'max_error_rate=0.2' to set error tolerance at 20%
    
    system2(cut_path, args = c(R1.flags, "-n", 2,
                               "-o", output_file,
                               input_file))
  }
}

#### Change configurations below ####
setwd("~/path_to_your_working_directory/remove_primers/") # Change to your working directory path
input_files <- c("db_nt_ribA.fasta", "db_nt_tycC.fasta", "db_nt_yvoA.fasta", "db_nt_glnG.fasta")
forward_primers <- c("CTGCCCTCGYTGGGCTA", "TGTACACCGARCAGGTCGAG", "CCGAGAGCGGCATGATCGA","GTGATGTCGGCCTAYACCG")
reverse_primers <- c("GACGATGATCGCSACCTGG", "TCTTGGCGTTGTGACGGATATC", "CAGGCARCGCATCGCCA","GCCACCAGYTCCTTGCC")
#Cutadapt executable path
cut_path <- "path_to_cutadapt/cutadapt.exe" 

#### Run function to remove primers ####
remove_primers(input_files, forward_primers, reverse_primers)

