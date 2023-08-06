#### __Step-1__: Configurations ####

## A. Modify paths to desired directory ##

# Change Locus name
LOCUS_HiSST <- "LOCUS" # locus name (e.g., "bssA")

# Path to directory
Directory_path <- "/path/to/Directory_path/" # "path/to/Directory_path/"
# Path to corresponding locus working directory
path <- setwd(paste(Directory_path,LOCUS_HiSST, sep = "")) # "path/to/Directory_path/LOCUS_HiSST"
# Change path_Illu.fasta to corresponding Illumina fasta files directory
path_Illu.fasta <- paste0(path,"/Illu_fasta") # "path/to/Directory_path/LOCUS_HiSST/fasta/files"
list.files(path_Illu.fasta)

## B. Modify the script configuration for the studied locus ##

# Change Locus Primers
FWD <-  "FORWARDPRIMER" # (e.g., "CGCAGTTTCTCAACGCYATCG" for bssA)
REV <-  "REVERSEPRIMER" # (e.g., "CGAATGGCCGTTGGATTCGATC" for bssA)

# Change 'Cutadapt' executable path
cutadapt <- "cutadapt"
system2(cutadapt, args = "--version")

# maxLen : (Optional). Default Inf (no maximum). Remove reads with length greater than maxLen. maxLen is enforced before trimming and truncation.
maxLen = Inf
# minLen : (Optional). Default 20. Remove reads with length less than minLen. minLen is enforced after trimming and truncation.
minLen = 20

## C. Optional Steps (3-5) : database "db_nt_LOCUS" and file 'SST_LOCUS_no_primers.fasta' are required ##

# Path to the BLASTn executable
blastn_path <- "/path/to/BLASTn.exe/" # Change path (e.g. "ncbi-blast-2.14.0+/bin/blastn")
system2(blastn_path, args = "-h")

# Pathogen names
Genus_name <- "GENUS" # e.g. "Serratia"
Species_name <- "SPECIES" # e.g. "marcescens"

# Data base of the HiSST locus for BLASTn (Step-3)
database_HiSST <- paste("db_nt_",LOCUS_HiSST,"/db_nt_",LOCUS_HiSST, sep = "")
# Data base of the HiSST locus, if needed to create "Samples_and_ASV-ST.txt" (Step-5)
SST_db_whithout_primers <- paste("SST",LOCUS_HiSST,"no_primers.fasta", sep = "_")


#Load Function Fun.HiSST.Dada2()
source("/path/to/FunHiSSTDada2.R")

#### ____ RUN FUNCTION: 'Fun.HiSST.Dada2' ____ ####
Fun.HiSST.Dada2(LOCUS_HiSST, Directory_path, path, path_Illu.fasta,
                FWD, REV, Amplicon_size, cutadapt,
                Genus_name = Genus_name, Species_name = Species_name,
                blastn_path = blastn_path, database_HiSST = database_HiSST,
                SST_db_whithout_primers = SST_db_whithout_primers,
                Plot_Quality_Profile = FALSE) # (Optional) Plot Quality_Profile ? Default = FALSE


