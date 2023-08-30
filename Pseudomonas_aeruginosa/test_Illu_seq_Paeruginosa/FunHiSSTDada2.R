Fun.HiSST.Dada2 <- function(LOCUS_HiSST, Directory_path, path, path_Illu.fasta, 
                            FWD, REV, Amplicon_size, cutadapt,
                            Genus_name = NULL, Species_name = NULL, 
                            blastn_path = NULL, database_HiSST = NULL,
                            SST_db_whithout_primers = NULL,
                            Plot_Quality_Profile = FALSE) {
  library(dada2)
  library(ShortRead)
  library(Biostrings)
  path
  # Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
  fnFs <- sort(list.files(path_Illu.fasta, pattern="_R1.fastq.gz", full.names = TRUE))
  fnRs <- sort(list.files(path_Illu.fasta, pattern="_R2.fastq.gz", full.names = TRUE))
  
  # Extract sample names, assuming filenames have format: SAMPLENAME_RXX.fastq
  sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
  
  #### Identify Primers ####
  
  allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  FWD.orients
  
  #"pre-filter" the sequences just to remove those with Ns
  fnFs.filtN <- file.path(path_Illu.fasta, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
  fnRs.filtN <- file.path(path_Illu.fasta, "filtN", basename(fnRs))
  filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, minLen = minLen, maxLen = maxLen, multithread = TRUE)
  
  #count the number of times the primers appear
  primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
  
  #### Step1_Remove Primers ####
  
  #create output filenames for the cutadapt-ed files, and define the parameters we are going to give the cutadapt command.
  path.cut <- file.path(path_Illu.fasta, "cutadapt")
  if(!dir.exists(path.cut)) dir.create(path.cut)
  fnFs.cut <- file.path(path.cut, basename(fnFs))
  fnRs.cut <- file.path(path.cut, basename(fnRs))
  
  FWD.RC <- dada2:::rc(FWD)
  REV.RC <- dada2:::rc(REV)
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  R1.flags <- paste("-g", FWD, "-a", REV.RC) 
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  R2.flags <- paste("-G", REV, "-A", FWD.RC) 
  # Run Cutadapt
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
  
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

 
  # Forward and reverse fastq filenames have the format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
  cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
  cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))
  
  # Extract sample names, assuming filenames have index names S501 to S522:
  
  # Generate the sequence of numbers from 501 to 522
  index_number <- seq(501, 522)
  # Create the sequence of strings by concatenating "S" and the numbers
  ID_index_primers <- paste0("S", index_number, ".")
  # Perform replacement using regular expressions
  cutFs.name <- gsub(paste(ID_index_primers,collapse="|"), "xyz.", cutFs)
  
  get.sample.name <- function(fname) strsplit(basename(fname), "xyz.")[[1]][2]
  sample.names <- unname(sapply(cutFs.name, get.sample.name))
  head(sample.names)
  
  # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
  sample.names <- sapply(strsplit(basename(sample.names), "_R"), `[`, 1)
  
  
  #### Step2_dada2 ####
  # Plot quality profiles
  
  if (Plot_Quality_Profile) {
    print(plotQualityProfile(cutFs[1]))
    print(plotQualityProfile(cutRs[1]))
  }
  
  # User choice for truncLen
  truncLen_choice <- readline("Enter truncLen values for Forward and Reverse reads - e.g., '220 150' : ")
  truncLen <- as.numeric(strsplit(truncLen_choice, " ")[[1]])
  
  # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(path.cut, "filtered", basename(cutFs))
  filtRs <- file.path(path.cut, "filtered", basename(cutRs))
  
  out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen = truncLen,
                       maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                       compress = TRUE, multithread = TRUE)
  
  head(out)
  
  
  # Learn the Error Rates
  # Tool to visualize the frequency of error rate as a function of quality score.
  # Necessary for the algortith - see the paper
  errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST = 20)
  errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST = 20)
  
  print(plotErrors(errF, nominalQ=TRUE))
  
  
  # Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding 
  # “abundance” equal to the number of reads with that unique sequence. Dereplication substantially reduces 
  # computation time by eliminating redundant comparisons.
  
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  
  
  #We are now ready to apply the core sample inference algorithm to the dereplicated data
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  
  #dadaFs[[1]]
  # dadaRs[[1]]
  
  
  # Merging paired ends
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  # Inspect the merger data.frame from the first sample
  # head(mergers[[1]])
  
  # We can now construct an amplicon sequence variant table (ASV) table
  seqtab <- makeSequenceTable(mergers)
  dim(seqtab)
  # Inspect distribution of sequence lengths
  print("Inspect distribution of sequence lengths :")
  print(table(nchar(getSequences(seqtab))))
  
  # User choice for truncLen
  Amplicon_size_choice <- readline("Enter the expected amplicon size whithout primers - e.g., ´269 270', |or| e.g., '269' : ")
  Amplicon_size <- as.numeric(strsplit(Amplicon_size_choice, " ")[[1]])
  
  #Optional : remove unexpected sequences regarding amplicon size
  seqtab.filt <- seqtab[,nchar(as.character(colnames(seqtab)))== Amplicon_size] #Change value by the expected amplicon size whithout primers
  dim(seqtab.filt)
  table(nchar(getSequences(seqtab.filt)))
  
  # Remove chimera
  seqtab.nochim <- removeBimeraDenovo(seqtab.filt, method="consensus", multithread=TRUE, verbose=TRUE)
  dim(seqtab.nochim)
  sum(seqtab.nochim)/sum(seqtab)
  
  # Track reads through the pipeline
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  print(head(track))
  write.csv(track, file = "Summary.csv")
  
  ## making and writing out standard output files:
  # giving our seq headers more manageable names (ASV_1, ASV_2...)
  asv_seqs <- colnames(seqtab.nochim)
  asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
  
  for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV",LOCUS_HiSST, i, sep="_")
  }
  
  # Write the fasta file
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  write(asv_fasta, paste("no_trim", LOCUS_HiSST, "ASVs.fasta", sep = "_"))
  
  # Write the count table
  asv_tab <- t(seqtab.nochim)
  row.names(asv_tab) <- sub(">", "", asv_headers)
  write.table(asv_tab, paste("no_trim", LOCUS_HiSST, "ASVs_counts.txt", sep = "_"), sep = "\t", quote = FALSE)
  
  # Check if optional arguments are provided
  if (!is.null(Genus_name) &&!is.null(Species_name) && !is.null(blastn_path) && !is.null(database_HiSST)){
    
    #### __Step-3__: BLASTn DNA sequence filtering to remove non-specific ASVs ####
    
    fasta_file <-  paste("no_trim", LOCUS_HiSST, "ASVs.fasta", sep="_")  # Import FASTA file
    #database_HiSST <- paste("db_nt_",LOCUS_HiSST,"/db_nt_",LOCUS_HiSST, sep = "")  # Data base of the HiSST locus for BLASTn
    output_blast <- paste(LOCUS_HiSST,"_blastn_HiSST_ASV.txt", sep = "") # BLASTn results
    
    results_blastn <- system2(blastn_path,args = c("-query",fasta_file,
                                                   "-db",database_HiSST,
                                                   "-outfmt", 6, 
                                                   "-out", output_blast,
                                                   "-perc_identity", 85, # minimum percentage identity for taxonomy assignation
                                                   "-max_target_seqs", 5))
    
    
    blastn_HiSST_ASV <- read.table(output_blast)
    
    # Create a taxon table only for ASV specifics to the pathogen studied #
    asv_tax_BLAST <- data.frame(ASV = sub(">", "", asv_headers),Taxonomy = NA, Percent_identity = NA)
    for (ASV_ID in unique(blastn_HiSST_ASV$V1)) {
      i <- which(blastn_HiSST_ASV$V1 == ASV_ID)
      max <- max(blastn_HiSST_ASV$V3[i[1]:i[length(i)]])
      asv_tax_BLAST.row <- which(asv_tax_BLAST$ASV == ASV_ID)
      asv_tax_BLAST$Percent_identity[asv_tax_BLAST.row] <- max
      asv_tax_BLAST$Taxonomy[asv_tax_BLAST.row] <- paste(blastn_HiSST_ASV$V2[i[1]:i[length(i)]][blastn_HiSST_ASV$V3[i[1]:i[length(i)]] == max], 
                                                         collapse = ", ")
    }
    write.table(asv_tax_BLAST, paste(LOCUS_HiSST,"ASVs_taxonomy.txt",sep = "_"), sep="\t", quote=F)
    
    
    ## keep the DNA sequences if at least one of the two highest scores corresponds to "Genus_Species_name" or if all three scores correspond to "Genus_sp"
    
    # Get the corresponding sequence names
    selected_names <- c()
    
    for (ASV_ID in unique(blastn_HiSST_ASV$V1)) {
      row_nb <- which(blastn_HiSST_ASV$V1 == ASV_ID)
      max_v3 <- max(blastn_HiSST_ASV$V3[row_nb])
      max_v3_row <- row_nb[which(blastn_HiSST_ASV$V3[row_nb]== max_v3)]
      
      if (any(grepl(paste(Genus_name,Species_name,"SST", sep = "_"), blastn_HiSST_ASV$V2[row_nb[1]:(row_nb[1]+2)])) && # Check if the 3 first Species are Genus_Species_SST, and :
          grepl(paste(Genus_name,"sp", sep = "_"), blastn_HiSST_ASV$V2[row_nb[1]]) | grepl(paste(Genus_name,Species_name,sep = "_"), blastn_HiSST_ASV$V2[row_nb[1]]) || # Check if the first Species is Genus_sp or Genus_Species
          any(grepl(paste(Genus_name, Species_name, "SST", sep = "_"), blastn_HiSST_ASV$V2[max_v3_row]))) { # OR, Check if the maximum percent_identity of the ASV include Genus_Species_SST
        selected_names <- c(selected_names, blastn_HiSST_ASV$V1[row_nb[1]]) # If the conditions are met, then keep the ASV
      }
    }
    
    
    # Read the DNAStringSet object
    ASVs_fasta_file <- readDNAStringSet(fasta_file)
    
    # Filter the DNA sequences based on selected names
    filtered_sequences <- ASVs_fasta_file[names(ASVs_fasta_file) %in% selected_names]
    
    # Write the filtered sequences to a new fasta file
    writeXStringSet(filtered_sequences, file = paste(LOCUS_HiSST,"filtered_ASV.fasta", sep = "_"))
    
    
    #### __Step-4__: Create "LOCUS_ASV_samples.txt" table, used for Jaccard dendrograms ####
    
    seqtab.nochim.trim <- seqtab.nochim[,colnames(seqtab.nochim) %in% as.character(filtered_sequences)]
    ASV.sample <- seqtab.nochim.trim
    colnames(ASV.sample) <- names(filtered_sequences)
    
    #Export results in a matrix of samples (rows) and AVS counts (columns), in descending order.
    ASV.sample.sort <- ASV.sample
    for (i in ncol(ASV.sample.sort):1) {
      ASV.sample.sort <- ASV.sample.sort[order(-ASV.sample.sort[,i]),]
    }
    ASV.sample.sort <- data.frame(Samples = row.names(ASV.sample.sort),ASV.sample.sort,row.names = NULL)
    
    write.table(ASV.sample.sort, paste(LOCUS_HiSST,"ASV_samples.txt",sep = "_"), sep="\t", quote=F)
    
    #### __Step-5__: Create "Samples_and_ASV-ST.txt" if needed for HiSST data base update, in 'HiSST-Assignation.R' script ####
    # Check if optional arguments are provided
    if (!is.null(SST_db_whithout_primers)){
      
      taxa_assign <- assignSpecies(seqtab.nochim.trim, SST_db_whithout_primers, allowMultiple = TRUE)#Taxonomic assignment to the ST level by exact matching
      colnames(taxa_assign) <- c("Species", "ST")
      
      ASV.sample.mat <- seqtab.nochim.trim
      colnames(ASV.sample.mat) <- names(filtered_sequences)
      name.ST.ASV <- data.frame(ST = taxa_assign[,2], ASV = row.names(taxa_assign), row.names = NULL)
      
      # Change column names without any ST-ID by ASV names 
      name.STASV <- name.ST.ASV
      for (i in 1:nrow(name.STASV)) {
        ifelse(is.na(name.STASV[i,1]),name.STASV[i,1]<-name.STASV[i,2],name.STASV[i,1])
      }
      colnames(ASV.sample.mat) <- name.STASV$ST
      
      # Compute a data frame with Sample names, and ST-ID or ASV-ID
      df <- data.frame(Samples=rownames(ASV.sample.mat), ST=nrow(ASV.sample.mat))
      for (j in 1:nrow(ASV.sample.mat)) {
        max(ASV.sample.mat[j,]) -> x
        for (i in 1:ncol(ASV.sample.mat)) {
          ifelse(ASV.sample.mat[j,i] == x,
                 z <- colnames(ASV.sample.mat)[i],z <- "NA")
          ifelse(z == "NA",NA, df$ST[j] <- z)
        }
      }
      
      write.table(df, paste(LOCUS_HiSST,"Samples_and_ASV-ST.txt",sep = "_"), sep="\t", quote=F, row.names = F)
    }
  }
}