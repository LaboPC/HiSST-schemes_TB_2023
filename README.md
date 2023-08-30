# HiSST schemes: databases and R scripts
HiSST scheme databases and R-scripts used for opportunistic pathogens pure culture and environmental DNA survey, using the Illumina next-generation sequencing (NGS) method.

### HiSST-dada2 pipeline : Bioinformatical pipeline for HiSST analysis:
Originate from [DADA2 pipeline](https://benjjneb.github.io/dada2/index.html) and was adapted for HiSST schemes analysis, including additional steps to enhance the accuracy of taxonomic assignment. 
The entire process can be executed with the R script [Script_RUN_FunHiSSTDada2.R](https://github.com/LaboPC/HiSST-schemes_TB_2023/blob/main/Script_RUN_FunHiSSTDada2.R), using a unified function called [FunHiSSTDada2.R](https://github.com/LaboPC/HiSST-schemes_TB_2023/blob/main/FunHiSSTDada2.R).

The pipeline involved these steps: (i) Processed raw sequencing reads using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) to remove primers, (ii) Applied default [DADA2](https://github.com/benjjneb/dada2) settings for error correction, denoising, and paired-end merging, (iii) Filtered DNA sequences with BLASTn to eliminate non-specific ASVs, (iv) Generated "LOCUS_ASV_samples.txt" for Jaccard dendrograms, and (v) Created "Samples_and_ASV-ST.txt" if HiSST database update required.

# HiSST scheme for _Serratia marcescens_ :

The databases for each locus are available in the [Serratia_marcescens](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main/Serratia_marcescens) directory file.

- R-scripts used to develop, test and validate the HiSST scheme for _Serratia marcescens_ are available at https://github.com/TBourd/R_scripts_for_HiSST_scheme

- HiSST application example: R-scripts and data used for _Serratia marcescens_ colonizations survey, using the HiSST scheme, are available at https://github.com/TBourd/R_scripts_HiSST_SM-colonizations (see article: Bourdin et al., 2023)

###### If you use this HiSST scheme in your research, please cite :
> Bourdin T, Monnier A, Benoit MÈ, Bédard E, Prévost M, Quach C, Déziel E, Constant P. A High-Throughput Short Sequence Typing Scheme for Serratia marcescens Pure Culture and Environmental DNA. Appl Environ Microbiol. 2021. DOI: https://doi.org/10.1128/AEM.01399-21


# HiSST scheme for _Pseudomonas aeruginosa_ :

The [Pseudomonas_aeruginosa](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main/Pseudomonas_aeruginosa) directory contains locus databases, comprising:

- "Pseudomonas_DataBase HiSST.xlsx": HiSST scheme database for the short sequence types (SSTs) of _btuB_, _bvgS_, _pheT_ and _sdaA_ loci.
- "version_HiSSTprofiles-Pa.txt": Table mapping HiSST profiles to their SSTs.
- The [dbHiSST_Pseudomonas_aeruginosa](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main/Pseudomonas_aeruginosa/dbHiSST_Pseudomonas_aeruginosa) folder:
  - Locus databases used in the HiSST-dada2 pipeline (see [above](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main#bioinformatical-pipeline-for-hisst-analysis)).
  - Five R Scripts to create _P. aeruginosa_ HiSST database in five steps (plus an optional script to remove primers if necessary).
  - [Example_Create_Paeruginosa_database](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main/Pseudomonas_aeruginosa/dbHiSST_Pseudomonas_aeruginosa/Example_Create_Paeruginosa_database) folder with initial locus database creation steps.
- The [test_Illu_seq_Paeruginosa](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main/Pseudomonas_aeruginosa/test_Illu_seq_Paeruginosa) folder: HiSST-dada2 pipeline example using real Illumina Miseq data for each HiSST locus.

###### If you use this HiSST scheme in your research, please cite :
> Bourdin T, Benoit M-È, Bédard E, Prévost M, Quach C, Déziel E, Constant P. _Submitted for publication_. High-Throughput Short Sequence Typing Schemes for _Pseudomonas aeruginosa_ and _Stenotrophomonas maltophilia_ pure culture and environmental DNA.


# HiSST scheme for _Stenotrophomonas maltophilia_ :

The [Stenotrophomonas_maltophilia](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main/Stenotrophomonas_maltophilia) directory contains locus databases, comprising:

- "Stenotrophomonas_DataBase HiSST.xlsx": HiSST scheme database for the short sequence types (SSTs) of _yvoA_, _glnG_, _tycC_ and _ribA_ loci.
- "version_HiSST-ST_database.txt": Table mapping HiSST profiles to their SSTs.
- The [dbHiSST_Stenotrophomonas_maltophilia](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main/Stenotrophomonas_maltophilia/dbHiSST_Stenotrophomonas_maltophilia) folder:
  - Locus databases used in the HiSST-dada2 pipeline (see [above](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main#bioinformatical-pipeline-for-hisst-analysis)).
  - Five R Scripts to create _S. maltophilia_ HiSST database in five steps (plus an optional script to remove primers if necessary).
  - [Example_Create_Smaltophilia_database](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main/Stenotrophomonas_maltophilia/dbHiSST_Stenotrophomonas_maltophilia/Example_Create_Smaltophilia_database) folder with initial locus database creation steps.
- The [test_Illu_seq_Smaltophilia](https://github.com/LaboPC/HiSST-schemes_TB_2023/tree/main/Stenotrophomonas_maltophilia/test_Illu_seq_Smaltophilia) folder: HiSST-dada2 pipeline example using real Illumina Miseq data for each HiSST locus.

###### If you use this HiSST scheme in your research, please cite :
> Bourdin T, Benoit M-È, Bédard E, Prévost M, Quach C, Déziel E, Constant P. _Submitted for publication_. High-Throughput Short Sequence Typing Schemes for _Pseudomonas aeruginosa_ and _Stenotrophomonas maltophilia_ pure culture and environmental DNA.

 _______________________________________________________

# Publications

 
- Bourdin T, Monnier A, Benoit MÈ, Bédard E, Prévost M, Quach C, Déziel E, Constant P. 2021. A High-Throughput Short Sequence Typing Scheme for _Serratia marcescens_ Pure Culture and Environmental DNA. Appl Environ Microbiol. DOI: https://doi.org/10.1128/AEM.01399-21

- Bourdin T, Benoit M-È, Monnier A, Bédard E, Prévost M, Charron D, Audy N, Gravel S, Sicard M, Quach C, Déziel E, Constant P. 2023. _Serratia marcescens_ Colonization in a Neonatal Intensive Care Unit Has Multiple Sources, with Sink Drains as a Major Reservoir. Applied and Environmental Microbiology 89:e00105-23. DOI: https://doi.org/10.1128/aem.00105-23

- Bourdin T, Benoit M-È, Bédard E, Prévost M, Quach C, Déziel E, Constant P. _Submitted for publication_. High-Throughput Short Sequence Typing Schemes for _Pseudomonas aeruginosa_ and _Stenotrophomonas maltophilia_ pure culture and environmental DNA.
