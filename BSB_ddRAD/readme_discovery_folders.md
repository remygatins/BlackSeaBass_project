## /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/

### Files
- all raw FASTQ files from sequencing facility  (moved to folder `raw`)
- both original PERL scripts from sequencing facility:
  - gbstrim.pl - trims off padding sequences & Illumina adapter  
  - resync.pl - resync R1 and R2 files post-trimming
- scripts edited by Thais: 
  - gbstrimedited.pl: edited by Thais to include more combinations of padding sequences (see notebook for full explanation).
  - gbstrimeditedkeepadapter: edited by Thais, same as gbstrimedited.pl plus has dummy sequences in place of the real Illumina adapter sequences as a trick to keep the Illumina adapters.
  

### Folders
- raw: all raw FASTQ files from sequencing facility 
- trim_padding: files trimmed using gbstrimedited.pl, i.e., padding sequences are cut off, and Illumina adapters are cut off.
- trim_padding_keepadapter: files trimmed using gbstrimeditedkeepadapter.pl, i.e., padding sequences are cut off and Illumina adapters are retained.


## /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding

### Files:  
- gbstrimedited.pl

### Folders:
R1_trimmed: 
 - gbstrimeditedkeepadapter.pl
 - gbstrimedited.pl input & output files: all R1 raw input FASTQ.GZ files plus 4 sets of files per R1 raw file (trim.fastq; trim.fastq.no3cutsite; trim.fastq.nocutsite; trim.fastq.nocutsite.adapters)
 - kept_seqR1.txt: contains a list of number of R1 sequences kept (seqs that will move further down the pipeline) per individual sequenced.

R2_trimmed: 
 - same as in R1_trimmed, but for R2 files.

synced_renamed:
 - all R1 and R2 trimmed files (trim.fastq).   
 - all renamed to match dDocent naming convention & resynced files (.F.fq.gz & .R.fq.gz)


## /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding_keepadapter

### Files:
- gbstrimedited.pl (edited to include more padding sequences - a fix to the original gbstrim.pl from the sequencing facility)
- gbstrimeditedkeepadapter.pl (same as gbstrimedited.pl plus a dummy adapter replaces the real Illumina adapter in order to keep the Illumina adapter)
- resync.pl (from sequencing facility, this script is used to resync R1 and R2 files following trimming of padding sequences)

### Folders:
R1_trimmed_keepadapter: 
 - gbstrimeditedkeepadapter.pl
 - gbstrimeditedkeepadapter.pl input & output files: all R1 raw input FASTQ.GZ files plus 4 sets of files per R1 raw file (trim.fastq; trim.fastq.no3cutsite; trim.fastq.nocutsite; trim.fastq.nocutsite.adapters)
 - kept_seqR1_keepadapter.txt: contains a list of number of R1 sequences kept (seqs that will move further down the pipeline) per individual sequenced.

R2_trimmed_keepadapter: 
 - same as in R1_trimmed_keepadapter, but for R2 files.

synced_renamed_keepadapter: 
 - all R1 and R2 trimmed files (trim.fastq).   
 - all renamed to match dDocent naming convention & resynced files (.F.fq.gz & .R.fq.gz)
 - resync.pl

ddocent_subset01: a set of 20 files for 10 individuals and 2 populations 
  - NC_242.F.fq  NC_242.R.fq  NC_243.F.fq  NC_243.R.fq  NC_244.F.fq  NC_244.R.fq  NC_245.F.fq  NC_245.R.fq  NC_246.F.fq  NC_246.R.fq  NJ_106.F.fq  NJ_106.R.fq  NJ_108.F.fq  NJ_108.R.fq  NJ_109.F.fq  NJ_109.R.fq  NJ_112.F.fq  NJ_112.R.fq  NJ_113.F.fq  NJ_113.R.fq

ddocent_subset02: a set of 40 files for 20 individuals and 2 populations 
 - NC_233.F.fq  NC_234.F.fq  NC_235.F.fq  NC_237.F.fq  NC_238.F.fq  NC_239.F.fq  NC_240.F.fq  NC_241.F.fq  NJ_114.F.fq  NJ_118.F.fq  NJ_119.F.fq  NJ_121.F.fq  NJ_122.F.fq  NJ_124.F.fq  NJ_128.F.fq  NJ_129.F.fq  NJ_130.F.fq  NJ_131.F.fq  NJ_132.F.fq  NJ_133.F.fq
NC_233.R.fq  NC_234.R.fq  NC_235.R.fq  NC_237.R.fq  NC_238.R.fq  NC_239.R.fq  NC_240.R.fq  NC_241.R.fq  NJ_114.R.fq  NJ_118.R.fq  NJ_119.R.fq  NJ_121.R.fq  NJ_122.R.fq  NJ_124.R.fq  NJ_128.R.fq  NJ_129.R.fq  NJ_130.R.fq  NJ_131.R.fq  NJ_132.R.fq  NJ_133.R.fq)

