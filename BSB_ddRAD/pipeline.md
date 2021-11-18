## 1. Get the files

- download FASTQ files from sequencing facility using:

```wget -nH -np -N -r --cut-dirs 2 --no-check-certificate --user lotterhos --password Ride-Stir-Secretary-Size-7 https://umgcdownload.msi.umn.edu/lotterhos/210713_A00223_0599_BH7YTHDRXY/Lotterhos_Project_001```

- check download was sucessfull using md5sum according to https://drk-lo.github.io/lotterhoslabprotocols/code_checksums/ - PASS

## 2. Enzyme pair used here: BamHI and Mspl

BamHI cut: sticky ends with 5' end overhangs of G/GATCC and palindromic complementary CCTAG/G

Mspl cut: sticky ends with 5' end overhangs of C/CGG and palindromic complementary GGC/C  

Barcode sequences: see project summary file.

## 3. name convention for dDocent (this works but we're not using those; renaming to convention was integrated into step 5).

file names provided by sequencing facility are in the following format (example: Cs_MA_298_S1_R1_001.fastq.gz): 
      - Cs_XX_YYY_S1_R1_001.fastq.gz (forward)
      - Cs_XX_YYY_S1_R2_001.fastq.gz (reverse)

dDocent requires the following name convention:
- for raw sequences:

    - Pop1_Sample1.F.fq.gz (raw, forward)
    - Pop1_Sample1.R.fq.gz (raw, reverse)

- for trimmed sequences:

    - Pop1_001.R1.fq.gz (trimmed, forward)
    - Pop1_001.R2.fq.gz (trimmed, reverse)


~~I believe our sequences are already trimmed - I see the barcode assigned to each sequence, but I don't see the overhang or the padding sequence; additionally, the sequence files are named ...R1... and ...R2... which makes me think they are trimmed~~ SEQUENCES ARE NOT TRIMMED, they contain both the padding sequences and the overhangs. See 4.

~~Assuming we have the trimmed sequences,~~ we need to do the following to comply with dDocent name convention:

1) replace **R1_001.fastq.gz** with **.R1.fq.gz** ```sed 's/_R1_001.fastq.gz/.R1.fq.gz/g'```
2) eliminate the underscores ```sed 's/_//g'```
      - the following line does both:

```for f in *_R1_001.fastq.gz ; do mv $f `echo "$f" | sed 's/_R1_001.fastq.gz/.R1.fq.gz/g' | sed 's/_//g'` ; done```

      
   - input: Cs_MA_298_S1_R1_001.fastq.gz
   - output: CsMA298S1.R1.fq.gz
      
3) add the prefix PopXX_ (where XX is the region or pop ID; this is the only underscore allowed in the file name, between popID and sampleID)
      - use the following two lines:
      
```prefix=PopXX_```

```for filename in *.fq.gz; do mv "$filename" "$prefix$filename" ; done```
    
  - input: CsMA298S1.R1.fq.gz
  - output: PopMA_CsMA298S1.R1.fq.gz


Meeting with Katie: file names should be MA_298.R1.fq.gz


## 4. Trim padding sequences using cutadapt

Emailed the facility and learned that our data still have the padding sequences and overhangs. Use their code to trim the data, the following steps:

i) Install cutadapter with conda (Python 3.6 or up is a dependency, already installed).
https://cutadapt.readthedocs.io/en/stable/installation.html

on OOD Discovery:

- logged in to home, where I had Miniconda3-latest-Linux-x86_64.sh already installed.
- ran `conda update conda` since I hadn't used it in a while
- then followed the steps in the link above:
     - conda create -n cutadaptenv cutadapt (environment location: /home/tbittar/miniconda3/envs/cutadaptenv) - this only needs to be done once.
     - conda activate cutadaptenv (activating the newly crated env) - this needs to be done each time we want to use the environment.
     - cutadapt --version (just to check it worked)

ii) Download the repository from https://bitbucket.org/jgarbe/gbstrim/src/master/ - only done once, must be inside the working folder.

- used Globus to transfer the repo to my working folder on Discovery (working folder = /work/lotterhos/NOAA...)

iii) ~~use their script:~~ DON'T - keep reading.

R1 reads:

`gbstrim.pl --enzyme1 ApeKI --enzyme2 ApeKI --fastqfile sample1_R1.fastq.gz --read R1 --outputfile sample1.trim.fastq --verbose --threads 24 --minlength 50`

NOTES on this line of code: 
  - add `perl` before the command; 
  - change the enzymes to what we used (mspi & bamhi); 
  - change input and output file names to match ours.

The script 'worked' as the padding sequences were removed and all retained sequences start with the overhang (see screen shots below), HOWEVER:
  - Approx. 40% of the sequences from R1 were discarded.
  - Approx. 99% of the sequences from R2 were discarded- see hurdle.


Padding sequences (green, variable lengths), R1 overhang (yellow, CGG):
![ddRAD_R1_padding_example](https://user-images.githubusercontent.com/52291277/138742166-121f016b-6cb7-402e-9182-d77c39dbd8db.png)

After padding trimmed out, R1 overhangs remain:
![ddRAD_R1_padding_cut](https://user-images.githubusercontent.com/52291277/138744664-ecaf076b-573d-4d08-8982-b47d613381b1.png)




-----------------------------------------------

**20211020 - hurdle**

  - An average of 42% of my sequences in R1 files are being discarded. Too high, but still leaving a lot to work with. 
  - However, ~99% of sequences in R2 files are being discarded, althought didn't run all of the files yet. Emailed facility for help. 
  - In the meantime, I looked into some of the raw sequences. I found that most sequences do have the padding sequence but they also have A SINGLE extra base before the padding sequence begins; here's what I'm looking at ('extra' base is in bold; middle sequence is the padding; crossed out is the overhang for BamHI -GATCC):
            
**C** ACTCTG ~~GATCC~~  
**G** TTCGACAT ~~GATCC~~  
**T** AG ~~GATCC~~  
**C** AAGT ~~GATCC~~  
**T** ACGAA ~~GATCC~~  
**C** ACTCTG ~~GATCC~~  
**C** ~~GATCC~~  
**T** CGATGTGCT ~~GATCC~~  
**T** CGATGTGCT ~~GATCC~~  
**T** AG ~~GATCC~~  
**T** CGATGTGCT ~~GATCC~~  
**C** ~~GATCC~~  
**G** TCA ~~GATCC~~  
**G** TCA ~~GATCC~~  
**C** ~~GATCC~~  
**T** ACGAA ~~GATCC~~  
**C** ~~GATCC~~  
**G** TTCGACAT ~~GATCC~~  
**T** AG ~~GATCC~~  
**C** ~~GATCC~~  
**T** ACGAA ~~GATCC~~  
**G** TTCGACAT ~~GATCC~~  
**A** GTACGGT ~~GATCC~~  

- So, I edited gbstrim.pl to include all combinations of each base A, C, T, G followed by each of the padding sequences corresponding to BamHI R2 (gatc_r2 in the script). All possible padding sequences in the edited script are: $gatc_r2 = ",G,AG,TCA,AAGT,ACGAA,ACTCTG,GTACGGT,TTCGACAT,CGATGTGCT,A,C,T,AAG,CAG,GAG,TAG,ATCA,CTCA,GTCA,TTCA,AAAGT,CAAGT,GAAGT,TAAGT,AACGAA,CACGAA,GACGAA,TACGAA,AACTCTG,CACTCTG,GACTCTG,TACTCTG,AGTACGGT,CGTACGGT,GGTACGGT,TGTACGGT,ATTCGACAT,CTTCGACAT,GTTCGACAT,TTTCGACAT,ACGATGTGCT,CCGATGTGCT,GCGATGTGCT,TCGATGTGCT";
- These edits were saved as gbstrimedited.pl
- After running the edited script, I'm getting ~17-19% discarded sequences on a few files ran manually, so will move on with the edited script.
------------------------------------------------------

iv) So re-run both R1 and R2 files using the edited script (edited in my local computer and transferred to the working folder on Discovery through Globus); the lines of code below work to run in batch **assuming raw (.fastq.gz) R1 and R2 files are in separate folders** (all R1 files plus the gbstrimedited.pl script file are in the R1_trimmed folder and the code is run out of that folder and all R2 files plus the gbstrimedited.pl script file are in the R2_trimmed folder and the code is run a second time out of that folder).

`for i in *.fastq.gz; do perl gbstrimedited.pl --enzyme1 mspi --enzyme2 bamhi --fastqfile "$i" --read R1 --outputfile "${i%%.*}".trim.fastq --verbose --threads 24 --minlength 50; done`

`for i in *.fastq.gz; do perl gbstrimedited.pl --enzyme1 mspi --enzyme2 bamhi --fastqfile "$i" --read R2 --outputfile "${i%%.*}".trim.fastq --verbose --threads 24 --minlength 50; done`

v) Next, we count the number of kept sequences from each R1 file:

`for i in *.trim.fastq; do grep '^CGG' "$i" | wc -l; done > kept_seqsR1.txt`

and each R2 file:

`for i in *trim.fastq; do grep '^GATCC' "$i" | wc -l; done > kept_seqsR2.txt`
  

Number of kept sequences was used to calculate the number and percentage of discarded sequences.
  - R1 files: on average 42% were discarded when running the UNedited script; dropped to 13.5% (9%-46%) once the edited/corrected script was used. 
  - R2 files: on average 16% (8%-49%) were discarded once the edited/corrected script was used.
  - see plots below (% discarded vs # total raw reads); most samples had <20% removed samples, which is considered 'OK' ('higher than 20% is unusual').
      

![pct_removed_after_trim](https://user-images.githubusercontent.com/52291277/139146836-1ab68222-1b61-40cd-8346-6298cb1c1356.png)


## 5) Resync trimmed sequences & rename to comply with ddocent naming convension

using the script resync.pl fromthe facility to 're-pair' the pair-end files (R1 and R2) after trimming off the padding sequences.


usage for a single pair of files:

`resync.pl sample1_R1.trim.fastq sample1_R2.trim.fastq sample1_R1.trim.sync.fastq sample1_R2.trim.sync.fastq`

for loop to batch process all 118 pairs (written by Katie):

```
for i in *_R1_001.*; 
do echo "$i"; 
j=`echo "$i" | sed -r 's/'R1'/'R2'/g'`
i2=`echo "$i" | sed -r 's/'Cs_'/''/g'`
i3=`echo "$i2" | sed -r 's/'_S[0-9]+_R1_001'/'.R1'/g'`
j3=`echo "$i3" | sed -r 's/'R1'/'R2'/g'`
perl resync.pl "$i" "$j" "$i3" "$j3"
done
```

 - input: Cs_MA_298_S1_R1_001.trim.fastq
 - output: MA_298.R1.trim.fastq

realized these files are not trimmed for the overhang so the extension we want is .F.fg.gz (forward reads, R1 files) and .R.fg.gz (reverse reads, R2 files), so:

R1 files:

`for f in *R1.trim.fastq; do mv "$f" "${f%.R1.trim.fastq}.F.fq"; done`

R2 files:

`for f in *.R2.trim.fastq; do mv "$f" "${f%.R2.trim.fastq}.R.fq"; done` 

then zip all files in the folder:

`gzip -r folderID`

- input: MA_298.R1.trim.fastq and MA_298.R2.trim.fastq
- output: MA_298.F.fq.gz and MA_298.R.fq.gz

## 6) dDocent

Using a small subset of individuals, **21 MA individuals**, in the folder 'ddocent_subset', I will run ddocent (using module load method).  
      
cd to working folder `cd /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/ddocent_subset`  
load module `module load lotterhos/2020-08-24`  
activate dDocent `source activate ddocent2`  
start running dDocent, type `dDocent` - this initiates the interactive version of dDocent
      
**2021-01-08 - interactive settings:**

Confirm # of individuals - yes  
Choose # of processors - 3  
Limit memory use - 10  
Quality trim? yes  
Perform assembly? yes  
What type of assembly? PE  
New c-parameter? no  
Map reads? yes  
Adjust -A -B -O; new parameters? no  
Use FreeBayes to call SNPs? yes  
Enter email address to get a message when done - will this work on Discovery?

Interactive takes about 20min. 
datacutoff01 =4
datacutoff02 =3

### Run took about ~6h (for only 21 ind!!!) but there was some kind of error - here is a screenshot of the whole run with the error message at the very bottom:


![ddocent_error](https://user-images.githubusercontent.com/52291277/141155240-a4e421a8-b3ba-4f32-a400-3af4526a32f0.png)
                                                                                                                                                                             !
                                                                                                              
vcfcombine: /shared/centos7/salmon/1.1.0/lib/libm.so.6: version `GLIBC_2.15' not found (required by vcfcombine)

All output files are created, but they are empty.

### seff to see time and memory usage:

> (base) [tbittar@d3037 trim_padding]$ seff 21857509  
Job ID: 21857509  
Cluster: discovery  
User/Group: tbittar/users  
State: TIMEOUT (exit code 0)  
Cores: 1  
CPU Utilized: 05:41:24  
CPU Efficiency: 94.74% of 06:00:22 core-walltime  
Job Wall-clock time: 06:00:22  
Memory Utilized: 3.68 GB  
Memory Efficiency: 188.64% of 1.95 GB


notes from Remy:
run a fastqc - we found out adapters are gone.
will try to remove overhang sites with gbstrim.pl & run ddocent with trimmed files .R1 & .R2 


