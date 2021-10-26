## 1. Get the files

- download FASTQ files from sequencing facility using:

```wget -nH -np -N -r --cut-dirs 2 --no-check-certificate --user lotterhos --password Ride-Stir-Secretary-Size-7 https://umgcdownload.msi.umn.edu/lotterhos/210713_A00223_0599_BH7YTHDRXY/Lotterhos_Project_001```

- check download was sucessfull using md5sum according to https://drk-lo.github.io/lotterhoslabprotocols/code_checksums/ - PASS

## 2. Enzyme pair used here: BamHI and Mspl

BamHI cut: sticky ends with 5' end overhangs of G/GATCC and palindromic complementary CCTAG/G

Mspl cut: sticky ends with 5' end overhangs of C/CGG and palindromic complementary GGC/C  

Barcode sequences: see project summary file.

## 3. name convention for dDocent

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

iii) use their script: 

`perl gbstrim.pl --enzyme1 mspi --enzyme2 bamhi --fastqfile sample1_R1.fastq.gz --read R1 --outputfile sample1.trim.fastq --verbose --threads 24 --minlength 50`

NOTES on this line: add `perl` before the command; change the enzymes to what we used (mspi & bamhi); change input and output file names to match ours.

iv) batch trim all 118 **R1** files:

`for i in *.fastq.gz; do perl gbstrim.pl --enzyme1 mspi --enzyme2 bamhi --fastqfile "$i" --read R1 --outputfile "${i%%.*}".trim.fastq --verbose --threads 24 --minlength 50; done`

Padding sequences (green, variable lengths), R1 overhang (yellow, CGG):
![ddRAD_R1_padding_example](https://user-images.githubusercontent.com/52291277/138742166-121f016b-6cb7-402e-9182-d77c39dbd8db.png)

After padding trimmed out, R1 overhangs remain:
![ddRAD_R1_padding_cut](https://user-images.githubusercontent.com/52291277/138744664-ecaf076b-573d-4d08-8982-b47d613381b1.png)



and 118 **R2** files - DO NOT USE THE SCRIPT/LINE BELOW, see next hurdle.

`for i in *.fastq.gz; do perl gbstrim.pl --enzyme1 mspi --enzyme2 bamhi --fastqfile "$i" --read R2 --outputfile "${i%%.*}".trim.fastq --verbose --threads 24 --minlength 50; done`

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
- After running the edited script, I'm getting ~17-19% discarded sequences on a few files ran manually, so will move on with the edited script.
------------------------------------------------------

So for the 118 **R2** files we are using the edited script (edited in my local computer and transferred to the working folder on Discovery through Globus).

`for i in *.fastq.gz; do perl gbstrimedited.pl --enzyme1 mspi --enzyme2 bamhi --fastqfile "$i" --read R2 --outputfile "${i%%.*}".trim.fastq --verbose --threads 24 --minlength 50; done`


v) Next, we count the number of kept sequences from each R1 file:

`for i in *.trim.fastq; do grep '^CGG' "$i" | wc -l; done > kept_seqsR1.txt`

and each R2 file:

`for i in *trim.fastq; do grep '^GATCC' "$i" | wc -l; done > kept_seqsR2.txt`
  

Number of kept sequences was used to calculate the number and percentage of discarded sequences.
  - R1 files: on average 42% were discarded; percentage of discarded sequences a little too high; it correlates well with "percentage of adapter" metrics. Percentage of adapter is defined as "The cumulative proportion of each sample in which sequencing adapter sequences have been seen at each position. Once an adapter sequence has been seen in a read it is counted as being present right through to the end of the read so the percentages only increase as the read length goes on. It is common to see significant adapter sequence content at the ends of reads in short insert libraries (16s/18s, small RNA, amplicon). This metric is calculated by FastQC." (see plot below).
  - R2 files: on average 16% were discarded - this is low enough, I don't think we need to investigate.
      
![pct_discarded_R1_correl](https://user-images.githubusercontent.com/52291277/138901755-3611d40e-64e0-4ddb-9567-22cf66522749.png)


