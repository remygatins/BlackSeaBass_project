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
     - conda create -n cutadaptenv cutadapt (environment location: /home/tbittar/miniconda3/envs/cutadaptenv)
     - conda activate cutadaptenv (activating the newly crated env)
     - cutadapt --version (just to check it worked)

ii) Download the repository from https://bitbucket.org/jgarbe/gbstrim/src/master/ 

- used Globus to transfer the repo to my working folder on Discovery (working folder = /work/lotterhos/NOAA...)

iii) use their command: 

`perl gbstrim.pl --enzyme1 mspi --enzyme2 bamhi --fastqfile sample1_R1.fastq.gz --read R1 --outputfile sample1.trim.fastq --verbose --threads 24 --minlength 50`

NOTES on this line: add `perl` before the command; change the enzymes to what we used (mspi & bamhi); change input and output file names to match ours.

iv) batch trim all 118 files:

`for i in *.fastq.gz; do perl gbstrim.pl --enzyme1 mspi --enzyme2 bamhi --fastqfile "$i" --read R1 --outputfile "${i%%.*}".trim.fastq --verbose --threads 24 --minlength 50; done`




           

