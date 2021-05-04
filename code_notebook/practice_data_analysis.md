# Practice ddRAD data analysis using Ddocent on Discovery cluster 

Black Sea Bass practice data using dataset published by Zhao et al Marine Biotechnology 2018 doi 10.1007/s10126-017-9786-0

# 1. Get the data from NCBI

You will need: Accession List and SRA toolkit package

## Get the accession list.

*Option 1.*

Go to ncbi.nlm.nih.gov/bioproject/

In the search bar (top of page) enter the accession number, e.g., PRJNA356786

In the Project Data table, click in the number of links e.g., 221

On top of the page, there is a box “View Results as an Expanded Interactive Table using the Run Selector” – click the link “Send Results to run Selector”.

*Option 2.*

Go to ncbi.nlm.nih.gov/Traces/study/?

Type the accession number on the accession search box and select ‘Runs’.

**Both options 1 and 2 will land at the SRA Run Selector page.**

Under “Found (221) Items” (orange ribbon), select the files you want (or select all, by clicking on the check mark - first column in Found 221 Items).

Under “Select” (blue ribbon), Click on Accession List for Total (if you want all files) or selected (if you want only the selected files under “Found Items”. 

> This will create and download the accession list file, a text file that contains a list of file names that correspond to the SRA files from NCBI that will be downloaded. Save this file in the location from which you are running the SRA Toolkit. You can also download the metadata corresponding to those files.

**Now you have the accession list, you will need to use the SRA toolkit package to download the files.**

## Get the SRA toolkit package.

*Option 1.*

Use your local computer, desktop or laptop.

Download and install the SRA Toolkit for your OS following instructions at: github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

Make sure to install the latest release!

Make sure to use the correct release number and platform, e.g., in step 2 of the Mac OSX install, the correct command would be tar -vxzf sratoolkit.2.10.8.mac64.tar.gz (for the current release on July 2020 for Mac).

Configure SRA toolkit, following instructions at: github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration

Test that the toolkit is functional, following instructions at: github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

*Option 2.*

Use Discovery.

Login to discovery. SRA toolkit is an available module in the Lotterhos cluster.

*Here, I downloaded, installed and configured the SRA toolkit on my desktop and used it to download the SRA files to my desktop, then I uploaded the files to Discovery through Globus (see below). But the  SRA toolkit is available on Discovery so should be able to get it from NCBI directly to Discovery folder.*

# 2. Get the data to Discovery using Globus

If you used option 1 above and your data are in your local computer, or you are getting a hard drive with your data from the sequencing facility, you will need to transfer your data files to the Discovery folder.

You will need: a Discovery account, a Globus account and an endpoint computer.

## Login to Discovery.

In Terminal, login to discovery:

`ssh tbittar@login.discovery.neu.edu`

then enter your password

## Login to Globus.

Open your browser and go to Globus.org, then login or authenticate:

tbittar

Password

> This will get you to “File Manager” in Globus. 

## Get your SRA files into your Discovery folder using Globus.

In File Manager, top right corner, select the 2-panel icon. This is so you can open 2 collections (= folders), 1 to determine the origin folder and 1 to determine the destination folder.

On the first panel, in the box Collection search box, find northeastern#discovery (you might have to login/authenticate), then on the Path box type /~/ or select from you bookmarks northeastern#discovery /~/ (your home folder) – this is your **destination folder**.

On the second panel, in the search box, find your endpoint computer and on the Path box, go to the folder where your files are – this is the **origin folder**.

With the origin and destination folders open in each panel, drag and drop the files or folder that you want to transfer.

**Tip:** Files that need to be kept on backup (such as the original SRA and FASTQ files) should be stored in your home directory. When working on those files, it is best to use the scratch folder than the home folder on Discovery because you have more space there to work when you are creating a lot of large intermediate files. However, the scratch folder is not backed up, so transfer the output files or any other files you want to keep safe to the home directory on Discovery.

**Tip:** Look at “Activity” (left menu) or hit “refresh list” (on the panel corresponding to destination folder) to check progress of transfer.

**Now that all SRA files are in the destination folder, you can close the origin panel, and you will use the SRA toolkit to extract the FASTQ files.**

# 3. Extracting the FASTQ files

You will need: SRA toolkit installed and configured - should already be ready to go if you used it to download the SRA files.

*Option 1.*

You can do that in your local desk- or laptop but this option is not great because these files are large and you will need to work with them on Discovery anyways to proceed with the pipeline. So if possible, option 2 is best – do it directly on Discovery.

*Option 2.*

Using the SRA toolkit in the Lotterhos node.

In terminal, login if you have not yet. 

Check which packages are available in the Lotterhos module: 

`module show lotterhos/2020-07-21`

> This will show you a list of available packages

then load the one you need to use:

`module load sratoolkit/2.10.8`

If the toolkit has not been configured, you will get a message in terminal to do so. Follow the instructions to configure (they are the same instructions as the one you followed to configure the toolkit on your local computer) – github.com/ncbi/sra-tools/wiki/03.-Qiocl-Toolkit-Configuration

Once configured, go to terminal and type: 

`Fasterq-dump –split-files filename.sra`

To extract one file. Or to extract all sra files, type:

`Fasterq-dump –split-files \*.sra`

> This will create two FASTQ files for each SRA file (because the sequencing was pair-ended, so you get the forward and the reverse reads in two separate FASTQ files, labeled FASTQ1 and FASTQ2). For each file name, you will end up with 3 files: the 'original' SRA file, one .sra_1.fastq and one .sra_2.fastq file.

**Now that we have the FASTQ files, we will start using the pipeline described in Ddocent to process the data.**

As mentioned before, it is best to work from the scratch folder. So copy the fastq files from the home folder to the scratch folder. In Globus > File Manager, open two panels again:

one for the home folder where the fastq files are located: click the search box and select northeastern#discovery (you might need to re-authenticate) and in the path box type /~/ (the home folder), then navigate to where fastq files are; 

and one for your scratch folder: click the search box and select northeastern#discovery (you might need to re-authenticate) and in the path box, type /scratch/tbittar/ (the scratch folder), then naviate to where you will copy the fastq files to.

Drag and drop the fastq files from home to scratch.

**Tip**: Bookmark these folders to make it easier to find them in the future in Globus > Bookmarks list (left side panel).

# 4. Use dDocent pipeline on Discovery

> Note: Here we are assuming that the barcoding and adapters are already stripped off for this specific practice dataset. Since it was a ddRAD library prep, we should be left with the DNA sequence of interest and the overhang sequence corresponding to the recognition sequence where the enzymes attached. In this study, the two enzymes they used were EcoRI and MspI. “EcoRI is a restriction enzyme that creates 4 nucleotide sticky ends with 5' end overhangs of AATT. The nucleic acid recognition sequence where the enzyme cuts is G/AATTC, which has a palindromic, complementary sequence of CTTAA/G. The / in the sequence indicates which phosphodiester bond the enzyme will break in the DNA molecule.” MspI overhang sequence is C/CGG and the palindromic complementary sequence is GGC/C. 

> Note: Actually, rather than assuming, this can be checked by looking at (some of) the sequence files. All the random sra_1.fastq files that I looked at start with the overhang sequence AATTC - meaning they correspond to the forward reads, cut with EcoRI - and the sra_2.fastq files start with the overhang CGG - meaning they correspond to the reverse reads, cut with MspI.

> Note: In the Lotterhos lab module, ddocent-2.7.8 was added to the conda environment “lotterhos-py38” within the anaconda3/L2020-03 module.

## The files need to be in the following format: PopID_filename.F.fq.gz and PopID_filename.R.fq.gz

# 4.1 Rename files to comply with the above naming convention required:

> My fastq files came out of NCBI in the following format: SRRxxxxxxxx.sra_1.fastq (forward reads) and SRRxxxxxxxx.sra_2.fastq (reverse reads).

> They need to be in the following format: PopID_SRRxxxxxxxx.F.fq.gz (forward reads) and PopID_SRRxxxxxxxx.R.fq.gz (reverse reads).

> Note: We do not have information on the Pop ID in this study so we will call them all Pop1. 

> Note: gz means the file is compressed.

You will need: request resources in Discovery and a code to batch-rename the files.

## Request Discovery resources

In terminal, login to Discovery and navigate to where your working fastq files are:

`cd /scratch/tbittar/`

Type one fo the the following codes to get access to resources: 

`srun -p debug -N1 --pty /bin/bash`

> This will give you 20 mins to work. Retype when needed. Use this to check your code is working on 1-2 files.

`srun -p lotterhos -N 1 --pty /bin/bash`

> This will give you 24h to work. Use this if 20 min is not enough to rename and compress all files. 

*To use Lotterhos Lab resources, you need to ask Dr. L to be added as a user ahead of time. To check if you are a user type `groups $USER`. If you are successfully added to the group, you will see:

tbittar : users lotterhos 

## Script to batch-rename files (I could only come up with a two-step code (two for-loops) to batch-rename the files):

**1st step** - this will replace the extension .sra_1.fastq with an .F and the extension .sra_2.fastq with an .R in all files in the folder:

for f in *.sra_1.fastq; do 

mv "$f" "${f%.sra_1.fastq}.F"; 

done

for f in *.sra_2.fastq; do 

mv "$f" "${f%.sra_2.fastq}.R"; 

done

**2nd step** - this will add Pop1_ as a prefix to all file names; compress the file with gzip; and add the .fq as the extension in all files in the folder:

prefix=Pop1_

for name in SRR*; do

gzip < "$name" > "$prefix$name.fq.gz" 

done 

**Tip:** Keep the scratch folder open in Globus and use Refresh List to double-check that the file names are changing as expected and that the new files are smaller in size (because they are now compressed).

# 4.2 Load the module where dDocent is nested.

Login to Discovery and request resources using the lotterhos partition.

`srun -p lotterhos -N 1 --pty /bin/bash`

Load the module where Ddocent is nested and activate dDocent:

`module load lotterhos/2020-08-24`

`source activate ddocent2`

Start dDocent, type:

`dDocent`

**The interactive part:**

Confirm # of individuals - yes/no

Choose # of processors - 3

Limit memory use - 1

Quality trim? yes

Perform assembly? yes

What type of assembly? PE

New c-parameter? no

Map reads? yes

Adjust -A -B -O; new parameters? no

Use FreeBayes to call SNPs? yes

Enter email address to get a message when done - will this work on Discovery?

If all works, you get a plot of Unique Sequences vs Coverage; based on this graph, choose the data cutoff.

You get another plot of Unique Sequences vs Individuals, choose the data cutoff.

From here on, the pipeline should run on its own.

# Attempts at running the pipeline.

---
## 16-Sept-2020

**Attempt #1 - all settings as listed above under "4.2 the interactive part":**

Running dDocent on all 56 files (28 individuals) for which I have fastq files.
* cutoff chosen for unique seq vs coverage = 6
* cutoff chosen for unique seq vs #individuals = 6

ERROR: EXCEEDED JOB MEMORY LIMIT AT BWA TO MAP READS. 
 
**Attempt #2 - same as above, but changed "Limit memory use" to 10Gb:**

SAME ERROR: EXCEEDED JOB MEMORY LIMIT AT BWA TO MAP READS. 

**Attempt #3 - same as above, but changed "Limit memory use" to 0GB (based on dDocent user guide):**

SAME ERROR: EXCEEDED JOB MEMORY LIMIT AT BWA TO MAP READS. 

---
## 17-Sept-2020

Meeting with Katie - use seff to check memory usage and run less files to test the pipeline.

Running `seff JOBID` on yesterday's jobs tells me the job is using ~2Gb of memory.

**Attempt #4, same settings as above, except:**
* 10 files/5 individuals
* Choose # of processors - 5
* Limit memory use - 10Gb
* cutoff chosen for unique seq vs coverage = 6
* cutoff chosen for unique seq vs #individuals = 3


SAME ERROR: EXCEEDED JOB MEMORY LIMIT AT BWA TO MAP READS. 
Seff: tells me the job is using ~2Gb of memory.

**Attempt #5, same settings as above, except:**
* 4 files/2 individuals
* Choose # of processors - 3
* Limit memory use - 0Gb
* cutoff chosen for unique seq vs coverage = 6
* cutoff chosen for unique seq vs #individuals = 2

THIS WORKS and the pipeline runs all the way.
Seff: tells me memory utilized was 0Gb

---
## 18-Sept-2020

From here on, I'm incrementally increasing the number of files/individuals to see if/when I get the memory error.

**Attempt #6, same settings as above, except:**
* 8 files/4 individuals
* Choose # of processors - 3
* Limit memory use - 0Gb
* cutoff chosen for unique seq vs coverage = 6
* cutoff chosen for unique seq vs #individuals = 3

THIS WORKS and the pipeline runs all the way.
Seff: tells me memory utilized was 0Gb

**Attempt #7, same settings as above, except:**
* 16 files/8 individuals
* Choose # of processors - 3
* Limit memory use - 0Gb
* cutoff chosen for unique seq vs coverage = 6
* cutoff chosen for unique seq vs #individuals = 5

THIS WORKS and the pipeline runs all the way.
Seff: tells me memory utilized was 0Gb

**Attempt #8, same settings as above, except:**
* 32 files/16 individuals
* Choose # of processors - 3
* Limit memory use - 0Gb
* cutoff chosen for unique seq vs coverage = 6
* cutoff chosen for unique seq vs #individuals = 5

ERROR: EXCEEDED JOB MEMORY LIMIT AT BWA TO MAP READS.
Seff: tells me memory utilized was 2.17Gb

---
## 21-Sept-2020

**Attempt #9, same settings as above, except:**
* 20 files/10 individuals
* Choose # of processors - 3
* Limit memory use - 0Gb
* cutoff chosen for unique seq vs coverage = 6
* cutoff chosen for unique seq vs #individuals = 5

ERROR: EXCEEDED JOB MEMORY LIMIT AT BWA TO MAP READS.
Seff: tells me memory utilized was 1.97Gb out of 1.95Gb

**Attempt #10, same settings as above, except:**
* 18 files/9 individuals
* Choose # of processors - 3
* Limit memory use - 0Gb
* cutoff chosen for unique seq vs coverage = 6
* cutoff chosen for unique seq vs #individuals = 5

THIS WORKS and the pipeline runs all the way.
Seff: tells me memory utilized was 0Gb.

---

**Attempts 5-7 worked and the pipeline ran all the way through!!!!!!!!!!!!! :) - but with limited number of files**

Based on the 10 attempts listed above, it looks like the memory available in partition can handle 18 files (9 individuals) 
* how much data does that correspond to?

**There is a 2Gb memory limit in the partition, - or maybe it is a glitch and we should be able to request more memory?**

**Next steps: create a settings file for the interactive part of the pipeline and try to submit as a job so more memory is available.**

---
