# Discovery cluster 

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

In Terminal, login to discovery:

ssh tbittar@login.discovery.neu.edu 

then enter your password

Get your SRA files into your Discovery folder using Globus.

Open your browser and go to Globus.org, then login or authenticate:

Tbittar

Password

> This will get you to “File Manager” in Globus. 

In File Manager, top right corner, select the 2-panel icon. This is so you can open 2 collections (= folders), 1 to determine the origin folder and 1 to determine the destination folder.

On the first panel, in the box Collection search box, find northeastern#discovery (you might have to login/authenticate), then on the Path box type or select /scratch/tbittar (or your own scratch folder) – this is your **destination folder**.

On the second panel, in the search box, find your endpoint computer and on the Path box, go to the folder where your files are – this is the **origin folder**.

With the origin and destination folders open in each panel, drag and drop the files or folder that you want to transfer.

**Tip:** it is best to use the scratch folder than the home folder on Discovery because you have more space there to work when you are creating a lot of large intermediate files. However, the scratch folder is not backed up, so transfer the output files or any other files you want to keep safe to the home directory on Discovery.

**Tip:** Look at “Activity” (left menu) or hit “refresh list” (on the panel corresponding to destination folder) to check progress of transfer.

**Now that all SRA files are at the destination folder, you can close the origin panel, and you will use the SRA toolkit to extract the FASTQ files.**

# 3. Extracting the FASTQ files

*Option 1.*

You can do that in your local desk- or laptop but this option is not great because these files are large and you will need to work with them on Discovery anyways to proceed with the pipeline. So if possible, option 2 is best – do it directly on Discovery.

*Option 2.*

Using the SRA toolkit in the Lotterhos node.

In terminal, login if you have not yet. 

Check which packages are available in the Lotterhos module: 

module show lotterhos/2020-07-21 

> This will show you a list of available packages

then load the one you need to use:

module load sratoolkit/2.10.8

If the toolkit has not been configured, you will get a message in terminal to do so. Follow the instructions to configure (they are the same instructions as the one you followed to configure the toolkit on your local computer) – github.com/ncbi/sra-tools/wiki/03.-Qiocl-Toolkit-Configuration

Once configured, go to terminal and type: 

Fasterq-dump –split-files filename.sra

To extract one file. Or to extract all sra files, type:

Fasterq-dump –split-files \*.sra

> This will create two FASTQ files for each SRA file (because the sequencing was pair-ended, so you get the forward and the reverse reads in two separate FASTQ files, labeled FASTQ1 and FASTQ2).

**Now that we have the FASTQ files, we will start using the pipeline described in Ddocent to process the data.**

# 4. Use Ddocent pipeline on discovery

In the Lotterhos lab module, ddocent-2.7.8 was added to the conda environment “lotterhos-py38” within the anaconda3/L2020-03 module.

In terminal, login to Discovery if you haven’t yet.

Load the module where Ddocent is nested:

Module load anaconda3/L2020-03

> Here we are assuming that the barcoding and adapters are already stripped off for this specific practice dataset. Since it was a ddRAD library prep, we should be left with the DNA sequence of interest and the overhang sequence corresponding to the recognition sequence where the enzymes attached. In this study, the two enzymes they used were EcoRI and MspI. “EcoRI is a restriction enzyme that creates 4 nucleotide sticky ends with 5' end overhangs of AATT. The nucleic acid recognition sequence where the enzyme cuts is G/AATTC, which has a palindromic, complementary sequence of CTTAA/G. The / in the sequence indicates which phosphodiester bond the enzyme will break in the DNA molecule.” MspI overhang sequence is C/CGG and the palindromic complementary sequence is GGC/C.
