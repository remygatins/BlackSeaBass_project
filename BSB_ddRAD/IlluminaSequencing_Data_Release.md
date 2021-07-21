## 21-Jul-2021

Hello all,

The GBS fastq data for project Lotterhos_Project_001 (GBS: 122 Black Sea Bass DNA samples - BamHI + MspI) is now available.  Please find the NovaSeq summary and BasicQC report attached.

Notes:
1) Created 118 dual-indexed GBS libraries using enzyme combination BamHI + MspI.
2) Combined all libraries into a single pool and sequenced on a lane of a NovaSeq SP 2x150-bp run.
3) Generated ≥ 375M pass filter reads for the lane.
4) The majority of expected barcodes and samples detected.  There are a handful of libraries with low/no coverage and the library creation process appears to have failed and is likely related to low input mass.
5) Mean quality scores ≥Q30 for all libraries.
6) Configuration of data reads and a script for trimming can be found at the following link: https://bitbucket.org/jgarbe/gbstrim

Data access:
You have 30 days to download your data from our secure download website:

Username: lotterhos
Password: Sent in a separate email
URL: https://umgcdownload.msi.umn.edu/lotterhos/210713_A00223_0599_BH7YTHDRXY/Lotterhos_Project_001

Web browser download (Windows, Mac, and Linux):
-Open the link above in a web browser
-Enter your username and password at the prompt
-Click on the folder links until you find the individual fastq files, and click on each file to download it to your computer

Command-line download (Linux, MSI systems)
-Run the following command, replacing PASSWORD with your password:
wget -nH -np -N -r --cut-dirs 2 --no-check-certificate --user lotterhos --password PASSWORD https://umgcdownload.msi.umn.edu/lotterhos/210713_A00223_0599_BH7YTHDRXY/Lotterhos_Project_001

Globus download (Windows, Mac, and Linux)
-We can share the data with you through Globus, which is the best method for downloading large datasets. Let us know what your Globus email address is, or contact us for information about this data transfer method.

Data download help

Data Analysis
You will find the following analysis results in the project release directory:
illumina-basicQC report

file:///C:/Users/arons/OneDrive%20-%20Northeastern%20University(1)/1.NEU/SeaBass/BSB_ddRAD/Lotterhos_Project_001_BasicQC_Report.html

Data analysis help

Please let us know if you have any problems accessing the data.
---------------------------
Aaron Becker
next-gen@umn.edu
612.624.7427
---------------------------
