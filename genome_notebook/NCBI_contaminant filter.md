'Contamination.txt' from Genbank after submitting 'C_striata_01.fasta' which is also 'tom-rev2777-mb-hirise-nxqj6__03-05-2022__hic_output.fasta' in my files 

```{}
SUBID     	BioProject	BioSample	Organism
--------------------------------------------------------
SUB12274227	PRJNA900280	SAMN31685414	Centropristis striata

[] We ran your sequences through our Contamination Screen. The screen found 
contigs that need to be trimmed and/or excluded. The results are in the 
Contamination.txt file posted in your submission on the WGS submission portal 
https://submit.ncbi.nlm.nih.gov/subs/genome/.  

GenBank staff will automatically remove contaminants that are found to be 
the entire sequence or at the end of a sequence, and will post the reports 
and edited fasta file to the submission portal. Note that internal contamination 
will not be automatically removed since the sequence may be misassembled and 
therefore should be split at the contamination and resubmitted as separate sequences.
In addition, we do not automatically remove mitochondrial sequences in 
eukaryotic submissions. 

If you selected the submission portal option "Do not automatically trim or 
remove sequences identified as contamination" then you will need 
to adjust the sequences appropriately and then resubmit your sequences. 
After you remove the contamination, trim any Ns at the ends of the sequence 
and remove any sequences that are shorter than 200 nt and not part of a 
multi-component scaffold.

Note that mismatches between the name of the adaptor/primer identified in the screen 
and the sequencing technology used to generate the sequencing data should not be used 
to discount the validity of the screen results as the adaptors/primers of many 
different sequencing platforms share sequence similarity.


Adaptor:
[] Some of the sequences hit primers or adaptors used in Illumina, 
454, or other sequencing strategies or platforms.  Adaptor at the 
end of a sequence should be removed. However, if adaptors are 
present within sequences then you should strongly consider 
splitting the sequences at the adaptor because the adaptor sequence 
could have been the region of overlap, causing a misassembly.


Skipped 429 same as before; no new sequences to screen.
Note: 26 sequences with runs of Ns 10 bp or longer (or those longer that 20 MB) were split before comparing.
6 sequences with locations to mask/trim
(7 split spans with locations to mask/trim)

Trim:
Sequence name, length, span(s), apparent source
ScOpHrR_19;HRSCAF=179	38455635	29300242..29300329	adaptor:NGB00972.1-not_cleaned
ScOpHrR_1;HRSCAF=3	37471414	1820840..1820927	adaptor:NGB00972.1-not_cleaned
ScOpHrR_20;HRSCAF=182	41338807	16160790..16160872,31829880..31829924	adaptor:NGB00972.1-not_cleaned
ScOpHrR_26;HRSCAF=235	35856960	29908858..29908891	adaptor:NGB00972.1-not_cleaned
ScOpHrR_3;HRSCAF=56	38835999	30900595..30900682	adaptor:NGB00972.1-not_cleaned
ScOpHrR_5;HRSCAF=79	43113489	1836749..1836833	adaptor:NGB00972.1-not_cleaned
```


# Remove Mitochondrial genome from assembly

Downloaded the [Centropristis striata mitochondrial genome](https://www.ncbi.nlm.nih.gov/nuccore/MH645336.1/) from NCBI and stored as `C_striata_mito.fasta` in `/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome`

run interactive node
  `srun -p lotterhos -N 1 --pty /bin/bash`
  
load minimap
```{}
module load minimap2/2.17

minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]

minimap2 -a C_striata_mito.fasta C_striata_01.fasta > C_striata_mito_map.sam
```

get a summary text 

    cat C_striata_mito_map.sam | grep ^@* | less -S > C_striata_mito_map_summary.txt



# Remove adapter contaminants 

|Sequence name| length| span(s)| apparent source| bp's to remove|
|-------------|-------|--------|----------------|---------------|
|ScOpHrR_19;HRSCAF=179|	38455635	|29300242..29300329|	adaptor:NGB00972.1-not_cleaned|87 |
|ScOpHrR_1;HRSCAF=3|	37471414	|1820840..1820927|	adaptor:NGB00972.1-not_cleaned|87|
|ScOpHrR_20;HRSCAF=182|	41338807	|16160790..16160872,31829880..31829924|	adaptor:NGB00972.1-not_cleaned|82,44|
|ScOpHrR_26;HRSCAF=235|	35856960	|29908858..29908891|	adaptor:NGB00972.1-not_cleaned|33|
|ScOpHrR_3;HRSCAF=56|	38835999	|30900595..30900682	|adaptor:NGB00972.1-not_cleaned|87|
|ScOpHrR_5;HRSCAF=79|	43113489	|1836749..1836833	|adaptor:NGB00972.1-not_cleaned|84|


`cat C_striata_01.fasta | grep ScOpHrR_1 -n`

1:>ScOpHrR_1;HRSCAF=3

cat C_striata_01.fasta | sed -n '1p'

Produce a file with sequence length

```
cat C_striata_01.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > C_striata_seq_length.txt
```
output
```
ScOpHrR_1;HRSCAF=3      37471414
ScOpHrR_2;HRSCAF=47     309108
ScOpHrR_3;HRSCAF=56     38835999
ScOpHrR_4;HRSCAF=74     660109
ScOpHrR_5;HRSCAF=79     43113489
ScOpHrR_6;HRSCAF=101    32175207
ScOpHrR_7;HRSCAF=104    32032266
ScOpHrR_8;HRSCAF=108    431186
ScOpHrR_9;HRSCAF=120    43614
ScOpHrR_10;HRSCAF=121   36376519
ScOpHrR_11;HRSCAF=135   1084579
ScOpHrR_12;HRSCAF=138   30413
ScOpHrR_13;HRSCAF=141   42410642
ScOpHrR_14;HRSCAF=148   129042
ScOpHrR_15;HRSCAF=151   46365755
ScOpHrR_16;HRSCAF=162   19756
ScOpHrR_17;HRSCAF=164   36354236
ScOpHrR_18;HRSCAF=168   45437359
...
```
length's match those reported from ncbi 

Now create a test fasta file

    cat C_striata_01.fasta | sed -n '1,+2p' > test_1.fasta
    
get line number for second sequence

    cat C_striata_01.fasta | grep ScOpHrR_2 -n
    
output

    374717:>ScOpHrR_2;HRSCAF=47
    
output sequence 2 plus  2 lines 

    cat C_striata_01.fasta | sed -n '374717,+2p' > test_2.fasta
    
concatenate file

    cat test_1.fasta test_2.fasta > test.fasta

    cat test.fasta
```
>ScOpHrR_1;HRSCAF=3
ACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC
TAACCCTAACCCTACCTACCCTAACCCTAACCCTAACCCTAACCCACCCTAACCCTAACCCTAACCCTAACCCTACCTACCCACCCTAACCCTAACCCTA
>ScOpHrR_2;HRSCAF=47
TTTCAAAGTTCTGGAAAAGTCCCATGTTGCATTGCAACACTCCCACATGTATTCTCATTTTGTGTCTTTTTGTGTCTTTTTTTAGTCATTTTCATTTTGT
GTCTTTTTTTGGTCATTTTTTATCTTTTTTGGTCATTTTGTGTCTTTTTTAGTCATTTTGCATCTTTTTTTGGTCATTTTGTGTCTTTATGTATCTTTTT
```

Create a bed file with coordinants of contaminants


    cat contaminants.bed
```
ScOpHrR_19;HRSCAF=179	29300242	29300329
ScOpHrR_1;HRSCAF=3	1820840	1820927
ScOpHrR_20;HRSCAF=182	16160790	16160872
ScOpHrR_20;HRSCAF=182	31829880	31829924
ScOpHrR_26;HRSCAF=235	29908858	29908891
ScOpHrR_3;HRSCAF=56	30900595	30900682
ScOpHrR_5;HRSCAF=79	1836749	1836833
```
obtain contaminant sequences 
```
bedtools getfasta -fi C_striata_01.fasta -bed contaminants.bed -fo contaminants.out
```

    cat contaminants.out
    
```unix
>ScOpHrR_19;HRSCAF=179:29300242-29300329
CTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGATTCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT
>ScOpHrR_1;HRSCAF=3:1820840-1820927
CTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGATTCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT
>ScOpHrR_20;HRSCAF=182:16160790-16160872
TCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGATTCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT
>ScOpHrR_20;HRSCAF=182:31829880-31829924
CTCTCTCTTTTCCTCCTCCTCCGTTTGTTGTTGTTGAGAGAGAT
>ScOpHrR_26;HRSCAF=235:29908858-29908891
TCTCTCTCAACAACACAACGGAGGAGGAGGAAA
>ScOpHrR_3;HRSCAF=56:30900595-30900682
TCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGA
>ScOpHrR_5;HRSCAF=79:1836749-1836833
TCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAATCTCTCTCTTTTCCTCTCCTCCGTTGTTGTTGTTGAGAGAGA
```

This is how it looks when ligned up with eachother

![image](https://user-images.githubusercontent.com/26288352/234093752-17c7216d-d5cf-45f4-a222-2bd2ee6037de.png)

I added `contaminants_extra.bed` to see bases before and after

```
bedtools getfasta -fi C_striata_01.fasta -bed contaminants_extra.bed -fo contaminants_extra.out
```
![image](https://user-images.githubusercontent.com/26288352/234101304-6fa03e6a-e70f-49b4-a3ae-95cad5e71c35.png)

In red are bases identified as contaminants. Marked in green are differences between the sequences.

using Geneious Prime

<img width="930" alt="image" src="https://user-images.githubusercontent.com/26288352/234396497-adf44f68-6a92-4e6e-9227-292adbfa9304.png">


Map contaminant sequences to genome assembly to check if they map elsewhere in the genome and make sure we will not remove extra sequences.

run interactive node
  `srun -p short -N 1 --pty /bin/bash`
    
    minimap2 -a C_striata_01.fasta contaminants.fasta > contaminants_C_striata_01.sam
    

In Geneious I manually replaced contaminants with N's

```
>ScOpHrR_1;HRSCAF=3:1820840-1820927	
CTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGATTCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT
Geneious: 1820841-1820927- Edited
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

>ScOpHrR_20;HRSCAF=182:16160790-16160872		
TCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGATTCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT
Geneious: 16160791-16160872
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

>ScOpHrR_20;HRSCAF=182:31829880-31829924		
CTCTCTCTTTTCCTCCTCCTCCGTTTGTTGTTGTTGAGAGAGAT
Geneious: 31829881-31829924
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

>ScOpHrR_26;HRSCAF=235:29908858-29908891		
TCTCTCTCAACAACACAACGGAGGAGGAGGAAA
Geneious: 29908859-29908891
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

>ScOpHrR_3;HRSCAF=56:30900595-30900682			
TCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGA
Geneious: 30900596-30900682
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

>ScOpHrR_5;HRSCAF=79:1836749-1836833			
TCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAATCTCTCTCTTTTCCTCTCCTCCGTTGTTGTTGTTGAGAGAGA
Geneious: 1836750-1836833
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
````
*sequences searched in Geneious matched bases starting one base pair after. 

I saved and exported the contaminant free assembly as `C_striata_v1.2.fasta`


    stats for C_striata_01.fasta
    
```
sum = 926021963, n = 92, ave = 10065456.12, largest = 46365755
N50 = 39023748, n = 11
N60 = 37471414, n = 14
N70 = 36432694, n = 16
N80 = 35856960, n = 19
N90 = 34200995, n = 21
N100 = 8430, n = 92
N_count = 32476
Gaps = 334
```

    assembly-stats C_striata_v1.2.fasta
    
```
stats for C_striata_v1.2.fasta
sum = 926021963, n = 92, ave = 10065456.12, largest = 46365755
N50 = 39023748, n = 11
N60 = 37471414, n = 14
N70 = 36432694, n = 16
N80 = 35856960, n = 19
N90 = 34200995, n = 21
N100 = 8430, n = 92
N_count = 32980
Gaps = 341
```

Uploaded to genbank but it detedcted another contaminant

```bash
Skipped 436 same as before; no new sequences to screen.
Note: 26 sequences with runs of Ns 10 bp or longer (or those longer that 20 MB) were split before comparing.
1 sequence with locations to mask/trim
(1 split span with locations to mask/trim)

Trim:
Sequence name, length, span(s), apparent source
ScOpHrR_20;HRSCAF=182	41338807	31829834..31829880	adaptor:NGB00972.1-not_cleaned
```
I converted those to N's with Geneious and saved as file C_striata_v1.3.fasta




## Order Scaffolds by length and rename 

I would now like to order my contigs from largest to smallest and rename contigs to a sequential order.

Right now my fasta file looks like this:

`cat C_striata_v1.3.fasta | grep ^">"| head`

```bash
>ScOpHrR_1;HRSCAF=3
>ScOpHrR_2;HRSCAF=47
>ScOpHrR_3;HRSCAF=56
>ScOpHrR_4;HRSCAF=74
>ScOpHrR_5;HRSCAF=79
>ScOpHrR_6;HRSCAF=101
>ScOpHrR_7;HRSCAF=104
>ScOpHrR_8;HRSCAF=108
>ScOpHrR_9;HRSCAF=120
>ScOpHrR_10;HRSCAF=121
```

I will use `seqkit` to order contigs by length and rename

```bash
conda activate seqkit

#order contigs by length
seqkit sort --by-length --reverse C_striata_v1.3.fasta > C_striata_v1.3_ordered.fasta

#rename contigs
seqkit replace --pattern '.+' --replacement 'Scaffold_{nr}' C_striata_v1.3_ordered.fasta > C_striata_v1.3_ordered_renamed.fasta
```


- where `{nr}` is the contig of record
- `.+` means any character and more than one instance

After ordering
`cat C_striata_v1.3_ordered.fasta | grep ^">"| head`

```bash
>ScOpHrR_15;HRSCAF=151
>ScOpHrR_37;HRSCAF=290
>ScOpHrR_18;HRSCAF=168
>ScOpHrR_5;HRSCAF=79
>ScOpHrR_13;HRSCAF=141
>ScOpHrR_33;HRSCAF=273
>ScOpHrR_20;HRSCAF=182
>ScOpHrR_35;HRSCAF=285
>ScOpHrR_23;HRSCAF=205
>ScOpHrR_56;HRSCAF=349
```

Final output
`cat C_striata_v1.3_ordered_renamed.fasta | grep ^">"| head`

```bash
>Scaffold_1
>Scaffold_2
>Scaffold_3
>Scaffold_4
>Scaffold_5
>Scaffold_6
>Scaffold_7
>Scaffold_8
>Scaffold_9
>Scaffold_10
```

The final genome assembly `C_striata_v1.3_ordered_renamed.fasta` will be renamed to `C_striata_v2.fasta`

`cp C_striata_v1.3_ordered_renamed.fasta C_striata_v2.fasta`

