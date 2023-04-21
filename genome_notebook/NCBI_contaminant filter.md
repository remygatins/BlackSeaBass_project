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


# Remove contaminants and Mitochondrial genome

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





