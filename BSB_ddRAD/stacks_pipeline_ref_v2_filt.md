# BSB RADseq STACKS pipeline with an updated reference genome `C_striata_v2.fasta` and new filtered reads

### Paths:

Working directory: `/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2_filt`
or using the shortcut: `/home/r.gatins/BSB_ddRAD/stacks_ref_v2_filt`

Reference genome:   `/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_v2.fasta`

Demultiplexed sequences: `/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed`


## Set up directory
Within the working directory lets create the following subdirectories

```
mkdir 00.popmap 01.jobs 02.Trimmomatic_filtering 03.Assembly 04.bam_alignements 05.Stacks 06.Population_genetics
```

## Demultiplex
Samples have already been demultiplexed previously (see Thais'[pipeline.md](https://github.com/thais-neu/BlackSeaBass_project/blob/master/BSB_ddRAD/pipeline.md))

Demultiplexed sequences: `/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed`

## Trim

Remove the first 5 bp and 3 bp from the forward and reverse reads to remove the enzyme cut site, paired-end reads.

```



