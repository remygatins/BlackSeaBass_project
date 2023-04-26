# SNMF on LEA to get structure plot


## Convert VCF to Ped file
run interactive mode

            srun -p lotterhos -N 1 --pty /bin/bash

1. Create chromosome map

```bash
module load miniconda3
conda activate plink
conda activate bcftools

bcftools view -H populations.snps.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > populations.snps.vcf.chrom-map.txt
```

2. Make ped file using this chromosome map

```bash
module load vcftools

vcftools --vcf populations.snps.vcf --out BSB_r0.8_maf0.01 --plink --chrom-map populations.snps.vcf.chrom-map.txt
```

```bash
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf populations.snps.vcf
	--chrom-map populations.snps.vcf.chrom-map.txt
	--out BSB_r0.8_maf0.01
	--plink

After filtering, kept 117 out of 117 Individuals
Writing PLINK PED and MAP files ...
	Read 39 chromosome mapping file entries.
Done.
After filtering, kept 24896 out of a possible 24896 Sites
Run Time = 1.00 seconds
```
This should have created a .ped and .map file. Now download this to your local computer or wherever you will run R from

## Run SNMF on LEA (input file = .ped)

In R, install LEA and convert from a .ped file to a geno file to run SNMF (We will first need to convert to a lfmm file and then to a geno file)

http://membres-timc.imag.fr/Olivier.Francois/LEA/tutorial.html

