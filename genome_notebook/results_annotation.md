
# Methods
Repeat families found in the genome assemblies of Centropristis striata were identified de novo and classified using the software package RepeatModeler (version 2.0.1). RepeatModeler depends on the programs RECON (version  1.08) and RepeatScout (version 1.0.6) for the de novo identification of repeats within the genome. The custom repeat library obtained from RepeatModeler were used to discover, identify and mask the repeats in the assembly file using RepeatMasker (Version 4.1.0). Coding sequences from Amphiprion ocellaris, Epinephelus lanceolatus, Gadus morhua, Morone saxatilis and Plectropomus leopardus were used to train the initial ab initio model for Centropristis striata using the AUGUSTUS software (version 2.5.5). Six rounds of prediction optimisation were done with the software package provided by AUGUSTUS. The same coding sequences were also used to train a separate ab initio model for Centropristis striata using SNAP (version 2006-07-28). RNAseq reads were mapped onto the genome using the STAR aligner software (version 2.7) and intron hints generated with the bam2hints tools within the AUGUSTUS software. MAKER, SNAP and AUGUSTUS (with intron-exon boundary hints provided from RNA-Seq) were then used to predict for genes in the repeat-masked reference genome. To help guide the prediction process, Swiss-Prot peptide sequences from the UniProt database were downloaded and used in conjunction with the protein sequences from Amphiprion ocellaris, Epinephelus lanceolatus, Gadus morhua, Morone saxatilis and Plectropomus leopardus to generate peptide evidence in the Maker pipeline. Only genes that were predicted by both SNAP and AUGUSTUS softwares were retained in the final gene sets. To help assess the quality of the gene prediction, AED scores were generated for each of the predicted genes as part of the MAKER pipeline. Genes were further characterised for their putative function by performing a BLAST search of the peptide sequences against the UniProt database. tRNA were predicted using the software tRNAscan-SE (version 2.05).


# Delivery
The delivery package for your genome annotation contains:  
- repeat soft-masked genome assembly file 
- genome annotation file in GFF3 format 
- genes CDS sequences in fasta format 
- genes peptide sequence in fasta format 
- bam alignment files of your transcript dataset (if provided to us previously) 

![image](https://user-images.githubusercontent.com/26288352/195677522-d9054a2b-c658-4fdc-92ae-8f27feaecee4.png)
![image](https://user-images.githubusercontent.com/26288352/195677610-b699e203-df07-4412-ae4c-b73d3daebbb7.png)
![image](https://user-images.githubusercontent.com/26288352/195677677-e25b6a94-f418-492b-a89c-8ab7b9ab92f7.png)

