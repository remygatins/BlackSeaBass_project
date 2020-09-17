#! /bin/batch

for f in *.sra_1.fastq; do mv "$f" "${f%.sra_1.fastq}.F"; done

for f in *.sra_2.fastq; do mv "$f" "${f%.sra_2.fastq}.R"; done

prefix=Pop1_

for name in SRR*; do

gzip < "$name" > "$prefix$name.fq.gz" 

done
