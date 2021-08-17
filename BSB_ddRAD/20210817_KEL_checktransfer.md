
Previously Thais downloaded the files to Discovery with `wget`, but she had trouble accessing them

Katie went in to check the files were properly downloaded
```
cd /work/lotterhos/ADD PATH
md5sum Cs* > md5sum_new.txt # takes a few minutes
awk '{gsub("/panfs/roc/umgc/illumina_analysis/210713_A00223_0599_BH7YTHDRXY-analysis/demultiplex_20210720-12-02-46/demultiplex/Lotterhos_Project_001/","",$2)}1' md5.txt > md5_edit.txt # remove file path from original file
comm -3 md5_edit.txt md5sum_new.txt # this didn't work because the files weren't in the same order
```
