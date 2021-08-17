
Previously Thais downloaded the files to Discovery with `wget`, but she had trouble accessing them

Katie went in to check the files were properly downloaded
```
cd /work/lotterhos/ADD PATH
md5sum Cs* > md5sum_new.txt # takes a few minutes
comm -3 md5.txt md5sum_new.txt
```
