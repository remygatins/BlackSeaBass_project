See [R Discovery tutorial](https://rc-docs.northeastern.edu/en/latest/software/r.html)

1. Through terminal or iTerm connect to Discovery
2. Load the R version you need `module load R/4.0.3`
3. Create a new directory for your R project by typing, `mkdir /home/r.gatins/R`
4. Open the R interface `R`
5. Install packrat `install.packages("packrat")`
    - during the installation, you will be prompted to install in a local directory, as you cannot install as root 
7. Initialize the environment where R packages will be written to:
    `packrat::init("/home/r.gatins/R")`\
7. You can then install R packages that you need. For example, to install a package called adegenet, type the following:
```bash
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("adegenet")
```

```bash
BiocManager::install("vcfR")
```
install all the other packages needed within your packrat environment.

In R to call on all the packages you need to initialize the environment where the R packages have been written to:
  `packrat::init("/home/r.gatins/R")`

