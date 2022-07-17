See [R Discovery tutorial](https://rc-docs.northeastern.edu/en/latest/software/r.html)

# Install R packages and dependencies using [Packrat](https://rstudio.github.io/packrat/)


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
install all the other packages needed within your packrat environment.
```bash
BiocManager::install("vcfR")
```

In R, to call on all the packages you need to initialize the environment where the R packages have been written to:
  `packrat::init("/home/r.gatins/R")`
Now just call on each package as you would normally `library("<package>")`

or 

```R
libraries_needed <- c("vcfR", "ape","adegenet", "RColorBrewer", "poppr", "radiator", "dplyr", "igraph",
                      "ggplot2", "reshape2")

for (i in 1:length(libraries_needed)){
  library( libraries_needed[i], character.only = TRUE)
}
```
