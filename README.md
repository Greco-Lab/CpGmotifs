# CpGmotifs

A user-friendly graphical interface that allows searching and visualising motifs linked to CpG sites from differential methylation analyses.

## Running the CpGmotifs Docker image (suggested)

If you don't have docker installed on your system you can install it by following the instructions at  https://www.docker.com/get-docker.

The FunMappOne docker image is available at https://hub.docker.com/r/grecolab/cpgmotifs


## Using CpGmotifs source from GitHub

### Linux system library dependencies

```BASH
     pandoc
     pandoc-citeproc
     libexpat1-dev
     libcairo2-dev
     libxt-dev
     libssl-dev
     libssh2-1-dev
     libssl1.0.0
     libcurl4-openssl-dev
     libxml2-dev
     ghostscript
```

### Install R dependencies

```R
#Universal Bioconductor package installation function
  install.bioc <- function(pkg){
    vers <- getRversion()
    if (vers >= "3.6"){
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    }else{
      if (!requireNamespace("BiocInstaller", quietly = TRUE)){
        source("https://bioconductor.org/biocLite.R")
        biocLite(pkg, suppressUpdates=TRUE)
      }else{
        BiocInstaller::biocLite(pkg, suppressUpdates=TRUE)
      }
    }
  }

#Install Bioconductor dependencies
bioc_pkgs <- c('minfi','BSgenome.Hsapiens.UCSC.hg19','IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b2.hg19', 'Biostrings')
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
if(length(bioc_pkgs.inst)>0){
  print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
  for(pkg in bioc_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."))
    install.bioc(pkg)
    print("Installed!!!")
  }
}

#Install CRAN dependencies
cran_pkgs <- c('curl','RCurl','shiny', 'shinyjs', 'shinydashboard', 'readr', 'DT', 'tibble', 'gplots',
                              'dendextend', 'foreach', 'doParallel', 'XML', 'BiocManager')

cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(cran_pkgs.inst)>0){
  print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
  for(pkg in cran_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."))
    install.packages(pkg, repo="http://cran.rstudio.org", dependencies=TRUE)
    print("Installed!!!")
  }
}
```

### Run CpGmotifs From GitHub
```R
# Load 'shiny' library
library(shiny)
library(shinyjs)
# run on the host port 8787 (or whaterver port you want to map on your system)
runGitHub("CpGmotifs", "Greco-Lab")
```

### or from a local copy
```R
  # Clone the git repository
  git clone https://github.com/Greco-Lab/FunMappOne FunMappOne
  # Start R session and run by using runApp()
  library(shiny)
  library(shinyjs)
  # run on the host port 8787 (or whaterver port you want to map on your system)
  runApp("./CpGmotifs/")
```
