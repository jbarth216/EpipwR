# EpipwR

EpipwR is an R package for conducting power analysis for EWAS with continuous or binary outcomes. It is based on empirical EWAS to generate realistic data sets while leveraging classic statistical theory to cut down on computation time. Our corresponding paper can be found [here](https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf150/8173951). Users accustomed to pwrEWAS should first apply the `pwrE_to_EpipwR` function to put pwrEWAS parameterization into EpipwR format.

# Installation

EpipwR currently requires R version 4.4.0 or later.

```         
#Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EpipwR")

#github
#install.packages("devtools")
library(devtools)
install_github("jbarth216/EpipwR.data")
install_github("jbarth216/EpipwR", build_vignettes=TRUE)
```

# Usage

A detailed EpipwR workflow can be found in the EpipwR vignette:

```         
browseVignettes("EpipwR")
```
