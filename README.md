# rcell2
functions to analyze Cell-ID data in a tidy framework


# Installation

## dependencies

Install the [```tidyverse```][1] meta-package (and use it, you'll not regret it) and [```devtools```][2]. In addition, install [```EBImage```][3] package (required to look at cells) by copying and running the following script. 



```
install.packages(c("tidyverse", "devtools"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EBImage")
```

## installing the package

Install using devtools, from github repo:

```
devtools::install_github("gerbeldo/rcell2")

```


[1]:https://www.tidyverse.org/
[2]:https://github.com/r-lib/devtools
[3]:https://bioconductor.org/packages/release/bioc/html/EBImage.html
