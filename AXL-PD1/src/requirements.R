required_packages <- scan(text="
ape
ggplot2
ggrepel
dplyr
clipr
stats
pheatmap
limma
stringr
stringi
edgeR
fgsea
tidyr
plyr
DOSE
Seurat
Matrix
openxlsx
randomcoloR
VennDiagram
xfun
data.table
readr
ggplotify
ggpubr
tximport
SeuratWrappers
vegan
timechange
tidyverse
forcats
magrittr
matrixStats
gplots
cluster
NMF
tsne
RColorBrewer
remotes
SeuratObject
locfit
roxygen2
miceadds
",
what = "character")



bioconductor_packages <- scan(text="
limma
DESeq2
edgeR
batchelor
fgsea
DOSE
GSVA
AUCell
TCGAbiolinks
preprocessCore
cBioPortalData
AnVIL
sva",
what = "character")


install.packages(required_packages, repos='https://cran.uni-muenster.de/');

n_warnings <- length(warnings())
if (n_warnings) {
    print("Warnings:")
    print(n_warnings)
    print(warnings())
    stop(n_warnings)
} 


remotes::install_github('satijalab/seurat-wrappers')




if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(bioconductor_packages, force = T)

n_warnings <- n_warnings - length(warnings())
if (n_warnings) {
    print("Warnings:")
    print(n_warnings)
    print(warnings())
    stop(n_warnings)
}  
