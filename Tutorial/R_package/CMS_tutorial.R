library(Biobase)
library(CMScaller)
par(mfrow=c(1,2))

### CMS prediction of TCGA primary colorectal cancers
res <- CMScaller(crcTCGAsubset, RNAseq=TRUE, doPlot=TRUE)
head(res)

### Camera Gene Set Analysis with CMS informative gene sets
cam <- CMSgsa(emat=crcTCGAsubset, class=res$prediction, RNAseq=TRUE)
head(cam$CMS4)

### limma differential gene expression analysis and visualization
deg <- subDEG(emat=crcTCGAsubset, class=res$prediction, doVoom=TRUE)
subVolcano(deg, geneID="symbol")