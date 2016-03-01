library(phangorn)
library(multicore)
setwd("~/Desktop/Cichlid UCE's [SPRING]/fatlip/snphylo-trees")

phylip <- read.phyDat("fatlip.run2.phylip.txt", format="phylip", type="DNA")
newick <- read.tree("master-fatlip.snphylo.bs.tre")
fit <- pml(newick, phylip)

set.seed(18451)
bs <- bootstrap.pml(fit, bs = 1000, optNni=TRUE, multicore=T)
sample <- bs[floor(runif(100, min=1, max=1001))]

write.tree(sample,file="100-random-bs-trees.tre")