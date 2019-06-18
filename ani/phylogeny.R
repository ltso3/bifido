library(ape)
library(phyclust, quiet = TRUE)

m <- as.matrix(read.table("~/Desktop/repos/bifido/ani/ANIm_output/ANIm_percentage_identity.tsv", head=T, row.names=1)) 
map <- read.table(file = '~/Desktop/repos/bifido/ani/mapping.tsv', sep = '\t', header = TRUE)
labels <- as.vector(map[["X0"]])
rownames(m) <- labels
colnames(m) <- labels

arbol <- nj(as.dist(m)) 
png('rplot.png', width = 32400, height = 3400)
plot(arbol, main = "Neighbor Joining")
dev.off()

phylip <- lower.tri(m, diag = FALSE)

#------------------------------------------------------------------------------------------

n <- as.matrix(read.table("~/Desktop/repos/bifido/ani/ani_ada/ANIm_output/ANIm_percentage_identity.tab", head=T, row.names = 1))
ones <- matrix(1, 252, 252)

dm <- ones - n
tree <- nj(as.dist(dm))
png('rplot.png', width = 4000, height = 4400)
plot(tree, main="Neighbor Joining")
dev.off()

root <- root(tree, c("breve_GCF_000220135.1_ASM22013v1", "breve_GCF_000213865.1_ASM21386v1", 
                     "breve_GCF_000158015.1_ASM15801v1"))
png('rootplot.png', width = 4000, height = 4400)
plot(root, main="Neighbor Joining")
dev.off()