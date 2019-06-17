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
