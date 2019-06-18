library(ape)
library(ggplot2)
library(svglite)

n <- as.matrix(read.table("~/Desktop/repos/bifido/ani/ANIm_output/ANIm_percentage_identity.tab", 
                          head=T, row.names = 1))
ones <- matrix(1, 252, 252)

dm <- ones - n
tree <- nj(as.dist(dm))
png('tree.png', width = 4000, height = 4400)
plot(tree, main="Neighbor Joining")
dev.off()

root <- root(tree, c("breve_GCF_000220135.1_ASM22013v1", "breve_GCF_000213865.1_ASM21386v1", 
                     "breve_GCF_000158015.1_ASM15801v1"))

tipcol <- rep('black', length(root$tip.label))
categs <- c("longum", "blongum", "infantis", "suis", "breve")
colors <- c("dodgerblue2", "black", "firebrick2", "purple", "magenta")
for(i in 1:length(categs)) {
  tipcol[grep(categs[i], root$tip.label)] <- colors[i]
}

png('root.png', width = 4000, height = 4400)
plot(root, main="Neighbor Joining", tip.color=tipcol)
ggsave(filename='root.svg')
dev.off()

write.tree(root, file="tree.tre")
