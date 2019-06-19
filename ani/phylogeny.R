library(ape)
library(ggplot2)
library(Cairo)
library(RColorBrewer)
library(gplots)
library(phytools)

# read in identity matrix from pyani
n <- as.matrix(read.table("~/Desktop/repos/bifido/ani/ANIm_output/ANIm_percentage_identity.tab", 
                          head=T, row.names = 1))
ones <- matrix(1, 252, 252) # matrix as ones for element subtraction

dm <- ones - n # subtract to create distance matrix
tree <- nj(as.dist(dm)) # create a tree from the distance matrix

# plot tree
png('tree.png', width = 4000, height = 4400)
plot(tree, main="Neighbor Joining")
dev.off()

# set breve genomes as the root
root <- root(tree, c("breve_GCF_000220135.1_ASM22013v1", "breve_GCF_000213865.1_ASM21386v1", 
                     "breve_GCF_000158015.1_ASM15801v1"))

# color branch labels by species/subspecies
tipcol <- rep('black', length(root$tip.label))
categs <- c("longum", "blongum", "infantis", "suis", "breve")
colors <- c("dodgerblue2", "black", "firebrick2", "purple", "magenta")
for(i in 1:length(categs)) {
  tipcol[grep(categs[i], root$tip.label)] <- colors[i]
}

# plot rooted tree
png('root.png', width = 4000, height = 4400)
plot(root, main="Neighbor Joining", tip.color=tipcol)
dev.off()

svg('root.svg', width=44, height=54)
plot(root, main="Neighbor Joining", tip.color=tipcol)
dev.off()

# write rooted tree to a tree file
write.tree(root, file="tree.tre")

# read in heatmap csv
df <- as.matrix(read.table("~/Desktop/repos/bifido/scripts/blast/output/heatmap.csv", 
                           header=T, row.names = 1, sep=','))

# create heatmap
png('heatmap.png', width=800)
heatmap.2(df)
dev.off()

# pair phylogenetic tree with heatmap
phylo.heatmap(root, df)
