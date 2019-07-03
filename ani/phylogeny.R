library(ape)
library(ggplot2)
library(Cairo)
library(RColorBrewer)
library(gplots)
library(phytools)
library(phangorn)

# read in identity matrix from pyani
n <- as.matrix(read.table("~/Desktop/repos/bifido/ani/ANIm_output/ANIm_percentage_identity.tab", 
                          head=T, row.names = 1))
ones <- matrix(1, 274, 274) # matrix as ones for element subtraction

dm <- ones - n # subtract to create distance matrix
tree <- nj(as.dist(dm)) # create a tree from the distance matrix

# plot tree
png('tree.png', width = 4000, height = 4400)
plot(tree, main="Neighbor Joining")
dev.off()

# reduce number of breve genomes to be monophyletic
breve <- read.csv('~/Desktop/repos/bifido/scripts/blast/output/breve.csv')
breves <- as.character(breve[['Unnamed..0']])
keep <- c("breve_GCF_000220135.1_ASM22013v1", "breve_GCF_000213865.1_ASM21386v1", 
          "breve_GCF_000158015.1_ASM15801v1")
breves <- setdiff(breves, keep)
tree_red <- drop.tip(tree, breves)

png('tree_red.png', width = 4000, height = 4400)
plot(tree_red, main="Neighbor Joining")
dev.off()

# set three breve genomes as root
root <- root(tree_red, keep)

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
write.tree(root, file="root.tre")

# read in heatmap csv
df <- as.matrix(read.table("~/Desktop/repos/bifido/scripts/blast/output/heatmap.csv", 
                           header=T, row.names = 1, sep=','))

# create heatmap
png('tree_heatmap.png', width=800)
heatmap.2(df)
dev.off()

# drop all but three of the breve genomes in df
root_df <- df[!rownames(df) %in% breves,]

# create heatmap
png('root_heatmap.png', width=800)
heatmap.2(root_df)
dev.off()

# pair phylogenetic tree with heatmap
png('tree_combo.png', width=4000, height=4400)
phylo.heatmap(tree, df, colors=colorRampPalette(c("dodgerblue3", "white"))(200)[200:1])
dev.off()

# pair phylogenetic tree with heatmap
png('root_combo.png', width=4000, height=4400)
phylo.heatmap(root, root_df, colors=colorRampPalette(c("dodgerblue3", "white"))(200)[200:1])
dev.off()

# -------------------------------------------------------------------------------------------------

# read in blon protein genes from clustal omega
tree <- read.tree("~/Desktop/repos/bifido/scripts/blast/protein_blons/blon_2361.tree")
species <- tree$tip.label
split_spec <- sapply(strsplit(species, "!"), function(x) x[2], simplify=TRUE)

# color branch labels by unique names (needs some specification for final figure)
# add specific colors for infantis and longum and suis
# if infantis in the string replace with infantis, same with all subspecies
split1 <- replace(split_spec, grepl("infantis",split_spec), "infantis")
split2 <- replace(split1, grepl("subsp.-longum",split1), "subsp.-longum")

tipcol <- rep('black', length(tree$tip.label))
colors <- rainbow(length(unique(split2)))
for(i in 1:length(unique(split2))) {
  tipcol[grep(unique(split2)[i], tree$tip.label)] <- colors[i]
}

png('blon_2361.png', width = 2500, height = 2000)
plot(tree, main="Neighbor Joining", tip.color=tipcol)
add.scale.bar(cex = 4, font = 2, col = "red")
layout(1)
dev.off()

# -------------------------------------------------------------------------------------------------
# tree from concatenated infantis cluster
tree <- read.tree("~/Desktop/repos/bifido/scripts/blast/cluster_genomes.tree")
png('cluster_genomes.png', width = 2500, height = 2000)
plot(tree, main="Neighbor Joining")
dev.off()

