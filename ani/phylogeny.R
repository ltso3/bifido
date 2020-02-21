library(ape)
library(ggplot2)
library(Cairo)
library(RColorBrewer)
library(gplots)
library(phytools)
library(phangorn)

# #read in identity matrix from pyani
# n <- as.matrix(read.table("~/Desktop/repos/bifido/ani/ANIm_output/ANIm_percentage_identity.tab",
#                           head=T, row.names = 1))
# ones <- matrix(1, 274, 274) # matrix as ones for element subtraction
# 
# dm <- ones - n # subtract to create distance matrix
# tree <- nj(as.dist(dm)) # create a tree from the distance matrix

# plot tree
# png('tree.png', width = 4000, height = 4400)
# plot(tree, main="Neighbor Joining")
# dev.off()

# reduce number of breve genomes to be monophyletic
# breve <- read.csv('~/Desktop/repos/bifido/scripts/blast/output/breve.csv')
# breves <- as.character(breve[['Unnamed..0']])
# keep <- c("breve_GCF_000220135.1_ASM22013v1", "breve_GCF_000213865.1_ASM21386v1",
#           "breve_GCF_000158015.1_ASM15801v1")
# breves <- setdiff(breves, keep)
# tree_red <- drop.tip(tree, breves)
# 
# png('tree_red.png', width = 4000, height = 4400)
# plot(tree_red, main="Neighbor Joining")
# dev.off()

# set three breve genomes as root
# root <- root(tree_red, keep)

# color branch labels by species/subspecies
# tipcol <- rep('black', length(root$tip.label))
# categs <- c("longum", "blongum", "infantis", "suis", "breve")
# colors <- c("dodgerblue2", "black", "firebrick2", "purple", "magenta")
# for(i in 1:length(categs)) {
#   tipcol[grep(categs[i], root$tip.label)] <- colors[i]
# }

# plot rooted tree
# pdf('root.pdf', width = 100, height = 150)
# plot(root, main="Neighbor Joining", tip.color=tipcol)
# dev.off()
# 
# svg('root.svg', width=44, height=54)
# plot(root, main="Neighbor Joining", tip.color=tipcol)
# dev.off()

# write rooted tree to a tree file
# write.tree(root, file="root.tre")

# read in heatmap csv
# df <- as.matrix(read.table("~/Desktop/repos/bifido/scripts/blast/output/heatmap.csv",
#                           header=T, row.names = 1, sep=','))

# create heatmap
# png('tree_heatmap.png', width=800)
# heatmap.2(df)
# dev.off()

# drop all but three of the breve genomes in df
# root_df <- df[!rownames(df) %in% breves,]

# create heatmap
# png('root_heatmap.png', width=800)
# heatmap.2(root_df)
# dev.off()

# pair phylogenetic tree with heatmap
# png('tree_combo.png', width=4000, height=4400)
# phylo.heatmap(tree, df, colors=colorRampPalette(c("dodgerblue3", "white"))(200)[200:1])
# dev.off()

# pair phylogenetic tree with heatmap
# png('root_combo.png', width=4000, height=4400)
# phylo.heatmap(root, root_df, colors=colorRampPalette(c("dodgerblue3", "white"))(200)[200:1],
#              fsize=c(1,4,4), grid=TRUE, ftype="off", leg.txt="Num. of BLAST hits")
# dev.off()
 
# # -------------------------------------------------------------------------------------------------
# 
# # read in blon protein genes from clustal omega
# tree <- read.tree("~/Desktop/repos/bifido/scripts/blast/protein_blons/blon_2361.tree")
# species <- tree$tip.label
# split_spec <- sapply(strsplit(species, "!"), function(x) x[2], simplify=TRUE)
# 
# # color branch labels by unique names (needs some specification for final figure)
# # add specific colors for infantis and longum and suis
# # if infantis in the string replace with infantis, same with all subspecies
# split1 <- replace(split_spec, grepl("infantis",split_spec), "infantis")
# split2 <- replace(split1, grepl("subsp.-longum",split1), "subsp.-longum")
# 
# tipcol <- rep('black', length(tree$tip.label))
# colors <- rainbow(length(sort(unique(split2))))
# for(i in 1:length(unique(split2))) {
#   if(sort(unique(split2))[i] != "") {
#     tipcol[grep(sort(unique(split2))[i], tree$tip.label)] <- colors[i]
#   }
# }
# 
# png('blon_2361.png', width = 2500, height = 2000)
# plot(tree, main="Neighbor Joining", tip.color=tipcol)
# add.scale.bar(cex = 4, font = 2, col = "red")
# layout(1)
# dev.off()
# 
# # -------------------------------------------------------------------------------------------------
# # tree from concatenated infantis cluster
# tree <- read.tree("~/Desktop/repos/bifido/scripts/blast/cluster_genomes.tree")
# png('cluster_genomes.png', width = 2500, height = 2000)
# plot(tree, main="Neighbor Joining")
# dev.off()

# tree <- read.tree("~/Downloads/blon_2331.tree")
# plot(tree, main="Neighbor Joining")
# -------------------------------------------------------------------------------------------------
# creating maximum likelihood trees from protein_blons
# files <- list.files(path="~/Desktop/repos/bifido/scripts/blast/protein_blons", full.names = TRUE)
# blons <- list.files(path="~/Desktop/repos/bifido/scripts/blast/protein_blons", full.names = FALSE)
# likes <- numeric(length(files))
# blon_names <- numeric(length(files))
# 
# for(i in 1:length(files)) {
#   sequences <- read.FASTA(files[i], type="AA")
#   aligned <- muscle(sequences)
#   dat <- as.phyDat(aligned, type="AA")
#   dm = dist.ml(dat)
#   tree = NJ(dm)
#   
#   species <- tree$tip.label
#   split_spec <- sapply(strsplit(species, "!"), function(x) x[2], simplify=TRUE)
#   
#   # color branch labels by unique names (needs some specification for final figure)
#   # add specific colors for infantis and longum and suis
#   # if infantis in the string replace with infantis, same with all subspecies
#   split1 <- replace(split_spec, grepl("infantis",split_spec), "infantis")
#   split2 <- replace(split1, grepl("subsp.-longum",split1), "subsp.-longum")
#   
#   tipcol <- rep('black', length(tree$tip.label))
#   colors <- topo.colors(length(sort(unique(split2))))
#   for(i in 1:length(unique(split2))) {
#     if(sort(unique(split2))[i] != "") {
#       tipcol[grep(sort(unique(split2))[i], tree$tip.label)] <- colors[i]
#     }
#   }
#   
#   # png("~/Downloads/test.png", height = 1000, width = 1000)
#   fitStart = pml(tree, dat, k=4, inv=.2)
#   fit = optim.pml(fitStart, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
#   # root_name = "blon_2361-WP_094369814.1-ABC-transporter-ATP-binding-protein-!Romboutsia-weinsteinii!"
#   # plot(ladderize(root(fit$tree, match(root_name,fit$tree$tip.label))),
#   # plot(ladderize(fit$tree),type="phylogram", show.node.label=T,underscore=T,cex=0.8,no.margin=T, tip.color=tipcol)
#   # add.scale.bar(cex = 2, font = 2, col = "red")
#   # dev.off()
#   likes[i] <- fit["logLik"]
#   
#   blon <- strsplit(blons[i], ".faa")[[1]]
#   blon_names[i] <- blon
#   png(sprintf("%s.png", blon), height=3000, width=3000)
#   # need to add boostrapping step from NCBI paper
#   BS <- bootstrap.pml(fit, optNni=TRUE)
#   plotBS(fit$tree, BS, type = "phylogram", tip.color=tipcol)
#   add.scale.bar(cex = 2, font = 2, col = "red")
#   dev.off()
#   
#   write.tree(BS, file=sprintf("%s.tree", blon))
# }
# 
# df <- data.frame(blon_names, likes, stringsAsFactors=FALSE)
# write.csv(file="tree_logs.csv", x=df)

# -------------------------------------------------------------------------------------------------

# sequences <- read.FASTA("~/Desktop/repos/bifido/figure2/core_gene_genomes.fna", type="DNA")
# aligned <- muscle(sequences)
dat <- read.phyDat("~/Desktop/repos/bifido/figure2/core_alignment.phylip", type="DNA")
dm = dist.ml(dat)
tree = NJ(dm)
tree <- root(tree, c("breve_1", "breve_2", "breve_3", "breve_4"))

# color branch labels by species/subspecies
tipcol <- rep('black', length(tree$tip.label))
categs <- c("lon", "blo", "inf", "suis", "breve")
colors <- c("dodgerblue2", "black", "firebrick2", "purple", "magenta")
for(i in 1:length(categs)) {
  tipcol[grep(categs[i], tree$tip.label)] <- colors[i]
}

png('~/Downloads/nj_core.png', width = 4000, height = 4400)
plot(tree, main="Neighbor Joining", tip.color=tipcol)
dev.off()

# read in heatmap csv
df <- as.matrix(read.table("~/Desktop/repos/bifido/figure2/heatmap_all_species.csv",
                           header=T, row.names = 1, sep=','))

# create heatmap
png('~/Downloads/tree_heatmap.png', width=800)
heatmap.2(df)
dev.off()

# pair phylogenetic tree with heatmap
png('~/Downloads/tree_combo.png', width=4000, height=4400)
phylo.heatmap(tree, df, colors=colorRampPalette(c("dodgerblue3", "white"))(200)[200:1],
              fsize=c(1,4,4), grid=TRUE, ftype="off", leg.txt="Num. of BLAST hits")
dev.off()


png("~/Downloads/ml_core.png", width = 4000, height = 4400)
fitStart = pml(tree, dat, k=4, inv=.2)
fit = optim.pml(fitStart, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
# root_name = "blon_2361-WP_094369814.1-ABC-transporter-ATP-binding-protein-!Romboutsia-weinsteinii!"
# plot(ladderize(root(fit$tree, match(root_name,fit$tree$tip.label))),
# plot(ladderize(fit$tree),type="phylogram", show.node.label=T,underscore=T,cex=0.8,no.margin=T, tip.color=tipcol)
# add.scale.bar(cex = 2, font = 2, col = "red")
# dev.off()
likes[i] <- fit["logLik"]

# blon <- strsplit(blons[i], ".faa")[[1]]
# blon_names[i] <- blon
# png(sprintf("%s.png", blon), height=3000, width=3000)
# need to add boostrapping step from NCBI paper
BS <- bootstrap.pml(fit, optNni=TRUE)
plotBS(fit$tree, BS, type = "phylogram", tip.color=tipcol)
add.scale.bar(cex = 2, font = 2, col = "red")
dev.off()

png('~/Downloads/ml_tree_combo.png', width=4000, height=4400)
phylo.heatmap(fit$tree, df, colors=colorRampPalette(c("dodgerblue3", "white"))(200)[200:1],
              fsize=c(1,4,4), grid=TRUE, ftype="off", leg.txt="Num. of BLAST hits")
dev.off()


