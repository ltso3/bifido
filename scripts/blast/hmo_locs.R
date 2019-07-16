library(karyoploteR)

# read in every cytobands file and save plot
files <- list.files(path="~/Desktop/repos/bifido/scripts/blast/cytobands", full.names = TRUE)
fs <- list.files(path="~/Desktop/repos/bifido/scripts/blast/cytobands", full.names = FALSE)
lens <- c(2832748, 2828958, 2564809, 2780321, 2789037, 2879623, 2793888, 2791569, 2821883, 
          2787528, 2791524, 2786838, 2578115, 2749833, 2643621, 2640127, 2579732, 2794568,
          2812033, 2579388, 2788991, 2602591, 2832748)
df <- data.frame(files, lens)

for(i in 1:length(files)) {
  par(bg=NA)
  len <- df$lens[df$files == files[i]]
  png(sprintf("%s.png", strsplit(fs[i], ".txt")), width = 1000)
  custom.genome <- toGRanges(data.frame(chr=c("infantis"), start=c(1), end=c(len)))
  custom.cytobands <- toGRanges(files[i])
  genes <- read.table(file = files[i], sep = '\t', header = TRUE)
  markers <- data.frame(chr=genes$chr, pos=genes$start, labels=genes$name)
  kp <- plotKaryotype(plot.type=2, genome = custom.genome, cytobands = custom.cytobands,
                      chromosomes="infantis")
  kpAddBaseNumbers(kp, tick.dist = len/10, tick.len = 10, tick.col="red", cex=1,
                   minor.tick.dist = len/100, minor.tick.len = 5, minor.tick.col = "gray")
  kpPlotMarkers(kp, chr=markers$chr, x=markers$pos, labels=markers$labels, 
                marker.parts = c(1, 100, 1), label.dist = 0.01, max.iter=1000)
  dev.off()
}

#kpAddCytobandLabels(kp, force.all=TRUE, srt=90, col="orange", cex=0.5)
#kpAddBaseNumbers(kp)
