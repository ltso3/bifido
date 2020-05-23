library(ggplot2)
library(tidyr)
library(plyr)

df <- read.csv('~/Downloads/bigger.csv')

df <- df[df$bf != "Unreported",]
png("~/Downloads/boxplot_small.png", width=600, height = 400)
p <- ggplot(df, aes(x=bf, y=abundance)) + 
  geom_boxplot() + xlab("Breastfeeding Status") + ylab("Bifidobacterium Abundance") +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15))
plot(p)
tapply(df$abundance, df$bf, mean)
ggsave('~/Downloads/boxplot_small.png', antialias="none", dpi="print")

pairwise.wilcox.test(df$abundance, df$bf,
                     p.adjust.method = "BH")


# create PCR barplot
df <- read.csv("~/Downloads/table5.csv")
png("~/Downloads/barplot.png", width = 800)
df$taxa <- factor(df$taxa, levels = df$taxa)
p <- ggplot(df, aes(x=taxa, y=Frequency)) + 
  geom_bar(stat="identity") + theme(axis.text.x=element_text(size=25)) +
  theme(axis.text.y=element_text(size=50)) +
  geom_text(aes(label=Frequency), vjust=-0.25, size=30)
p
ggsave('~/Downloads/barplot.png', antialias="none", dpi="print", width=44.4, height=26.7)
dev.off()
