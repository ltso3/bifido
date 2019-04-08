library(reshape)

# read in data
species1 <- read.csv("tax_profiles/species1.csv")
species3 <- read.csv("tax_profiles/species3.csv")
species4 <- read.csv("tax_profiles/species4.csv")

# isolate b. longum
longum1 <- species1[species1$X.SampleID == "Bifidobacterium_longum",]
longum3 <- species3[species3$X.SampleID == "Bifidobacterium_longum",]
longum4 <- species4[species4$X.SampleID == "Bifidobacterium_longum",]

# concatenate all b. longum
total <- cbind(longum1, longum3)
total <- cbind(total, longum4)

# transpose data frame for reading
total <- melt(total, id="X.SampleID")
total <- subset(total, select = c(variable, value))
names(total) <- c("sample", "b. longum")

# filter data frame by counts
five <- total[total$`b. longum` > 5, ]
ten <- total[total$`b. longum` > 10, ]
twentyfive <- total[total$`b. longum` > 25, ]

# write files to csvs
write.csv(total, file = "total_longum.csv")
write.csv(five, file="five_longum.csv")
write.csv(ten, file="ten_longum.csv")
write.csv(twentyfive, file="twentyfive_longum.csv")
