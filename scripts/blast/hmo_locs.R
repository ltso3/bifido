library(karyoploteR)

custom.genome <- toGRanges(data.frame(chr=c("infantis"), start=c(1), end=c(2791569)))
custom.genome <- toGRanges(data.frame(chr=c("infantis"), start=c(1), end=c(300000)))

kp <- plotKaryotype(genome = custom.genome)
kp <- plotKaryotype(plot.type=1, main="GCF_000825065.1_BIC1401212621a.V1", genome = custom.genome)
kpDataBackground(kp, data.panel = 1)
kpText(kp, chr="infantis", x=49691.5, y=1, labels="blon_2331", data.panel = 1)
kpText(kp, chr="infantis", x=51438, y=0.96, labels="blon_2332", data.panel = 1)
kpText(kp, chr="infantis", x=728.5, y=1, labels="blon_2333", data.panel = 1)
kpText(kp, chr="infantis", x=1888.4, y=0.96, labels="blon_2334", data.panel = 1)
kpText(kp, chr="infantis", x=4775, y=0.92, labels="blon_2335", data.panel = 1)
kpText(kp, chr="infantis", x=6699, y=0.88, labels="blon_2336", data.panel = 1)
kpText(kp, chr="infantis", x=7657.5, y=0.84, labels="blon_2337", data.panel = 1)
kpText(kp, chr="infantis", x=8361, y=0.80, labels="blon_2338", data.panel = 1)
kpText(kp, chr="infantis", x=9348.5, y=0.76, labels="blon_2339", data.panel = 1)
kpText(kp, chr="infantis", x=10402.5, y=0.72, labels="blon_2340", data.panel = 1)
kpText(kp, chr="infantis", x=11521.5, y=0.68, labels="blon_2341", data.panel = 1)
kpText(kp, chr="infantis", x=520, y=0.64, labels="blon_2342", data.panel = 1)
kpText(kp, chr="infantis", x=1496.5, y=0.60, labels="blon_2343", data.panel = 1)
kpText(kp, chr="infantis", x=2865, y=0.56, labels="blon_2345", data.panel = 1)
kpText(kp, chr="infantis", x=520, y=0.52, labels="blon_2346", data.panel = 1)
kpText(kp, chr="infantis", x=1496.5, y=0.48, labels="blon_2347", data.panel = 1)
kpText(kp, chr="infantis", x=781.5, y=0.44, labels="blon_2348", data.panel = 1)
kpText(kp, chr="infantis", x=2289, y=0.40, labels="blon_2349", data.panel = 1)
kpText(kp, chr="infantis", x=3525, y=0.36, labels="blon_2350", data.panel = 1)
kpText(kp, chr="infantis", x=4981, y=0.32, labels="blon_2351", data.panel = 1)
kpText(kp, chr="infantis", x=916, y=0.28, labels="blon_2352", data.panel = 1)
kpText(kp, chr="infantis", x=913, y=0.24, labels="blon_2353", data.panel = 1)
kpText(kp, chr="infantis", x=1927, y=0.20, labels="blon_2354", data.panel = 1)
kpText(kp, chr="infantis", x=2962.5, y=0.16, labels="blon_2355", data.panel = 1)
kpText(kp, chr="infantis", x=4893.5, y=0.12, labels="blon_2356", data.panel = 1)
kpText(kp, chr="infantis", x=6338, y=0.08, labels="blon_2357", data.panel = 1)
kpText(kp, chr="infantis", x=7338, y=0.04, labels="blon_2358", data.panel = 1)
kpText(kp, chr="infantis", x=8823.5, y=0.03, labels="blon_2359", data.panel = 1)
kpText(kp, chr="infantis", x=9876, y=0.02, labels="blon_2360", data.panel = 1)
kpText(kp, chr="infantis", x=10755.5, y=0.01, labels="blon_2361", data.panel = 1)

# read in every cytobands file and save plot
files <- list.files(path="~/Desktop/repos/bifido/scripts/blast/cytobands", full.names = TRUE)

for(i in 1:length(files)) {
  custom.cytobands <- toGRanges(files[i])
  par(bg=NA)
  kp <- plotKaryotype(plot.type=1, genome = custom.genome, cytobands = custom.cytobands)
}


