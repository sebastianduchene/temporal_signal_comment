library(NELSI)
library(phangorn)
setwd("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/measurable_evolution/")

# For a TB-like example, use a 4.3Mb and a subst rate of 1e-7 subs/site/year, meaning one substitution every 2.3 years?

time_to_subst <- 1 / (1e-7 * 4.3e6)


trees <- read.tree(text = write.tree(read.nexus("measurable_evolution.trees"))) # Just some wrangling because the trees from the read.nexus have a single attribute with all the tip names


for(i in 1:length(trees)){
    dates <- round(allnode.times(trees[[i]], tipsonly = T, reverse = F), 2)
    temp_tree <- trees[[i]]
    temp_tree$tip.label <- paste(temp_tree$tip.label, dates, sep = "_")
    trees[[i]] <- temp_tree
}

write.tree(trees, file = "measurable_evolution_dated.trees")

# Simulate evolutionary rates and sequence alignment for each tree:
# just simulate 10 alignments for now:

for(i in 1:10){ #length(trees)){
    branch_rates <- rlnorm(length(trees[[i]]$edge.length), meanlog = log(1e-7), sdlog = 0.5)
    phylogram <- trees[[i]]
    phylogram$edge.length <- phylogram$edge.length * branch_rates
    site_rates <- phangorn::discrete.gamma(0.1, 4)
    print(paste("simulating alignment", i, sep = " "))
    segment1 <- simSeq(phylogram, l = 4300000/4, rate = site_rates[1])
    segment2 <- simSeq(phylogram, l = 4300000/4, rate = site_rates[2])
    segment3 <- simSeq(phylogram, l = 4300000/4, rate = site_rates[3])
    segment4 <- simSeq(phylogram, l = 4300000/4, rate = site_rates[4])
    simulated_genome <- as.DNAbin(c(segment1, segment2, segment3, segment4))
    write.dna(simulated_genome, file = paste("measurable_evolution_sim_", i, ".fasta", sep = ""), 
        format = "fasta", nbcol = -1, colsep = "")
}
