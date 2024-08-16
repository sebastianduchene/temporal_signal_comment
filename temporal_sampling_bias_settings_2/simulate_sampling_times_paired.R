# Script to sample biased and uniform sampling times. I am basing this from HBV virus data
# sampled over about 12K years
# dsDNA genome of about 3.2 Kb
# clock rate of about 1e-5 subs/site/year
library(NELSI)
library(phangorn)

setwd("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/temporal_sampling_bias_settings_2/")

complete_trees <- list()
balanced_trees <- list()
unbalanced_trees <- list()

for(i in 1:100){
  # The sampling times will be exponentially distributed with a mean of 4000.
  # This means that on average the difference in sampling times is about 8,500, but can be about twice as high.
  sampling_times <- sort(round(rexp(5, 1/4000), 2))

  sampled_data <- rep(sampling_times, 100) # simulate 100 samples per time point
  
  sampled_data <- sampled_data - min(sampled_data)
  sampling_summary <- table(sampled_data)
  
  times_for_xml <- paste(names(sampling_summary), collapse = " ")
  samples_for_xml <- paste(as.numeric(sampling_summary), collapse = " ")
  
  xml_input <- readLines("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/temporal_sampling_bias/coalescent_heterochronous.xml")
  xml_input <- gsub("INPUT_TIMES", times_for_xml, xml_input)
  xml_input <- gsub("INPUT_NUM_SAMPLES", samples_for_xml, xml_input)
  writeLines(xml_input, con = "temp_out.xml")
  
  system("~/phyloApps/beast2.7.6/bin/beast -overwrite temp_out.xml")
  
  tree_temp <- read.nexus("coal_trees.tree")
  tree_temp$tip.label <- paste(tree_temp$tip.label, round(allnode.times(tree_temp, tipsonly = T, reverse = F), 3), sep = "_")
  tree_temp$tip.label
  #balanced sampling

  dates <- gsub(".+_", "", tree_temp$tip.label)
  dates_unique <- sort(as.numeric(unique(dates)), decreasing = T)
  balanced_samples <- c()
  for(s in dates_unique){
    balanced_samples <- c(balanced_samples, sample(which(dates == s), 20))
  }

  unbalanced_samples <- c()
  num_samples <- c(1, 1, 3, 5, 90)
  for(s in 1:length(dates_unique)){
    unbalanced_samples <- c(unbalanced_samples, sample(which(dates == dates_unique[s]), size = num_samples[s]))
  }
  
  
  balanced_trees[[i]] <- keep.tip(tree_temp, balanced_samples)
  unbalanced_trees[[i]] <- keep.tip(tree_temp, unbalanced_samples)
  complete_trees[[i]] <- tree_temp
  
  branch_rates <- rlnorm(nrow(tree_temp$edge), meanlog = log(1e-5) - 0.5^2 / 2, sdlog = 0.5)
  phylogram_temp <- tree_temp
  phylogram_temp$edge.length <- branch_rates * phylogram_temp$edge.length
  
  site_rates <- phangorn::discrete.gamma(0.1, 4)
  segment1 <- simSeq(phylogram_temp, l = 800, rate = site_rates[1])
  segment2 <- simSeq(phylogram_temp, l = 800, rate = site_rates[2])
  segment3 <- simSeq(phylogram_temp, l = 800, rate = site_rates[3])
  segment4 <- simSeq(phylogram_temp, l = 800, rate = site_rates[4])
  complete_sampling_aln <- as.DNAbin(c(segment1, segment2, segment3, segment4))
  balanced_sampling_aln <- complete_sampling_aln[balanced_trees[[i]]$tip.label, ]
  unbalanced_sampling_aln <- complete_sampling_aln[unbalanced_trees[[i]]$tip.label, ]
  write.dna(complete_sampling_aln, file = paste0("complete_sampling_aln_", i, ".fasta"), format = "fasta", nbcol = -1, colsep = "")
  write.dna(balanced_sampling_aln, file = paste0("balanced_sampling_aln_", i, ".fasta"), format = "fasta", nbcol = -1, colsep = "")
  write.dna(unbalanced_sampling_aln, file = paste0("unbalanced_sampling_aln_", i, ".fasta"), format = "fasta", nbcol = -1, colsep = "")
  #write.dna(simulated_genome, file = paste0("balanced_sampling_genome_", i, ".fasta"), format = "fasta", nbcol = -1, colsep = "")
  #plot(tree_temp, show.tip.label = F)
}

class(simulated_trees) <- "multiPhylo"

write.tree(complete_trees, file = "complete_sampling_trees.trees")
write.tree(balanced_trees, file = "balanced_sampling_trees.trees")
write.tree(unbalanced_trees, file = "unbalanced_sampling_trees.trees")
