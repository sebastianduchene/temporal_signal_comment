library(RColorBrewer)
library(NELSI)

setwd("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/temporal_sampling_bias/")

dir(pattern = "log")

balanced <- dir(pattern = "^balanced.+log")
balanced
unbalanced <- dir(pattern = "^unbalanced.+log")
unbalanced
# Compare estimates of the clock rate because it is the same across simulations. 

get_ci_tree_height <- function(log_file){
  posterior <- tail(read.table(log_file, header = TRUE), 2000)
  tree_height <- posterior$tree.height
  ci <- quantile(tree_height, c(0.025, 0.975))
  return(c(mean(tree_height), ci))
}

balanced_trees <- read.tree("balanced_sampling_trees.trees")

balanced_stats <- matrix(NA, nrow = length(balanced), ncol = 4)
for(i in 1:length(balanced)){
  corresponding_tree <- balanced_trees[[i]]
  balanced_stats[i,] <- c(balanced[i], get_ci_tree_height(balanced[i]) - max(allnode.times(corresponding_tree)))
}


head(balanced_stats)


#####

#balanced_stats <- balanced_stats[order(balanced_stats[,2]),]
colnames(balanced_stats) <- c("log_file", "mean", "lower", "upper")
balanced_stats[, 2:4] <- log10(as.numeric(balanced_stats[, 2:4]))
balanced_stats <- as.data.frame(balanced_stats)
balanced_stats$mean <- as.numeric(as.character(balanced_stats$mean))
balanced_stats$lower <- as.numeric(as.character(balanced_stats$lower))
balanced_stats$upper <- as.numeric(as.character(balanced_stats$upper))
head(balanced_stats)

balanced_stats <- balanced_stats[order(balanced_stats$mean),]

rate_quantiles <- quantile(rlnorm(1000, meanlog = log(1e-5), sdlog = 0.1), c(0.025, 0.975))
rate_quantiles

# True value is 1e-5 subs/site/year
plot(0, 0, ylim = log10(c(7e-6, 3e-5)), xlim = c(-10, 110), 
     xlab = "", xaxt = "n", yaxt = "n", 
     ylab = "Clock rate 10^-5(subs/site/year)")
axis(2, at = log10(c(7.5e-6, 1e-5, 1.25e-5, 1.5e-5, 2e-5)), 
     labels = c("0.75", "1.0", "1.25", "1.5", "2"), las = 2)
for(i in 1:nrow(balanced_stats)){
  lines(rep(i, 2), balanced_stats[i, 3:4], col = "darkblue")
  points(i, balanced_stats[i, 2], pch = 20, col = "darkblue")
}
lines(c(-15, 115), log10(c(1e-5, 1e-5)), lty = 2)
lines(c(-15, 115), rep(log10(rate_quantiles[1]), 2), lty = 3)
lines(c(-15, 115), rep(log10(rate_quantiles[2]), 2), lty = 3)

##############################
##############################
##############################

unbalanced_stats <- matrix(NA, nrow = length(unbalanced), ncol = 4)
for(i in 1:length(unbalanced)){
  unbalanced_stats[i,] <- c(balanced[i], get_ci_clock_rate(unbalanced[i]))
}

#balanced_stats <- balanced_stats[order(balanced_stats[,2]),]
colnames(unbalanced_stats) <- c("log_file", "mean", "lower", "upper")
unbalanced_stats[, 2:4] <- log10(as.numeric(unbalanced_stats[, 2:4]))
unbalanced_stats <- as.data.frame(unbalanced_stats)
unbalanced_stats$mean <- as.numeric(as.character(unbalanced_stats$mean))
unbalanced_stats$lower <- as.numeric(as.character(unbalanced_stats$lower))
unbalanced_stats$upper <- as.numeric(as.character(unbalanced_stats$upper))
head(unbalanced_stats)

unbalanced_stats <- unbalanced_stats[order(unbalanced_stats$mean),]

rate_quantiles <- quantile(rlnorm(1000, meanlog = log(1e-5), sdlog = 0.1), c(0.025, 0.975))
rate_quantiles

# True value is 1e-5 subs/site/year
#plot(0, 0, ylim = log10(c(7.25e-6, 2e-5)), xlim = c(-10, 110), 
#     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#axis(2, at = log10(c(7.5e-6, 1e-5, 1.25e-5, 1.5e-5, 2e-5)), 
#     labels = c("0.75", "1.0", "1.25", "1.5", "2"), las = 2)
for(i in 1:nrow(unbalanced_stats)){
  lines(rep(i-0.5, 2), unbalanced_stats[i, 3:4], col = "darkorange")
  points(i-0.5, unbalanced_stats[i, 2], pch = 20, col = "darkorange")
}

balanced_stats$hpd_width <- abs((balanced_stats$upper - balanced_stats$lower) / balanced_stats$mean)
unbalanced_stats$hpd_width <- abs((unbalanced_stats$upper - unbalanced_stats$lower) / unbalanced_stats$mean)

hist(balanced_stats$hpd_width, breaks = 20, col = rgb(0, 0, 1, 0.1), 
     xlab = "HPD width", main = "", border = "darkblue")
hist(unbalanced_stats$hpd_width, breaks = 100, col = rgb(0.7, 0.5, 0, 0.1), border = "darkorange", add = TRUE)

