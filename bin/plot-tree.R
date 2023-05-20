#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load required packages
library(ape)
library(adephylo)

# Read the tree file
tree <- read.tree(args[1])

# define year-month
yearmonth = as.character(args[2])

# define reference sequence id
ref_id <- as.character(args[3])

# root the tree, if it isn't already
if (is.rooted(tree)==F){
  tree <- root(tree, ref_id, resolve.root = T)
}

# fixing non-unique node labels
node_support <- tree$node.label
tree$node.label <- 1:length(node_support)

# Create a plot of the tree
pdf("tree_plot.pdf", width = 8.5, height = 10)
plot(tree, show.node.label = TRUE)
dev.off()