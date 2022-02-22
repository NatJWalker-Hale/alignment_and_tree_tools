setwd("") # dir for files

# require("ape")
require("ggplot2")
require("ggseqlogo")

# aa1 <- read.table("tmp_DODAa1_683.txt")
# 
# aa2 <- read.table("tmp_DODAa2_683.txt")
# 
# ggseqlogo(aa1, method="prob", col_scheme="clustalx")
# ggseqlogo(aa2, method="prob", col_scheme="clustalx")

# read in JSDs to scale

df <- read.table("", header=T) # tsv file of sizes to scale to here

# read in CSV formatted probs from fastas of sites and seqs of interest
# produced by prep_props_for_custom_height_logo.py
# format for custom heights ggseqlogo

aa1 <- read.csv("", header = FALSE) # probs file for alignment 1 here
# format for input to ggseqlogo
aa1 <- t(aa1) 
dimnames(aa1)[[1]] <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
                        "P", "S", "T", "W", "Y", "V")
aa1 <- sweep(aa1, 2, df$size, FUN="*")

# make custom col scheme
# my version of CLUSTALX colours

cs1 <- make_col_scheme(chars = dimnames(aa1)[[1]], cols = c("#2a7fff", "#ff0000", "#3db83d", "#a93aa3", "#f08080",
                                                            "#3db83d", "#a93aa3", "#f18237", "#2ad4ff", "#899bd6",
                                                            "#899bd6", "#ff0000", "#899bd6", "#899bd6", "#ffcc00",
                                                            "#3db83d", "#3db83d", "#2a7fff", "#2ad4ff", "#2a7fff"))

ggseqlogo(aa1, method="custom", col_scheme=cs1)

# same as above but for alignment 2 

aa2 <- read.csv("", header = FALSE) # probs file for alignment 2 here
aa2 <- t(aa2)
dimnames(aa2)[[1]] <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
                        "P", "S", "T", "W", "Y", "V")
aa2 <- sweep(aa2, 2, df$size, FUN="*")

ggseqlogo(aa2, method="custom", col_scheme=cs1)
