#!/usr/local/bin/Rscript

args <- commandArgs(TRUE)


library("mbGraphics")

meta.frame <- read.table("mb_graphics_data.txt", skip=1, header=TRUE, stringsAsFactors=FALSE)
data.types <- strsplit(readLines("mb_graphics_data.txt")[1], '\t')[[1]]
adiv(meta.frame, data.types)

meta.frame <- read.table("mb_graphics_data.txt", skip=1, header=TRUE, stringsAsFactors=FALSE)
data.types <- strsplit(readLines("mb_graphics_data.txt")[1], '\t')[[1]]
beta.frame <- read.table("beta_data.out", header=TRUE, stringsAsFactors=FALSE)
bdiv(meta.frame, data.types, beta.frame)

meta.frame <- read.table("mb_graphics_data.txt", skip=1, header=TRUE, stringsAsFactors=FALSE)
data.types <- strsplit(readLines("mb_graphics_data.txt")[1], '\t')[[1]]
numseq(meta.frame, data.types)

meta.frame <- read.table("mb_graphics_data.txt", skip=1, header=TRUE, stringsAsFactors=FALSE)
data.types <- strsplit(readLines("mb_graphics_data.txt")[1], '\t')[[1]]
tax.key <- read.table(args[1], header = TRUE, stringsAsFactors = FALSE)
tax.data <- read.table(args[2], header = TRUE, stringsAsFactors = FALSE)
phylo(meta.frame, data.types, tax.key, tax.data, args[3])
