#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Histogram of distances from fragment start/end to nearest HindIII cut site")
parser$add_argument("-f", help="distance file (ex. hindIII_distances)", type="character", dest="f", required=TRUE)
parser$add_argument("-w", help="wildcard for sample (ex. HBCRACKHiC-K562-HN-TD-R1_ACAGTG_L008)", type="character", dest="w", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()

distances <- scan(paste(wkdir, args$f, sep="/"))

pdf(paste(wkdir, "/plots/", args$w, "_distance_hist.pdf", sep=""), width=8, height=5)
par(mar=c(5, 5, 4, 2) + 0.1)
hist(distances, breaks=seq(0,max(distances)), xlim=c(0, 100),  main=args$w,
	ylab="# of mapped fragments", xlab="bp distance to nearest HindIII site (AAGCTT)", cex.lab=1.5, cex.main=1.5)
box(bty="l")
dev.off()