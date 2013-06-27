#!/usr/bin/Rscript

# #############################################################################
# gencon: analyze genomic context for a set of genomic ranges
# Author: Jeffrey Bhasin <jeffb@case.edu>
# Created: 2013-06-27
#
# See readme.md for more information
#
# #############################################################################

# =============================================================================
# Packages and Globals
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(genomesyncr))
# =============================================================================

# =============================================================================
# Local Functions

# =============================================================================


# =============================================================================
# Analysis Code

# -----------------------------------------------------------------------------
# Parse command line options
option_list <- list(
    make_option(c("-i", "--input-bed"), action="store", help="BED format file of input ranges [required]"),
    make_option(c("--genome"), action="store", help="UCSC Genome build of coordinates in input BED file (e.g. hg19, mm10, etc) [required]"),
    make_option(c("-u", "--ucsc-path"), action="store", help="Path to local mirror of UCSC files as created by genomesyncr [required]"),
    #make_option(c("-o", "--output-prefix"), action="store", help="Base file name of output files (.wide, .long, and .sum extensions will be added) [optional, default: input file name]"),
    #make_option(c("-s", "--strand-col"), action="store", help="Take into account strandedness when computing overlaps, and specify which column number of the input BED has strand information [optional]"),
    make_option(c("-r", "--overlap-funcs"), action="store", help="Use a different function file to control what overlaps are output [optional, default: overlap.func.R]")
)

opt <- parse_args(OptionParser(option_list=option_list))

#genome <- "hg19"
#ucsc.path <- "~/Documents/myucsc/"
#input.bed <- "ignore/test.bed"
#output.prefix <- input.bed
#strand.col <- ""
#overlap.funcs <- "overlap.func.R"

genome <- opt$genome
ucsc.path <- opt$u
input.bed <- opt$i
output.prefix <- input.bed
#strand.col <- ""
overlap.funcs <- "overlap.func.R"

debug <- function()
{
    genome <- "hg19"
    ucsc.path <- "~/Documents/myucsc/"
    input.bed <- "ignore/test.bed"
    output.prefix <- input.bed
    strand.col <- ""
    overlap.funcs <- "overlap.func.R"
}

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Load input BED and turn into GRanges object

input <- read.table(file=input.bed, sep="\t", header=TRUE)
input.gr <- GRanges(seqnames=input[,1],ranges=IRanges(start=input[,2]+1,end=input[,3]))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Run overlaps

# load in the overlap function code (this can be swapped to try different annotation schemes)

source(overlap.funcs)

# call the main function to produce overlap output list
over.out <- doOverlaps(input.gr, genome, ucsc.path)

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Save output files

#over.out$wide
write.csv(over.out$wide, file=paste(output.prefix,".wide.csv",sep=""), row.names=FALSE)

#over.out$long
#write.csv(over.out$long, file=paste(output.prefix,".long.csv",sep=""), row.names=FALSE)

#over.out$summary
write.csv(over.out$summary, file=paste(output.prefix,".summary.csv",sep=""), row.names=FALSE)

# -----------------------------------------------------------------------------
# =============================================================================