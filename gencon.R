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
    make_option(c("-g", "--genome"), action="store", help="UCSC Genome build of coordinates in input BED file (e.g. hg19, mm10, etc) [required]"),
    make_option(c("-o", "--output-prefix"), action="store", help="Base file name of output files (.wide, .long, and .sum extensions will be added) [optional, default: input file name]"),
    make_option(c("-s", "--strand-col"), action="store", help="Take into account strandedness when computing overlaps, and specify which column number of the input BED has strand information [optional]")
)

opt <- parse_args(OptionParser(option_list=option_list))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Run overlaps

# -----------------------------------------------------------------------------
# =============================================================================