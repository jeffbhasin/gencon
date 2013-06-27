# gencon: analyze genomic context for a set of genomic ranges

## Basics
* Given a set of genomic ranges in a BED file, gencon outputs an annotated BED file with a column for each annotation category indicating if there is an overlap, a long format table where each overlap event is its own row (with full metadata from the associated annotation), and a summary table producing total counts of all input ranges that overlap with each of the annotation categories
* Implemented as an R script that takes command line options
* Uses the genomesyncr package (http://github.com/bluecranium/genomesyncr) to fetch and update annotation files from the UCSC Genome Browser tables
* Works on Linux, Unix, and Mac OS X (no Windows version due to rsync dependency of genomesyncr)

## Installation
First install the genomesyncr dependency. From an R console:

	library(devtools)
	options(repos=c("http://cran.rstudio.com","http://www.bioconductor.org/packages/release/bioc"))
	install_github("genomesyncr", username="bluecranium")

Clone this repo and run the script from a BASH terminal to output usage information:

	./gencon --help

## Usage
* See the "--help" option

## TODO
* Design a configuration file to specify the contexts and annotations to use, rather than this being hard coded into the R script