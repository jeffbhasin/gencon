# gencon: analyze genomic context for a set of genomic ranges

Warning: Under active development, all features may not yet perform as documented.

## Basics
* Given a set of genomic ranges in a BED+header file, gencon outputs an annotated BED file with a column for each annotation category indicating if there is an overlap, a long format table where each overlap event is its own row (with full metadata from the associated annotation), and a summary table producing total counts of all input ranges that overlap with each of the annotation categories
* Implemented as an R script that takes command line options
* Uses the genomesyncr package (http://github.com/bluecranium/genomesyncr) to fetch and update annotation files from the UCSC Genome Browser tables
* Works on Linux, Unix, and Mac OS X (no Windows version due to rsync dependency of genomesyncr)

## Installation
First install the genomesyncr dependency. From an R console:

	library(devtools)
	options(repos=c("http://cran.rstudio.com","http://www.bioconductor.org/packages/release/bioc"))
	install_github("genomesyncr", username="bluecranium")

Clone this repo and run the script from a BASH terminal to output usage information:

	git clone git://github.com/bluecranium/gencon.git
	cd gencon
	./gencon --help

## Usage
	# Example:
	./gencon -i regions.bed --genome hg19 -u ~/Documents/myucsc/

* See also the "--help" option
* BED files must contain a header row (so that the output has a name for any custom fields you add) and at least 3 columns (chromosome, start, end), be tab-delimited, have chromosome numbers prefixed with "chr" (e.g. "chr1", "chr2", "chr3"), and use zero-based, half-open coordinates (see http://genome.ucsc.edu/FAQ/FAQformat#format1 and https://code.google.com/p/bedtools/wiki/FAQ#What_does_zero-based,_half-open_mean? for more information on the BED format)
* All columns after column 3 will be retained in the "wide" output, so they can contain any arbitrary metadata you'd like to keep attached with your ranges in the output
* By default, strandedness is not taken into account. Use the "-s" option to specify a column in the input file that has strand information (coded using "+" or "-") and then the overlaps will be done strand-aware.

## TODO
* Design a configuration file to specify the contexts and annotations to use, rather than this being hard coded into the R script (right now this is done by putting all overlap code into a separate R file so this function file can be swappable to change the details of what overlaps and contexts are output)