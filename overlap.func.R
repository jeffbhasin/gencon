# #############################################################################
# overlap.func.R
# Author: Jeffrey Bhasin <jeffb@case.edu>
# Created: 2013-06-27
#
# This is a function file for gencon which performs a specific annotation.
# These function files can be swapped out to interchange annotation types.
# There must be a "doOverlaps()" function that returns a list.
# The list must contain 3 dataframes: wide, long, summary.
# #############################################################################

# =============================================================================
# Packages and Globals
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
nCores <- 2
registerDoMC(nCores)
# =============================================================================

# =============================================================================
# Main Function
doOverlaps <- function(input.gr, genome, ucsc.path)
{
	# Sync/download annotation
	syncTable(genome=genome, table="knownGene", local=ucsc.path)
	syncTable(genome=genome, table="knownCanonical", local=ucsc.path)
	syncTable(genome=genome, table="kgXref", local=ucsc.path)
	#syncTable(genome=genome, table="cytoBand", local=ucsc.path)

	# Create the "ann" big table that joins in gene names
	ann <- readUCSCAnnotation(genome, ucsc.path)
	#ann.ex <- parseExons(ann)

	# Make the wide format
	genic <- getGenicOverlap(input.gr,ann)
	genic[genic>0] <- "yes"
	genic[genic==0] <- "no"

	genic.genes <- getGenicOverlapGenes(input.gr,ann)

	upstream <- getUpstreamOverlap(input.gr,ann,1000,500)
	upstream[upstream>0] <- "yes"
	upstream[upstream==0] <- "no"

	downstream <- getDownstreamOverlap(input.gr,ann,500,1000)
	downstream[downstream>0] <- "yes"
	downstream[downstream==0] <- "no"

	utr <- get3primeUTROverlap(input.gr,ann)
	utr[utr>0] <- "yes"
	utr[utr==0] <- "no"

	link <- getBrowserURLs(input.gr, genome)

	wide <- cbind(chr=as.character(seqnames(input.gr)), start=start(input.gr), end=end(input.gr),genic,upstream,downstream,utr,genic.genes, link)

	# Make the long format
	long <- data.frame()

	# Make the summary
	types <- c("genic", "utr", "upstream", "downstream")

	count_genic <- length(genic[genic=="yes"])
	count_utr <- length(utr[utr=="yes"])
	count_upstream <- length(upstream[upstream=="yes"])
	count_downstream <- length(downstream[downstream=="yes"])
	counts <- c(count_genic, count_utr, count_upstream, count_downstream)

	percents <- round((counts/length(input.gr))*100, digits=2)

	summary <- data.frame(type=types, count=counts, percent=percents)

	# Return the list object
	final.list <- list(wide=wide, long=long, summary=summary)
	final.list
}
# =============================================================================

# =============================================================================
# Local Functions

getBrowserURLs <- function(input.gr, genome)
{
	# Give links to UCSC at the position of this region
	#http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr22%3A1-500
	paste("http://genome.ucsc.edu/cgi-bin/hgTracks?db=", genome, "&position=", as.vector(seqnames(input.gr)), "%3A", start(input.gr), "-", end(input.gr),sep="")
}

readUCSCAnnotation <- function(genome, ucsc.path)
{
	#load downloaded tables
	knownGene <- loadTable(genome=genome, table="knownGene", local=ucsc.path)
	knownCanonical <- loadTable(genome=genome, table="knownCanonical", local=ucsc.path)
	kgXref <- loadTable(genome=genome, table="kgXref", local=ucsc.path)
	#cytoBand <- loadTable(genome=genome, table="cytoBand", local=ucsc.path)

	#join in gene symbols
	kg.sub <- knownGene
	names(kg.sub)[1] <- "name" 
	kgXref.sub <- kgXref
	names(kgXref.sub)[1] <- "name"
	kg.ann.1 <- join(kg.sub,kgXref.sub,by="name",type="left")

	#join in canonical annotation (largest coding seq from nearest cluster)
	kgCanon.sub <- data.frame(name=knownCanonical$transcript,canonical="1",stringsAsFactors=FALSE)
	kg.ann.2 <- join(kg.ann.1,kgCanon.sub,by="name",type="left")
	kg.ann <- kg.ann.2
	kg.ann[is.na(kg.ann$canonical),]$canonical <- 0

	#join in cytological bands based on gene position
	#kg.ranges <- GRanges(seqnames=kg.ann$chrom,ranges=IRanges(start=kg.ann$txStart,end=kg.ann$txEnd),strand=kg.ann$strand,id=kg.ann$name)
	#cb.ranges <- with(cytoBand, GRanges(seqnames=chrom,ranges=IRanges(start=chromStart,end=chromEnd),band=name,gstain=gieStain))
	
	#band genes start in
	#kg.ranges$CytoBandStartIndex <- findOverlaps(kg.ranges,cb.ranges,select="first")
	#band genes end in (in case they span more than one)
	#kg.ranges$CytoBandEndIndex <- findOverlaps(kg.ranges,cb.ranges,select="last")
	#total bands spanned
	#kg.ranges$CytoBandsSpan <- countOverlaps(kg.ranges,cb.ranges)
	#create DF from GR metadata
	#cyto.map <- data.frame(name=kg.ranges$id,count=kg.ranges$CytoBandsSpan,start=kg.ranges$CytoBandStartIndex,end=kg.ranges$CytoBandEndIndex,stringsAsFactors=FALSE)

	#cyto.map$startname <- "NA"
	#cyto.map[is.na(cyto.map$start)==FALSE,]$startname <- cb.ranges[cyto.map[is.na(cyto.map$start)==FALSE,]$start]$band
	#cyto.map$endname <- "NA"
	#cyto.map[is.na(cyto.map$start)==FALSE,]$endname <- cb.ranges[cyto.map[is.na(cyto.map$start)==FALSE,]$end]$band
	#cyto.map <- data.frame(name=cyto.map$name,cytoBandSpan=cyto.map$count,cytoBandStart=cyto.map$startname,cytoBandEnd=cyto.map$endname,stringsAsFactors=FALSE)

	#join cyto map data back into main annotation by id
	#kg.ann <- join(kg.ann,cyto.map,by="name",type="left")

	#add 1-based start coordinate column to remind viewer about this aspect of UCSC data
	kg.ann$txStart.1based <- kg.ann$txStart + 1
	kg.ann$cdsStart.1based <- kg.ann$cdsStart + 1

	kg.ann
}
parseExons <- function(ann)
{
	#split by chrs to make it go faster
	chrs <- unique(ann$chrom)

	#my.ann <- ann[ann$chrom=="chr1",]	

	# parse out exon lists to give regions list where each individual exon is a range

	#ex <- foreach(j=1:length(chrs))
	ex <- foreach(j=1:length(chrs),.verbose=TRUE,.combine="rbind") %dopar%
	{
		print(paste("Parsing exons for chr ",chrs[j],sep=""))
		my.ann <- ann[ann$chrom==chrs[j],]
		foreach(i=1:nrow(my.ann),.verbose=FALSE,.combine="rbind") %do%
		{
			chr <- my.ann[i,]$chrom
			starts <- as.numeric(unlist(strsplit(my.ann[i,]$exonStarts,",")))
			# correct for UCSC's 0-based system
			starts <- starts + 1
			ends <- as.numeric(unlist(strsplit(my.ann[i,]$exonEnds,",")))

			data.frame(chr=chr,start=starts,end=ends)
		}
	}


	ex
}
getGenicOverlap <- function(regions.ranges, ann)
{
	ann.ranges <- with(ann, GRanges(seqnames=chrom,ranges=IRanges(start=txStart.1based,end=txEnd)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}

getGenicOverlapGenes <- function(regions.ranges, ann)
{
	ann.ranges <- with(ann, GRanges(seqnames=chrom,ranges=IRanges(start=txStart.1based,end=txEnd)))
	overlap <- findOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0

	overlap <- as.data.frame(overlap)
	overlap$name <- ann[overlap$subjectHits,]$geneSymbol

	out <- foreach(i=1:length(regions.ranges),.verbose=FALSE,.combine="c") %do%
	{
		hits <- overlap[overlap$queryHits==i,]
		genes <- unique(hits$name)
		paste(genes,collapse=", ")
	}

	out
}
getUpstreamOverlap <- function(regions.ranges, ann, before=1000, after=500)
{
	# add offsets, accounting for strandedness
	ups <- with(ann,data.frame(chr=chrom,start=txStart.1based,end=txEnd,strand=strand))
	ups$start.us <- NA
	ups$end.us <- NA

	ups[ups$strand=="+",]$start.us <- ups[ups$strand=="+",]$start - before
	ups[ups$strand=="+",]$end.us <- ups[ups$strand=="+",]$start + after


	ups[ups$strand=="-",]$start.us <- ups[ups$strand=="-",]$end - after
	ups[ups$strand=="-",]$end.us <- ups[ups$strand=="-",]$end + before

	ann.ranges <- with(ups, GRanges(seqnames=chr,ranges=IRanges(start=start.us,end=end.us)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
getDownstreamOverlap <- function(regions.ranges, ann, before=500, after=1000)
{
	# add offsets, accounting for strandedness
	downs <- with(ann,data.frame(chr=chrom,start=txStart.1based,end=txEnd,strand=strand))
	downs$start.ds <- NA
	downs$end.ds <- NA

	downs[downs$strand=="+",]$start.ds <- downs[downs$strand=="+",]$end - before
	downs[downs$strand=="+",]$end.ds <- downs[downs$strand=="+",]$end + after


	downs[downs$strand=="-",]$start.ds <- downs[downs$strand=="-",]$start - after
	downs[downs$strand=="-",]$end.ds <- downs[downs$strand=="-",]$start + before

	ann.ranges <- with(downs, GRanges(seqnames=chr,ranges=IRanges(start=start.ds,end=end.ds)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
get3primeUTROverlap <- function(regions.ranges, ann)
{
	# filter for genes with a 3' UTR, accounting for strandedness
	utr <- with(ann,data.frame(chr=chrom, txStart=txStart.1based, txEnd=txEnd, strand=strand, cdsStart=cdsStart.1based, cdsEnd=cdsEnd))

	# filter non-coding transcripts which UCSC codes as cdsStart==cdsEnd
	utr <- utr[utr$cdsStart!=(utr$cdsEnd+1),]

	# filter out if cdsEnd == txEnd for (+) strand
	utr.p <- utr[(utr$strand=="+")&(utr$cdsEnd!=utr$txEnd),]

	# filter out if cdsStart == txStart for (-) strand because these are really the ends
	utr.m <- utr[(utr$strand=="-")&(utr$cdsStart!=utr$txStart),]

	# build list of utr regions
	utr.p$start.utr <- utr.p$cdsEnd
	utr.p$end.utr <- utr.p$txEnd

	utr.m$start.utr <- utr.m$txStart
	utr.m$end.utr <- utr.m$cdsStart

	reg <- rbind(utr.p,utr.m)

	ann.ranges <- with(reg, GRanges(seqnames=chr,ranges=IRanges(start=start.utr,end=end.utr)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
getExonOverlap <- function(regions.ranges, ann.ex)
{
	ann.ranges <- with(ann.ex, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
# =============================================================================