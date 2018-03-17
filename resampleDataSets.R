NREPS <- 30 
SAMPLESIZE <- 300
ODIR <- 'resamples/'


#~ MK2008_2015-02-12

require(ape)
require(lubridate)
d <- read.dna('data/Makona_1610_cds_ig.fas' , format = 'fasta' )

aln2sampleTimes <- function( aln){
	nms <- rownames(aln)
	dtstrs <- sapply( strsplit( nms, split='_'), 'tail',1)
	dts <- strptime( dtstrs, format = '%Y-%m-%d')
	data.frame( sampleTimes=setNames( decimal_date(dts), nms ) )
}
 
for (i in 1:NREPS){
	aln <- d[ sample(1:nrow(d), size=SAMPLESIZE, replace=FALSE) , ]
	sts <- aln2sampleTimes( aln )
	write.csv( sts, file = paste0(ODIR, 'sts-', i, '.csv') )
	write.dna( aln, file = paste0(ODIR, 'aln-', i, '.fas') )
}
