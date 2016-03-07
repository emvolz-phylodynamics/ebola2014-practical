NREPS <- 1e2 
ODIR <- 'bootstrap/'

# NOTE will only save data for bootstrap trees where likelihood can be evaluated at starting conditions for optimisation

require(phydynR)
require(phangorn)
d <- read.phyDat( 'ebov_panafr_082014.fasta' , format = 'fasta')
bsseqs <- bootstrap.phyDat( d,bs = NREPS, FUN=function(x) as.DNAbin(x) )

# make lsd 
lsd.start.conditions <- function( intree, sampleTimes, sequenceLength, minEdgeLength = 1e-3)
{
	print('NOTE: this version of lsd.initial.conditions is designed to work with v0.2 of lsd. May give erroneous or unexpected results with different versions.')
	DIR <- tempdir() # NOTE difficult to find an option that works on HPC
	lsdfn <- tempfile(pattern = "lsd", tmpdir = DIR, fileext = ".lsd")
	intreefn <- tempfile(pattern = "intree", tmpdir = DIR, fileext = ".nwk")
	stfn <- tempfile(pattern = "st", tmpdir = DIR, fileext = ".tsv")
	lsdtreefn  <- paste(sep='.', lsdfn, 'date', 'newick' )
	lsdtree_evo_fn  <- paste(sep='.', lsdfn,  'newick' )
	
	write.tree( intree, file = intreefn )
	n <- length(intree$tip.label)
	write.table( file = stfn,  sampleTimes, row.names=TRUE, col.names=c( n) , quote=FALSE)
	
	LSDCOMMAND <- paste(sep=' ', 'lsd', '-i',  intreefn,  '-d',  stfn,  '-o', lsdfn,  '-c -r a -s ', sequenceLength )
	
	o <- system('lsd -h')
	if (o != 0){
		stop('There is a problem with the least-squares-dating installation. Returning NULL. ')
	}
	system(LSDCOMMAND , wait = TRUE)
	lsdrate <- scan( lsdfn, what = character(0))[101] ##TODO should do a more accurate scan to extract rate
	lsdrate <- as.numeric( substr( lsdrate, 1, nchar(lsdrate)-1) )
	
	lsdtree <- read.tree(lsdtreefn)
	#lsdtree$edge.length <- pmax( minEdgeLength, lsdtree$edge.length )
	bdt <- DatedTree( lsdtree, sampleTimes , tol = Inf, minEdgeLength = minEdgeLength)
	
	lsdtree_evo <- read.tree( lsdtree_evo_fn ) 
	
	# internal node order is precalc in DatedTree
	#ino <- n + sort(bdt$heights[(n+1):(n+bdt$Nnode)], index.return=TRUE)$ix
	ino <- bdt$eventIndicatorNode[ bdt$events==1 ] 
	nh  <- bdt$heights
	
	print(lsdrate)
	list( bdt,  nh , ino, lsdrate, lsdfn, lsdtree_evo )
}


sampleDates_table <- read.csv( 'ebov_sampleDates.csv' ) 
sampleDates <- setNames( sampleDates_table$year, sampleDates_table[,1]) # we will also need this in vector format


# start conditions for optimisation
t0 <- 2014
I0 <- .1 # initial conditions 
births <- c( I = 'parms$beta * I' )
deaths <- c( I = 'parms$gamma * I' )
GAMMA <- 365 / 15 
theta_start <- c( beta = 2 * GAMMA, gamma = GAMMA) # starting conditions for maximum likelihood
ebov_model1 <- build.demographic.process(births=births
  , deaths = deaths
  , parameterNames=c('beta', 'gamma') 
  , rcpp=FALSE # specifies that equations are R-code as opposed to C-code
)

for ( imsa in 1:length( bsseqs )){
	ebov_algn <- bsseqs[[imsa]] 
	D <- dist.dna( ebov_algn, model = 'F84', pairwise.deletion=TRUE, as.matrix=TRUE)
	ebov_nj <- nj( D )
	o_lsd <- lsd.start.conditions( ebov_nj, sampleDates, sequenceLength = 18961)
	lsdtree <- o_lsd[[1]]
	toremove <- lsdtree$tip.label[ !grepl( 'SierraLeone_G',  lsdtree$tip.label) ]
	sl_lsdtree <- drop.tip ( lsdtree, toremove ) 
	sl_lsdtree$edge.length <- pmax (0, sl_lsdtree$edge.length )
	sl_sampleTimes <- sampleDates[ sl_lsdtree$tip.label ]
	sl_lsdtree <- DatedTree( sl_lsdtree, sampleTimes = sl_sampleTimes, tol = Inf, minEdgeLength = 1e-3 )
	
	ll <- colik( sl_lsdtree
	  , theta_start 
	  , ebov_model1
	  , x0 = c( I = I0 )
	  , t0 = t0
	)
	print( paste(date(), ll ))
	
	if (!is.infinite(ll)){
		write.dna( bsseqs[[imsa]], paste(sep='', ODIR, 'ebov_panafr_bs', imsa, '.interleaved' ) )
		write.tree( sl_lsdtree, file = paste(sep='', ODIR, 'sl_lsdtree', imsa, '.nwk') )
		write.tree( lsdtree, file = paste(sep='', ODIR, 'ebov_lsdtree', imsa, '.nwk') )
	}
}
