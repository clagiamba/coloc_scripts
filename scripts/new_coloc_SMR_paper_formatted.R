args <- commandArgs(trailingOnly=TRUE)
biom.fname=args[1]
eqtl.fname=args[2]
outfolder = args[3]
prefix = args[4]


source("/sc/orga/projects/epigenAD/coloc/coloc2_gitrepo/coloc_scripts/scripts/functions_coloc_pipeline.R")
source("/sc/orga/projects/epigenAD/coloc/coloc2_gitrepo/coloc_scripts/scripts/functions_coloc_likelihood_summary_integrated.R")
#source("/sc/orga/projects/epigenAD/coloc/coloc2_gitrepo/coloc_scripts/scripts/functions_coloc_likelihood_summary_integrated_add_beta.R")

# filter by maf and imp quality when doing colocalization
# liver eQTL is messed up: data= fread(fname, select = colsAll[[1]], col.names=names(colsAll[[1]]), colClasses=c(SNPID="character", ProbeID="character", BETA="numeric", PVAL="numeric", CHR="character", POS="numeric", F="numeric", SE="numeric", N="numeric"))
# data= fread(fname, stringsAsFactors=FALSE)
#biom.df = formatColoc(fname = biom.fname, type=type, N=NA, Ncases=NA, info_filter=0, maf_filter=0, fread=T, eqtl=FALSE)
#eqtl.df = formatColoc(fname = eqtl.fname, type="quant", N=NA, Ncases=NA, info_filter=0, maf_filter=0, fread=T, eqtl=TRUE)

 library(data.table)
 biom.df = data.frame(fread(biom.fname, header=T))
 eqtl.df = data.frame(fread(eqtl.fname, header=T))

 p12=1e-6 # this is for coloc with set priors
 plot = FALSE
 useBETA=TRUE
 save.coloc.output=FALSE

 # estimate per locus likelihoods (for now also outputs the ppa using set priors
 res1 = coloc.eqtl.biom(eqtl.df=eqtl.df, biom.df=biom.df, p12=p12, useBETA=useBETA, outfolder=outfolder, prefix=prefix, plot=plot, save.coloc.output=save.coloc.output, match_snpid=FALSE)
