args <- commandArgs(trailingOnly=TRUE)
biom.fname=args[1]
eqtl.fname=args[2]
outfolder = args[3]
prefix = args[4]


source("/sc/orga/projects/epigenAD/coloc/coloc2_gitrepo/coloc_scripts/scripts/functions_coloc_likelihood_summary_integrated.R")

 # optimization and posteriors using these priors estimates
 outfname = paste(outfolder, prefix, '_summary.tab', sep='')
 if (!(file.exists(outfname)) & !all(file.exists(paste( gsub("_summary.tab", "", outfname), 1:22, "summary.tab", sep="_")))) {
     write(outfname, file = "/sc/orga/projects/psychgen/resources/COLOC2/temp_results/results/missingAllChr.tab", append = TRUE, sep = "\t")
     stop(outfname, " does not exist")
  }

# if (file.exists(outfname)) {
#    res.all = read.table(outfname, header=T)
# }
if (all(file.exists(paste( gsub("_summary.tab", "", outfname), 1:22, "summary.tab", sep="_")))) {
     library(plyr)
     res.all = ldply(lapply(paste( gsub("_summary.tab", "", outfname), 1:22, "summary.tab", sep="_"), function(i){read.table(i, header=TRUE)}), data.frame)
     table(res.all$Chr);
     write.table(res.all, file=outfname, row.names = FALSE, quote = FALSE, col.names = TRUE)
  } 

if (file.exists(outfname) & !all(file.exists(paste( gsub("_summary.tab", "", outfname), 1:22, "summary.tab", sep="_")))) {
    res.all = read.table(outfname, header=T)
}

# this is only for now to test different models
 models = grep(".lkl", names(res.all), value=T)
 models = models[-grep("set.priors", models)]
 for (i in 1:length(models)) {
   col.with.lkl = models[i] # e.g. "coloc.old.pval.lkl"
   # split the lklds
   # out <- strsplit(as.character(res.all$coloc.old.pval.lkl),',')
   out <- strsplit(as.character(res.all[,col.with.lkl]),',')
# coloc.old.var.lkl
# coloc.supplied.var.lkl
# coloc.var.Neff.lkl
   t= data.frame(res.all, do.call(rbind, out))
   names(t)[(ncol(t)-4):ncol(t)] = c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf") 
   res2 = est_lkl(t, outfolder=outfolder, bootstrap=F, no_bootstraps=1000)

   res2 = addGeneNames(res2, biomart=FALSE, geneFileNames = "/sc/orga/projects/psychgen/resources/COLOC2/data/ENSEMBL_v70_TO_HGNC.tsv")

   optim.res =  paste(outfolder, 'maximization_results.txt', sep='')
   optim.res.new =  paste(outfolder, 'maximization_results_', col.with.lkl, '.txt', sep='')
   file.rename(optim.res, optim.res.new)

   outfname.new = paste(outfname, col.with.lkl, sep="_")
   #res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
 write.table(x =  res2 , file = outfname.new, row.names = FALSE, quote = FALSE, sep = '\t')
}

