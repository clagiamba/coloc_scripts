args <- commandArgs(trailingOnly=TRUE)
biom.fname=args[1]
eqtl.fname=args[2]
outfolder = args[3]
prefix = args[4]


library(data.table)
# Add  N MAF BETA SE P

 # optimization and posteriors using these priors estimates
 outfname = paste(outfolder, prefix, '_summary.tab', sep='')

if (!file.exists(outfname)) {
     stop(outfname, " does not exist")
}

   res.all = data.frame(fread(outfname, header=T, stringsAsFactors=F))
   biom.df = data.frame(fread(biom.fname, header=T, stringsAsFactors=F))
   eqtl.df = data.frame(fread(eqtl.fname, header=T, stringsAsFactors=F))

   eqtl.df$chrpos= paste(eqtl.df$CHR, eqtl.df$POS, sep=":")
   biom.df$chrpos= paste(biom.df$CHR, biom.df$POS, sep=":")

   res.all$N.snp.biom = biom.df[match(res.all$snp.biom, biom.df$chrpos), "N"]
   res.all$N.snp.eqtl = eqtl.df[match(res.all$snp.eqtl, eqtl.df$chrpos), "N"]
   if ("Ncases" %in% names(biom.df)) {
      res.all$Ncases.snp.biom = biom.df[match(res.all$snp.biom, biom.df$chrpos), "Ncases"]
   }

   res.all$BETA.snp.biom = biom.df[match(res.all$snp.biom, biom.df$chrpos), "BETA"]
   res.all$BETA.snp.eqtl = eqtl.df[match(res.all$snp.eqtl, eqtl.df$chrpos), "BETA"]
   res.all$SE.snp.biom = biom.df[match(res.all$snp.biom, biom.df$chrpos), "SE"]
   res.all$SE.snp.eqtl = eqtl.df[match(res.all$snp.eqtl, eqtl.df$chrpos), "SE"]


# MAF: take it from the eqtl first, if no MAF then go to biom
maf.eqtl = ifelse(any(c("MAF","F") %in% names(eqtl.df)), TRUE, FALSE)
maf.biom = ifelse(any(c("MAF","F") %in% names(biom.df)), TRUE, FALSE)
if (!maf.eqtl & !maf.biom) {
  message("There is no MAF information in neither datasets, must use external")
  }

if (maf.eqtl) {
   if ("F" %in% names(eqtl.df) & !("MAF" %in% names(eqtl.df))){
     eqtl.df$MAF = ifelse(eqtl.df$F<0.5, eqtl.df$F, 1-eqtl.df$F)
     res.all$MAF.snp.eqtl = eqtl.df[match(res.all$snp.eqtl, eqtl.df$chrpos), "MAF"]
     }
    }
if (maf.biom) {
   if ("F" %in% names(biom.df) & !("MAF" %in% names(biom.df))){
     biom.df$MAF = ifelse(biom.df$F<0.5, biom.df$F, 1-biom.df$F)
     res.all$MAF.snp.biom = biom.df[match(res.all$snp.biom, biom.df$chrpos), "MAF"]
     }
   }

if( !maf.eqtl & maf.biom ){
     res.all$MAF.snp.eqtl = biom.df[match(res.all$snp.eqtl, biom.df$chrpos), "MAF"]
   }
if( !maf.biom & maf.eqtl ){
     res.all$MAF.snp.biom = biom.df[match(res.all$snp.eqtl, eqtl.df$chrpos), "MAF"]
}

indx <- apply(res.all[,c("N.snp.eqtl", "N.snp.biom", "BETA.snp.biom", "BETA.snp.eqtl", "SE.snp.biom", "SE.snp.eqtl", "MAF.snp.eqtl", "MAF.snp.biom")], 2, function(x) any(is.na(x) | is.infinite(x)))

 if (any(indx)) {
    message("Some values are missing")
    print(colnames[indx])
}


   outfname.new = paste(outfolder, prefix, '_summary_addinfo.tab', sep='')
   #res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
   write.table(x =  res.all , file = outfname.new, row.names = FALSE, quote = FALSE, sep = '\t')

   message("File saved in ", outfname.new)

