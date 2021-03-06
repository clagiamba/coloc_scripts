match_snpid = F
submit = TRUE
#submit=F
GWAS_DIR = "/sc/orga/projects/psychgen/resources/COLOC2/files/GWAS/"
EQTL_DIR = "/sc/orga/projects/psychgen/resources/COLOC2/files/eQTLs2/"
#*#EQTL_DIR = "/sc/orga/projects/epigenAD/coloc/coloc_STARNET/eQTLs/"
OUT_DIR = "/sc/orga/projects/psychgen/resources/COLOC2/temp_results"
#*#OUT_DIR = "/sc/orga/projects/epigenAD/coloc/coloc_STARNET/"
## OUT_DIR = "/sc/orga/projects/epigenAD/coloc/coloc_SCZ/"

# case1:
# eqtl.fname = "/sc/orga/projects/epigenAD/coloc/data/AJAFFEmeQTL/rawData/meQTLs_adults_cis_allPvalues_formatted.txt"
# biom.fname="/sc/orga/projects/psychgen/resources/COLOC2/files/GWAS/SCZ_formatted"
# case2:
# eqtl.fname = "/sc/orga/projects/psychgen/resources/COLOC2/files/eQTLs2/DLPFC_formatted"
# biom.fname= "/sc/orga/projects/epigenAD/coloc/data/AJAFFEmeQTL/rawData/meQTLs_adults_cis_allPvalues_formatted.txt"
eqtl_files= list.files(EQTL_DIR, full.names=T, pattern="_formatted$")
biom_files = list.files(GWAS_DIR, full.names=T, pattern="_formatted$")

table_pairs = expand.grid(biom_files, eqtl_files)
names(table_pairs)=c("biom_files", "eqtl_files")

outmain = paste(OUT_DIR, "/results/", sep="")
      if (!file.exists (outmain)) dir.create(outmain, recursive=TRUE)
      if (!file.exists (paste(OUT_DIR, "/scripts/", sep="")))  dir.create(paste(OUT_DIR, "/scripts/", sep=""))
      if (!file.exists (paste(OUT_DIR, "/log/", sep="")))  dir.create(paste(OUT_DIR, "/log/", sep=""))

main_script = '/sc/orga/projects/epigenAD/coloc/coloc2_gitrepo/coloc_scripts/scripts/new_coloc_SMR_paper_formatted.R'

# do not run eRNA
table_pairs = table_pairs[!grepl("CMC_eGenes_1Mb|CMC_eRNA_1Mb", table_pairs$eqtl_files),]

for (i in 1:nrow(table_pairs)) {

biom.fname = as.character(table_pairs$biom_files[i])
eqtl.fname = as.character(table_pairs$eqtl_files[i])

biom.name = gsub("_formatted", "", basename(biom.fname))
eqtl.name = gsub("_formatted|_Analysis_cis-eQTLs.coloc.txt_formatted", "", basename(eqtl.fname))
#eqtl.name = as.character(lapply(strsplit(as.character(eqtl.fname), "/", fixed=TRUE), "[", 9))
message("Using biom ", biom.name, " and eQTL ", eqtl.name)
prefix = paste(biom.name, eqtl.name, sep="_")
outfolder = paste(outmain, prefix, "/", sep="")

# alloc
# premium
scriptname=paste(OUT_DIR, "/scripts/Submit_main_script_", prefix, ".sh", sep="")

     write(file=scriptname, paste("#!/bin/bash
          #BSUB -J coloc_", prefix,"
          #BSUB -q alloc
          #BSUB -P acc_psychgen 
          #BSUB -n 20
          #BSUB -R span[hosts=1]
          #BSUB -R 'rusage[mem=5000]'
          #BSUB -W 60:00 
          #BSUB -L /bin/bash
          #BSUB -oo ", OUT_DIR, "/log/", prefix, ".out
          #BSUB -eo ", OUT_DIR, "/log/", prefix, ".err", sep=""), append=F)
          write(file=scriptname, paste("Rscript", main_script, biom.fname, eqtl.fname, outfolder, prefix,sep=" "), append=T)

          message("Submit script: ", scriptname)
          if (submit) {
           system(paste("bsub < ", scriptname, sep=""))
          }
}
 
