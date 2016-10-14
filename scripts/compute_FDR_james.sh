cat /sc/orga/projects/psychgen/resources/COLOC2/temp_results/results/SCZ_DLPFC/SCZ_DLPFC_summary.tab_coloc.supplied.var.sdY1.lkl | awk -v OFS='\t' '{if(NR == 1) {x ="PEP"}else{x =1-$32}{print x"\t"$0}}' | sort -k 1,1g | awk '{if(NR < 2){print "FDR""\t"$0}else{x+=$1; print x/NR"\t"$0}}' > /sc/orga/projects/psychgen/resources/COLOC2/temp_results/results/SCZ_DLPFC/try

# dim(x[x$FDR<0.05,])
