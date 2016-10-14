for dir in $(ls -d *); do the_file=${dir}/${dir}_12_summary.tab;    if [[ ! -f $the_file ]]; then       printf '%s found in %s:\n' "$dir"; fi; done > missing

x=read.table("missing", header=F)
x$eQTL = as.character(lapply(strsplit(as.character(x$V1), "_\\s*(?=[^_]+$)", perl=TRUE), "[", 2))
x$biom = as.character(lapply(strsplit(as.character(x$V1), "_\\s*(?=[^_]+$)", perl=TRUE), "[", 1))


