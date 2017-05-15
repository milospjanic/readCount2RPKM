#!/bin/bash

echo "Proccessing input mastertable from file" $1

cp $1 mastertable

cut -f1 mastertable > gene.name

###write R script for ID conversion, needs biomaRt

echo "#!/usr/bin/Rscript
library(biomaRt)
listMarts(host=\"grch37.ensembl.org\")

ensembl = useMart(\"ENSEMBL_MART_ENSEMBL\",dataset=\"mmusculus_gene_ensembl\", host=\"grch37.ensembl.org\")

id_merge_ens_length = getBM(attributes=c(\"ensembl_gene_id\",\"transcript_length\"),mart=ensembl)

write.table(id_merge_ens_length, file=\"id_merge_ens_length.txt\", sep = \"\t\", quote =F, col.names=F, row.names=F)

" > script.r

#run R script

chmod 755 script.r
./script.r

echo "

# Read first file and put first column into hash table h1 and second column into hash h2 with keys from first column 

NR==FNR { h1[\$1] = \$1; h2[\$1] = \$0; next; }

# Read $1 from second, third file etc, compare first column with hash h1, if positive add to hash column 2

BEGIN{OFS=\"\t\";}  NF{ if(\$1 in h1) h2[\$1] = h2[\$1] OFS \$2;}

# Read $1 and if it is empty append 0 to hash h2

NF{ if(\$1 in h1 && \$2==\"\") h2[\$1] = h2[\$1] OFS 0}

# At the end calculate lenght of the hash h2 if it is >ARGC-2 print, this eliminates truncated rows that couldnt be found in all files

END { for(k in h2)
  if(split(h2[k], a) > ARGC-2)
    print k OFS h2[k]
}
" > fileMulti2TableMod1.all.col.file1.awk

#remove duplicates from biomart output

sort -k1,1 -k2,2nr id_merge_ens_length.txt | awk '!a[$1]++' > id_merge_ens_length.unique.txt

#merge lengths to mastertable

awk -f fileMulti2TableMod1.all.col.file1.awk mastertable id_merge_ens_length.unique.txt | cut -f2- > mastertable_length

echo "Calculating RPKM values"

echo "#!/usr/bin/Rscript
x<-read.delim(\"mastertable_length\", header=F, na.strings = \"NA\")
x<-na.omit(x)
x\$V1<-as.character(x\$V1)
idx <- sapply(x, is.factor)
x[idx]<-lapply(x[idx], function(x) as.numeric(as.character(x)))
x[,2:(ncol(x)-1)]<-x[,2:(ncol(x)-1)]/x\$V10*1000
x[2:(ncol(x)-1)]<-sweep(x[2:(ncol(x)-1)],2,colSums(x[2:(ncol(x)-1)]),\`/\`)
x[2:(ncol(x)-1)]<-x[2:(ncol(x)-1)]*1000000
write.table(file=\"mastertable.RPKM\", x[1:(ncol(x)-1)], sep=\"\t\", quote=F, col.names = F, row.names = F)
" > script.2.r

echo "Writing RPKM values in mastertable.RPKM"

chmod 755 script.2.r
./script.2.r

rm script.2.r
rm id_merge_ens_length.unique.txt
rm fileMulti2TableMod1.all.col.file1.awk
rm id_merge_ens_length.txt
rm gene.name
rm script.r
