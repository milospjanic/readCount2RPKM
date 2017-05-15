# readCount2RPKM


readCount2RPKM is a combined bash/awk/R script for the conversion of mastertable of per-gene read counts from an RNA-Seq experiment to the mastertable of per-gene RPKM values. It connects to biomaRt to collect transcript lengths and selects the longest transcript to calculate RPKM values.
readCount2RPKM uses mastertable of per-gene read counts as input and allocates the longest transcript length to each gene, by writing and running an R script to connect to biomaRt to obtain transcript lengths. Then it creates extended mastertable by assigning maximum transcript lenght to each gene in the mastertable using custom made awk script fileMulti2TableMod1.all.col.file1.awk. It then writes and R script to convertread count to RPKM values with the following R code:

```R
x[,2:(ncol(x)-1)]<-x[,2:(ncol(x)-1)]/x\$V10*1000
x[2:(ncol(x)-1)]<-sweep(x[2:(ncol(x)-1)],2,colSums(x[2:(ncol(x)-1)]),\`/\`)
x[2:(ncol(x)-1)]<-x[2:(ncol(x)-1)]*1000000
```

RPKMs are calculated using the formula:

</pre>
reads per kb = (read count / gene length) * 1000
reads per kb per million mapped reads (RPKM) = (reads per kb / per sample total read count) * 1,000,000
<pre>

# Usage


readCount2RPKM will 
