# readCount2RPKM


readCount2RPKM is a combined bash/awk/R script for the conversion of mastertable of per-gene read counts from an RNA-Seq experiment to the mastertable of per-gene RPKM values. It connects to biomaRt to collect transcript lengths and selects the longest transcript to calculate RPKM values.

readCount2RPKM will 
