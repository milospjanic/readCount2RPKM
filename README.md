# readCount2RPKM


readCount2RPKM is a combined bash/awk/R script for the conversion of mastertable of per-gene read counts from an RNA-Seq experiment to the mastertable of per-gene RPKM values. It connects to biomaRt to collect transcript lengths and selects the longest transcript to calculate RPKM values.
readCount2RPKM uses mastertable of per-gene read counts as input and allocates the longest transcript length to each gene, by writing and running an R script to connect to biomaRt to obtain transcript lengths. Then it creates extended mastertable by assigning maximum transcript lenght to each gene in the mastertable, by writing custom made awk script fileMulti2TableMod1.all.col.file1.awk. It then writes and R script to convertread count to RPKM values with the following R code:

```R
x[,2:(ncol(x)-1)]<-x[,2:(ncol(x)-1)]/x\$V10*1000
x[2:(ncol(x)-1)]<-sweep(x[2:(ncol(x)-1)],2,colSums(x[2:(ncol(x)-1)]),\`/\`)
x[2:(ncol(x)-1)]<-x[2:(ncol(x)-1)]*1000000
```

RPKMs are calculated using the formula:

<pre>
reads per kb = (read count / gene length) * 1000
reads per kb per million mapped reads (RPKM) = (reads per kb / per sample total read count) * 1,000,000
</pre>

# Dependencies

readCount2RPKM requires awk, R and Rscript to run intermediary awk and R scripts, and R library biomaRt installed in R. 

# Usage

readCount2RPKM.mm.sh requires mastertable of mouse genes in format:

<pre>
ENSMUSG00000089707	22	0	0	0	0	0	10	0
ENSMUSG00000073747	1	0	0	0	1	0	8	0
ENSMUSG00000035000	1291	672	1219	1300	1308	683	563	2192
ENSMUSG00000104919	0	0	0	0	0	0	0	0
ENSMUSG00000088677	2	0	0	0	2	0	0	0
ENSMUSG00000103889	0	0	0	0	0	0	0	0
ENSMUSG00000081037	4	0	0	0	2	0	0	0
ENSMUSG00000075164	0	0	0	0	0	0	0	0
ENSMUSG00000030510	22	0	4	6	8	4	2	0
</pre>

Example run:

<pre>
./readCount2RPKM mastertable.txt 
Proccessing input mastertable from file mastertable.txt
Loading required package: methods
               biomart            version
1 ENSEMBL_MART_ENSEMBL      Ensembl Genes
2     ENSEMBL_MART_SNP  Ensembl Variation
3 ENSEMBL_MART_FUNCGEN Ensembl Regulation
4    ENSEMBL_MART_VEGA               Vega
Calculating RPKM values
Writing RPKM values in mastertable.RPKM
</pre>

Example output:

<pre>
head mastertable.RPKM 
ENSMUSG00000053783	0	0	0	0	0	0	0	0
ENSMUSG00000078607	14.3315290267182	5.27148998470254	9.47657223005189	14.1782105971938	25.9570809377262	12.5352263065328	3.49215569538464	8.12296582986335
ENSMUSG00000021900	65.6250210929769	22.948792387045	60.228922019752	38.7788513986779	52.1120026703478	57.297708533105	19.665464723166	62.9006962103502
ENSMUSG00000021901	63.5325181369401	29.4644784214023	66.5475638175278	26.4470153993863	29.9656600887317	19.5248137747344	9.16367902027884	29.6891300712868
ENSMUSG00000081820	0	0	0	0	0	0	0.111648246991865	0
ENSMUSG00000021902	5.34723506150612	5.33759154872005	9.07451378908043	8.48829107065946	9.13791262030437	7.06216010001703	2.64534314243628	3.86463200913539
ENSMUSG00000021903	15.8545088541307	9.74915150330365	22.1922681668968	0.63286285488678	0.530720326421349	3.22314506289442	3.27454460225776	11.8394474075886
ENSMUSG00000081821	0	0	0	0	0	0	0	0
ENSMUSG00000019710	26.4100432213795	21.5133687933026	30.3889090935622	30.3215943263374	28.0822793064558	19.990967660423	28.4085170598329	24.6477111636308
ENSMUSG00000080500	0	0	0	0	0	0	0	0
</pre>
