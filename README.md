# skrTools
> a little program

# install
```
git clone https://github.com/sharkLoc/skrTools.git
cd skrTools && make 
```

# usage
```
Program: skr 

Usage: skr <command> [options]

	fq2fa      translate fastq file to fasta
	fqstat     summary statistics of fastq file
	mergeVcf   merge vcf files from list
	statVcf    summary statistics of vcf file
	makewind   make bed from a list file
```



## skr fqstat -i xx1.fq.gz -I xx2.fq.gz
#### output:
```
Iterm	reads_1.fq	reads_2.fq
read average length:	150	150
read GC content(%):	48.42	48.48
total read Count:	34946389	34946389
total base Count:	5241958350	5241958350

base A Count:	1352284833(25.80%)	1342903044(25.62%)
base C Count:	1270459966(24.24%)	1246706604(23.78%)
base G Count:	1267522866(24.18%)	1294357728(24.69%)
base T Count:	1351401800(25.78%)	1357986115(25.91%)
base N Count:	288885(0.01%)	4859(0.00%)

Number of base calls with quality value of 20 or higher (Q20+) (%)	5113248711(97.54%)	5092440219(97.15%)
Number of base calls with quality value of 30 or higher (Q30+) (%)	4886887711(93.23%)	4832524601(92.19%)
```
