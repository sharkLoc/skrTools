#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "skr.h"

static int usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: skr (名字我瞎编的,嘿嘿！)\n");
	fprintf(stderr, "Usage: skr <command> [options]\n\n");
	fprintf(stderr, "\tmergeVcf\tmerge vcf files from list\n");
	fprintf(stderr, "\tfq2fa\ttranslate fastq file to fasta\n");
	fprintf(stderr, "\n");
	return 1;
}


int main(int argc, char *argv[])
{
	if(argc<2) return usage();
	if (strcmp(argv[1], "mergeVcf") == 0)
	{
		mergeVcfs(argc, argv);	
	}
	if(strcmp(argv[1],"fq2fa") == 0)
	{
		fq2fa(argc, argv);
	}
	exit(0);
}
