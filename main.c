#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "skr.h"

static int usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: skr \n");
	fprintf(stderr, "Usage: skr <command> [options]\n\n");
	fprintf(stderr, "\tfq2fa      translate fastq file to fasta\n");
	fprintf(stderr, "\tfqstat     summary statistics of fastq file\n");
	fprintf(stderr, "\tmergeVcf   merge vcf files from list\n");
	fprintf(stderr, "\tstatVcf    summary statistics of vcf file\n");
	fprintf(stderr, "\tmakewind   make bed from a list file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "\n");
	return 1;
}


int main(int argc, char *argv[])
{
	if(argc<2) return usage();
	if (strcmp(argv[1], "mergeVcf") == 0){
		mergeVcfs(argc, argv);	
	}
	else if(strcmp(argv[1],"fq2fa") == 0){
		fq2fa(argc, argv);
	}
	else if(strcmp(argv[1],"fqstat") == 0){
		fqstat(argc, argv);	
	}
	else if(strcmp(argv[1],"statVcf")==0){
		statVcfs(argc,argv);
	}
	else if(strcmp(argv[1],"makewind")==0){
		makewind(argc,argv);	
	}
	else{
		fprintf(stderr,"%s: invalid option -- '%s'\n",argv[0],argv[1]);	
	}

	exit(0);
}
