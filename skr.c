#include <zlib.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>

#include "skr.h"

char *readline(gzFile file)
{
	size_t baselen = 256;
	char *line = (char *)calloc(baselen,sizeof(char));
	if(!line) exit(1);

	int ch;
	size_t index = 0;
	while((ch=gzgetc(file)) != -1 && ch != 10)
	{
		line[index] = ch;
		index++;
		if(index == baselen)
		{
			baselen += 128;
			line=(char *)realloc(line,baselen);
			if(!line)
			{
				free(line);
				exit(1);
			}
		}
	}
	line[index] = '\0';  // tail '\0' 
	if(ch == -1)  return NULL;// end of file
	return line;
}

void mergeVcfs(int argc, char *argv[])
{
	gzFile fp, fo; int r;
	while((r=getopt(argc,argv,"l:o:h")) != -1)
	{
		switch(r)
		{
			case 'l': fp = gzopen(optarg,"r"); break;
			case 'o': fo = gzopen(optarg,"w"); break;
			case 'h': 
			case '?':
				fprintf(stderr, "mergeVcf    merge vcf files from list\n\n");
				fprintf(stderr, "            -l :<char>    vcf list file, one file per line and no blank line\n");
				fprintf(stderr, "            -o :<char>    merged vcf[.gz] file name\n");
				fprintf(stderr, "            -h :<char>    show this help\n\n");
				break;
		}
	}
	if(argc<=2)
	{
		fprintf(stderr, "mergeVcf    merge vcf files from list\n\n");
		fprintf(stderr, "            -l :<char>    vcf list file, one file per line and no blank line\n");
		fprintf(stderr, "            -o :<char>    merged vcf[.gz] file name\n");
		fprintf(stderr, "            -h :<char>    show this help\n\n");
	}
	if(!fp || !fo) exit(1);

	char *list=NULL; 
	int i=0;  unsigned long long count=0;
	while((list=readline(fp)) != NULL)
	{
		i++;
		gzFile vcf = gzopen(list,"r");
		if(!vcf) exit(1);
		char *line=NULL;
		
		while((line=readline(vcf)) != NULL)
		{
			if(1==1)
			{ 
				if(line[0] != '#') count++;
				gzprintf(fo,"%s\n",line);
			}
			else
			{
				if(line[0]=='#')
				{
					continue;
				}
				count++;
				gzprintf(fo,"%s\n",line);
			}
			free(line);
		}
		gzclose(vcf);
		free(list);
	}
	gzclose(fp);
	gzclose(fo);
	fprintf(stderr,"note: total %d files merged, and %llu SNPs!\n",i,count);
	return ;
}

void fq2fa(int argc, char *argv[])
{
	int r,flag=1; gzFile fp,fo;
	while((r=getopt(argc,argv,"i:o:h"))!= -1)
	{
		switch(r)
		{
			case 'i': fp = gzopen(optarg,"r"); break;
			case 'o': fo = gzopen(optarg,"w"); break;
			case 'h':
			case '?':
				fprintf(stderr, "fq2fq       merge vcf files from list\n\n");
				fprintf(stderr, "            -i :<char>    fastq file\n");
				fprintf(stderr, "            -o :<char>    fasta file[.gz]\n");
				fprintf(stderr, "            -h :<char>    show this help\n\n");

				break;
		}
	}
	if(argc<=2)
	{
		fprintf(stderr, "fq2fq       merge vcf files from list\n\n");
		fprintf(stderr, "            -i :<char>    fastq file\n");
		fprintf(stderr, "            -o :<char>    fasta file[.gz]\n");
		fprintf(stderr, "            -h :<char>    show this help\n\n");
	}
	if(!fp || !fo) exit(1);

	while(flag)
	{
		for(int i=0; i<4; i++)
		{
			char *read=readline(fp);
			if(read)
			{
				if(i==0)
				{
					memcpy(read,">",1);
					gzprintf(fo,"%s\n",read);
				}
				if(i==1)
				{
					gzprintf(fo,"%s\n",read);
				}
			}
			else
			{
				flag=0;
			}
			free(read);
		}
	}
	gzclose(fp);
	gzclose(fo);
	return ;
}

void statVcfs(int argc, char *argv[])
{
	int r; gzFile fp; FILE *fo;
	unsigned long long count = 0;
	while((r=getopt(argc,argv,"i:o:h"))!= -1)
	{
		switch(r)
		{
			case 'i': fp = gzopen(optarg,"r"); break;
			case 'o': fo = fopen(optarg,"w"); break;
			case 'h':
			case '?':
				fprintf(stderr, "statVcf    state vcf file\n\n");
				fprintf(stderr, "            -i :<char>    input vcf file\n");
				fprintf(stderr, "            -o :<char>    output state file name\n");
				fprintf(stderr, "            -h :<char>    show this help\n\n");
				break;
		}
	}
	if(argc<=2)
	{
		fprintf(stderr, "statVcf    state vcf file\n\n");
		fprintf(stderr, "            -i :<char>    input vcf file\n");
		fprintf(stderr, "            -o :<char>    output state file name\n");
		fprintf(stderr, "            -h :<char>    show this help\n\n");
	}
	if(!fp || !fo) exit(1);
	struct GT {
		unsigned long long AT; unsigned long long AG; unsigned long long AC;
		unsigned long long TA; unsigned long long TG; unsigned long long TC;
		unsigned long long GA; unsigned long long GT; unsigned long long GC;
		unsigned long long CA; unsigned long long CT; unsigned long long CG;
	};
	struct GT geno ={0,0,0,0,0,0,0,0,0,0,0,0};
	char *tmp,*line=NULL; char ref,alt;
	while((line=readline(fp)) != NULL)
	{
		if(line[0]=='#')
		{
			continue;
		}		
		else
		{
			count++;
			int i = 0;
			tmp=strtok(line,"\t");
			while(tmp!= NULL && i<=4) 
			{
				i++;
				tmp=strtok(NULL,"\t");
				if(i==3) ref=*tmp;
				if(i==4) alt=*tmp;
			}
			if(ref =='A' && alt == 'T')	geno.AT++;
			if(ref =='A' && alt == 'G') geno.AG++;
			if(ref =='A' && alt == 'C') geno.AC++;
			if(ref =='T' && alt == 'A') geno.TA++;
			if(ref =='T' && alt == 'G') geno.TG++;
			if(ref =='T' && alt == 'C') geno.TC++;
			if(ref =='G' && alt == 'A') geno.GA++;
			if(ref =='G' && alt == 'T') geno.GT++;
			if(ref =='G' && alt == 'C') geno.GC++;
			if(ref =='C' && alt == 'A') geno.CA++;
			if(ref =='C' && alt == 'T') geno.CT++;
			if(ref =='C' && alt == 'G') geno.CG++;

		}
		free(line);
	}
	unsigned long long transversion = geno.GT + geno.GC + geno.CA + geno.CG + \
									  geno.AT + geno.AC + geno.TA + geno.TG;
	unsigned long long transition = geno.GA + geno.CT + geno.AG + geno.TC;
	fprintf(fo,"SNP AT: %llu\n",geno.AT); fprintf(fo,"SNP AG: %llu\n",geno.AG); fprintf(fo,"SNP AC: %llu\n",geno.AC);
	fprintf(fo,"SNP TA: %llu\n",geno.TA); fprintf(fo,"SNP TG: %llu\n",geno.TG); fprintf(fo,"SNP TC: %llu\n",geno.TC);
	fprintf(fo,"SNP GA: %llu\n",geno.GA); fprintf(fo,"SNP GT: %llu\n",geno.GT); fprintf(fo,"SNP GC: %llu\n",geno.GC);
	fprintf(fo,"SNP CA: %llu\n",geno.CA); fprintf(fo,"SNP CT: %llu\n",geno.CT); fprintf(fo,"SNP CG: %llu\n",geno.CG);
	fprintf(fo,"SNP transition: %llu\n",transition);
	fprintf(fo,"SNP transversion: %llu\n",transversion);
	fprintf(fo,"SNP count: %llu\n",count);

	gzclose(fp);
	fclose(fo);
	return ;
}
