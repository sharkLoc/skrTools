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
				fprintf(stderr, "            -l :<char>    vcf list file, one file per line\n");
				fprintf(stderr, "            -o :<char>    merged vcf[.gz] file name\n");
				fprintf(stderr, "            -h :<char>    show this help\n\n");
				break;
		}
	}
	if(argc<=2)
	{
		fprintf(stderr, "mergeVcf    merge vcf files from list\n\n");
		fprintf(stderr, "            -l :<char>    vcf list file, one file per line\n");
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
