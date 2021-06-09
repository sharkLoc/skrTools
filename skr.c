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
		exit(1);
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
			if(i==1)
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
		fprintf(stderr,"file %s done !\n",list);
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
		exit(1);
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
		exit(1);
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
			if(ref =='A' && alt == 'T') geno.AT++;
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

void makewind(int argc, char *argv[])
{
	gzFile fp; int r;  
	long int w=10000,i,tmp;
	while((r=getopt(argc,argv,"i:w:h")) != -1)
	{
		switch(r)
		{
			case 'i': fp =gzopen(optarg,"r"); break;
			case 'w': w = atol(optarg); break;
			case 'h':
			case '?':
				fprintf(stderr, "makewind	make bed region from chrlen list\n\n");
				fprintf(stderr, "           -i  :<char> chrlen file, two columns per line are split by tabs and no blank line\n");
				fprintf(stderr, "                       example: chr1\t15000000\n");
				fprintf(stderr, "                                chr2\t12000000\n\n");
				fprintf(stderr, "           -w  :<int>	window length, default[%ld]\n",w);
				fprintf(stderr, "           -h  :<char>	show this help\n\n");
				break;
		}
	}
	if(argc<=2)
	{
		fprintf(stderr, "makewind   make bed region from chrlen list\n\n");
		fprintf(stderr, "           -i  :<char> chrlen file, two columns per line are split by tabs and no blank line\n");
		fprintf(stderr, "                       example: chr1\t15000000\n");
		fprintf(stderr, "                                chr2\t12000000\n\n");
		fprintf(stderr, "           -w  :<int>  window length, default[%ld]\n",w);
		fprintf(stderr, "           -h  :<char> show this help\n\n");
		exit(1);
	}
	if(!fp) exit(1);

	char *line, *col1, *col2 = NULL;
	while((line=readline(fp)) != NULL)
	{
		
		col1=strtok(line,"\t");
		col2=strtok(NULL,"\t");
		for(i=1; i<=atol(col2); i+=w)
		{
			tmp = i + w - 1;
			if(tmp>atol(col2))
				tmp = atol(col2);
			fprintf(stdout,"%s\t%ld\t%ld\n",col1,i,tmp);
		}
		
		free(line);
	}
	gzclose(fp);
	return ;	
}

uint64_t qualityN(char *line,size_t q)
{
	int i=0; uint64_t qCount=0;
	while(line[i]!='\0')
	{
		if(line[i]-33>=q) qCount++;
		i++;
	}
	return qCount;
}

void freeRead(rinfo read)
{
	for(int i=0; i<4; i++)
	{
		free(read.line[i]);
		read.line[i]=NULL;
	}
	return ;
}

void fqstat(int argc, char *argv[])
{
	gzFile fp, fd;
	int result;
	while((result = getopt(argc, argv, "i:I:h")) != -1)
	{
		switch(result)
		{
			case 'i': fp = gzopen(optarg,"r"); break;
			case 'I': fd = gzopen(optarg,"r"); break;
			case 'h':
			case '?':
				fprintf(stderr, "fqstat     summary statistics of PE fastq file\n\n");
				fprintf(stderr, "           -i  :<char> clean fastq1[.gz] file name\n");
				fprintf(stderr, "           -I  :<char> clean fastq2[.gz] file name\n");
				fprintf(stderr, "           -h  :<char> show this help\n\n");			
				break;
		}
	}
	if(argc<=2)
	{
		fprintf(stderr, "fqstat     summary statistics of vcf file\n\n");
		fprintf(stderr, "           -i  :<char> clean fastq1[.gz] file name\n");
		fprintf(stderr, "           -I  :<char> clean fastq2[.gz] file name\n");
		fprintf(stderr, "           -h  :<char> show this help\n\n");
		exit(1);
	}

	if(!fp || !fd)  exit(1);
	stat stat_1={0.,0.,0,0,0,0,0,0,0,0,0}; 
	stat stat_2={0.,0.,0,0,0,0,0,0,0,0,0}; 
	int flag=1;
	while(flag)
	{
		rinfo read1={{NULL,NULL,NULL,NULL},0};
		rinfo read2={{NULL,NULL,NULL,NULL},0};
		for(int i=0; i<4; i++)
		{
			char *r1=readline(fp);
			char *r2=readline(fd);
			if(r1 && r2)
			{
				size_t len1 = strlen(r1);
				size_t len2 = strlen(r2);
				read1.line[i]=(char *)malloc(len1+1); 
				read2.line[i]=(char *)malloc(len2+1);
				if(!read1.line[i] || !read2.line[i]) exit(1);
				strncpy(read1.line[i],r1,len1+1);
				strncpy(read2.line[i],r2,len2+1);
				if(i==1)
				{
					read1.len=len1; 
					read2.len=len2;
				}
			}
			else
			{
				flag=0;
			}
			free(r1);
			free(r2);
		}
	
		if(flag)
		{
			for(int i=0; i<4; i++)
			{
				if(i==1)
				{
					int k=0; stat_1.readCount++;
					while(read1.line[i][k] != '\0')
					{
						if(read1.line[i][k] == 'A') stat_1.aCount++;
						if(read1.line[i][k] == 'T') stat_1.tCount++;
						if(read1.line[i][k] == 'G') stat_1.gCount++;
						if(read1.line[i][k] == 'C') stat_1.cCount++;
						if(read1.line[i][k] == 'N') stat_1.nCount++;
						k++; 
						stat_1.baseCount++;
					}
					int s=0; stat_2.readCount++;
					while(read2.line[i][s] != '\0')
					{
						if(read2.line[i][s] == 'A') stat_2.aCount++;
						if(read2.line[i][s] == 'T') stat_2.tCount++;
						if(read2.line[i][s] == 'G') stat_2.gCount++;
						if(read2.line[i][s] == 'C') stat_2.cCount++;
						if(read2.line[i][s] == 'N') stat_2.nCount++;
						s++; 
						stat_2.baseCount++;
					}
				}
				if(i==3)
				{
					stat_1.q20 += qualityN(read1.line[i],Q20);
					stat_1.q30 += qualityN(read1.line[i],Q30);
					stat_2.q20 += qualityN(read2.line[i],Q20);
					stat_2.q30 += qualityN(read2.line[i],Q30);
				}
			}
			freeRead(read1); 
			freeRead(read2);
		}
	}
	gzclose(fp); 
	gzclose(fd);

	stat_1.averageLen = (float)stat_1.baseCount / stat_1.readCount;
	stat_2.averageLen = (float)stat_2.baseCount / stat_2.readCount;
	stat_1.gc = (float)(stat_1.gCount + stat_1.cCount) / stat_1.baseCount * 100;
	stat_2.gc = (float)(stat_2.gCount + stat_2.cCount) / stat_2.baseCount * 100;

	fprintf(stdout,"Iterm\treads_1.fq\treads_2.fq\n");
	fprintf(stdout,"read average length:\t%d\t%d\n",(int)stat_1.averageLen,(int)stat_2.averageLen);
	fprintf(stdout,"read GC content(%%):\t%.2f\t%.2f\n",stat_1.gc,stat_2.gc);
	fprintf(stdout,"total read Count:\t%lu\t%lu\n",stat_1.readCount,stat_2.readCount);
	fprintf(stdout,"total base Count:\t%lu\t%lu\n",stat_1.baseCount,stat_2.baseCount);
	
	fprintf(stdout,"\n");
	fprintf(stdout,"base A Count:\t%lu(%.2f%%)\t%lu(%.2f%%)\n", \
			stat_1.aCount,(float)stat_1.aCount/stat_1.baseCount*100.,stat_2.aCount,(float)stat_2.aCount/stat_2.baseCount*100.);

	fprintf(stdout,"base C Count:\t%lu(%.2f%%)\t%lu(%.2f%%)\n", \
			stat_1.cCount,(float)stat_1.cCount/stat_1.baseCount*100.,stat_2.cCount,(float)stat_2.cCount/stat_2.baseCount*100.);

	fprintf(stdout,"base G Count:\t%lu(%.2f%%)\t%lu(%.2f%%)\n", \
			stat_1.gCount,(float)stat_1.gCount/stat_1.baseCount*100.,stat_2.gCount,(float)stat_2.gCount/stat_2.baseCount*100.);

	fprintf(stdout,"base T Count:\t%lu(%.2f%%)\t%lu(%.2f%%)\n", \
			stat_1.tCount,(float)stat_1.tCount/stat_1.baseCount*100.,stat_2.tCount,(float)stat_2.tCount/stat_2.baseCount*100.); 

	fprintf(stdout,"base N Count:\t%lu(%.2f%%)\t%lu(%.2f%%)\n", \
			stat_1.nCount,(float)stat_1.nCount/stat_1.baseCount*100.,stat_2.nCount,(float)stat_2.nCount/stat_2.baseCount*100.);

	fprintf(stdout,"\n");
	fprintf(stdout,"Number of base calls with quality value of 20 or higher (Q20+) (%%)\t%lu(%.2f%%)\t%lu(%.2f%%)\n",stat_1.q20,(float)stat_1.q20/stat_1.baseCount*100.,stat_2.q20,(float)stat_2.q20/stat_2.baseCount*100.);
	fprintf(stdout,"Number of base calls with quality value of 30 or higher (Q30+) (%%)\t%lu(%.2f%%)\t%lu(%.2f%%)\n",stat_1.q30,(float)stat_1.q30/stat_1.baseCount*100.,stat_2.q30,(float)stat_2.q30/stat_2.baseCount*100.);

	return ;
}
