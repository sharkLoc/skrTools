#ifndef _H_SEQ_
#define _H_SEQ_

#include <zlib.h>
#include <stdint.h>

#define Q20 20
#define Q30 30

typedef struct Stat {
	float averageLen,gc;
	uint64_t readCount,baseCount,nCount,aCount,tCount,gCount,cCount,q20,q30;
} stat;

typedef struct readInfo {
	char *line[4];
	size_t len;
} rinfo;


char *readline(gzFile file);

void mergeVcfs(int argc, char *argv[]);

void fq2fa(int argc, char *argv[]);

uint64_t qualityN(char *line,size_t q);

void freeRead(rinfo read);

void fqstat(int argc, char *argv[]);

void statVcfs(int argc, char *argv[]);

void makewind(int argc, char *argv[]);

#endif
