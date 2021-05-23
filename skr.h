#ifndef _H_SEQ
#define _H_SEQ

#include <zlib.h>

char *readline(gzFile file);

void mergeVcfs(int argc, char *argv[]);

void fq2fa(int argc, char *argv[]);

#endif
