CC=gcc
OBJS=main.o skr.o 
CFLAGS+=-c -Wall -g -lz -std=c99

skr:$(OBJS)
	$(CC) $^ -lz -O3 -o $@

%.o:%.c
	$(CC) $^ $(CFLAGS) -o $@

clean:
	$(RM) *.o skr -r
