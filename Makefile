CC=gcc
CFLAGS= -DBSD=0 -g -ansi -lm -O -Wall -I..

PROGS=test3

all: $(PROGS)

.c:
	$(CC) $(CFLAGS) -o $* $*.c

$(PROGS): 

clean:
	rm -rf $(PROGS) *.o




