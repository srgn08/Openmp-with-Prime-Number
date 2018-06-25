#
#
#

CC = gcc
CFLAGS = -fopenmp -O
DEBUGFLAGS = -fopenmp -g

all: asal asal_debug
debug: asal_debug

.PHONY: clean

clean:
	rm -rf asal asal_debug *~ *.o

asal: asal.c
	$(CC) $(CFLAGS) $< -o $@

asal_debug: asal.c
	$(CC) $(DEBUGFLAGS) $< -o $@

