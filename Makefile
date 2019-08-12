CC=gcc
CFLAGS=-std=c99 -O2

OBJECTS = flac

all: $(OBJECTS)

flac: main.c hierarch.c paf.c
	${CC} ${CFLAGS} -o flac main.c hierarch.c paf.c -lz

clean:
	rm -f flac
