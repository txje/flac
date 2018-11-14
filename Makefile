CC=gcc
CFLAGS=-std=c99 -O2

OBJECTS = clust

all: $(OBJECTS)

clust: main.c hierarch.c paf.c
	${CC} ${CFLAGS} -o clust main.c hierarch.c paf.c -lz

agglomerate: agglomerate.c
	${CC} ${CFLAGS} -o agglomerate agglomerate.c -lm

clean:
	rm -f agglomerate clust
