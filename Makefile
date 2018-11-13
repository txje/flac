OBJECTS = clust

all: $(OBJECTS)

clust: main.c hierarch.c paf.c
	gcc -std=c99 -o clust main.c hierarch.c paf.c

agglomerate: agglomerate.c
	gcc -std=c99 -o agglomerate agglomerate.c -lm

clean:
	rm -f agglomerate clust
