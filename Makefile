OBJECTS = clust

all: $(OBJECTS)

clust: main.c hierarch.c
	gcc -std=c99 -o clust main.c hierarch.c

agglomerate: agglomerate.c
	gcc -std=c99 -o agglomerate agglomerate.c -lm

clean:
	rm -f agglomerate clust
