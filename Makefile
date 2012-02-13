genetrack: genetrack.c gopt.o
	gcc gopt.o genetrack.c -o genetrack

gopt.o: gopt.c gopt.h
	gcc -c gopt.c

