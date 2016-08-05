CC=xlc++_r

exec: main.o
	$(CC) -O3 -qsmp=omp -o exec.x main.o

main.o: main.c
	$(CC) -O3 -qsmp=omp -c main.c

clean:
	rm *.o *.x
