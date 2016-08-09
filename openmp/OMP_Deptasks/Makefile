CC=icpc

exec: deptask.o
	$(CC) -O3 -qopenmp -o exec.x deptask.o

deptask.o: deptask.c
	$(CC) -O3 -qopenmp -c deptask.c

clean:
	rm *.o *.x
