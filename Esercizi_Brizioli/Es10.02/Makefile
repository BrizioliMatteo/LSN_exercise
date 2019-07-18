CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

Gas.exe : Gas.o random.o Maps.o Chrom.o
	$(CC) random.o Maps.o Chrom.o Gas.o -o Gas.exe

Gas.o : Gas.cpp Gas.h
	$(CC) -c Gas.cpp -o Gas.o $(CFLAGS)

Chrom.o : Chrom.cpp Chrom.h 
	$(CC) -c Chrom.cpp -o Chrom.o $(CFLAGS)

Maps.o : Maps.cpp Maps.h 
	$(CC) -c Maps.cpp -o Maps.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Gas.exe
