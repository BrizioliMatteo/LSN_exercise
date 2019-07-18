CC = g++
CFLAGS = -Wall -O3 --std=c++11

Hydrogen_atom.exe : Hydrogen_atom.o random.o
	$(CC) random.o Hydrogen_atom.o -o Hydrogen_atom.exe
Hydrogen_atom.o : Hydrogen_atom.cpp Hydrogen_atom.h
	$(CC) -c Hydrogen_atom.cpp -o Hydrogen_atom.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Hydrogen_atom.exe
