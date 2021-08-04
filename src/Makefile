CC = gcc
OBJECTS = main.o sort.o print.o locationstatistics.o hashmap.o mergesort.o
CFLAGS = -g -c -Wall

sort : $(OBJECTS)
	$(CC) -o ../sort $(OBJECTS)

main.o : main.c main.h
	$(CC) $(CFLAGS) main.c
sort.o : sort.c sort.h
	$(CC) $(CFLAGS) sort.c
print.o : print.c print.h
	$(CC) $(CFLAGS) print.c
locationstatistics.o : locationstatistics.c locationstatistics.h
	$(CC) $(CFLAGS) locationstatistics.c
mergesort.o : mergesort.c mergesort.h
	$(CC) $(CFLAGS) mergesort.c
hashmap.o : hashmap.c hashmap.h
	$(CC) $(CFLAGS) hashmap.c

.PHONY : clean
clean :
	rm *.o
