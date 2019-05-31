CC = gcc
CFLAGS = -O2 -Wall

MT = mt19937ar-cok

all: main xcorr

xcorr: xcorr.o
	$(CC) $(CFLAGS) -o $@ $@.o -lm -lgd
main: main.o $(MT).o
	$(CC) $(CFLAGS) -o $@ $@.o $(MT).o -lm
clean:
	rm -f *.o main xcorr
