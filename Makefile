CC = gcc
CFLAGS= -Wall -O3 -ffast-math -std=gnu17
LIBS = -lusb-1.0 -lrtlsdr -lpthread -lfftw3f -lcurl -lm

OBJS = rtlsdr_wsprd.o wsprd.o wsprsim_utils.o wsprd_utils.o tab.o fano.o nhash.o

TARGETS = rtlsdr_wsprd

.PHONY: all clean

all: $(TARGETS)

%.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

rtlsdr_wsprd: $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(TARGETS) wspr_wisdom.dat hashtable.txt