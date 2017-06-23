#CFLAGS= -Wall -O3 -ffast-math -std=gnu99

CC = gcc
CFLAGS= -Wall -O2 -std=gnu99
LDFLAGS = -L/usr/lib
LIBS = -lusb-1.0 -lrtlsdr -lpthread -lfftw3f -lcurl -lm

OBJS = rtlsdr_wsprd.o wsprd.o wsprsim_utils.o wsprd_utils.o tab.o fano.o nhash.o

%.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

rtlsdr_wsprd: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o rtlsdr_wsprd wspr_wisdom.dat hashtable.txt