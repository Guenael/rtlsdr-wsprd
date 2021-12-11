CC = clang
CFLAGS= -O3 -std=gnu17 -Wall # -fsanitize=address
LIBS = -lusb-1.0 -lrtlsdr -lpthread -lfftw3f -lcurl -lm # -fsanitize=address

OBJS = rtlsdr_wsprd.o wsprd/wsprd.o wsprd/wsprsim_utils.o wsprd/wsprd_utils.o wsprd/tab.o wsprd/fano.o wsprd/nhash.o

TARGETS = rtlsdr_wsprd

.PHONY: all clean

all: $(TARGETS)

%.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

rtlsdr_wsprd: $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o wsprd/*.o $(TARGETS) fftw_wisdom.dat hashtable.txt selftest.iq

install:
	install rtlsdr_wsprd /usr/local/lib/rtlsdr_wsprd
