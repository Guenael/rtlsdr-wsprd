CC = clang
CFLAGS = -O3 -std=gnu17
LIBS = -lusb-1.0 -lrtlsdr -lpthread -lfftw3f -lcurl -lm

# Note
#   gcc is a bit faster that clang on this app
#   for dbg: -Wall -fsanitize=address

ifeq ($(findstring armv6, $(shell uname -m)), armv6)
	# Broadcom BCM2835 SoC with 700 MHz 32-bit ARM 1176JZF-S (ARMv6 arch)
	#   Used in Raspberry Pi1 (A,A+,B,B+), Pi-Zero, Pi-Zero W, Pi-Compute-Module1
	EXTRA_OPTS = -DRPI1 --target=arm-linux-gnueabihf -mcpu=arm1176jzf-s -mfloat-abi=hard
endif
ifeq ($(findstring armv7, $(shell uname -m)), armv7)
	# Broadcom BCM2836 SoC with 900 MHz 32-bit quad-core ARM Cortex-A7 (ARMv7 arch)
    #   Used in Raspberry Pi 2 Model B
	# Broadcom BCM2837 SoC with 1.2 GHz 64-bit quad-core ARM Cortex-A53 (ARMv8 arch)
    #   Used in Raspberry Pi 3 Model B +later models of the Raspberry Pi 2 Model B, 
    #   and Pi-Compute-Module3.
	EXTRA_OPTS = -DRPI23
endif

OBJS = rtlsdr_wsprd.o wsprd/wsprd.o wsprd/wsprsim_utils.o wsprd/wsprd_utils.o wsprd/tab.o wsprd/fano.o wsprd/nhash.o

TARGETS = rtlsdr_wsprd

.PHONY: all clean

all: $(TARGETS)

%.o: %.c
	${CC} ${CFLAGS} $(EXTRA_OPTS) -c $< -o $@

rtlsdr_wsprd: $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o wsprd/*.o $(TARGETS) fftw_wisdom.dat hashtable.txt selftest.iq

install:
	install rtlsdr_wsprd /usr/local/bin/rtlsdr_wsprd
