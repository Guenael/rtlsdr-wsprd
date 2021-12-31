# -------------------
# The build container
# -------------------
FROM debian:bullseye-slim AS build

RUN apt-get update && \
  apt-get -y --no-install-recommends install \
    build-essential \
    clang \
    cmake \
    libcurl4-openssl-dev \
    libfftw3-dev \
    libusb-1.0-0-dev \
    pkg-config \
    unzip && \
  rm -rf /var/lib/apt/lists/*

ADD https://github.com/steve-m/librtlsdr/archive/master.zip /root/librtlsdr-master.zip
RUN unzip /root/librtlsdr-master.zip -d /root && \
  rm /root/librtlsdr-master.zip && \
  cd /root/librtlsdr-master && \
  mkdir -p /root/librtlsdr-master/build && \
  cd /root/librtlsdr-master/build && \
  cmake -Wno-dev ../ && \
  make && \
  make install && \
  rm -rf /root/librtlsdr-master

COPY . /root/rtlsdr-wsprd

RUN cd /root/rtlsdr-wsprd && \
  make && \
  make install

# -------------------------
# The application container
# -------------------------
FROM debian:bullseye-slim

RUN apt-get update && \
  apt-get -y --no-install-recommends install \
   	libcurl4 \
    libfftw3-single3 \
    usbutils && \
  rm -rf /var/lib/apt/lists/*

COPY --from=build /usr/local/lib/librtlsdr.so.0 /usr/local/lib/librtlsdr.so.0
COPY --from=build /usr/local/bin/rtlsdr_wsprd /usr/local/bin/rtlsdr_wsprd
RUN ldconfig

ENTRYPOINT ["/usr/local/bin/rtlsdr_wsprd"]
