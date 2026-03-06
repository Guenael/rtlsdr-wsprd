# ----------------------
# --- STAGE 1: Build ---
# ----------------------

# Vulnerability check: https://hub.docker.com/layers/library/debian/trixie-slim
FROM debian:trixie-slim AS build

# Install build dependencies
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -qq update \
 && apt-get -q -y install \
    build-essential \
    git \
    clang \
    cmake \
    libcurl4-openssl-dev \
    libfftw3-dev \
    libusb-1.0-0-dev \
    pkg-config \
 && apt-get clean all \
 && rm -rf /var/lib/apt/lists/*

# Copy rtlsdr-wsprd source code
COPY . /rtlsdr-wsprd

# Build rtl-sdr and rtlsdr-wsprd
RUN git clone https://github.com/osmocom/rtl-sdr \
 && cd rtl-sdr \
 && mkdir -p make \
 && cd make \
 && cmake -Wno-dev .. \
 && make \
 && make install \
 && ldconfig \
 && cd /rtlsdr-wsprd \
 && make \
 && make test \
 && make install


# ----------------------------
# --- STAGE 2: Application ---
# ----------------------------

FROM debian:trixie-slim

ARG HOST_UID=1000
ARG HOST_GID=1000
ARG LOCAL_USER=app
ARG LOCAL_GROUP=users

# Only copy the necessary files from the build stage
COPY --from=build /usr/local/lib/librtlsdr.so.0 /usr/local/lib/librtlsdr.so.0
COPY --from=build /usr/local/bin/rtlsdr_wsprd /usr/local/bin/rtlsdr_wsprd

RUN apt-get -qq update \
 && apt-get -q -y install \
    ca-certificates \
    libcurl4t64 \
    libfftw3-single3 \
    usbutils \
  && apt-get clean all \
  && rm -rf /var/lib/apt/lists/* \
  && ldconfig

RUN groupmod -g $HOST_GID $LOCAL_GROUP \
 && useradd -u $HOST_UID -g $LOCAL_GROUP -m $LOCAL_USER

WORKDIR /home/$LOCAL_USER

# Drop privileges to the non-root user
USER $LOCAL_USER

ENTRYPOINT ["/usr/local/bin/rtlsdr_wsprd"]
