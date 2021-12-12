FROM gitpod/workspace-full

RUN sudo apt-get update \
 && sudo apt-get install -y \
         valgrind libcmocka-dev cmocka-doc libcmocka0 \
         build-essential clang cmake \
         libfftw3-dev libusb-1.0-0-dev librtlsdr-dev libcurl4-gnutls-dev \
 && sudo rm -rf /var/lib/apt/lists/*
