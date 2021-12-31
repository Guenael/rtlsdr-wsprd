# rtlsdr-wsprd -- WSPR daemon for RTL receivers

![rtlsdr-wsprd](art/rtlsdr-wsprd-web.jpg)

![Project Status](https://img.shields.io/badge/status-OK-green)
![Workflow Status](https://img.shields.io/github/workflow/status/Guenael/rtlsdr-wsprd/CI)
![Last commit](https://img.shields.io/github/last-commit/Guenael/rtlsdr-wsprd)
![Commit activity](https://img.shields.io/github/commit-activity/m/Guenael/rtlsdr-wsprd)
![LGTM Alerts](https://img.shields.io/lgtm/alerts/github/Guenael/rtlsdr-wsprd)
![LGTM Grade](https://img.shields.io/lgtm/grade/cpp/github/Guenael/rtlsdr-wsprd)

## TL;DR

This project aim at decoding [WSPR](https://en.wikipedia.org/wiki/WSPR_(amateur_radio_software)) signals using an [RTL device](https://osmocom.org/projects/rtl-sdr/wiki/Rtl-sdr), usually connected to a [Raspberry Pi](https://www.raspberrypi.org/).
To install and use your dongle on a Raspberry Pi with a Raspberry Pi OS, follow these steps:

```bash
echo "== Install dependencies"
sudo apt-get update && sudo apt-get -y install build-essential clang cmake libfftw3-dev libusb-1.0-0-dev libcurl4-gnutls-dev help2man ntp git

echo "== Install rtl-sdr library (on RPi, don't use your distro package)"
git clone https://github.com/osmocom/rtl-sdr
cd rtl-sdr
mkdir -p make
cd make
cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr -DDETACH_KERNEL_DRIVER=ON -Wno-dev ..
make
sudo make install
cd ../..

echo "== Install rtlsdr-wsprd"
git clone https://github.com/Guenael/rtlsdr-wsprd
cd rtlsdr-wsprd
make
sudo make install

echo "== Start/test rtlsdr-wsprd"
rtlsdr_wsprd -f 2m -c A1XYZ -l AB12cd -g 29
```

## Overview

This non-interactive application allows automatic reporting of WSPR spots on [WSPRnet](https://wsprnet.org). The initial idea was to allow a small computer like a Raspberry Pi and a RTL-SDR device to send WSPR reports for [VHF/UHF](https://en.wikipedia.org/wiki/Amateur_radio_frequency_allocations#Very_high_frequencies_and_ultra_high_frequencies) bands. This kind of lightweight setup could run continuously without maintenance and help to get additional propagation reports. The code is massively based on Steven Franke ([K9AN](https://github.com/k9an)) implementation of Joe Taylor ([K1JT](https://en.wikipedia.org/wiki/Joseph_Hooton_Taylor_Jr.)) publication and work.

This application written in C does:

- A time alignment (2 mins, required NTPd to run on the OS)
- Start the reception using the RTL lib
- Decimate the IQ data (2.4Msps to 375 sps)
- Decode WSPR signals
- Push any spots on WSPRnet
- Repeat, again and again...

## Installation

  1. Install a Linux compatible distro on your device.

     For Raspberry Pi, you can download official images [here](https://www.raspberrypi.com/software/operating-systems/).

  2. It's a good practice to update your OS. With Pi OS, run this command as usual:
     ```bash
     sudo apt-get update && sudo apt-get upgrade
     ```

  3. Install dependencies & useful tools (for example, [NTP](https://en.wikipedia.org/wiki/Network_Time_Protocol) for time synchronization). Example with a Debian based OS, like Rasbian, or Raspberry Pi OS:
     ```bash
     sudo apt-get update && sudo apt-get -y install build-essential clang cmake libfftw3-dev libusb-1.0-0-dev libcurl4-gnutls-dev help2man ntp git
     ```

  4. Install `rtl-sdr` library manually. **Do not use the `librtlsdr-dev` package on Raspberry PiOS**. There is a know bug with this lib and rtlsdr_wsprd will not be able to get enough samples (don't decode anything & 100% CPU pattern).
     ```bash
     git clone https://github.com/osmocom/rtl-sdr
     cd rtl-sdr
     mkdir -p make
     cd make
     cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr -DDETACH_KERNEL_DRIVER=ON -Wno-dev ..
     make
     sudo make install
     cd ../..
     ```
  Note: You may have to re-plug you dongle if it was already connected, or play with `udev` if not automatically detected.

  5. Clone this repository:
     ```bash
     git clone https://github.com/Guenael/rtlsdr-wsprd
     ```

  6. Build the application:
     ```bash
     cd rtlsdr-wsprd
     make
     sudo make install
     ```

  7. Finally, start the application with the right parameters/options for you (frequency, callsign, locator etc... Fake example below):
     ```bash
     rtlsdr_wsprd -f 2m -c A1XYZ -l AB12cd -g 29
     ```

## Container Image

As an alternative to the above steps, a pre-built container image containing rtlsdr-wsprd is available for use with [Docker](https://www.docker.com/) or [Podman](https://podman.io/).

The RTL DVB kernel modules must first be blacklisted on the host running the container. RTL-SDR itself is not required on the host running the container. This can be permanently accomplished using the following commands:

```bash
echo 'blacklist dvb_usb_rtl28xxu' | sudo tee /etc/modprobe.d/blacklist-dvb_usb_rtl28xxu.conf
sudo modprobe -r dvb_usb_rtl28xxu
```

If the `modprobe -r` command errors, a reboot is recommended to unload the module.

You can then start the container with the right parameters/options for you (frequency, callsign, locator etc... Fake example below):

```bash
docker run --rm -it --pull=always --device=/dev/bus/usb ghcr.io/guenael/rtlsdr-wsprd:latest -f 2m -c A1XYZ -l AB12cd -g 29
```

## Tips (for your Raspberry Pi and SDR dongles)

  - Use ferrite bead on the USB cable to limit the QRN
  - Use an external clean power supply
  - Cut off the display (could help to reduce QRN)
    ```bash
    /opt/vc/bin/tvservice -o
    ```
  - Remove unused modules (for example, /etc/modules: #snd-bcm2835)
  - Use an enclosure, and ground it

## Crystal stability

Most of RTL dongles use a cheap crystal, and frequency drift can effect the decoding & performance. The use of no-name RTL dongle for VHF/UHF bands usually require crystal modification, for a better one. External clock could be also used, like GPSDO or rubidium reference clock, aligned on 28.8MHz.

Some manufacturers integrate a 0.5ppm TCXO. It's the best second option, after an external clock. Based on my personal experience:

- NooElec NESDR SMART : Works fine out of the box
- RTL-SDR Blog 1PPM TCXO : Works with some drift, require additional mass, or a better enclosure
- Other no-name like : RT820, E4000, FC0012, FC0013, can work, but require modification and usually drift a lot

## Performance & hardware tests

Some performance tests using:
- Raspbian GNU/Linux 11 (bullseye) for Raspberry Pi devices
- rtlsdr-wsprd version 0.4.2
- Build with `clang -O3 -std=gnu17`

| Hardware      | Supported          | RX Load | Decode burst |
| ------------- | ------------------ | ------- | ------------ |
| RPi-1         | :heavy_check_mark: | 23.2%   | 8.4s         |
| RPi-2         | :heavy_check_mark: | 13.5%   | 4.1s         |
| RPi-3         | :heavy_check_mark: | 10.9%   | 2.1s         |
| RPi-4         | :heavy_check_mark: |  5.8%   | 1.1s         |
| PC (i7-5820K) | :heavy_check_mark: |  1.7%   | 0.5s         |
