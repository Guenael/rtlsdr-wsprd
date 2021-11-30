# rtlsdr-wsprd -- WSPR daemon for RTL receivers

This non-interactive application allows automatic reporting of WSPR spots on WSPRnet. The idea is to allow the use of small computer like RaspberryPi or Beaglebone boards, with a simple daemon. This kind of very lightweight setup could run continuously without maintenance and help to increase the WSPR network. The code is massively based on Steven Franke (K9AN) implementation and Joe Taylor (K1JT) work. This code was originally written for AirSpy receiver.

## WARNING -- Crystal stability

Most of RTL dongles use a cheap crystal, and frequency drift can effect WSPR decoding. The use of no-name RTL dongle for VHF/UHF bands usually require crystal modification, for a better one. External clock could be also used, like GPSDO or rubidium reference clock, aligned on 28.8MHz. 
In some case, it's possible to use the factory crystal (usually HC49, through hole), using a good thermal isolation. I successfully used two devices with no modification, but it's tricky, easy to miss the window, and RTL devices do not allow for fine frequency tuning.

For now, a good option is to buy an RTL device designed for SDR applications and integrating a 0.5ppm TCXO. After many tests, I would recommend this version :

NooElec NESDR SMArt - Premium RTL-SDR w/ Aluminum Enclosure, 0.5PPM TCXO
https://www.nooelec.com/store/nesdr-smart.html

## Basically, this application

- Perform a time alignment (2 mins)
- Start the reception using the RTL lib
- Decimate the IQ data (2.4Msps to 375 sps)
- Decode WSPR signal
- Push the spots on WSPRnet
- Repeat...

## Installation
  1. Install a Linux compatible distro on your device (ex. Raspbian for RaspberryPi)
  1. Install dependencies & useful tools (for example, NTP for time synchronization)
     ```
     sudo apt-get install build-essential cmake libfftw3-dev libusb-1.0-0-dev curl libcurl4-gnutls-dev ntp
     ```
  1. Install rtlsdr library: https://github.com/steve-m/librtlsdr
  1. Install rtlsdr-wsprd (this app): https://github.com/Guenael/rtlsdr-wsprd
  1. Enjoy it with ```./rtlsdr_wsprd <your options>```

## Tips (for Raspberry Pi)
  - Use ferrite bead to limit the interferences
  - Cut off the display (could help to reduce QRN)
    ```
    /opt/vc/bin/tvservice -o
    ``` 
  - Remove unused modules (for example, /etc/modules: #snd-bcm2835)

## RTL devices & tests
  - NooElec NESDR SMAr : Works fine out of the box
  - RTL-SDR Blog 1PPM TCXO : Works with some drift, require additional mass, or additional enclosure
  - Other no-name like : RT820, E4000, FC0012, FC0013, can work, but require modification and drift a lot

## Raspberry devices & tests
  - RaspberryPi 2, 15% of one core @1GHz (rx & decimation), and a burst during 10s on the second core), using Raspbian GNU Linux (Wheezy, 2015-02-16)
  - RaspberryPi 1, 23% @700MHz (rx & decimation), and a burst during 30s on a second thread, using Raspbian GNU Linux (Wheezy, 2015-02-16)
  - RaspberryPi 3 : TODO

## Rasberian to use with a Raspberry PI
I noticed some disconnection problems with USB port while using Raspbian Jessie (2016-05-27-raspbian-jessie-lite.img). I rolled back to Raspbian Wheezy (2015-02-16-raspbian-wheezy.img) and the problems solved by themselves magically... For now, I have not investigated this issue, but if you experience some "Caught signal 11" error message, it could be this same problem. For now, I would recommend Raspbian Wheezy v2015-02-16.
