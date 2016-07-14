# rtlsdr-wsprd -- WSPR daemon for RTL receivers

This non-interactive application allows automatic reporting of WSPR spots on WSPRnet. The idea is to allow the use of small computer like RasberryPi or Beaglebone boards, with a simple deamon. This kind of very lightweight setup could run continuously without maintenance and help to increase the WSPR network. The code is massively based on Steven Franke (K9AN) implementation and Joe Taylor (K1JT) work. This code was originally written for AirSpy receiver.

<h3>WARNING -- Crystal stability</h3>
Most of RTL dongles use a cheap crystal, and frequency drift can avoid WSPR decoding. The use of no-name RTL dongle for VHF/UHF bands usually require crystal modification, for a better one. External clock could be also used, like GPSDO or rubidium reference clock, aligned on 28.8MHz. 
In some case, it's possible to use the factory crystal (usually HC49, through hole), using a good thermal isolation. I successfully used two device with no modification, but it's tricky, easy to miss the window, and RTL devices does not allows fine frequency tuning.
For now, a good option is to buy a RTL device designed for SDR applications and integrating a 0.5ppm TCXO. After many tests, I would recommend this version :
NooElec NESDR SMArt - Premium RTL-SDR w/ Aluminum Enclosure, 0.5PPM TCXO
https://www.nooelec.com/store/nesdr-smart.html

<h3>Basically, this application :</h3>
- Perform a time alignment (2 mins)
- Start the reception using the RTL lib
- Decimate the IQ data (2.4Msps to 375 sps)
- Decode WSPR signal
- Push the spots on WSPRnet
- Repeat...

<h3>Howto :</h3>
1. Install a Linux compatible disto on your device (ex. Raspbian for RaspberryPi)
2. Install dependencies & useful tools (ex. ntp for time synchronization)
   ex: sudo apt-get install build-essential cmake libfftw3-dev libusb-1.0-0-dev curl libcurl4-gnutls-dev ntp 
3. Install rtlsdr library : https://github.com/steve-m/librtlsdr
4. Install rtlsdr-wsprd (this app) : https://github.com/Guenael/rtlsdr-wsprd
5. Enjoy it with ./rtlsdr_wsprd <your options>

<h3>Tips (for Raspberry Pi):</h3>
- Use ferrite bead to limit the interferences
- Cut off the display (could help to reduce QRN) : /opt/vc/bin/tvservice -o 
- Remove unused modules (ex: /etc/modules: #snd-bcm2835)

<h3>RTL devices & tests</h3>
- NooElec NESDR SMAr : Works fine out of the box
- RTL-SDR Blog 1PPM TCXO : Works with some drift, require additional mass, or additional enclosure
- Other no-name like : RT820, E4000, FC0012, FC0013, can work, but require modification and drift a lot

<h3>Raspberry devices & tests</h3>
- RaspberryPi 2, 15% of one core @1GHz (rx & decimation), and a burst during 10s on the second core), using Raspbian GNU Linux (Wheezy, 2015-02-16)
- RaspberryPi 1, 23% @700MHz (rx & decimation), and a burst during 30s on a second thread, using Raspbian GNU Linux (Wheezy, 2015-02-16)
- RaspberryPi 3 : TODO

<h3>Rasberian to use with a Raspberry PI</h3>
I noticed some disconnection problems with USB port while using Raspbian Jessie (2016-05-27-raspbian-jessie-lite.img). I rollbacked to Raspbian Wheezy (2015-02-16-raspbian-wheezy.img) and the problems solved by themself magically... For now, I have not investigated this issue, but if you experience some "Caught signal 11" error message, it could be this same problem. For now, I would recommend Raspbian Wheezy v2015-02-16.

