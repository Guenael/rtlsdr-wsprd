# rtlsdr-wsprd -- WSPR daemon for RTL receivers

This non-interactive application allows automatic reporting of WSPR spots on WSPRnet. The idea is to allow the use of small computer like RasberryPi or Beaglebone boards, with a simple deamon. This kind of very lightweight setup could run continuously without maintenance and help to increase the WSPR network. The code is massively based on Steven Franke (K9AN) implementation and Joe Taylor (K1JT) work. This code was originally written for AirSpy receiver.

<h3>WARNING -- Crystal stability</h3>
Most of RTL dongles use a cheap crystal, and frequency drift can avoid WSPR decoding. The use of RTL dongle for VHF/UHF band usually require a crystal modification, for a better one. External clock could be also used, like GPSDO or rubidium reference clock aligned on 28.8MHz. In some case, it's possible to use the factory crystal, using a good thermal isolation. In my opinion, AirSpy receiver is a best and easier option. In addition, it includes an LNA, a metal enclosure and gives better results.

<h3>Basically, this application :</h3>
- Perform a time alignment (2 mins)
- Start the reception using the RTL lib
- Decimate the IQ data (2.4Msps to 375 sps)
- Decode WSPR signal
- Push the spots on WSPRnet
- Repeat...

*Tested with Raspbian GNU Linux, using a RaspberryPi 2*
(20% of one core @1GHz (rx & decimation), and a burst during 10s on the second core)

<h3>Howto :</h3>
1. Install a Linux compatible disto on your device (ex. Raspbian for RaspberryPi)
2. Install dependencies & useful tools (ex. ntp for time synchronization)
   ex: sudo apt-get install build-essential cmake libfftw3-dev libusb-1.0-0-dev curl libcurl4-gnutls-dev ntp 
3. Install rtlsdr library : https://github.com/steve-m/librtlsdr
4. Install rtlsdr-wsprd (this app) : http://github.com/Guenael/rtlsdr-wsrd
5. Enjoy it with ./rtlsdr_wsprd <your options>

<h3>Tips (for Raspberry Pi):</h3>
- Use ferrite bead to limit the interferences
- Cut off the display (could help to reduce QRN) : /opt/vc/bin/tvservice -o 
- Remove unused modules (ex: /etc/modules: #snd-bcm2835)
