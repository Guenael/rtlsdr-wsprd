---
name: Bug report
about: ''
title: ''
labels: 'bug'
assignees: Guenael
---

## Bug description
A clear and concise description of what the bug is.

## Setup & hardware used
- Version of `rtlsdr-wsprd` & source used (git pull or release archive for example)
- Version of `librtlsdr` intalled (with `ls -la /usr/lib/arm-linux-gnueabihf/librtlsdr.so.0`)
- Computer hardware (ex. RaspberryPi 3, Laptop with Intel processor, etc...)
- OS/Distro version (ex. `Raspbian GNU/Linux 10 (buster)`, or use `cat /etc/os-release`)
- Kernel used (with `uname -a`)
- Reference of your RTLsdr dongle in use (& info reported by `rtlsdr-wsprd` at start)

## Options used
Please provide the command line used with `rtlsdr-wsprd` and copy the first bloc that show the init. process.

## Samples (optionnal)
`rtlsdr-wsprd` allows to save samples with `-w <file prefix>` option, for debugging purpose. You can open these files with an editor like Audacity and use the spectrogram to inspect your signal. Zipping and adding these files could help to solve some issues.

## Additional context
Add any other context about the problem here.
