# SASIR-library

Requirements:
- Cmake
- c/c++ compiler
- cfitsio

To install:

$ mkdir build
$ cd build
$ cmake ..
$ make
$ make install

Then to finalize the installation, follow instruction or run the install.sh script if you have root permission on your system.

This will create a library file:   libCEA_comp_sens.a


For building the altered version of LWimager with SASIR compressed sensing imaging code, you need to locate two files:

- libCEA_comp_sens.a (we just compiled and packaged it)
- CEA_comp_sens.h

The are copied in the SASIR-library/library folder.