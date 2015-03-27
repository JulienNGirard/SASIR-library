#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "GlobalInc.h"
#include "FloatTrans.h"
#include "writefits3d.h"
#include "CEA_comp_sens.h"

static float_array data, mask;		// Treated as 2-D arrays

void read_fltarr(char *filename, float_array& fa) {
  fltarray data;

  fits_read_fltarr(filename, data);
  if ( data.nx() != data.ny() ) {
    fprintf(stderr,"FITS images must be square\n");
    exit(1);
  }
  fa.alloc(data.nx());
  for (int i=0; i<data.nx(); ++i) for (int j=0; j<data.ny(); ++j) fa[i][j] = data(i,j);
}

void write_fltarr(char *filename, float_array& fa) {
  fltarray data(fa.size(),fa.size());
   for (int i=0; i<data.nx(); ++i) for (int j=0; j<data.ny(); ++j) data(i,j) = fa[i][j];
  writefltarr(filename,data);
}

int main(int argc, char *argv[]) {
  const int FISTA_ALGORITHM = 2;
  const int MCA_ALGORITHM = 3;
  float_array dirty_image, next_image, difference;
  bool positivity = true;
  int algorithm = FISTA_ALGORITHM;

  if ( argc != 4 ) {
    fprintf(stderr,"Usage: %s <input FITS image> <input FITS FFT mask> <output FITS image>\n",argv[0]);
    fprintf(stderr,"\nThe program runs CS using the library that is also used by the LOFAR imagers.\n");
    fprintf(stderr,"Its purpose is to test the CS algorithms used by the imagers, but in a simpler way,\n");
    fprintf(stderr,"using a test image and a random mask like mr_uv_inpainting.\n");
    fprintf(stderr,"Include .fits extension on file names. CS parameters are hardwired,\n");
    fprintf(stderr,"change them in the code in main.cc\n");
    return 1;
  }
  read_fltarr(argv[1],data);
  read_fltarr(argv[2],mask);
  
 
  dirty_image.alloc(data.size());  next_image.alloc(data.size()); difference.alloc(data.size());

  // Create the dirty image, used to prime the compressed sensing. 
  fourier_transform_and_mask(data, mask, dirty_image);      

  // The first wavelet reconstruction is generated from dirty_image and placed in next_image. 
  compressed_sensing_begin(dirty_image, next_image, true, 200, true, 1, 0.01, 1e-6, 7, algorithm, 3, false, 3, false);

  do
    // Take the original image and the next wavelet image and calculate the difference between their 
    // masked Fourier transforms. This implements |y-Hx|
   
    fft_difference(data, next_image, difference, mask, positivity);
  
    // Pass the difference back to the compressed sensing, to get the next wavelet image, or stop
  while ( compressed_sensing_next(difference, next_image, algorithm) );

  write_fltarr(argv[3],next_image);
  return 0;
}


