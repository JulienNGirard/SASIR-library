#include <complex>
using namespace std;
// Simple classes to allow shipping of data from one side to the other (Imager to CS). 2-D arrays implemented in 1-D. All inline.


class complex_array {
private:
  complex<float> *data;
  int dimension;
public:
  complex_array() { dimension = 0; data = NULL; }
  complex_array(int dim) { alloc(dim);  }
  ~complex_array() { delete[] data; }
  void alloc(int dim) { dimension = dim; data = new complex<float>[dim*dim]; }
  complex<float> *operator[](int i) { return data+i*dimension; }
  int size() { return dimension; }
};

class float_array {
private:
  float *data;
  int dimension;
public:
  float_array() { dimension = 0; data = NULL; }
  float_array(int dim) { alloc(dim); }
  ~float_array() { delete[] data; }
  void alloc(int dim) { dimension = dim; data = new float[dim*dim]; }
  float *operator[](int i) { return data+i*dimension; }
  int empty() { return (dimension == 0); }
  int size() { return dimension; }
};


// The interface to the Compressed Sensing library -----------------------------------------------


extern void fourier_transform_and_mask(float_array& image, float_array& mask, float_array& dirty_image);

extern void fft_difference(float_array& original, float_array& current, float_array& difference, float_array& mask, bool pos);

    // Iterator routines for AWImager -----------------------------------------------------------------------------------------------
extern void compressed_sensing_begin(float_array& dirty_image, float_array& image, bool verbose, 
		int Niter, bool Positivity, float mu, float noise, float TolVar, int NbrScale, 
		int Type_algo, int Type_domain, bool no_coarse, double ksigma);


extern bool compressed_sensing_next(float_array& diff, float_array& image, int algorithm);

extern void compressed_sensing_load(float s);

extern void compressed_sensing_zero_fine();









