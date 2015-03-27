
#ifndef _MCA_LOFAR_H
#define _MCA_LOFAR_H

#include <fstream>
#include "GlobalInc.h"
#include "FloatTrans.h"
#include "writefits3d.h"

extern Bool Verbose;
extern float correlation(fltarray&, fltarray&);

class MCA_params
{
public:
	void reset();			// Initialise the class with default parameters
	char NameOut[256];
	int MaxNiter;			// Maximum number of iterations
	float TolVar;			// Stopping criterium : ||xt+1 - xt||_infty < TolVar
	float Mu;				// Gradient regularization parameter
	bool Positivity;		// The reconstruction must be >0
	bool No_coarse;			// if true : kills the coarse scale before apaplying fista, then thresholds it
	float ksigma;           // threshold level = k * SigmaNoise.
	double SigmaNoise;      // Noise standard deviation:  lambda = ksigma * SigmaNoise
};

void MCA_params::reset()
{
	MaxNiter=10;
	TolVar = 1e-4;
	Mu=1;
	Positivity = false;
	No_coarse = false;
	ksigma=3.;
	SigmaNoise=0.;
}


class MCA_LOFAR : public MCA_params
{
public:
	MCA_LOFAR(FloatTrans *domain);
	~MCA_LOFAR(){};
	
	FloatTrans* Domain;
	
	void begin(fltarray &image0);
	// b : observed signal, in degraded domain
	// z : recovered signal, in degraded domain
	// _degrade : degradation operator. 1bool to apply the degradation, 2bool to transpose
    
    bool next(fltarray &difference, fltarray& image);
	
 
    float lvl;
    bool done, allocD;
    float speed, old_speed, top_level;
    int i;
    float *wavelet_coeff;
    float *CS;
    fltarray z, y;

};

MCA_LOFAR::MCA_LOFAR(FloatTrans *domain)
{
        reset();
        Domain = domain;
}

void MCA_LOFAR::begin(fltarray &image0)
{
  float array_max;
	
	
    i = 0; done = false; speed = 0; old_speed = 0; allocD = true;

  
	//cout << "FISTA params " <<  P.MaxNiter << " " << P.Threshold << " " << P.TolVar << " " <<  P.Mu << " " << P.Fast << " " << P.Decreasing << " " << P.Positivity << " " << P.No_coarse << endl;
	
// Initial point
  

	if(No_coarse)   // turn the FFT back into an image, then wavelet space, then fiddle with wavelets, then back to Fourier space
	{
		Domain->transform(image0, wavelet_coeff, allocD); allocD=false; // allocate wavelet_coeff and localTB
		Domain->substract_coarse_scale(wavelet_coeff,CS,true);
		Domain->recons(wavelet_coeff, image0);// allocate wavelet_coeff and localTB
	}

	cerr<<"##########\nBegin ";
	cerr<<"MCA"; if ( Domain->UseMad ) cerr<<" with MAD, ksigma=" << ksigma << ",";
	cerr<<" with mu="<<Mu<<", and "<<MaxNiter<<" iterations.\n##########"<<endl;

       
	// z holds the current image
	z = image0; 
	
        Domain->transform(z,wavelet_coeff,allocD);   // just to allocate
 
        // Find wavelet max for top threshold
        array_max = wavelet_coeff[0];
        for (int i=0; i<Domain->size_transform()-z.nx()*z.ny(); ++i)  if ( wavelet_coeff[i] > array_max )  array_max = wavelet_coeff[i]; 
        top_level = 0.99*array_max;
        cout << "Calculated top level: " << array_max  << endl << endl;   
}



bool MCA_LOFAR::next(fltarray &difference, fltarray& image) {

// MCA algorithm.
	if( i<MaxNiter && !done )
	{
	  // Threshold update
	  lvl = 1-i*1.0/float(MaxNiter-1);			// 1..0

	  if ( access("original.fits",F_OK) != -1 ) {
              fltarray orig; fits_read_fltarr("original.fits",orig);
              cout << "Corr=" << correlation(orig,z) << endl;
          }
		
	// Evaluate the evolution. Compare z to y, which holds the previous image
	  y = y - z; // error between the last two estimates
	  speed = abs(y.maxfabs()/z.maxfabs());
	  done = (speed < TolVar) && (i>0) && (speed<old_speed);
	  old_speed=speed;
	  y = z; // Save the current solution image zreal to yreal
	
	// Gradient step
	  if ( Domain->UseMad ) 
	  {
	    for(int k=0;k<z.n_elem();k++) difference(k) *= Mu;
	    Domain->transform(difference, wavelet_coeff);
	    Domain->mad_calculation(wavelet_coeff);
	    for(int k=0;k<z.n_elem();k++) difference(k) /= Mu;
	  }
				
		 
	  for(int k=0;k<z.n_elem();k++) z(k) = z(k) + Mu*difference(k);

	
	// Proximal operator : hard thresholding using wavelets
	  Domain->transform(z, wavelet_coeff); 
	  if ( Domain->UseMad )
	  {
           Domain->mad_normalize(wavelet_coeff);
           Domain->hard_threshold(wavelet_coeff,ksigma+lvl,No_coarse);
           Domain->mad_unnormalize(wavelet_coeff);
           Domain->adjoint_recons(wavelet_coeff, z);
	   if(Verbose) cerr<<" Step "<<i<<", lvl=MD use "<<ksigma+lvl<<", || (z - zt)/z ||_infty = "<<abs(speed)<<", TolVar="<<TolVar<<endl; 
	  }
	  else
	  {
	   lvl = (pow((double)lvl,3)+lvl/25.)/1.04;
	   lvl = SigmaNoise+lvl*(top_level-SigmaNoise);    
           Domain->normalize(wavelet_coeff);
           Domain->hard_threshold(wavelet_coeff,lvl, No_coarse);
           Domain->unnormalize(wavelet_coeff);
           Domain->recons(wavelet_coeff, z);
	   if(Verbose) cerr<<" Step "<<i<<", lvl="<<lvl<<", || (z - zt)/z ||_infty = "<<abs(speed)<<", TolVar="<<TolVar<<endl; 
	  }
		
	  if(Positivity) for(int k=0;k<z.n_elem();k++) z(k) = z(k) *(z(k)>0);	
	  image = z;
        
	  i++;
        
	  return true;
	} else {
	  if(No_coarse) {
	    Domain->transform(z,wavelet_coeff);
	    Domain->normalize(wavelet_coeff);
	    Domain->unnormalize(wavelet_coeff);
  	    Domain->add_coarse_scale(wavelet_coeff,CS);
	    Domain->recons(wavelet_coeff,z);
	  }
	
	  
	  image = z;
	  cerr<<"##########\nEnd CS.\n##########"<<endl;
	  return false;
	}
};
#endif


