
#ifndef _FISTA_LOFAR_H
#define _FISTA_LOFAR_H

#include <fstream>
#include "GlobalInc.h"
#include "FloatTrans.h"
#include "writefits3d.h"

extern Bool Verbose;
extern float correlation(fltarray&, fltarray&);

class Fista_params
{
public:
        void reset();                   // Initialise the class with default parameters
        char NameOut[256];
        int MaxNiter;                   // Maximum number of iterations
        float TolVar;                   // Stopping criterium : ||xt+1 - xt||_infty < TolVar
        float Mu;                               // Gradient regularization parameter
        bool Fast;                              // (F)ISTA
        bool Decreasing;                // if true : linearily decreasing threshold
        bool Positivity;                // The reconstruction must be >0
        bool No_coarse;                 // if true : kills the coarse scale before apaplying fista, then thresholds it
    float ksigma;           // threshold level = k * SigmaNoise.
    double SigmaNoise;      // Noise standard deviation:  lambda = ksigma * SigmaNoise
};
void Fista_params::reset()
{
        MaxNiter=10;
        TolVar = 1e-4;
        Mu=1;
        Fast = false;
        Decreasing = false;
        Positivity = false;
        No_coarse = false;
    ksigma=3.;
    SigmaNoise=0.;
}

class Fista_LOFAR : public Fista_params
{
public:
	Fista_LOFAR(FloatTrans *domain);
	~Fista_LOFAR(){};
	
	FloatTrans* Domain;
	
	void begin(fltarray &image0);
	// b : observed signal, in degraded domain
	// z : recovered signal, in degraded domain
	// _degrade : degradation operator. 1bool to apply the degradation, 2bool to transpose
	bool next(fltarray &difference, fltarray& image);
	// b : observed signal, in degraded domain
	// z : recovered signal, in degraded domain
	// _degrade : degradation operator. 1bool to apply the degradation, 2bool to transpose
	// _degrade(z,*,true) is assumed real (imaginary part ignored)
	
	void load(float s);	
	void zero_fine();
	fltarray z;
	float *x, *xold, *xtmp; // current estimate, previous estimate, and temp variable

// Other variables
	float tk, told;
	int i, n; bool done;
	float speed; float old_speed;
	fltarray y;   // previous reconstruction
	bool allocD; // Domain allocation
	float* CS;
};


Fista_LOFAR::Fista_LOFAR(FloatTrans *domain)
{
	reset(); 
	Domain = domain;
}


void Fista_LOFAR::begin(fltarray &image0)
{
    TolVar = TolVar;  tk = 1;
    i = 0; done = false; speed = 0; old_speed = 0; allocD = true;
	
//cerr<<"z nx ny nl nc "<<z.nx()<<","<<z.ny()<<","<<z.nl()<<","<<z.nc()<<endl;
//cerr<<"b nx ny nl nc "<<b.nx()<<","<<b.ny()<<","<<b.nl()<<","<<b.nc()<<endl;
// Initial point z=b

	if(No_coarse) // Substract the coarse scale
	{
		Domain->transform(image0, xtmp, allocD); allocD=false;// allocate xtmp and localTB
		Domain->substract_coarse_scale(xtmp,CS,true);
		Domain->recons(xtmp, image0);// allocate xtmp and localTB
	}

	
	z = image0;
	
	Domain->transform(z, xtmp, allocD);// allocate xtmp and localTB
	n=Domain->size_transform();

	x = new float[n]; for(int k=0;k<n;k++) x[k] = xtmp[k];
	xold = new float[n];

	cerr<<"##########\nBegin ";
	if(Fast) cerr<<"FISTA";
	else cerr<<"ISTA";
	if ( Domain->UseMad ) cerr << " with MAD, ksigma=" << ksigma << ",";
	cerr<<" with mu="<<Mu<<", and "<<MaxNiter<<" iterations.\n##########"<<endl;
}



bool Fista_LOFAR::next(fltarray &difference, fltarray& image) {	
	
// Fista algorithm. 
	if ( i<MaxNiter && !done  )
	{
//cerr<<"z i nx ny nl nc "<<i<<":"<<z.nx()<<","<<z.ny()<<","<<z.nl()<<","<<z.nc()<<endl;
	// Threshold update
		//if ( access("original.fits",F_OK) != -1 ) {
                //  fltarray orig; fits_read_fltarr("original.fits",orig);
                //  cout << "Corr=" << correlation(orig,z) << endl;
                //}
		
	// Evaluate the evolution
		y = y-z; // error between the last two estimates
		if(  z.maxfabs() == 0 ) speed = 0;
		else speed = abs(y.maxfabs()/z.maxfabs()); 
		done = (speed < TolVar) && (i>0) && (speed<old_speed);
		old_speed=speed;
		y=z; // Save the new solution
		
	// Gradient step
		Domain->transform(difference, xold);
		if ( Domain->UseMad ) {
		  for(int k=0;k<n;k++) xold[k] *= Mu;
		  Domain->mad_calculation(xold);
		  for(int k=0;k<n;k++) xold[k] /= Mu;
		}
		  Domain->transform(difference, xold);
		for(int k=0;k<n;k++) xtmp[k] = xtmp[k] + Mu*xold[k];
		
	// Save estimate
		if(Fast) for(int k=0;k<n;k++) xold[k] = x[k];

	// Proximal operator : soft thresholding
		for(int k=0;k<n;k++) x[k] = xtmp[k];
		if ( Domain->UseMad )
		{
		  Domain->mad_normalize(x);
		  Domain->soft_threshold(x,ksigma,No_coarse);    // threshold(SigmaNoise) = 0
		  Domain->mad_unnormalize(x);
		  if(Verbose) cerr<<" Step "<<i<<", lvl=MAD using "<<ksigma<<", || (z - zt)/z ||_infty = "<<speed<<", TolVar="<<TolVar<<endl;
		}
		else
		{
		  Domain->normalize(x);
		  Domain->soft_threshold(x,SigmaNoise,No_coarse);		// ksigma=0
		  Domain->unnormalize(x);
		  if(Verbose) cerr<<" Step "<<i<<", lvl="<<SigmaNoise<<", || (z - zt)/z ||_infty = "<<speed<<", TolVar="<<TolVar<<endl;
		}	
   
		
	// New point
		if(Fast) 
		{
			told=tk;
			tk = (1.+sqrt(1.+4.*told*told))/2.;
			for(int k=0;k<n;k++) xtmp[k] = x[k] + (told-1)/tk * (x[k]-xold[k]);
		}
		else
			for(int k=0;k<n;k++) xtmp[k] = x[k];
		
		++i;
		
		Domain->recons(x,z);
		if(Positivity) for(int k=0;k<z.n_elem();k++) z(k) = z(k) *(z(k)>0);
 		image = z;
		
		
		return true;
	} else {
	  if(No_coarse) {
	    Domain->add_coarse_scale(x,CS);
	    Domain->recons(x,z);
	  }
	  
          image = z;
	  cerr<<"##########\nEnd CS.\n##########"<<endl;
	  return false;
	}
}

void Fista_LOFAR::load(float scale) {
	cout << "Loading original over reconstruction\n";
	FILE *f = fopen("simulating/original.map","r");
	int dim;
	fscanf(f,"%d\n",&dim);
	fltarray in;
	in.alloc(dim,dim,0); float ff;
	for (int i=0; i<dim; ++i)
	for (int j=0; j<dim; ++j) { fscanf(f,"%f\n",&ff); in(i,j) = ff*scale; }
	fclose(f);
	Domain->transform(in,x); Domain->transform(in,xtmp);
}


void Fista_LOFAR::zero_fine() {
        cout << "Zero fine scale\n";
	for (int i=0; i<z.nx()*z.ny(); ++i)  { x[i] = xtmp[i] = 0; }
}

#endif


