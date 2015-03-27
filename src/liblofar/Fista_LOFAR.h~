
#ifndef _FISTA_LOFAR_H
#define _FISTA_LOFAR_H

#include <fstream>
#include "GlobalInc.h"
#include "FloatTrans.h"
#include "writefits3d.h"

extern Bool Verbose;
//extern float correlation(fltarray&, fltarray&);

class Fista_params
{
public:
        void reset();                   // Initialise the class with default parameters
        char NameOut[256], nameSCL[256]; // JG Mod
        int MaxNiter;                   // Maximum number of iterations
        float Threshold;				// Thresholding level
	float TolVar;                   // Stopping criterium : ||xt+1 - xt||_infty < TolVar
        float Mu;                               // Gradient regularization parameter
	bool Adj;
        bool Fast;                              // (F)ISTA
        bool Decreasing;                // if true : linearily decreasing threshold
        bool Positivity;                // The reconstruction must be >0
        bool No_coarse;                 // if true : kills the coarse scale before apaplying fista, then thresholds it
    float ksigma;           // threshold level = k * SigmaNoise.
    double SigmaNoise;      // Noise standard deviation:  lambda = ksigma * SigmaNoise
        bool AllOutput;
};
void Fista_params::reset()
{
  Adj= true;
	MaxNiter=10;
	Threshold = 1;
	TolVar = 1e-4;
	Mu=1;
	Fast = false;
	Decreasing = false;
	Positivity = false;
	No_coarse = false;
    ksigma=3.;
    SigmaNoise=0.;
	AllOutput=false;
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
	
//	void load(float s);	
	fltarray z;
	float *x, *xold, *xtmp; // current estimate, previous estimate, and temp variable
// JG
	fltarray SCLData;
	float *xs; // to store current scale

// Other variables
	float tk, told, lvl;
	int i, n, nscales; // JG
    bool done;
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
    TolVar = TolVar;
    tk = 1;
    lvl = Threshold;
    i = 0; done = false; speed = 0; old_speed = 0; allocD = true;
	
//cerr<<"z nx ny nl nc "<<z.nx()<<","<<z.ny()<<","<<z.nl()<<","<<z.nc()<<endl;
//cerr<<"b nx ny nl nc "<<b.nx()<<","<<b.ny()<<","<<b.nl()<<","<<b.nc()<<endl;
// Initial point z=b

	if(No_coarse) // Substract the coarse scale
	{

		Domain->transform(image0, xtmp, allocD); 
        allocD=false;// allocate xtmp and localTB
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
	if ( Domain->UseMad ) cout << " with MAD, ksigma=" << ksigma << ",";
	cout<<" with mu="<<Mu<<", and "<<MaxNiter<<" iterations.\n##########"<<endl;
}

bool Fista_LOFAR::next(fltarray &difference, fltarray& image) {	
	
// Fista algorithm. 
	if ( i<MaxNiter && !done  )
	{
cout<<"z i nx ny nl nc "<<i<<":"<<z.nx()<<","<<z.ny()<<","<<z.nl()<<","<<z.nc()<<endl;
	// Threshold update


        
		//if(Positivity) for(int k=0;k<z.n_elem();k++) z(k) = z(k) *(z(k)>0);
		// sprintf(filename,"%s_%05d.fits",NameOut,i);
        // writefltarr(filename, z);
		
	// Evaluate the evolution
		y = y-z; // error between the last two estimates
		if(  z.maxfabs() == 0 ) speed = 0;
		else speed = abs(y.maxfabs()/z.maxfabs()); 
		done = (speed < TolVar) && (i>0) && (speed<old_speed);
		old_speed=speed;
		cerr <<"before step" << endl;
		//if(Verbose) cerr<<" Step "<<i<<", lvl="<<lvl<<", || (z - zt)/z ||_infty = "<<speed<<", TolVar="<<TolVar<<endl;
		y=z; // Save the new solution
		cerr <<"step1: entering transform"<<endl;
	// Gradient step
		Domain->transform(difference, xold);
		cerr <<"step2"<<endl;
		if ( Domain->UseMad ) {
		for(int k=0;k<n;k++) xold[k] *= Mu;
		Domain->mad_calculation(xold);
		for(int k=0;k<n;k++) xold[k] /= Mu;
		}
		Domain->transform(difference, xold);
		for(int k=0;k<n;k++) xtmp[k] = xtmp[k] + Mu * xold[k];
		cerr <<"step3"<<endl;
	// Save estimate
		if(Fast) for(int k=0;k<n;k++) xold[k] = x[k];
		cerr <<"step4"<<endl;
	// Proximal operator : soft thresholding
		for(int k=0;k<n;k++) x[k] = xtmp[k];
		if ( Domain->UseMad == true)
		{
		cerr <<"Use MAD branch"<<endl;
		  double Threshold = ksigma + lvl;
		  Domain->mad_normalize(x);
		cerr<<"going out mad_normalize"<<endl;
		  Domain->soft_threshold(x,Threshold,No_coarse);    // threshold(SigmaNoise) = 0
		cerr<<"going out soft threshold"<<endl; 
		 Domain->mad_unnormalize(x);
		cerr<<"going out mad unormalize"<<endl;
		  if(Verbose) cerr<<" Step "<<i<<", lvl=MAD using "<<ksigma<<", || (z - zt)/z ||_infty = "<<speed<<", TolVar="<<TolVar<<endl;
		if(Verbose) cout<<" Step "<<i<<", lvl=MAD using "<<ksigma<<", || (z - zt)/z ||_infty = "<<speed<<", TolVar="<<TolVar<<endl;
}
		else
		{
		cerr<<"NO MAD branch"<<endl;
		  Domain->normalize(x);
		cerr<<"going out normalize"<<endl;
		  Domain->soft_threshold(x,lvl,No_coarse);		// ksigma=0
		cerr<<"going out softthreshold"<<endl;
		  Domain->unnormalize(x);
		cerr<<"going out unnormalize"<<endl;
		  if(Verbose) cerr<<" Step "<<i<<", lvl="<<SigmaNoise<<", || (z - zt)/z ||_infty = "<<speed<<", TolVar="<<TolVar<<endl;
		  if(Verbose) cout<<" Step "<<i<<", lvl="<<SigmaNoise<<", || (z - zt)/z ||_infty = "<<speed<<", TolVar="<<TolVar<<endl;
		
}	
   
		
	// New point
		if(Fast) 
		{
			told=tk;
			tk = (1.+sqrt(1.+4.*told*told))/2.;
			for(int k=0;k<n;k++) xtmp[k] = x[k] + (told-1)/tk * (x[k]-xold[k]);
		}
		else
		cerr<<"x to xtmp"<<endl;
			for(int k=0;k<n;k++) xtmp[k] = x[k];
		
		i++;
		
		//Domain->recons(x,z);
        if (Adj == True) 
        {
cerr<<"entering domain.recons Adj"<<endl;
           if(i!=0)Domain->adjoint_recons(x,z);
        }
        else
        {
//cerr<<"entering domain.recons"<<endl;
           if(i!=0) Domain->recons(x,z);
        }

		
		if(Positivity) for(int k=0;k<z.n_elem();k++) z(k) = z(k) *(z(k)>0);
//		cerr<<"after positivity"<<endl;
 		image = z;
//		cerr<<"after image=z"<<endl;
if (AllOutput) {
cerr << " Saving all scales ";
// INSPECTING SCALES
		xs=new float[n];// for(int k=0;k<n;k++) xs[k]=x[k]
		cerr <<" n="<<n<< endl;
		cerr <<" z.n_elem()="<<z.n_elem()<<endl;
		cerr <<" iter="<<i<<endl;
		nscales=n/z.n_elem();
		cerr <<" Nscales="<<nscales<<endl;

for (int iscale=0; iscale < nscales; iscale++)  // loop over the different scales 
	{
	cout <<"SAVING SCALE N#"<< iscale <<endl;
	for (int kk=0; kk<n; ++kk)
	{	
		if ( ( iscale*(z.n_elem()) <= kk && kk < (iscale+1)*(z.n_elem()) ) ) 
		{
			xs[kk]=x[kk]; // keep only one scale
		}
		else 
		{
			xs[kk]=0.; // put 0 to other scales
		}
	}
	SCLData=image;
//	SCLData.alloc(image.size(),image.size(),0);
	Domain->recons(xs,SCLData);
	sprintf(nameSCL,"iter%d-S%d",i,iscale);

	fits_write_fltarr(nameSCL,SCLData);
	}

// END INSPECTING SCALES
}
		
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


#endif


