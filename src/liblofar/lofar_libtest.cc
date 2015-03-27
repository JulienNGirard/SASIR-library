/******************************************************************************
**                   Copyright (C) 2010, 2012 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  11 Jan. 2010
**
**
**    Author: Hugh Garsden
**    Date: July 2012
**
**/


#include <math.h>
#include "IM_IO.h"
#include "SB_Filter_float.h"
#include "writefits3d.h"
#include "Fista_LOFAR.h"
#include "MCA_LOFAR.h"
#include "FCur_float.h"
#include "FFTN_2D.h"
#include "CEA_comp_sens.h"


// Transforms parameters
type_sb_filter SB_Filter = F_MALLAT_7_9;

// Algorithms
enum {ALG_ISTA, ALG_FISTA, ALG_MCA, ALG_NESTEROV} typedef type_algo;

// Domains
enum {DOM_DWT, DOM_UWT, DOM_IWT, DOM_FCT, DOM_TV} typedef type_domain;

const int NbrDir2d = 16;

Bool Verbose=False;
static cfarray cMask, cRecons, cData;
static FFTN_2D *DataFT2D = new FFTN_2D;


/***************/

static void degrade_fftmask2d(cfarray &x, bool degrade, bool reverse)
{
      // The FFT is normalized so that the sum of the norm()'s is the same before and after
	if(!reverse)
	{

		Icomplex_f xFrame;
		xFrame.alloc(x.buffer(), x.ny(),x.nx());
		DataFT2D->fftn2d(xFrame, False, true);

		
		if(degrade) { x *= cMask;}
	}
	else
	{
		if(degrade) { x *= cMask;}
		Icomplex_f xFrame;
		xFrame.alloc(x.buffer(), x.ny(),x.nx());
		DataFT2D->fftn2d(xFrame, True, true);
	}

}

static void float_to_flt(float_array& f1, fltarray& f2) {
  for (int i=0; i<f1.size(); ++i) for (int j=0; j<f1.size(); ++j) f2(i,j) = f1[i][j];
}
static void float_to_cf(float_array& f1, cfarray& f2) {
  for (int i=0; i<f1.size(); ++i) for (int j=0; j<f1.size(); ++j) f2(i,j) = complex_f(f1[i][j], 0);
}
static void flt_to_float(fltarray& f1, float_array& f2) {
  for (int i=0; i<f1.nx(); ++i) for (int j=0; j<f1.ny(); ++j) f2[i][j] = f1(i,j);
}
static void float_to_float(float_array& f1, float_array& f2) {
  for (int i=0; i<f1.size(); ++i) for (int j=0; j<f1.size(); ++j) f2[i][j] = f1[i][j];
}
static void complex_to_cf(complex_array& f1, cfarray& f2) {
  for (int i=0; i<f1.size(); ++i) for (int j=0; j<f1.size(); ++j) f2(i,j) = complex_f(f1[i][j].real(), f1[i][j].imag());
}
static void cf_to_complex(cfarray& f1, complex_array& f2) {
  for (int i=0; i<f1.nx(); ++i) for (int j=0; j<f1.ny(); ++j) f2[i][j] = complex<float>(f1(i,j).real(), f1(i,j).imag());
}
static void cf_to_float(cfarray& f1, float_array& f2) {
  for (int i=0; i<f1.nx(); ++i) for (int j=0; j<f1.ny(); ++j) f2[i][j] = f1(i,j).real();
}



float correlation(fltarray& a1, fltarray& a2) {
  float mean_a1=0.0;
  float mean_a2=0.0;
  float s1=0.0;
  float s2=0.0;
  float covar=0.0;

  for (int i=0; i<a1.n_elem(); ++i) {
    mean_a1 += a1(i); mean_a2 += a2(i);
  } 
  
  mean_a1 /= a1.n_elem(); mean_a2 /= a2.n_elem();
   for (int i=0; i<a1.n_elem(); ++i) {
    covar += (a1(i)-mean_a1)*(a2(i)-mean_a2);
    s1 += (a1(i)-mean_a1)*(a1(i)-mean_a1);
    s2 += (a2(i)-mean_a2)*(a2(i)-mean_a2);
  }

  if ( s1 > 0 && s2 > 0 ) return covar/sqrt(s1)/sqrt(s2);
  else return 0;
}


void fourier_transform_and_mask(float_array& image, float_array& mask, float_array& dirty_image) {
  cfarray image_fft;
  
  image_fft.alloc(image.size(),image.size(),0); cMask.alloc(image.size(),image.size(),0);
  
  float_to_cf(image, image_fft); float_to_cf(mask, cMask);
  degrade_fftmask2d(image_fft, true, false);

  degrade_fftmask2d(image_fft, false, true);
  cf_to_float(image_fft, dirty_image);

  cMask.free();
}

void fft_difference(float_array& original, float_array& current, float_array& difference, float_array& mask, bool pos) {
  cfarray im1, im2;

  im1.alloc(original.size(), original.size(), 0); im2.alloc(original.size(), original.size(), 0); cMask.alloc(original.size(), original.size(), 0);
  
  if ( pos )  for (int i=0; i<current.size(); ++i) for (int j=0; j<current.size(); ++j) current[i][j] *= (current[i][j]>0); 
  float_to_cf(original, im1); float_to_cf(current, im2); float_to_cf(mask,cMask);
 

  degrade_fftmask2d(im1, true, false); degrade_fftmask2d(im2, true, false);
    
  im2 = im1-im2; 
  degrade_fftmask2d(im2, false, true);   
  cf_to_float(im2, difference);
   
  cMask.free();
}    




// Iterator routines -----------------------------------------------------------------------------	
static Fista_LOFAR *FistaExecutor;		// Being in a library these don't work unless dynamic
static MCA_LOFAR *MCAExecutor;
void compressed_sensing_begin(float_array& dirty_image, float_array& image, bool verbose, 
		int Niter, bool Positivity, float mu, float SigmaNoise, float TolVar, int NbrScale, 
		int Type_algo, int Type_domain, bool no_coarse, double ksigma, bool AllOutput) {

  
  fltarray Data;
  
  Verbose = verbose?True:False;  // convert to internal T/F type
  

  // Print parameters
  cout << "Compressed Sensing -----------" << endl << endl << "PARAMETERS: " << endl << endl;

  Type_algo = (type_algo)(Type_algo-1); Type_domain = (type_domain)(Type_domain-1); 

  
  cout << "NbrScale = " <<  NbrScale    << endl;  

	// Filtering
  if ( Niter>1 ) cout << " Number of iterations : " <<  Niter    << endl;
  if(Positivity)	cout << " The reconstruction is positive" << endl;


  // Generic transform
  FloatTrans *Domain = NULL;

  // Specific transforms
  Ortho_2D_WT_float *Domain_dwt = NULL;
  PAVE_2D_WT_float *Domain_uwt = NULL;
  ATROUS_2D_WT_float *Domain_iwt = NULL;
  FCUR_float *Domain_fct = NULL;
  PAVE_2D_WT_float *Domain_tv = NULL;

  // Filters fir the wavelet transforms and TV
  FilterAnaSynt SelectFilter;
  if(Type_domain == DOM_TV)
	  SelectFilter.alloc(F_HAAR);
  else
	  SelectFilter.alloc(F_MALLAT_7_9);
  SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
  
  Data.alloc(image.size(), image.size(), 0);

  switch(Type_domain)
  {
	  case DOM_DWT:
		  Domain_dwt = new Ortho_2D_WT_float(*SB1D);
		  Domain_dwt->alloc(Data, NbrScale); 
		  Domain = Domain_dwt;
		  break;
	  case DOM_UWT:
		  Domain_uwt = new PAVE_2D_WT_float(*SB1D); 
		  Domain_uwt->alloc(Data, NbrScale);
		  Domain = Domain_uwt;
		  break;
	  case DOM_IWT:
		  Domain_iwt = new ATROUS_2D_WT_float();
		  Domain_iwt->alloc(Data, NbrScale);
		  Domain_iwt->AdjointRec = True;
		  //Domain_iwt->set_positive_coef(true);
		  Domain = Domain_iwt;
		  break;
	  case DOM_TV:
		  Domain_tv = new PAVE_2D_WT_float(*SB1D); 
		  Domain_tv->alloc(Data, 2);
		  Domain = Domain_tv;
		  break;
          case DOM_FCT:
                  Domain_fct = new FCUR_float();
                  Domain_fct->alloc(Data, NbrScale,  NbrDir2d);
                  Domain_fct->get_norm_coeff(3);
                  Domain = Domain_fct;
		  cout << "  Selected curvelets\n";
                  break;
  }

  if(Type_algo == ALG_ISTA || Type_algo == ALG_FISTA)
  {
  // Fista init
	  FistaExecutor = new Fista_LOFAR(Domain);
	  FistaExecutor->Fast = (Type_algo == ALG_FISTA);
	  FistaExecutor->SigmaNoise = SigmaNoise;   
	  FistaExecutor->MaxNiter = Niter;
	  FistaExecutor->No_coarse = no_coarse;
	  FistaExecutor->Positivity = Positivity;
	  if ( ksigma > 0 ) FistaExecutor->Domain->UseMad = true;
	  FistaExecutor->ksigma = ksigma;
	  if(mu>0) FistaExecutor->Mu = mu;
	  if(TolVar>0) FistaExecutor->TolVar = TolVar;
	  FistaExecutor->AllOutput = AllOutput;
	  strcpy(FistaExecutor->NameOut,"awimager");

	    // Load the image 
	  float_to_flt(dirty_image, Data); 
	  cerr << "Entering Fista.begin" << endl;
	  // Fista algorithm
	  FistaExecutor->begin(Data); 
	  cerr << "Going out from Fista.begin" << endl;
	  // Outgoing is the first reconstructed image
	  flt_to_float(Data,image);
  }

  else if(Type_algo == ALG_MCA)
  {
      MCAExecutor = new MCA_LOFAR(Domain);          
      MCAExecutor->SigmaNoise = SigmaNoise;
      MCAExecutor->MaxNiter = Niter;
      MCAExecutor->No_coarse = no_coarse;
      MCAExecutor->Positivity = Positivity;
      if(mu>0)MCAExecutor->Mu = mu;
      if ( ksigma > 0 ) MCAExecutor->Domain->UseMad = true;
      MCAExecutor->ksigma = ksigma;
      if(TolVar>0) MCAExecutor->TolVar = TolVar;
      strcpy(MCAExecutor->NameOut,"awimager");
	    // Load the image 
      float_to_flt(dirty_image, Data); 

	// Fista algorithm
      MCAExecutor->begin(Data); 
	  
	  // Outgoing is the first reconstructed image
      flt_to_float(Data,image);
    
  }
  else {
	  cerr<<"Unknown or unimplemented algorithm for AWImager"<<endl;
	  exit(0);
  }

} 

bool compressed_sensing_next(float_array& difference, float_array& image, int _Type_algo) { 
  fltarray diff, im;	
  bool keep_going;
  type_algo Type_algo;
   
  diff.alloc(difference.size(), difference.size(), 0); im.alloc(difference.size(), difference.size(), 0);
  
  // Incoming dirty map fft, which is the difference |y-Hx|
  float_to_flt(difference, diff);

  Type_algo = (type_algo)(_Type_algo-1); 
  cerr << "Entering Fista.next" << endl;
  if ( Type_algo == ALG_FISTA || Type_algo == ALG_ISTA )  keep_going = FistaExecutor->next(diff,im);
  else keep_going = MCAExecutor->next(diff,im);
  cerr << "Going out from Fista.next" <<endl;
  // Outgoing is the next reconstructed image
  flt_to_float(im, image);
 
  return keep_going;
  
}  

// JG
//void compressed_sensing_load(float s) { FistaExecutor->load(s); }
