
#define LOWE_COMPATIBLE

#include"mexutils.c"
#include<stdlib.h>
#include<math.h>

#ifdef WINDOWS
#include<string.h>
#ifndef __cplusplus
#define sqrtf(x)    ((float)sqrt((double)(x)))
#define powf(x,y)   ((float)pow((double)(x),(double)(y)))
#define fabsf(x)    ((float)fabs((double)(x)))
#define sinf(x)     ((float)sin((double)(x)))
#define cosf(x)     ((float)cos((double)(x)))
#define expf(x)     ((float)exp((double)(x)))
#define atan2f(x,y) ((float)atan2((double)(x),(double)(y)))
#endif
#else
#include<strings.h>
#endif

/* Altivec and Accelerate support.
 * Very crude at this time.
 */
#if defined( MACOSX ) && defined( __ALTIVEC__ )
#include<Accelerate/Accelerate.h>
typedef union 
{
  float x[4] ;
  vFloat vec ;
} float4 ;
#endif

#define greater(a,b) a > b
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


enum {SCALESPACE, NOSCALESPACE} ;

enum {PROP_MAGNIF=0,
      PROP_NBP,
      PROP_NBO,
      PROP_UNKNOWN} ;

char const * properties [4] = 
  { "Magnif",
    "NumSpatialBins",
    "NumOrientBins",
     0L
  } ;

/** Fast fmodf for 2*PI
 **/
/*inline*/
float fast_mod(float th)
{
  while(th < 0) th += 2*M_PI ;
  while(th > 2*M_PI) th -= 2*M_PI ;
  return th ;
}

/** Fast floor. Equivalent to (int) floor(x)
 **/
/*inline*/
int fast_floor(float x)
{
  return (int)( x - ((x>=0)?0:1) ) ; 
}

/** Normalizes in norm L_2 a descriptor.
 **/
void
normalize_histogram(float* L_begin, float* L_end)
{
  float* L_iter ;
  float norm=0.0 ;

  for(L_iter = L_begin; L_iter != L_end ; ++L_iter)
    norm += (*L_iter) * (*L_iter) ;

  norm = sqrtf(norm) ;
  /*  mexPrintf("%f\n",norm) ;*/

  for(L_iter = L_begin; L_iter != L_end ; ++L_iter)
    *L_iter /= norm ;
}

/** @brief MATLAB Driver.
 **/
void
mexFunction(int nout, mxArray *out[], 
            int nin, const mxArray *in[])
{
  int M,N,S=0,smin=0,K,num_levels=0 ;
  // const int* dimensions ;
  const size_t * dimensions ;
  const double* P_pt ;
  const double* G_pt ;
  float* descr_pt ;
  float* buffer_pt ;
  float sigma0 ;
  float magnif = 3.0f ; /* Spatial bin extension factor. */
  int NBP = 4 ;         /* Number of bins for one spatial direction (even). */
  int NBO = 8 ;         /* Number of bins for the ortientation. */
  int mode = NOSCALESPACE ;
  int buffer_size=0;

  enum {IN_G=0,IN_P,IN_SIGMA0,IN_S,IN_SMIN} ;
  enum {OUT_M=0,OUT_N} ;

  /* ------------------------------------------------------------------
  **                                                Check the arguments
  ** --------------------------------------------------------------- */ 
 
  if (nin < 3) {
    mexErrMsgTxt("At least three arguments are required") ;
  } else if (nout > 2) {
    mexErrMsgTxt("Too many output arguments.");
  }
		
  if( !uIsRealScalar(in[IN_SIGMA0]) ) {
    mexErrMsgTxt("SIGMA0 should be a real scalar") ;
  }
	
  if(!mxIsDouble(in[IN_G]) ||
     mxGetNumberOfDimensions(in[IN_G]) > 3) {
    mexErrMsgTxt("G should be a real matrix or 3-D array") ;
  }
  
  sigma0 = (float) *mxGetPr(in[IN_SIGMA0]) ;
  
  dimensions = mxGetDimensions(in[IN_G]) ;
  M = dimensions[0] ;
  N = dimensions[1] ;
  G_pt = mxGetPr(in[IN_G]) ;
  mexPrintf("\n M value (number of rows): %f, N value (number of columns) %f.\n", dimensions[0],dimensions[1]);
  P_pt = mxGetPr(in[IN_P]) ;	
  K = mxGetN(in[IN_P]) ;
  
  if( !uIsRealMatrix(in[IN_P],-1,-1)) {
    mexErrMsgTxt("P should be a real matrix") ;
  }

  if ( mxGetM(in[IN_P])  == 4) {
    /* Standard (scale space) mode */ 
    mode = SCALESPACE ;
    num_levels = dimensions[2] ;
    
    if(nin < 5) {
      mexErrMsgTxt("Five arguments are required in standard mode") ;
    }
    
    if( !uIsRealScalar(in[IN_S]) ) {
      mexErrMsgTxt("S should be a real scalar") ;
    }
    
    if( !uIsRealScalar(in[IN_SMIN]) ) {
      mexErrMsgTxt("SMIN should be a real scalar") ;
    }
    
    if( !uIsRealMatrix(in[IN_P],4,-1)) {
      mexErrMsgTxt("When the e mode P should be a 4xK matrix.") ;
    }
    
    S = (int)(*mxGetPr(in[IN_S])) ;
    smin = (int)(*mxGetPr(in[IN_SMIN])) ;
    
  } else if (  mxGetM(in[IN_P])  == 3 ) {
    mode = NOSCALESPACE ;
    num_levels = 1 ;
    S      = 1 ;
    smin   = 0 ;
  } else {
    mexErrMsgTxt("P should be either a 3xK or a 4xK matrix.") ;
  }

  /* Parse the property-value pairs */
  {
    char str [80] ;
    int arg = (mode == SCALESPACE) ? IN_SMIN + 1 : IN_SIGMA0 + 1 ;

    while(arg < nin) {
      int k ;

      if( !uIsString(in[arg],-1) ) {
        mexErrMsgTxt("The first argument in a property-value pair"
                     " should be a string") ;
      }
      mxGetString(in[arg], str, 80) ;

#ifdef WINDOWS      
      for(k = 0 ; properties[k] && strcmpi(str, properties[k]) ; ++k) ;
#else
      for(k = 0 ; properties[k] && strcasecmp(str, properties[k]) ; ++k) ;
#endif

      switch (k) {
      case PROP_NBP:
        if( !uIsRealScalar(in[arg+1]) ) {
          mexErrMsgTxt("'NumSpatialBins' should be real scalar") ;
        }
        NBP = (int) *mxGetPr(in[arg+1]) ;
        if( NBP <= 0 || (NBP & 0x1) ) {
          mexErrMsgTxt("'NumSpatialBins' must be positive and even") ;
        }
        break ;

      case PROP_NBO:
        if( !uIsRealScalar(in[arg+1]) ) {
          mexErrMsgTxt("'NumOrientBins' should be a real scalar") ;
        }
        NBO = (int) *mxGetPr(in[arg+1]) ;
        if( NBO <= 0 ) {
          mexErrMsgTxt("'NumOrientlBins' must be positive") ;
        }
        break ;

      case PROP_MAGNIF:
        if( !uIsRealScalar(in[arg+1]) ) {
          mexErrMsgTxt("'Magnif' should be a real scalar") ;
        }
        magnif = (float) *mxGetPr(in[arg+1]) ;
        if( magnif <= 0 ) {
          mexErrMsgTxt("'Magnif' must be positive") ;
        }
        break ;

      case PROP_UNKNOWN:
        mexErrMsgTxt("Property unknown.") ;
        break ;
      }
      arg += 2 ;
    }
  }
  
  /* -----------------------------------------------------------------
   *                                   Pre-compute gradient and angles
   * -------------------------------------------------------------- */
  /* Alloc two buffers and make sure their size is multiple of 128 for
   * better alignment (used also by the Altivec code below.)
   */
  buffer_size = (M*N*num_levels + 0x7f) & (~ 0x7f) ;
  buffer_pt = (float*) mxMalloc( sizeof(float) * 2 * buffer_size ) ;

  {
    /* Offsets to move in the scale space. */
    const int yo = 1 ;
    const int xo = M ;
    const int so = M*N ;
    int x,y,s ;

#define at(x,y) (*(pt + (x)*xo + (y)*yo))

    /* Compute the gradient */
    
      const double* pt = G_pt + s*so ;
      for(x = 1 ; x < N-1 ; ++x) {
        for(y = 1 ; y < M-1 ; ++y) {
          float Dx = 0.5 * ( at(x+1,y) - at(x-1,y) ) ;
          float Dy = 0.5 * ( at(x,y+1) - at(x,y-1) ) ;
          buffer_pt[(x*xo+y*yo+s*so) + 0          ] = Dx ;
          buffer_pt[(x*xo+y*yo+s*so) + buffer_size] = Dy ;
        }
      }
    
    

/* Restore pointer to the beginning of the descriptors. */

  {
    int k ;
    double* M_pt ;
    double* N_pt ;
    

    

        mwSize sz[2];
        sz[0] = M; // matlab is row first
        sz[1] = N;
        out[OUT_M] = mxCreateNumericArray( 2, sz, mxDOUBLE_CLASS, // create double array, you can change the type here
                                      mxREAL ); // create real 
        M_pt = mxGetPr(out[OUT_M]) ;
        
        out[OUT_N] = mxCreateNumericArray( 2, sz, mxDOUBLE_CLASS, // create double array, you can change the type here
                                      mxREAL ); // create real 
        N_pt = mxGetPr(out[OUT_N]) ;
        
        int x, y;
        int yo1 = 1;
        int xo1 = M;

            for(x=0; x < N; x++)
            {
                for(y = 0; y < M; y++)
                {
                    int index = x*xo1 + y*yo1; 
                    M_pt[index] = buffer_pt[index];
                }
            }
        
        

            for(x=0; x < N; x++)
            {
                for(y = 0; y < M; y++)
                {
                    int index = x*xo1 + y*yo1;
                    N_pt[index] = buffer_pt[index + buffer_size];
                }
            }

    }
  }
  mxFree(buffer_pt) ;
}
