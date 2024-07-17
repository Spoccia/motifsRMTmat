#include"mex.h"
#include"mexutils.c"
#include<stdlib.h>

/** Matlab driver.
 **/
#undef max
#undef min
#define greater(a,b) ((a) > (b)+threshold)
#define max(a,b)     (((a)>(b))?(a):(b))
#define min(a,b)     (((a)<(b))?(a):(b))

#define abs(a)       (((a)>0)?(a):(-a))


void
        mexFunction(int nout, mxArray *out[],
        int nin, const mxArray *in[])
{
    int MBoth, NBoth ;
    
    const double* FBoth_pt ;
    
    
    const double* H_pt ;
    const double* HT_pt ;
    
    int ndimsBoth ;
    
    int pdimsBoth = -1 ;
    
    int* offsetsBoth ;
    
    int* midxBoth ;
    
    double* neighborsBoth ;
    
    int nneighborsBoth ;
    
    int* dimsBoth ;
    
    
    enum {BothDiff = 0, THRESHOLD, H, HT, scaleDiff, P} ;
    enum {MAXIMA=0} ;
    double threshold = - mxGetInf() ;
    double val;
    int scaleD;
    
    /* ------------------------------------------------------------------
     *                                                Check the arguments
     * --------------------------------------------------------------- */
    if (nin < 1) {
        mexErrMsgTxt("At least one input argument is required.");
    }  else if (nout > 1) {
        mexErrMsgTxt("Too many output arguments");
    }
    
    if(nin > 1) {
        if(!uIsRealScalar(in[THRESHOLD])) {
            mexErrMsgTxt("THRESHOLD must be a real scalar.") ;
        }
        threshold = *mxGetPr(in[THRESHOLD]) ;
        scaleD = *mxGetPr(in[scaleDiff]) ;
        scaleD = scaleD-1;
    }
    
    /*if(nin > 4) {
        if(!uIsRealScalar(in[P]))
            mexErrMsgTxt("P must be a non-negative integer") ;
        pdims = (int) *mxGetPr(in[P])  ;
     * if(pdimsBoth < 0)
     * mexErrMsgTxt("PDepd must be a non-negative integer") ;
     * }*/
    
    ndimsBoth = mxGetNumberOfDimensions(in[BothDiff]) ;
    {
        /* We need to make a copy because in one special case (see below)
         * we need to adjust dims[].
         */
        int d ;
        const int* const_dimsBoth = (int*) mxGetDimensions(in[BothDiff]) ;
        
        dimsBoth = mxMalloc(sizeof(int)*ndimsBoth) ;
        for(d=0 ; d < ndimsBoth ; ++d) dimsBoth[d] = const_dimsBoth[d] ;
    }
    MBoth = dimsBoth[0] ;
    NBoth = dimsBoth[1] ;
    
    FBoth_pt = mxGetPr(in[BothDiff]) ;
    
    H_pt = mxGetPr(in[H]) ;
    HT_pt = mxGetPr(in[HT]) ;
    /*
     * If there are only two dimensions and if one is singleton, then
     * assume that a vector has been provided as input (and treat this
     * as a COLUMN matrix with p=1). We do this because Matlab does not
     * distinguish between vectors and 1xN or Mx1 matrices and because
     * the cases 1xN and Mx1 are trivial (the result is alway empty).
     */
    
    if((ndimsBoth == 2) && (pdimsBoth < 0) && (MBoth == 1 || NBoth == 1)) {
        pdimsBoth = 1 ;
        MBoth = (MBoth > NBoth)?MBoth:NBoth ;
        NBoth = 1 ;
        dimsBoth[0]=MBoth ;
        dimsBoth[1]=NBoth ;
    }
    
    /* search the local maxima along the first p dimensions only */
    if(pdimsBoth < 0){
        pdimsBoth = ndimsBoth ;
    }
    
    if(pdimsBoth > ndimsBoth) {
        mxFree(dimsBoth) ;
        mexErrMsgTxt("P must not be greater than the number of dimensions") ;
    }
    
    /* ------------------------------------------------------------------
     *                                                         Do the job
     * --------------------------------------------------------------- */
    {
        
        int maxima_sizeBoth = MBoth*NBoth ;
        int* maxima_startBoth = (int*) mxMalloc(sizeof(int) * maxima_sizeBoth) ;
        int* maxima_iteratorBoth = maxima_startBoth ;
        int* maxima_endBoth = maxima_startBoth + maxima_sizeBoth ;
        
        int tempHBoth;
        int iBoth,hBoth,oBoth,tBoth,kBoth;
        
        
        bool is_greaterBoth;
        int indBoth;
        const double* ptBoth = FBoth_pt ;
        
        /* Compute the offsets between dimensions -- Both. */
        offsetsBoth = (int*) mxMalloc(sizeof(int) * ndimsBoth) ;
        offsetsBoth[0] = 1 ;
        
        for(hBoth = 1  ; hBoth < ndimsBoth ; ++hBoth)
        {
            offsetsBoth[hBoth] = offsetsBoth[hBoth-1]*dimsBoth[hBoth-1] ;
        }
        
                /* Multi-index -- Both. */
        midxBoth = (int*) mxMalloc(sizeof(int) * ndimsBoth) ;
        for(hBoth = 0 ; hBoth < ndimsBoth ; ++hBoth)
        {
            midxBoth[hBoth] = 1 ;
        }
         
        nneighborsBoth = 1 ;
        oBoth=0 ;
        for(hBoth = 0; hBoth < pdimsBoth ; ++hBoth) {
            //nneighborsBoth *= 3 ;
            midxBoth[hBoth] = -1 ;
            oBoth -= offsetsBoth[hBoth] ;
        }
        nneighborsBoth = 26 ;
        
        /* Precompute offsets from offset(-1,...,-1,0,...0) to
         * offset(+1,...,+1,0,...,0). */
        iBoth = 0 ;
         /* Starts at the corner (1,1,...,1,0,0,...0) */
        for(hBoth = 0 ; hBoth < pdimsBoth ; ++hBoth) {
            midxBoth[hBoth] = 1 ;
            ptBoth += offsetsBoth[hBoth] ;
        }
        for(hBoth = pdimsBoth ; hBoth < ndimsBoth ; ++hBoth) {
            midxBoth[hBoth] = 0 ;
        }
        
        
        /* ---------------------------------------------------------------
         *                                                            Loop
         * ------------------------------------------------------------ */
        
        /*
         * If any dimension in the first P is less than 3 elements wide
         * then just return the empty matrix (if we proceed without doing
         * anything we break the carry reporting algorithm below).
         */
         for(hBoth=0 ; hBoth < pdimsBoth ; ++hBoth)
             if(dimsBoth[hBoth] < 3) goto end ;
        while(true) {
            neighborsBoth = (double*) mxMalloc( nneighborsBoth*sizeof(double) ) ;
            /* Propagate carry along multi index midx */
            hBoth = 0 ;
            /* skip first cube el. */
            while((midxBoth[hBoth]) >= midxBoth[hBoth] - 1) {
                ptBoth += 2*offsetsBoth[hBoth] ; /* skip first and last el. */
                midxBoth[hBoth] = 1 ;
                
                if(++hBoth >= pdimsBoth)
                    goto next_layer ;
                ++midxBoth[hBoth] ;
                
                
            }
            for(kBoth=0;kBoth<nneighborsBoth;kBoth++)
            {
                neighborsBoth[kBoth] = 0;
            }
            
            /*  Compute neighbor values*/
            
            indBoth = 0;
            for(tBoth=-1;tBoth<2;tBoth++){
                tBoth = tBoth;
                
                neighborsBoth[indBoth++] = (double)*(ptBoth+tBoth*offsetsBoth[2] +1);
                neighborsBoth[indBoth++] = (double) *(ptBoth+tBoth*offsetsBoth[2]-1);
                if(tBoth!=0)
                    neighborsBoth[indBoth++] = (double)*(ptBoth+tBoth*offsetsBoth[2] );
            }
            
            for(tBoth=-1;tBoth<2;tBoth++){
                tBoth = tBoth;
                neighborsBoth[indBoth] = 0;
                for(kBoth=0;kBoth<NBoth;kBoth++)
                {
                    kBoth = kBoth;
                    neighborsBoth[indBoth] += (double)*(H_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2]  );
                }
                indBoth++;
                neighborsBoth[indBoth] = 0;
                for(kBoth=0;kBoth<NBoth;kBoth++)
                {
                    kBoth = kBoth;
                    neighborsBoth[indBoth] += (double)*(H_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2] +1 );
                }
                indBoth++;
                neighborsBoth[indBoth] = 0;
                for(kBoth=0;kBoth<NBoth;kBoth++)
                {
                    kBoth = kBoth;
                    neighborsBoth[indBoth] += (double)*(H_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2] -1 ) ;
                }
                indBoth++;
            }
            
            for(tBoth=-1;tBoth<2;tBoth++){
                tBoth = tBoth;
                neighborsBoth[indBoth] = 0;
                for(kBoth=0;kBoth<NBoth;kBoth++)
                {
                    kBoth = kBoth;
                    neighborsBoth[indBoth] += (double)*(HT_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2] );
                }
                indBoth++;
                neighborsBoth[indBoth] = 0;
                for(kBoth=0;kBoth<NBoth;kBoth++)
                {
                    kBoth = kBoth;
                    neighborsBoth[indBoth] += (double)*(HT_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2] +1);
                }
                indBoth++;
                neighborsBoth[indBoth] = 0;
                for(kBoth=0;kBoth<NBoth;kBoth++)
                {
                    kBoth = kBoth;
                    neighborsBoth[indBoth] += (double)*(HT_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2] -1) ;
                }
                indBoth++;           
            }
            
            /*  Scan neighbors */
            {
                double vBoth = *ptBoth ;
                mexPrintf("%d  ", vBoth) ;
                // double v = max(tempV, vBoth);
                double v = vBoth;
                
                is_greaterBoth = (vBoth >= threshold) ;
                
                iBoth = 0  ;
                while(is_greaterBoth && iBoth < nneighborsBoth)
                    is_greaterBoth &= v >(neighborsBoth[iBoth++]) ;
                 
                /* Add the local maximum */
                if(is_greaterBoth) {
                    /* Need more space? */
                    if(maxima_iteratorBoth == maxima_endBoth) {
                        maxima_sizeBoth += MBoth*NBoth ;
                        maxima_startBoth = (int*) mxRealloc(maxima_startBoth,
                                maxima_sizeBoth*sizeof(int)) ;
                        maxima_endBoth = maxima_startBoth + maxima_sizeBoth ;
                        maxima_iteratorBoth = maxima_endBoth - MBoth*NBoth ;
                    }
                    
                    *maxima_iteratorBoth++ = ptBoth - FBoth_pt + 1 ;
                }
                
                /* Go to next element */
                ptBoth += 1 ;
                
                ++midxBoth[0] ;
                
                continue ;
                
                next_layer: ;
                if( hBoth >= ndimsBoth ){
                    goto end ;
                }
                 
                while((++midxBoth[hBoth]) >= dimsBoth[hBoth]) {
                    midxBoth[hBoth] = 0 ;
                    tempHBoth = ++hBoth;
                    if(tempHBoth >= ndimsBoth)
                    goto end ;
                }
            }
        }
        end:;
        /* Return. */
        {
            double* M_pt ;
            out[MAXIMA] = mxCreateDoubleMatrix
                    (1, (maxima_iteratorBoth-maxima_startBoth), mxREAL) ;
            maxima_endBoth = maxima_iteratorBoth ;
            maxima_iteratorBoth = maxima_startBoth ;
            M_pt = mxGetPr(out[MAXIMA]) ;
            while(maxima_iteratorBoth != maxima_endBoth) {
                *M_pt++ = *maxima_iteratorBoth++ ;
            }
        }
        
        /* Release space. */   
        mxFree(offsetsBoth) ;
        mxFree(neighborsBoth) ;
        mxFree(midxBoth) ;
        mxFree(maxima_startBoth) ;
    }
    mxFree(dimsBoth) ;
}
