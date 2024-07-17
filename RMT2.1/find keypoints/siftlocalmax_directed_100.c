#include"mex.h"
#include<mexutils.c>
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
    int MDepd, NDepd ;
    int MTime, NTime ;
    int MBoth, NBoth ;
    
    const double* FDepd_pt ;
    const double* FTime_pt ;
    const double* FBoth_pt ;
    
    
    const double* H_pt ;
    const double* HT_pt ;
    int ndimsDepd ;
    int ndimsTime ;
    int ndimsBoth ;
    
    int pdimsDepd = -1 ;
    int pdimsTime = -1 ;
    int pdimsBoth = -1 ;
    
    int* offsetsDepd ;
    int* offsetsTime ;
    int* offsetsBoth ;
    
    int* midxDepd ;
    int* midxTime ;
    int* midxBoth ;
    
    
    double* neighborsDepd ;
    double* neighborsTime ;
    double* neighborsBoth ;
    
    int nneighborsDepd ;
    int nneighborsTime ;
    int nneighborsBoth ;
    
    
    int* dimsDepd ;
    int* dimsTime ;
    int* dimsBoth ;
    
    
    enum {DepdDiff=0, TimeDiff, BothDiff, THRESHOLD, H, HT, scaleDiff, P} ;
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
     * if(pdimsDepd < 0)
     * mexErrMsgTxt("PDepd must be a non-negative integer") ;
     * }*/
    
    ndimsDepd = mxGetNumberOfDimensions(in[DepdDiff]) ;
    ndimsTime = mxGetNumberOfDimensions(in[TimeDiff]) ;
    ndimsBoth = mxGetNumberOfDimensions(in[BothDiff]) ;
//     ndimsDepd = ndimsDepd-1;
//     ndimsTime = ndimsTime-1;
//     ndimsBoth = ndimsBoth-1;
    
    
    {
        /* We need to make a copy because in one special case (see below)
         * we need to adjust dims[].
         */
        int d ;
        const int* const_dimsDepd = (int*) mxGetDimensions(in[DepdDiff]) ;
        const int* const_dimsTime = (int*) mxGetDimensions(in[TimeDiff]) ;
        const int* const_dimsBoth = (int*) mxGetDimensions(in[BothDiff]) ;
        
        
        dimsDepd = mxMalloc(sizeof(int)*ndimsDepd) ;
        dimsTime = mxMalloc(sizeof(int)*ndimsTime) ;
        dimsBoth = mxMalloc(sizeof(int)*ndimsBoth) ;
        for(d=0 ; d < ndimsDepd ; ++d) dimsDepd[d] = const_dimsDepd[d] ;
        for(d=0 ; d < ndimsTime ; ++d) dimsTime[d] = const_dimsTime[d] ;
        for(d=0 ; d < ndimsBoth ; ++d) dimsBoth[d] = const_dimsDepd[d] ;
    }
    MDepd = dimsDepd[0] ;
    NDepd = dimsDepd[1] ;
    
    
    MTime = dimsTime[0] ;
    NTime = dimsTime[1] ;
    
    
    MBoth = dimsBoth[0] ;
    NBoth = dimsBoth[1] ;
    
    FDepd_pt = mxGetPr(in[DepdDiff]) ;
    FTime_pt = mxGetPr(in[TimeDiff]) ;
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
    if((ndimsDepd == 2) && (pdimsDepd < 0) && (MDepd == 1 || NDepd == 1)) {
        pdimsDepd = 1 ;
        MDepd = (MDepd > NDepd)?MDepd:NDepd ;
        NDepd = 1 ;
        dimsDepd[0]=MDepd ;
        dimsDepd[1]=NDepd ;
    }
    
    
    if((ndimsTime == 2) && (pdimsTime < 0) && (MTime == 1 || NTime == 1)) {
        pdimsTime = 1 ;
        MTime = (MTime > NTime)?MTime:NTime ;
        NTime = 1 ;
        dimsTime[0]=MTime ;
        dimsTime[1]=NTime ;
    }
    
    
    
    if((ndimsBoth == 2) && (pdimsBoth < 0) && (MBoth == 1 || NBoth == 1)) {
        pdimsBoth = 1 ;
        MBoth = (MBoth > NBoth)?MBoth:NBoth ;
        NBoth = 1 ;
        dimsBoth[0]=MBoth ;
        dimsBoth[1]=NBoth ;
    }
    
    
    
    /* search the local maxima along the first p dimensions only */
    if(pdimsDepd < 0){
        pdimsDepd = ndimsDepd ;
    }
    if(pdimsTime < 0){
        pdimsTime = ndimsTime ;
    }
    if(pdimsBoth < 0){
        pdimsBoth = ndimsBoth ;
    }
    
    if(pdimsDepd > ndimsDepd) {
        mxFree(dimsDepd) ;
        mexErrMsgTxt("P must not be greater than the number of dimensions") ;
    }
    if(pdimsTime > ndimsTime) {
        mxFree(dimsTime) ;
        mexErrMsgTxt("P must not be greater than the number of dimensions") ;
    }
    if(pdimsBoth > ndimsBoth) {
        mxFree(dimsBoth) ;
        mexErrMsgTxt("P must not be greater than the number of dimensions") ;
    }
    
    /* ------------------------------------------------------------------
     *                                                         Do the job
     * --------------------------------------------------------------- */
    {
//         int maxima_sizeDepd = MDepd*NDepd ;
//         int* maxima_startDepd = (int*) mxMalloc(sizeof(int) * maxima_sizeDepd) ;
//         int* maxima_iteratorDepd = maxima_startDepd ;
//         int* maxima_endDepd = maxima_startDepd + maxima_sizeDepd ;
//         
//         
//         int maxima_sizeTime = MTime*NTime ;
//         int* maxima_startTime = (int*) mxMalloc(sizeof(int) * maxima_sizeTime) ;
//         int* maxima_iteratorTime = maxima_startTime ;
//         int* maxima_endTime = maxima_startTime + maxima_sizeTime ;
        
        
        
        int maxima_sizeBoth = MBoth*NBoth ;
        int* maxima_startBoth = (int*) mxMalloc(sizeof(int) * maxima_sizeBoth) ;
        int* maxima_iteratorBoth = maxima_startBoth ;
        int* maxima_endBoth = maxima_startBoth + maxima_sizeBoth ;
        
        
        
        
        int tempHDepd, tempHTime, tempHBoth;
        int iDepd,hDepd,oDepd,tDepd,kDepd;
        int iTime,hTime,oTime,tTime,kTime;
        int iBoth,hBoth,oBoth,tBoth,kBoth;
        
        
        bool is_greaterDepd, is_greaterTime, is_greaterBoth;
        int indDepd, indTime, indBoth;
        const double* ptDepd = FDepd_pt ;
        const double* ptTime = FTime_pt ;
        const double* ptBoth = FBoth_pt ;
        
        
        /* Compute the offsets between dimensions -- Depd. */
        offsetsDepd = (int*) mxMalloc(sizeof(int) * ndimsDepd) ;
        offsetsDepd[0] = 1 ;
        
        /* Compute the offsets between dimensions -- Time. */
        offsetsTime = (int*) mxMalloc(sizeof(int) * ndimsTime) ;
        offsetsTime[0] = 1 ;
        
        /* Compute the offsets between dimensions -- Both. */
        offsetsBoth = (int*) mxMalloc(sizeof(int) * ndimsBoth) ;
        offsetsBoth[0] = 1 ;
        
        //for(hTime = 1 ,hDepd = 1,hBoth = 1  ; hDepd < ndimsDepd,hTime < ndimsTime,hBoth < ndimsBoth ; ++hDepd,++hTime,++hBoth)
        for(hDepd = 1  ; hDepd < ndimsDepd ; ++hDepd)
        {
            offsetsDepd[hDepd] = offsetsDepd[hDepd-1]*dimsDepd[hDepd-1] ;
            offsetsTime[hDepd] = offsetsTime[hDepd-1]*dimsTime[hDepd-1] ;
            offsetsBoth[hDepd] = offsetsBoth[hDepd-1]*dimsBoth[hDepd-1] ;
        }
        
        
        
        
        
        
        /* Multi-index -- Depd. */
        midxDepd = (int*) mxMalloc(sizeof(int) * ndimsDepd) ;
                /* Multi-index -- Time. */
        midxTime = (int*) mxMalloc(sizeof(int) * ndimsTime) ;
                /* Multi-index -- Both. */
        midxBoth = (int*) mxMalloc(sizeof(int) * ndimsBoth) ;
        for(hDepd = 0 ; hDepd < ndimsDepd ; ++hDepd)
        {
            midxDepd[hDepd] = 1 ;
            midxTime[hDepd] = 1 ;
            midxBoth[hDepd] = 1 ;
        }
            
        
        
        /* Neighbors -- Depd. */
        nneighborsDepd = 1 ;
        nneighborsTime = 1 ;
        nneighborsBoth = 1 ;
        oDepd=0 ;
        oTime=0 ;
        oBoth=0 ;
        //for(hDepd = 0,hTime = 0,hBoth = 0 ; hDepd < pdimsDepd,hTime < pdimsTime, hBoth < pdimsBoth ; ++hDepd,++hTime, ++hBoth) {
        for(hDepd = 0; hDepd < pdimsDepd ; ++hDepd) {
            //nneighborsDepd *= 3 ;
            midxDepd[hDepd] = -1 ;
            oDepd -= offsetsDepd[hDepd] ;
            midxTime[hDepd] = -1 ;
            oTime -= offsetsTime[hDepd] ;
            midxBoth[hDepd] = -1 ;
            oBoth -= offsetsBoth[hDepd] ;
        }
        nneighborsDepd  = 26;
        nneighborsTime = 26 ;
        nneighborsBoth = 26 ;
        
        //mexPrintf("pdimsBoth :%d, neighbor :%d.. \n ", pdimsBoth, nneighborsBoth);
        //nneighborsBoth -= 1 ;
        
        
        
        
        /* Precompute offsets from offset(-1,...,-1,0,...0) to
         * offset(+1,...,+1,0,...,0). */
        iDepd = 0 ;
        iTime=0;
        iBoth = 0 ;
         /* Starts at the corner (1,1,...,1,0,0,...0) */
        for(hDepd = 0 ; hDepd < pdimsDepd ; ++hDepd) {
            midxDepd[hDepd] = 1 ;
            ptDepd += offsetsDepd[hDepd] ;
            
            midxTime[hDepd] = 1 ;
            ptTime += offsetsTime[hDepd] ;
            
            
            midxBoth[hDepd] = 1 ;
            ptBoth += offsetsBoth[hDepd] ;
        }
        for(hDepd = pdimsDepd ; hDepd < ndimsDepd ; ++hDepd) {
            midxDepd[hDepd] = 0 ;
            midxTime[hDepd] = 0 ;
            midxBoth[hDepd] = 0 ;
        }
        
//         /* Precompute offsets from offset(-1,...,-1,0,...0) to
//          * offset(+1,...,+1,0,...,0). */
//         
//          /* Starts at the corner (1,1,...,1,0,0,...0) */
//         
//         for(hTime = 0 ; hTime < pdimsTime ; ++hTime) {
//            
//         }
//         for(hTime = pdimsTime ; hTime < ndimsTime ; ++hTime) {
//             
//         }
//         
//         
//         /* Precompute offsets from offset(-1,...,-1,0,...0) to
//          * offset(+1,...,+1,0,...,0). */
//         
//          /* Starts at the corner (1,1,...,1,0,0,...0) */
//         for(hBoth = 0 ; hBoth < pdimsBoth ; ++hBoth) {
//             
//         }
//         for(hBoth = pdimsBoth ; hBoth < ndimsBoth ; ++hBoth) {
//             
//         }
        
        
        
        
        /*mexPrintf("pdimsDepd:%d, pdimsTime:%d, pdimsBoth: %d. \n ", pdimsDepd, pdimsTime, pdimsBoth) ;
        mexPrintf("offsetsDepd1 :%d, offsetsDepd2 :%d, offsetsDepd3 : %d. \n ", offsetsDepd[0], offsetsDepd[1], offsetsDepd[2]) ;
        mexPrintf("offsetsTime1 :%d, offsetsTime2 :%d, offsetsTime3 : %d. \n ", offsetsTime[0], offsetsTime[1], offsetsTime[2]) ;
        mexPrintf("offsetsBoth1 :%d, offsetsBoth2 :%d, offsetsBoth3 : %d. \n ", offsetsBoth[0], offsetsBoth[1], offsetsBoth[2]) ;
        
        
        mexPrintf("midxDepd1 :%d, midxDepd2 :%d, midxDepd3 : %d. \n ", midxDepd[0], midxDepd[1], midxDepd[2]) ;
        mexPrintf("midxTime1 :%d, midxTime2 :%d, midxTime3 : %d. \n ", midxTime[0], midxTime[1], midxTime[2]) ;
        mexPrintf("midxBoth1 :%d, midxBoth2 :%d, midxBoth3 : %d. \n ", midxBoth[0], midxBoth[1], midxBoth[2]) ;
        
        mexPrintf("scaleD:%d. \n ", scaleD) ;*/
        
        
        
        /* ---------------------------------------------------------------
         *                                                            Loop
         * ------------------------------------------------------------ */
        
        /*
         * If any dimension in the first P is less than 3 elements wide
         * then just return the empty matrix (if we proceed without doing
         * anything we break the carry reporting algorithm below).
         */
        
      
        
        
        for(hDepd=0 ; hDepd < pdimsDepd ; ++hDepd)
            if(dimsDepd[hDepd] < 3) goto end ;
//         for(hTime=0 ; hTime < pdimsTime ; ++hTime)
//             if(dimsTime[hTime] < 3) goto end ;
//         for(hBoth=0 ; hBoth < pdimsBoth ; ++hBoth)
//             if(dimsBoth[hBoth] < 3) goto end ;
        
        
        while(true) {
            
            // Depd
            neighborsDepd = (double*) mxMalloc( nneighborsDepd*sizeof(double) ) ;
            neighborsTime = (double*) mxMalloc( nneighborsTime*sizeof(double) ) ;
            neighborsBoth = (double*) mxMalloc( nneighborsBoth*sizeof(double) ) ;
            /* Propagate carry along multi index midx */
            hDepd = 0 ;
            hTime = 0 ;
            hBoth = 0 ;
            /* skip first cube el. */
            //ptDepd += MDepd*NDepd*scaleD ;
            while((midxDepd[hDepd]) >= dimsDepd[hDepd] - 1) {
                ptDepd += 2*offsetsDepd[hDepd] ; /* skip first and last el. */
                midxDepd[hDepd] = 1 ;
                
                ptTime += 2*offsetsTime[hDepd] ; /* skip first and last el. */
                midxTime[hDepd] = 1 ;
                
                
                ptBoth += 2*offsetsBoth[hDepd] ; /* skip first and last el. */
                midxBoth[hDepd] = 1 ;
                
                if(++hDepd >= pdimsDepd)
                    goto next_layer ;
                ++midxDepd[hDepd] ;
                ++midxTime[hDepd] ;
                ++midxBoth[hDepd] ;
//                 if(++hTime >= pdimsTime)
//                     goto next_layer ;
//                 if(++hBoth >= pdimsBoth)
//                     goto next_layer ;
                
                
            }
            for(kDepd=0;kDepd<nneighborsDepd;kDepd++)
            {
                neighborsDepd[kDepd] = 0;
                neighborsTime[kDepd] = 0;
                neighborsBoth[kDepd] = 0;
            }
            
                
            
            
            
            
            /*
             * mexPrintf("%f  \n", *pt) ;
             * for(h = 0 ; h < ndims ; ++h )
             * mexPrintf("%d  ", midx[h]) ;
             * mexPrintf(" -- %d -- pdims %d \n", pt - F_pt,pdims) ;
             */
            /*  Compute neighbor values*/
            
            //current scale -- Depd
//             mexPrintf("Current Depd. \n") ;
            indDepd = 0;
            indTime = 0;
            indBoth = 0;
            for(tDepd=-1;tDepd<2;tDepd++){
                tTime = tDepd;
                tBoth = tDepd;
                neighborsDepd[indDepd++] = (double)*(ptDepd+tDepd*offsetsDepd[2] +1);
                neighborsDepd[indDepd++] = (double) *(ptDepd+tDepd*offsetsDepd[2] -1);
                if(tDepd!=0)
                    neighborsDepd[indDepd++] = (double)*(ptDepd+tDepd*offsetsDepd[2]);
                
                neighborsTime[indTime++] = (double)*(ptTime+tTime*offsetsTime[2] +1);
                neighborsTime[indTime++] = (double) *(ptTime+tTime*offsetsTime[2] -1);
                if(tTime!=0)
                    neighborsTime[indTime++] = (double)*(ptTime+tTime*offsetsTime[2] );
                
                
                neighborsBoth[indBoth++] = (double)*(ptBoth+tBoth*offsetsBoth[2] +1);
                neighborsBoth[indBoth++] = (double) *(ptBoth+tBoth*offsetsBoth[2]-1);
                if(tBoth!=0)
                    neighborsBoth[indBoth++] = (double)*(ptBoth+tBoth*offsetsBoth[2] );
            }
            
            
            // previous scale -- Depd
//             mexPrintf("Previous Depd. \n") ;
            for(tDepd=-1;tDepd<2;tDepd++){
                tTime = tDepd;
                tBoth = tDepd;
                neighborsDepd[indDepd] = 0;
                neighborsTime[indTime] = 0;
                neighborsBoth[indBoth] = 0;
                for(kDepd=0;kDepd<NDepd;kDepd++)
                {
                    kTime = kDepd;
                    kBoth = kDepd;
                    neighborsDepd[indDepd] += (double)*(H_pt+NDepd*kDepd+midxDepd[1]) * *(ptDepd-offsetsDepd[1]*midxDepd[1]+offsetsDepd[1]*kDepd+tDepd*offsetsDepd[2]);
                    neighborsTime[indTime] += (double)*(H_pt+NTime*kTime+midxTime[1]) * *(ptTime-offsetsTime[1]*midxTime[1]+offsetsTime[1]*kTime+tTime*offsetsTime[2]  );
                    neighborsBoth[indBoth] += (double)*(H_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2]  );
                }
                indDepd++;
                indTime++;
                indBoth++;
                neighborsDepd[indDepd] = 0;
                neighborsTime[indTime] = 0;
                neighborsBoth[indBoth] = 0;
                for(kDepd=0;kDepd<NDepd;kDepd++)
                {
                    kTime = kDepd;
                    kBoth = kDepd;
                    neighborsDepd[indDepd] += (double)*(H_pt+NDepd*kDepd+midxDepd[1]) * *(ptDepd-offsetsDepd[1]*midxDepd[1]+offsetsDepd[1]*kDepd+tDepd*offsetsDepd[2] +1 );
                    neighborsTime[indTime] += (double)*(H_pt+NTime*kTime+midxTime[1]) * *(ptTime-offsetsTime[1]*midxTime[1]+offsetsTime[1]*kTime+tTime*offsetsTime[2] +1);
                    neighborsBoth[indBoth] += (double)*(H_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2] +1 );
                }
                indDepd++;
                indTime++;
                indBoth++;
                neighborsDepd[indDepd] = 0;
                neighborsTime[indTime] = 0;
                neighborsBoth[indBoth] = 0;
                for(kDepd=0;kDepd<NDepd;kDepd++)
                {
                    kTime = kDepd;
                    kBoth = kDepd;
                    neighborsDepd[indDepd] += (double)*(H_pt+NDepd*kDepd+midxDepd[1]) * *(ptDepd-offsetsDepd[1]*midxDepd[1]+offsetsDepd[1]*kDepd+tDepd*offsetsDepd[2] -1 ) ;
                    neighborsTime[indTime] += (double)*(H_pt+NTime*kTime+midxTime[1]) * *(ptTime-offsetsTime[1]*midxTime[1]+offsetsTime[1]*kTime+tTime*offsetsTime[2] -1) ;
                    neighborsBoth[indBoth] += (double)*(H_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2] -1 ) ;
                }
                indDepd++;
                indTime++;
                indBoth++;
            }
            
            
            
            
            // next scale -- Depd
//             mexPrintf("Next Depd. \n") ;
            for(tDepd=-1;tDepd<2;tDepd++){
                tTime = tDepd;
                tBoth = tDepd;
                neighborsDepd[indDepd] = 0;
                neighborsTime[indTime] = 0;
                neighborsBoth[indBoth] = 0;
                for(kDepd=0;kDepd<NDepd;kDepd++)
                {
                    kTime = kDepd;
                    kBoth = kDepd;
                    neighborsDepd[indDepd] += (double)*(HT_pt+NDepd*kDepd+midxDepd[1]) * *(ptDepd-offsetsDepd[1]*midxDepd[1]+offsetsDepd[1]*kDepd+tDepd*offsetsDepd[2] );
                    neighborsTime[indTime] += (double)*(HT_pt+NTime*kTime+midxTime[1]) * *(ptTime-offsetsTime[1]*midxTime[1]+offsetsTime[1]*kTime+tTime*offsetsTime[2] );
                    neighborsBoth[indBoth] += (double)*(HT_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2] );
                }
                indDepd++;
                indTime++;
                indBoth++;
                neighborsBoth[indBoth] = 0;
                neighborsTime[indTime] = 0;
                neighborsDepd[indDepd] = 0;
                for(kDepd=0;kDepd<NDepd;kDepd++)
                {
                    kTime = kDepd;
                    kBoth = kDepd;
                    neighborsDepd[indDepd] += (double)*(HT_pt+NDepd*kDepd+midxDepd[1]) * *(ptDepd-offsetsDepd[1]*midxDepd[1]+offsetsDepd[1]*kDepd+tDepd*offsetsDepd[2] +1);
                    neighborsTime[indTime] += (double)*(HT_pt+NTime*kTime+midxTime[1]) * *(ptTime-offsetsTime[1]*midxTime[1]+offsetsTime[1]*kTime+tTime*offsetsTime[2] +1);
                    neighborsBoth[indBoth] += (double)*(HT_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2] +1);
                }
                indDepd++;
                indTime++;
                indBoth++;
                neighborsBoth[indBoth] = 0;
                neighborsTime[indTime] = 0;
                neighborsDepd[indDepd] = 0;
                for(kDepd=0;kDepd<NDepd;kDepd++)
                {
                    kTime = kDepd;
                    kBoth = kDepd;
                    neighborsDepd[indDepd] += (double)*(HT_pt+NDepd*kDepd+midxDepd[1]) * *(ptDepd-offsetsDepd[1]*midxDepd[1]+offsetsDepd[1]*kDepd+tDepd*offsetsDepd[2] -1) ;
                    neighborsTime[indTime] += (double)*(HT_pt+NTime*kTime+midxTime[1]) * *(ptTime-offsetsTime[1]*midxTime[1]+offsetsTime[1]*kTime+tTime*offsetsTime[2] -1) ;
                    neighborsBoth[indBoth] += (double)*(HT_pt+NBoth*kBoth+midxBoth[1]) * *(ptBoth-offsetsBoth[1]*midxBoth[1]+offsetsBoth[1]*kBoth+tBoth*offsetsBoth[2] -1) ;
                }
                indDepd++;
                indTime++;
                indBoth++;           
            }
            
            /*  Scan neighbors */
            {
                double vDepd = *ptDepd ;
                double vTime = *ptTime ;
                double vBoth = *ptBoth ;
                mexPrintf("%d  ", vDepd) ;
                mexPrintf("%d  ", vTime) ;
                mexPrintf("%d  ", vBoth) ;
                double tempV = max(vDepd, vTime);
                double v = max(tempV, vBoth);
                
                // double v = vDepd;
                
                is_greaterDepd = (vDepd >= threshold) ;
                //is_greaterTime = true ;
                //is_greaterBoth = true ;
                is_greaterTime = (vTime >= threshold) ;
                is_greaterBoth = (vBoth >= threshold) ;
                iDepd = 0  ;
                while(is_greaterDepd && iDepd < nneighborsDepd)
                    is_greaterDepd &= v >(neighborsDepd[iDepd++]) ;
                
                
                ////////////////
                iTime = 0  ;
                while(is_greaterTime && iTime < nneighborsTime)
                  is_greaterTime &= v >(neighborsTime[iTime++]) ;
                
                
                iBoth = 0  ;
                while(is_greaterBoth && iBoth < nneighborsBoth)
                    is_greaterBoth &= v >(neighborsBoth[iBoth++]) ;
                ///////////////
                
                /*bool is_greaterDepd = (vDepd >= threshold) ;
                i = 0  ;
                while(is_greaterDepd && i < nneighborsDepd)
                    is_greaterDepd &= v > neighborsDepd[i++] ;
                
                
                
                bool is_greaterTime = (vTime >= threshold) ;
                i = 0  ;
                while(is_greaterTime && i < nneighborsTime)
                    is_greaterTime &= v > neighborsTime[i++] ;
                
                
                
                
                bool is_greaterBoth = (vBoth >= threshold) ;
                i = 0  ;
                while(is_greaterBoth && i < nneighborsBoth)
                    is_greaterBoth &= v > neighborsBoth[i++] ;*/
                
                
                
                
                
                /* Add the local maximum */
                if(is_greaterDepd && is_greaterTime && is_greaterBoth) {
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
                ptDepd += 1 ;
                ptTime += 1 ;
                ptBoth += 1 ;
                
                ++midxDepd[0] ;
                ++midxTime[0] ;
                ++midxBoth[0] ;
                
                continue ;
                
                next_layer: ;
                //if( hDepd >= ndimsDepd || hTime >= ndimsTime || hBoth >= ndimsBoth ){
                if( hDepd >= ndimsDepd ){
                    goto end ;
                }
                
                //while((++midxDepd[hDepd]) >= dimsDepd[hDepd] || (++midxTime[hTime]) >= dimsTime[hTime]|| (++midxBoth[hBoth]) >= dimsBoth[hBoth]) {
                while((++midxDepd[hDepd]) >= dimsDepd[hDepd]) {
                    midxDepd[hDepd] = 0 ;
                    midxTime[hDepd] = 0 ;
                    midxBoth[hDepd] = 0 ;
                    tempHDepd = ++hDepd;
//                     tempHTime = ++hTime;
//                     tempHBoth = ++hBoth;
                    //if(tempHDepd >= ndimsDepd || tempHTime >= ndimsTime || tempHBoth >= ndimsBoth)
                    if(tempHDepd >= ndimsDepd)
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
        /*{
            double* M_pt ;
            out[MAXIMA] = mxCreateDoubleMatrix
                    (1, maxima_iterator-maxima_start, mxREAL) ;
            maxima_end = maxima_iterator ;
            maxima_iterator = maxima_start ;
            M_pt = mxGetPr(out[MAXIMA]) ;
            while(maxima_iterator != maxima_end) {
                *M_pt++ = *maxima_iterator++ ;
            }
        }*/
        
        /* Release space. */
        mxFree(offsetsDepd) ;
        mxFree(neighborsDepd) ;
        mxFree(midxDepd) ;
//         mxFree(maxima_startDepd) ;
        
        mxFree(offsetsTime) ;
        mxFree(neighborsTime) ;
        mxFree(midxTime) ;
//         mxFree(maxima_startTime) ;
        
        
        
        mxFree(offsetsBoth) ;
        mxFree(neighborsBoth) ;
        mxFree(midxBoth) ;
        mxFree(maxima_startBoth) ;
    }
    mxFree(dimsDepd) ;
    mxFree(dimsTime) ;
    mxFree(dimsBoth) ;
}
