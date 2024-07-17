/* file:        siftlocalmax.c
** author:      Andrea Vedaldi
** description: Find local maximizer of multi-dimensional array.
**/

/* AUTORIGHTS
Copyright (c) 2006 The Regents of the University of California.
All Rights Reserved.

Created by Andrea Vedaldi
UCLA Vision Lab - Department of Computer Science

Permission to use, copy, modify, and distribute this software and its
documentation for educational, research and non-profit purposes,
without fee, and without a written agreement is hereby granted,
provided that the above copyright notice, this paragraph and the
following three paragraphs appear in all copies.

This software program and documentation are copyrighted by The Regents
of the University of California. The software program and
documentation are supplied "as is", without any accompanying services
from The Regents. The Regents does not warrant that the operation of
the program will be uninterrupted or error-free. The end-user
understands that the program was developed for research purposes and
is advised not to rely exclusively on the program for any reason.

This software embodies a method for which the following patent has
been issued: "Method and apparatus for identifying scale invariant
features in an image and use of same for locating an object in an
image," David G. Lowe, US Patent 6,711,293 (March 23,
2004). Provisional application filed March 8, 1999. Asignee: The
University of British Columbia.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND
ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE
MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include"mex.h"
#include<mexutils.c>
#include<stdlib.h>

/** Matlab driver.
 **/
#define greater(a,b) ((a) > (b)+threshold)

void
mexFunction(int nout, mxArray *out[],
            int nin, const mxArray *in[])
{
  int M, N, nM, nN ;
  const double* F_pt ;
  const double* N_pt ;
  int ndims ;
  int pdims = -1 ;
  int layerNo=0;
  int* offsets ;
  int* midx ;
  int* neighbors ;
  int nneighbors ;
  int* note;
  int* dims ;
  int* Ndims ;
  enum {F=0,THRESHOLD,Nb,P} ;
  enum {MAXIMA=0} ;
  double threshold = - mxGetInf() ;
  double threshold2 = 0.1;
  const int  *const_ndims;
  double ratio1, ratio2, nratio1, nratio2;

  /* ------------------------------------------------------------------
   *                                                Check the arguments
   * --------------------------------------------------------------- */
  if (nin < 1) {
    mexErrMsgTxt("At least one input argument is required.");
  } else if (nin > 4) {
    mexErrMsgTxt("At most three arguments are allowed.") ;
  } else if (nout > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

  /* The input must be a real matrix. */
  if (!mxIsDouble(in[F]) || mxIsComplex(in[F])) {
    mexErrMsgTxt("Input must be real matrix.");
  }
  if (!mxIsDouble(in[Nb]) || mxIsComplex(in[Nb])) {
    mexErrMsgTxt("Neighbor must be real matrix.");
  }
  if(nin > 1) {
    if(!uIsRealScalar(in[THRESHOLD])) {
      mexErrMsgTxt("THRESHOLD must be a real scalar.") ;
    }
    threshold = *mxGetPr(in[THRESHOLD]) ;
  }
  if(nin > 3) {
    if(!uIsRealScalar(in[P]))
      mexErrMsgTxt("P must be a non-negative integer") ;
    pdims = (int) *mxGetPr(in[P])  ;
    if(pdims < 0)
      mexErrMsgTxt("P must be a non-negative integer") ;
  }
  ndims = mxGetNumberOfDimensions(in[F]) ;
  {
    /* We need to make a copy because in one special case (see below)
       we need to adjust dims[].
    */
    int d ;
    const int* const_dims = (int*) mxGetDimensions(in[F]) ;
    dims = mxMalloc(sizeof(int)*ndims) ;
    for(d=0 ; d < ndims ; ++d) dims[d] = const_dims[d] ;
    const_ndims = mxGetDimensions(in[Nb]);
    Ndims = mxMalloc(sizeof(int)*2) ;
    Ndims[0] = *const_ndims; Ndims[1] = *(const_ndims+1);
  }
  M = dims[0] ;
  N = dims[1] ;
  nM = Ndims[0];
  nN = Ndims[1];
  F_pt = mxGetPr(in[F]) ;
  N_pt = mxGetPr(in[Nb]) ;
 
  //if(M!=Hdims[0] ||M!=Hdims[0] )
  //    mexErrMsgTxt("Dependency matrix not matching") ;
  /*
     If there are only two dimensions and if one is singleton, then
     assume that a vector has been provided as input (and treat this
     as a COLUMN matrix with p=1). We do this because Matlab does not
     distinguish between vectors and 1xN or Mx1 matrices and because
     the cases 1xN and Mx1 are trivial (the result is alway empty).
  */
  if((ndims == 2) && (pdims < 0) && (M == 1 || N == 1)) {
    pdims = 1 ;
    M = (M>N)?M:N ;
    N = 1 ;
    dims[0]=M ;
    dims[1]=N ;
  }

  /* search the local maxima along the first p dimensions only */
  if(pdims < 0)
    pdims = ndims ;

  if(pdims > ndims) {
    mxFree(dims) ;
    mexErrMsgTxt("P must not be greater than the number of dimensions") ;
  }

  /* ------------------------------------------------------------------
   *                                                         Do the job
   * --------------------------------------------------------------- */
  {
    int maxima_size = M*N ;
    int* maxima_start = (int*) mxMalloc(sizeof(int) * maxima_size) ;
    int* maxima_iterator = maxima_start ;
    int* maxima_end = maxima_start + maxima_size ;
    
    int maximadept_size = N;
    int* dept_start = (int*) mxMalloc(sizeof(int) * maximadept_size) ;
    int* dept_iterator = dept_start;
    int* dept_end = dept_start + maximadept_size;
    int dept_size=0;
    
    int i,j,h,o ;
    const double* pt = F_pt ;
    //const double* pt_start = F_pt; 
    const double* nt = N_pt ; 
    //mexPrintf("pt is =: %d \n", pt);
    //mexPrintf("pt_start is =: %d \n", F_pt);
    /* Compute the offsets between dimensions. */
    offsets = (int*) mxMalloc(sizeof(int) * ndims) ;
    offsets[0] = 1 ;
    for(h = 1 ; h < ndims ; ++h)
      offsets[h] = offsets[h-1]*dims[h-1] ;

    /* Multi-index. */
    midx = (int*) mxMalloc(sizeof(int) * ndims) ;
    for(h = 0 ; h < ndims ; ++h)
      midx[h] = 1 ;

    for(h = 0 ; h < pdims ; ++h) {
      //nneighbors *= 3 ;
      midx[h] = -1 ;
      //o -= offsets[h] ;
    }

    /* Starts at the corner (1,1,...,1,0,0,...0) */
    for(h = 0 ; h < pdims ; ++h) {
      midx[h] = 1 ;
      pt += offsets[h] ;
    }
    for(h = pdims ; h < ndims ; ++h) {
      midx[h] = 0 ;
    }

    /* ---------------------------------------------------------------
     *                                                            Loop
     * ------------------------------------------------------------ */

    /*
      If any dimension in the first P is less than 3 elements wide
      then just return the empty matrix (if we proceed without doing
      anything we break the carry reporting algorithm below).
    */
    for(h=0 ; h < pdims ; ++h)
      if(dims[h] < 3) goto end ;

    while(true) {

      /* Propagate carry along multi index midx */
      h = 0 ;
      while((midx[h]) >= dims[h] - 1) {
        pt += 2*offsets[h] ; /* skip first and last el. */
        midx[h] = 1 ;
        if(++h >= pdims)
          goto next_layer ;
        ++midx[h] ;
      }
      /*retrieve neighbor info*/
      //int nindex = ((pt - pt_start) % offset[1])*nM;
      //nindex *= nM;
      layerNo = ((pt-F_pt)%offsets[2])/offsets[1];
      //mexPrintf("%d  \n", layerNo) ;
      nneighbors = *(nt+layerNo);
      neighbors = (int*) mxMalloc(sizeof(int) * nneighbors) ;
      note = (int*) mxMalloc(sizeof(int) * 3) ;
      for(i = 0; i < nneighbors; i++)
      {neighbors[i] = *(nt+(i+3)*nM+layerNo);}
      note[0] = 0;
      note[1] = *(nt+nM+layerNo);
      note[2] = *(nt+2*nM+layerNo);
      /*
        for(h = 0 ; h < ndims ; ++h )
        mexPrintf("%d  ", midx[h]) ;
        mexPrintf(" -- %d -- pdims %d \n", pt - F_pt,pdims) ;
      */

      /*  Scan neighbors */
      {
        double v;
        bool is_greater = true; // dominate along time and equivalent along dependency
        bool is_greater2 = true; // dominate along time and dependency
        bool temp = false;
        bool temp2 = false;
        /*dominate through time*/
        for(j = 0; j < 3; j++){
            v = *(pt + note[j]);
            is_greater &= (v >= threshold) ;
            is_greater2 &= (v >= threshold) ;
            /*dominate through time*/
            is_greater &= (v > *(pt + note[j]+1))&&(v > *(pt + note[j]-1));
            is_greater2 &= (v > *(pt + note[j]+1))&&(v > *(pt + note[j]-1));
            if(j!=0)
                is_greater2 &= v > *(pt + note[j]);
            ratio1 = *(pt + note[j]+1)/v;
            ratio2 = *(pt + note[j]-1)/v;
            /*equivalent through neighbors in current scale*/
            i = 0;
            while((is_greater||is_greater2) && i < nneighbors)
            {//if(((pt + neighbors[i++]) >= F_pt) && ((pt + neighbors[i++]) <=F_pt+offsets[2]*dims[ndims]))
                is_greater2 &= v > *(pt + note[j] + neighbors[i]);
                nratio1 = *(pt + note[j] + neighbors[i]+1) / *(pt + note[j] + neighbors[i]);
                nratio2 = *(pt + note[j] + neighbors[i]-1) / *(pt + note[j] + neighbors[i]);
                temp2 = (fabs(nratio1-ratio1) < threshold2) && (fabs(nratio2-ratio2) < threshold2);
                if (neighbors[i] < 0 && temp2 == true)
                {temp = false; break;}
                else
                    temp |= temp2;
                i++;
            }
            is_greater &= temp;}
        /* Add the local maximum */
        if(is_greater || is_greater2) {
          /* Need more space? */
          if(maxima_iterator == maxima_end) {
            maxima_size += M*N ;
            maxima_start = (int*) mxRealloc(maxima_start,
                                            maxima_size*sizeof(int)) ;
            maxima_end = maxima_start + maxima_size ;
            maxima_iterator = maxima_end - M*N ;
          }

          *maxima_iterator++ = pt - F_pt + 1 ;
        }

        /* Go to next element */
        pt += 1 ;
        ++midx[0] ;
        mxFree(neighbors) ;
        continue ;

      next_layer: ;
        //pt_start = pt; 
        //layerNo++;
        //mexPrintf("next layer info %d\n",layerNo) ;
        if( h >= ndims )
          goto end ;

        while((++midx[h]) >= dims[h]) {
          midx[h] = 0 ;
          if(++h >= ndims)
            goto end ;
        }
      }
    }
  end:;
    /* Return. */
    {
      double* M_pt ;
      out[MAXIMA] = mxCreateDoubleMatrix
        (1, maxima_iterator-maxima_start, mxREAL) ;
      maxima_end = maxima_iterator ;
      maxima_iterator = maxima_start ;
      M_pt = mxGetPr(out[MAXIMA]) ;
      while(maxima_iterator != maxima_end) {
        *M_pt++ = *maxima_iterator++ ;
      }
    }

    /* Release space. */
    mxFree(offsets) ;
    mxFree(midx) ;
    mxFree(maxima_start) ;
    mxFree(Ndims);
    mxFree(dept_start);
  }
  mxFree(dims) ;
}
