/*
* Yongxiang Huang, last modification: 03/07/2008
* yongxianghuang@gmail.com
*
*/

/* This function is to estimate the joint probability density function of amplitude and frequency*/

#include <stdlib.h>
#include <stdio.h>
#include "mex.h"




/************************************************************************/
/*                                                                      */
/* MAIN FUNCTION                                                        */
/*                                                                      */
/************************************************************************/

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
  
    /* declarations */
    int i,j,nwf,nwa,k1,k2,Nf,Na,ntmp;
    
    double *instf,*insta,*wa,*wf; /*instantaneous frequency, instantaneous amplitude, resolution of amplitude resolution frequency*/
    double *zz;
  
    
/*     check input*/
    if (nrhs!=4)     mexErrMsgTxt("You have to give four parameters!");
    if (mxIsEmpty(prhs[0]))mexErrMsgTxt("Instantaneous frequency is empty!");
    if (mxIsEmpty(prhs[1]))mexErrMsgTxt("Instantaneous amplitude is empty!");
    if (mxIsEmpty(prhs[2]))mexErrMsgTxt("Frequency resolution is empty!");
    if (mxIsEmpty(prhs[3]))mexErrMsgTxt("Amplitude resolution is empty!");
    /* get input data */
    instf=mxGetPr(prhs[0]);
    insta=mxGetPr(prhs[1]);
    wf=mxGetPr(prhs[2]);
    wa=mxGetPr(prhs[3]);
  
    Nf=mxGetN(prhs[0]);
    ntmp=mxGetM(prhs[0]);
    if(ntmp>Nf)Nf=ntmp;
    
    Na=mxGetN(prhs[1]);
    ntmp=mxGetM(prhs[1]);
    if(ntmp>Na)Na=ntmp;
     
    nwf=mxGetN(prhs[2])-1;
    ntmp=mxGetM(prhs[2])-1;
    if(ntmp>nwf)nwf=ntmp;
    
    
    nwa=mxGetN(prhs[3])-1;
    ntmp=mxGetM(prhs[3])-1;
    if(ntmp>nwa)nwa=ntmp;
    
    if (Nf!=Na) mexErrMsgTxt("The length of instantaneous frequency and amplitude must be the same!");
    plhs[0]=mxCreateDoubleMatrix(nwa,nwf,mxREAL);
    zz=mxGetPr(plhs[0]);
    
    
  /*  for (i=0;i<n+1;i++) mexPrintf("%f\n",wa[i]);*/
    
    for (i=0;i<Nf;i++) /*bigest loop for each instantaneous frequency and amplitude*/
    {
        k1=0;
        for(j=0;j<nwa;j++) /*loop for amplitude*/
        {
            if(insta[i]<wa[0]){k1=-1;break;}
            if(insta[i]>=wa[j] && insta[i]<wa[j+1])
            {
                break;
            }
            k1=k1+1;
        
        }
        
        
        k2=0;
        for(j=0;j<nwf;j++) /*loop for instantaneous frequency*/
        {
            if(instf[i]<wf[0]){k2=-1;break;}
            if(instf[i]>=wf[j] && instf[i]<wf[j+1])
            {
                break;
            }
            k2=k2+1;
        
        }
        
        if (k1>0 && k1<nwa && k2>0 && k2<nwf)   {  zz[k2*nwa+k1]=zz[k2*nwa+k1]+1;}
    } /*finish boxcounting*/
    for(i=0;i<nwf;i++) for (j=0;j<nwa;j++) zz[i*nwa+j]=zz[i*nwa+j]/(wf[i+1]-wf[i])/(wa[j+1]-wa[j]);
}
