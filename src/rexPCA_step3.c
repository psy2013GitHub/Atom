// d3(:)*data(i,:);

// arguments
// 1, eigData: nTimePoints * nDim (svd)
// 2, mapData: nTimePoints * nVox 

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  

  double *eigData, *mapData, *result;
  int nTp, nDim, nVox;
  int i,j,k;
  //ouble tmp;
  
  // Input
  eigData=mxGetPr(prhs[0]); mapData=mxGetPr(prhs[1]);
  nTp=mxGetM(prhs[0]); nDim=mxGetN(prhs[0]); nVox=mxGetN(prhs[1]);
  if(nTp!=mxGetM(prhs[1])) 
  	mexErrMsgTxt("row dimension not match between eigData and mapData");

  // Output
  plhs[0]=mxCreateDoubleMatrix(nVox,nDim,mxREAL);
  result=mxGetPr(plhs[0]);

  // eigDatarow * mapDatarow
  for(i=0;i<nDim;++i)
  	for(j=0;j<nVox;++j)
  	{
  		*(result+i*nVox+j)=0; // init as 0
  		for(k=0;k<nTp;++k)
  			*(result+i*nVox+j)+=(*(eigData+i*nTp+k)) * (*(mapData+j*nTp+k));
  	}
}