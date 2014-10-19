
/*Attention!
 1, Not to Include Any NaN or Inf*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdlib.h" // rand, srand
#include "time.h"

double  REALMIN;

#define abs(a)   (a>0?a:(-1*a))
#define max(a,b) (a>b?a:b)
#define sign(x)  (x>0?x:(-1*x))

// fminsearch tolerence
#define MAXITER           100000
#define MAXFUNEVA  
#define TOLX
#define TOLF
// fminsearch default value
#define RHO               1 
#define CHI               2 
#define PSI               0.5 
#define SIGMA             0.5
#define USUAL_DELTA       0.05        // 5 percent deltas for non-zero terms
#define ZERO_TERM_DELTA   0.00025 // Even smaller delta for zero elements of x     
// fminsearch how value
#define INITAL_SIMPLEX    1
#define EXPAND            2
#define REFLECT           3
#define CONSTRACT_OUTSIDE 4
#define CONSTRACT_INSIDE  5
#define SHRINK            6

// icatb_myquafun/icatb_splinefun
#define QUAFUN            1
#define SPLINEFUN         2
#define PI                3.141592653589793


double meps(double num)
{
  int n;
  if(abs(num)<REALMIN)
    return pow(2.0,-1074);
  n=int(log(abs(num))/log(2));
  return pow(2.0,n)<=num?pow(2..0,n-52):pow(2.0,n-53); // if single in matlab, pow(2.0,n-)
}


void indvec(int *order, int length) //zero to n minus one
{
  int i;
  for(i=0;i<length;++i)
    *(order+i)=i;
}

void quickSort(double *data,int *rawOrder, int l, int r) //l,left; r,right
{
  double x;
  int i=l,j=r;
  x=*(data+i);
  if(i<j)
  {
    // back2front
    while(i<j&&*(data+j)>x)
      --j;
    if(i<j)
    {
      *(data+i)=*(data+j);
      *(rawOrder+(i++))=*(rawOder+j); // Anyway, dont forget to procede i.e. i++
    }
    // front2back
    while(i<j&&*(data+i)<x)
      ++i;
    if(i<j)
    {
      *(data+j)=*(data+i);
      *(rawOrder+(j++))=*(rawOder+i);
    }
    *(data+i)=x;
    quickSort(data,rawOrder,l,i-1);
    quickSort(data,rawOrder,j+1,r);
  }
}

double median(double *data)
{
  int length;
  length=sizeof(data)/sizeof(double);
  rawOrder=malloc(sizeof(data));
  for(i=0;i<length;++i)
    *(rawOrder+i)=i;
  quickSort(data,rawOrder,0,length-1);
  free(rawOrder);
  if(mod(length,2))
    return *(data+(length-1)/2)
  else
    return (*(data+(length-1)/2) + *(data+(length+1)/2))/2.0;
}

double icatb_quadObj(double *vTC, double *t,double *coff)
{ 
  int i=0, t2=0, t1=0, t0=0;
  int nTp;
  int val;
  nTp=sizeof(vTC)/sizeof(double);
  // smoothish quadratic fit
  for(i=0;i<nTp;++i)
  {
  	val=double(*(t+i));
  	diffval=*vTC - *coff * val * val + *(coff+1) *val + *(coff+2);
  	*err+=diffval*diffval;
  }
  return err;
}

double icatb_quadFit(doule *yfit,double *coff,int len,double TR)
{
  int i=0;
  int nTp;
  double tval;
  for(i=0;i<nTp;++i)
  {
    tval=i*TR;
    *(yfit+i)=0;
    *(yfit+i)=*(yfit+i) + *coff * tval * tval + *(coff+1) *tval + *(coff+2);
  }
}

double icatb_splObj(double *vTC, double *t,double *coff)
{
  int nTp, numP;
  int i;
  double val, diffval, err=0;
  double TR;

  nTp=sizeof(vTC)/sizeof(double);
  numP=nTp/30; //floor
  numP=(numP*30<nTp)?numP:(numP-1);
  TR=(*(t+nTp-1)-*(t))/(nTp-1);
  for(i=0;i<nTp;++i)
    {
      val=*(t+i);
      for(j=0;j<numP+1;++j) // 0, init diffval; 1~numP, spline fit;
        if(!j)
          diffval=*vTC - *(coff) * val -*(coff+1) * val*val;
        else
          diffval-=*(coff+j+2) * sin(2*PI*j*val/(nTp*TR)) + *(coff+j+2) * cos(2*PI*j*val/(nTp*TR));
      err+=diffval*diffval; 
    }
  return err;
}

void icatb_splFit(doule *yfit,double *coff,int len,double TR)
{
  int numP;
  numP=nTp/30; //floor
  numP=(numP*30<nTp)?numP:(numP-1);
  for(i=0;i<len;++i)
    {
     tval=i*TR;
     if(i==0)
      *(yfit+i)=*(coff)*tval;
     else if(i==1)
      *(yfit+i)+=*(coff+1)*tval*tval;
     else
      for(j=0;j<numP;++j)
        *(yfit+i)=*(coff+j+2) * sin(2*PI*j*val/(nTp*TR)) + *(coff+j+2) * cos(2*PI*j*val/(nTp*TR));
    }
}

void icatb_fun(double *xestimate, double *yestimate, double *vTC, double TR, int method)
{
  int nTp, numP, n;
	double *t;
	double err;
  double *x, *y;
	int i,iter,funeval,how;
  double *xbar, *xr, *xe, *xc;
  mxArray *mxv, *mxtmpv;
  double *v, *tmpv, *fv;
  double *rawOrder;

	nTp=sizeof(vTC)/sizeof(double);
  n=sizeof(xestimate)/sizeof(double);
	// init t
  t = malloc(nTp*double);
  for(i=0;i<nTp;++i)
    *(t+i)=i*TR;
  if(method==SPLINEFUN)
    numP=nt/30;
  //start point
  if(method==QUAFUN)
  {
    
  }
  else if(method==SPLINEFUN)
  {
    startpoint;
    srand((unsigned int)time(NULL));
    for(i=0;i<n;++i)
      *(xestimate)=double(rand());
  }

/*############### FMINSEARCH START ###############*/
  n=sizeof(xestimate)/sizeof(double);
  // Set up a simplex near the initial guess.
  mxv=mxCreateDoubleMatrix(n,n+1,mxReal);
  fv=malloc((n+1)*sizeof(double));
  x=malloc(n*sizeof(double));
  v=mxGetPr(mxv);
  for(i=0;i<n+1;++i) //column
    {
      for(j=0;j<n;++j)
        {
          val=*(x+j);
          if(i&&j==i-1) //exlude first column
            {
              if(val)
                val*=(1 + USUAL_DELTA);
              else
                val=ZERO_TERM_DELTA;
            }
            *(v+i*n+j)=val;
            *(x+j)=val;
        }
        y=icatb_quadObj(vTC,t,x);
        free(x);
        *(fv+i)=y;
    }
    // sort v & fv
    mxtmpv=mxDuplicateArray(v);
    tmpv=mxGetPr(mxtmpv);
    rawOrder=malloc((n+1)*sizeof(int));
    indvec(rawOrder, n+1);
    quickSort(fv,rawOrder,0,n);
    for(i=0;i<n+1;++i) //column
      {
        for(j=0;j<n;++j)
          *(v+i*n+j)=*(tmp_v+*(rawOrder+i)*n+j); // change column
      }
    mxDestroyArray(mxtmpv);
    how = INITAL_SIMPLEX;
    
    /*
    % Main algorithm: iterate until 
    % (a) the maximum coordinate difference between the current best point and the 
    % other points in the simplex is less than or equal to TOLX. Specifically,
    % until max(||v2-v1||,||v2-v1||,...,||v(n+1)-v1||) <= TOLX,
    % where ||.|| is the infinity-norm, and v1 holds the 
    % vertex with the current lowest value; AND
    % (b) the corresponding difference in function values is less than or equal
    % to TolFun. (Cannot use OR instead of AND.)
    % The iteration stops if the maximum number of iterations or function evaluations 
    % are exceeded
    */
    xbar=malloc(sizeof(double)*n);
    while(iter<maxIer)
      {
        // max(abs(fv(1)-fv(two2np1)))
        for(i=1;i<n+1;++i)
          if(abs(*fv-*(fv+i))>max(TOLF,10*eps(fv(1))))
            break;
        if(i>=n+1)
          {
            epsv=*v;
            for(j=1;j<n;++j)
              if(epsv<*(v+j))
                epsv=*(v+j);
            for(i=1;i<n+1;++i)
              for(j=0;j<n;++j)
                if(abs(*(v+i*n+j)-*(v+(i-1)*n+j))>max(TOLX,10*eps(fv(1))))
                  break;
          }
        if(i>=n+1)
          break;
        for(i=0;i<n;++i) //0~n-1
          for(j=0;j<n;++j)
              if(!j)
                *(xbar+j)=*(v+i*n+j)/n;
              else
                *(xbar+j)+=*(v+i*n+j)/n;
      
        // xr = (1 + RHO)*xbar - RHO*v(:,end);
        xr=malloc(sizeof(xbar));
        for(j=0;j<n;++j)
          *(xr+j)=*(xbar+j)*(1+RHO) - *(v+n*n+j)*RHO;
        fxr=icatb_quadObj(vTC, t, xr);
        if(fxr<*(fv))
         {
          xe=malloc(sizeof(xbar));
          for(j=0;j<n;++j)
              *(xe+j)=*(xbar+j)*(1+RHO*CHI) - *(v+n*n+j)*RHO*CHI;
          fxe=icatb_quadObj(vTC, t, xe);
          if(fxe<fxr)
           {
            for(j=0;j<n;++j)
              *(v+n*n+j)=*(xe+j);
            *(fv+n)=fxe;
            how=EXPAND; //how='expand'
           }
          else
           {
            for(j=0;j<n;++j)
              *(v+n*n+j)=*(xr+j);
            *(fv+n)=fxr;
            how=REFLECT; //how='reflect'
           }
          free(xe);
         }
        else
         {
          if(fxr<*(fv+n-1)) //fxr<fv(:,n)
           {
            for(j=0;j<n;++j)
              *(v+n*n+j)=*(xr+j)
            *(fv+n)=fxr;
            how=REFLECT; //how='reflect'
           }
          else
           {
            xc=malloc(sizeof(xbar));
            if(fxr<*(fv+n)) //if(fxr<fv(:,end));
            {
              for(j=0;j<n;++j)
                *(xc+j)=*(xbar+j)*(1+PSI*RHO) - *(v+n*n+j)*PSI*RHO;
              fxc=icatb_quadObj(vTC, t, xc);
              if(fxc<=fxr)
              {
                for(j=0;j<n;++j)
                  *(v+n*n+j)=*(xc+j);
                *(fv+n)=fxc;
                how=CONTRCT_OUTSIDE; //how='contract outside'
              }
              else
                how=SHRINK;
            }
            else
            {
              for(j=0;j<n;++j)
                *(xc+j)=*(xbar+j)*(1-PSI) - *(v+n*n+j)*PSI; //xcc = (1-PSI)*xbar + PSI*v(:,end);
              fxc=icatb_quadObj(vTC, t, xc);
              if(fxc<*(fv+n))
              {
                for(j=0;j<n;++j)
                  *(v+n*n+j)=*(xc+j);
                *(fv+n)=fxc;
                how=CONTRACT_INSIDE; //how='contract inside'                
              }
              else
                how=SHRINK; //how='shrink';
            }
            if(how==SHRINK)
              for(i=1;i<n+1;++i)
                {
                  for(j=0;j<n;++j)
                    {
                      *(v+i*n+j)=*(v+j) + (*(v+i*n+j) - *(v+j))*SIGMA; 
                      *(xc+j)=*(v+i*n+j);
                    }
                  *(fv+i)=icatb_quadObj(vTC, t, xc);
                }
            }
            free(xc);
           }
        indvec(rawOrder, n+1);
        quickSort(fv,rawOrder,0,n);
        mxtmpv=mxDuplicateArray(v);
        for(i=0;i<n+1;++i) //column
          {
           for(j=0;j<n;++j)
           *(v+i*n+j)=*(tmp_v + *(rawOrder+i)*n + j); // change column
          }
       }
    for(j=0;j<n;++j)
      *(xestimate+j)=*(v+j);
    *yestimate=*fv;
    //
    free(rawOrder);
    mxDestroyArray(mxv);
    mxDestroyArray(mxtmpv);
    free(fv);
    free(xbar);

/*############### FMINSEARCH  END  ###############*/

  free(t);
}






void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *data;
  int c=0;
  mxArray *mxout, *mxin, *mxstr;
  double *out;
  char *str;
  double *yfit;
  double median_tc,mad_res;
  double mad_res_array;
  struct fminestimate *estimate  

	data=mxGetPr(prhs[0]);
	nTp=mxGetM(data);
	nVox=mxGetN(data);
  
  mxout=mxCreateDoubleMatrix(1,1,mxREAL);
  out=mxGetPr(mxout);
  mxstr=mxCreateString("double");
  mexCallMatlab(1,mxout,0,mxstr,'REALMIN');
  REALMIN=*(mxGetPr(mxout));
  mxDestroyArray(mxout);
  mxDestroyArray(mxstr);
  
  yfit=malloc(sizeof(double)*nTp);
  tc=malloc(sizeof(double)*nTp);
  mad_res_array=malloc(sizeof(yfit));
  for(c=0;c<nVox;++c)
    {
      for(i=0;i<nTp;++i)
        *(tc+i)=*(data+nTp*c+i);

      ylerr=
      icatb_fun(xestimate, &yqerr, tc, TR);
      icatb_fun(xestimates,&yserr, tc, TR);
      //
      yerr=ylerr;
      miny=1;
      if(yerr>yqerr)
        miny=2;
      if(yerr>yserr)
        miny=3;
      //
      if(miny==1)
        yfit=[];
      else if(miny==2)
        icatb_quadFit(yfit,coff,nTp,TR);
      else if(miny==3)
        icatb_splFit(yfit,coff,nTp,TR);
      // residue
      for(i=0;i<nTp;++i)
        *(tc+i)-=*(yfit+i);
      quickSort(tc,rawOrder,0,nTp-1);
      median_tc=median(tc);
      for(i=0;i<nTp;++i)
        *(mad_res_array+i)=abs(*(tc+i)-median_tc);
      mad_res=median(mad_res_array);
      sigma=mad_res*pow(PI/2.0,0.5);
      for(i=0;i<nTp;++i)
        *(tc+i)=*(tc+i)/sigma;
      for(t=0;t<nTp;++t)
        if(abs(*(tc+t))>c1)
        {
          *(Data+nTp*c+t)=sign(*(tc+t)) * (c1+((c2-c1) * tanh((abs(*(tc+t)))-c1) / (c2-c1)));
          *(Data+nTp*c+t)=*(tc+t)*sigma + *(yfit+t);
        }
    }
    free(yfit);
    free(tc);
    free(mad_res_array); 
}
