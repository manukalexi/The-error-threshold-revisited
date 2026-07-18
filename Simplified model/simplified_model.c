#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cvode/cvode.h>                  
#include <nvector/nvector_serial.h>       
#include <cvode/cvode_ls.h>               
#include <sunlinsol/sunlinsol_dense.h>    
#include <sundials/sundials_matrix.h>     
#include <sundials/sundials_types.h>

#define L (25) //sequence length
#define N (3) //number of Hamming classes = number of equations, H0 (master), H1 and junk


double mu; //mutation rate
double p00, p11, p10, p01, p02, p12; //transition probabilities
double sigma, lambda; //fitness of the master and the 1-error mutant
double totprod; 


static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  double x  = NV_Ith_S(y,0);
  double yy  = NV_Ith_S(y,1);
  double z  = NV_Ith_S(y,2);
  double Phi = sigma * x + lambda * yy + z; //outflow
  
  NV_Ith_S(ydot,0) = p00 * sigma * x + p10 * lambda * yy - x * Phi; //x
  NV_Ith_S(ydot,1) = p01 * sigma * x + p11 * lambda * yy  - yy * Phi; //y
  NV_Ith_S(ydot,2) = p02 * sigma * x + p12 * lambda * yy + z - z * Phi; //z
  
  return 0;
}

int main(void)
{
  double t0 = 0.0;
  double tmax = 10000.0;
  double dt_out = 0.01;
  double reltol = 1e-14;
  double abstol = 1e-16;
  
  int i, negind;
  
  sigma = 10.0;
  lambda = 5.0;

  for(int nn = 0; nn <= 1000; nn++)  
  {
    mu = 0.1 * nn/1000;
    
    p00 = pow((1.0 - mu), L);
    p01 = L * mu * pow((1.0 - mu), (L-1));    
    p02 = 1.0 - p00 - p01;

    //CASE 1: original, including back mutations
    p10 = mu / 3.0 * pow((1.0 - mu), (L-1)); 

    //CASE 2: modified, excluding back mutations
    //p10 = 0.0; 

    p11 = (1.0 - mu / 3.0)*pow((1.0 - mu), (L - 1));
    p12 = 1.0 - p10 - p11;
    
    SUNContext sunctx;
    if (SUNContext_Create(0, &sunctx))
    {
      fprintf(stderr, "Error: SUNContext_Create failed\n");
      return 1;
    }
  
    N_Vector y = N_VNew_Serial(N, sunctx);
    if (y == NULL)
    {
      fprintf(stderr, "Error: could not allocate N_Vector y\n");
      SUNContext_Free(&sunctx);
      return 1;
    }
  
    //initially, only masters are present
    NV_Ith_S(y,0) = 1.0; //x 
    NV_Ith_S(y,1) = 0.0; //y
    NV_Ith_S(y,2) = 0.0; //z
    
    void *cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (cvode_mem == NULL)
    {
      fprintf(stderr, "Error: could not create CVode memory\n");
      N_VDestroy(y);
      SUNContext_Free(&sunctx);
      return 1;
    }
  
    int flag = CVodeInit(cvode_mem, rhs, t0, y);
    if (flag != CV_SUCCESS)
    {
      fprintf(stderr, "Error: CVodeInit failed (%d)\n", flag);
      CVodeFree(&cvode_mem);
      N_VDestroy(y);
      SUNContext_Free(&sunctx);
      return 1;
    }
  
    CVodeSStolerances(cvode_mem, reltol, abstol);
    CVodeSetUserData(cvode_mem, NULL);
  
    SUNMatrix A = SUNDenseMatrix(N, N, sunctx);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    if(LS == NULL)
    {
      fprintf(stderr, "Error: could not create dense linear solver\n");
      CVodeFree(&cvode_mem);
      N_VDestroy(y);
      SUNMatDestroy(A);
      SUNContext_Free(&sunctx);
      return 1;
    }
  
    int flag_ls = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (flag_ls != CV_SUCCESS)
    {
      fprintf(stderr, "Error: CVodeSetLinearSolver failed (%d)\n", flag_ls);
      CVodeFree(&cvode_mem);
      N_VDestroy(y);
      SUNLinSolFree(LS);
      SUNMatDestroy(A);
      SUNContext_Free(&sunctx);
      return 1;
    }
  
    double t = t0;
    double tout = dt_out;
  
    while (t < tmax)
    {
      flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      if (flag < 0)
      {
        fprintf(stderr, "Integration failed, flag = %d\n", flag);
        break;
      }
  
      negind = 0;
      for(i = 0; i <= N; i++)
      {
        if(NV_Ith_S(y,i) < 0.0)
        {
          NV_Ith_S(y,i) = 0.0;
          negind = 1;
        }
      }
  
      if(negind == 1)
        flag = CVodeReInit(cvode_mem, t, y);
  
      tout += dt_out;
    }
    
    printf("%lf\t%le\t%le\t%le\n", mu, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
    
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(y);
    SUNContext_Free(&sunctx);
    CVodeFree(&cvode_mem);
  }
  
  return 0;
}
