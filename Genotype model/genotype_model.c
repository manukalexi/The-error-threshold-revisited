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
#define sigma_master (10) //reproductive superiority of the master
#define n_scal (0.2) //fitness landscape geometry factor

double p[L+1][L+1], sigma[L+1], mu;
int N; //number of Hamming classes = number of equations, N = L+1

void seed(long seed);

long double factorial(int k){
  long double fac = 1.0;
  int i;

  if((k == 0) || (k==1)){
    return(1.0);
  }
    
  for(i = 2; i <= k; i++){
    fac *= i;
  }
  return(fac);  
}

long double binom(int n, int k){
  if(k>n){
    return(0.0);
  }
  return(factorial(n)/factorial(n-k)/factorial(k));
}

void fitness(double nn){
  int i, j;
  j = ceil((double)L/10);
  for(i = 0; i <= L; i++){
    if(i<=j){
      sigma[i] = ((sigma_master-1) * pow((1-pow(((double)i/(double)j), nn)),(1/nn))) + 1;
      if(sigma[i] < 1.0){
      sigma[i] = 1.0;
      } 
    }else{
      sigma[i] = 1.0;
    }
  }
}

//CASE 1: original, including back mutations
//calculating transition probabilities between Hamming classes  
void init_p(){
  int n, k, i;
  for(n = 0; n <= L; n++){
    for(k = 0; k <= L; k++){
      p[n][k] = 0.0;
      if(k >= n){
        for(i = 0; i <= n; i++)
        p[n][k] += binom(n,i)*pow(mu/3.0,i)*pow(1.0-(mu/3.0),n-i) * binom(L-n,k-n+i)*pow(mu,k-n+i)*pow(1.0-mu,L-k-i);
      }else{
        for(i = (n - k); i <= n; i++)
        p[n][k] += binom(n,i)*pow(mu/3.0,i)*pow(1.0-(mu/3.0),n-i) * binom(L-n,k-n+i)*pow(mu,k-n+i)*pow(1.0-mu,L-k-i);
      }
    }
  }
} 

//CASE 2: modified, excluding back mutations
//calculating transition probabilities between Hamming classes 
/*void init_p(){
  int n, k, i;
  for(n = 0; n <= L; n++){
    for(k = 0; k <= L; k++){
      p[n][k] = 0.0;
      if(k >= n){
        p[n][k] = pow(1.0-(mu/3.0),n) * binom(L-n,k-n) * pow(mu,k-n) * pow(1.0-mu,L-k);
      }
      if (n>k){
        p[n][k] = 0.0;
      }
    }
  }
}*/

static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data){
    int i,j;
    double Phi = 0.0; //outflow
    double psum = 0.0;

    for(i=0;i<N;i++){
      psum = 0.0;

      for(j=0;j<N;j++){
        psum += p[i][j];
      }

      Phi += sigma[i] * NV_Ith_S(y,i) * psum;
    }

    for(i=0;i<N;i++){ 

      NV_Ith_S(ydot,i) = 0;

      for(j=0;j<N;j++){
      NV_Ith_S(ydot,i) += sigma[j] * p[j][i] * NV_Ith_S(y,j);
      }
      
      NV_Ith_S(ydot,i) -= NV_Ith_S(y,i) * Phi;
    }

    return 0;
}

int main(void){
  double t0 = 0.0;
  double tmax = 10000.0;
  double dt_out = 0.01;
  double reltol = 1e-14;
  double abstol = 1e-16;
  double analx;
  
  int i, nmu, negind;  
  
  for(nmu = 0; nmu <=10000; nmu++){
    mu=0.0001*nmu;
  
    fitness(n_scal);
    N = L + 1; //number of equations
    init_p();
    
    SUNContext sunctx;
    if (SUNContext_Create(0, &sunctx)){
      fprintf(stderr, "Error: SUNContext_Create failed\n");
      return 1;
    }
  
    N_Vector y = N_VNew_Serial(N, sunctx);
    if (y == NULL){
      fprintf(stderr, "Error: could not allocate N_Vector y\n");
      SUNContext_Free(&sunctx);
      return 1;
    }
  
    NV_Ith_S(y,0) = 1.0; //initially, only masters are present
    for (int i = 1; i < N; i++){
      NV_Ith_S(y,i) = 0.0;
    }
    
    void *cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (cvode_mem == NULL){
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
    if(LS == NULL){
      fprintf(stderr, "Error: could not create dense linear solver\n");
      CVodeFree(&cvode_mem);
      N_VDestroy(y);
      SUNMatDestroy(A);
      SUNContext_Free(&sunctx);
      return 1;
    }
  
    int flag_ls = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (flag_ls != CV_SUCCESS){
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
  
    while (t < tmax){
      flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      if (flag < 0){
        fprintf(stderr, "Integration failed, flag = %d\n", flag);
        break;
      }
  
      negind = 0;
      for(i = 0; i < N; i++){
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
  
    printf("%e", mu);
    for(i = 0; i <= L; i++){
      printf("\t%e", NV_Ith_S(y,i));
      }
    printf("\n");
    
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(y);
    SUNContext_Free(&sunctx);
  }
  
  return 0;
}
