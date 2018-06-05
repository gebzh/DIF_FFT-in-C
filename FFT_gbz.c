/*=========================*function description*==============================
-the input file named  fft_in.txt stores  real numbers the input signal with 
the length no more than 8192.
-the output file named fft_out.txt stores the modulus of FFT of the signal 
padded with zeros that makes length of data to be power of 2.
=============================================================================*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#define PI 3.1415926535898

int N,v,j;
double x[8197][2][2], wn[8197][2];

void input();
void fft();
void output();
int bitReversal(int i);

int main(){
  input();
  fft();
  output();
}

void fft(){
  //DIF-FFT
  int st,j1,i,i1,t,m;
  double tmp,ac,bd,apbmcpd;
  //calculation of Wn
  for (i=0;i<N>>1;++i){
  	tmp = 2.0 * PI * i / N;
    wn[i][0] = cos(tmp); //real part of Wn
    wn[i][1] = -sin(tmp); //imag part of Wn
    //printf("%.4f + %.4f i\n", wn[i][0], wn[i][1]);
  } 
  //calculation of DFT i.e. sigma(i:0->N-1)[x(n)*(Wn^i*k)] 
  m = 1;
  for (st=N>>1;st>0;st>>=1){
    j1 = j; j = j^1; 
    for (i=0;i<N;i++){
   	  t = (st & i) > 0;
   	  i1 = i + st - t * st * 2;
      x[i][j][0] = x[i1][j1][0]+(1-2*t)*x[i][j1][0];
      x[i][j][1] = x[i1][j1][1]+(1-2*t)*x[i][j1][1];
      if (t) { //complex multiplication
        i1 = (i % st) * m;
      	//c = x[i][j][0]; d = x[i][j][1];
        //a = wn[i1][0]; b = wn[i1][1];
        ac = wn[i1][0] * x[i][j][0];
        bd = wn[i1][1] * x[i][j][1];
        apbmcpd = (wn[i1][0]+wn[i1][1])*(x[i][j][0]+x[i][j][1]);
      	x[i][j][0] = ac - bd;
		    x[i][j][1] = apbmcpd - ac - bd; 
      }
      //printf("[i,j]=%d %d x[i,j]=%.4f + %.4f i\n",i,j,x[i][j][0],x[i][j][1]);
    }
    m <<= 1;
  }
}

void input(){
  FILE *fp;
  N = 0;
  fp=fopen("fft_in.txt","r"); 
  if (fp == NULL) {
    printf("File cannot be opened./n");
    exit(1);		
  }
  while (!feof(fp)){
    fscanf(fp,"%lf", &x[N][0][0]);
    //fscanf(fp,"%lf", &x[N][0][1]);
    N++;
  }
  fclose(fp);
  if (N>1024) {
    printf("The number of data is out of range.\n");
    exit(1);
  } 
  v = (int)ceil(log(1.0*N)/log(2.0));
  N = pow( 2, v ); //makes the length to be power of 2
}

void output(){
  int i,k;
  double mag;
  FILE *fp;
  fp = fopen("fft_out.txt","w");
  for (i=0;i<N;i++){
  	k = bitReversal(i);
  	mag = sqrt( pow(x[k][j][0], 2) + pow(x[k][j][1], 2) );  
  	//printf("%.4f + %.4f i\n", x[k][j][0], x[k][j][1]);
    fprintf(fp,"%.4f ",mag);
  }
  fclose(fp);
}

int bitReversal(int i){
  int t = 0, n = 0;
  for (; i; i >>= 1){
  	n <<= 1;
  	n |= i & 1;
  	++t;
  }
  n = n << (v-t);
  return n;
}
