#define PI 3.141592653589793238462643383279502884L //From the bible, I mean  Knuth's "The art of Computer Programming"
#include <math.h>
#include <stdlib.h>
#include <jdh_math.h>
#include <stdio.h>
#include <string.h>
/**************************************************
***************************************************
***************** FFTs and DFTs *******************
***************************************************
**************************************************/


//1D FFT
void jdh_math_FFT (double* data, int N, int isign){
  int n, n2, nb, j, k, i0, i1;
  double wr, wi, wrk, wik;
  double d, dr, di, d0r, d0i, d1r, d1i;
  double *cp;

  j = 0;
  n2 = N / 2;
  for( k = 0; k < N; ++k ){
   	if( k < j ){
      		i0 = k << 1;
      		i1 = j << 1;
      		dr = data[i0];
      		di = data[i0+1];
      		data[i0] = data[i1];
      		data[i0+1] = data[i1+1];
      		data[i1] = dr;
      		data[i1+1] = di;
	}
    	n = N >> 1;
    	while( (n >= 2) && (j >= n) ){
    		j -= n;
	  	n = n >> 1;
  		}
	j += n;
	}

	for( n = 2; n <= N; n = n << 1 ){
    		wr = cos( (double)(2.0 * PI / (double) n ) );
    		wi = sin( (double)(2.0 * PI / (double) n ) );
    		if( isign == 1 ) wi = -wi;
    		cp = data;
    		nb = N / n;
    		n2 = n >> 1;
		for( j = 0; j < nb; ++j ){
			wrk = 1.0;
			wik = 0.0;
      			for( k = 0; k < n2; ++k ){
    				i0 = k << 1;
  				i1 = i0 + n;
				d0r = cp[i0];
				d0i = cp[i0+1];
				d1r = cp[i1];
     				d1i = cp[i1+1];
   				dr = wrk * d1r - wik * d1i;
    				di = wrk * d1i + wik * d1r;
       				cp[i0] = d0r + dr;
        			cp[i0+1] = d0i + di;
        			cp[i1] = d0r - dr;
        			cp[i1+1] = d0i - di;
        			d = wrk;
        			wrk = wr * wrk - wi * wik;
				wik = wr * wik + wi * d;
			}
			cp += n << 1;
		}
	}
	if(isign==-1){
		for(int i = 0; i<N; i++){
			data[2*i] /= N;
			data[2*i+1] /= N;
		}
	}
}

//Derivative Operation on Fourier Coefficients
void jdh_math_Fourier_Derivative(double* data, int n, double L){
	double tmp;
	double k;
	for(int i = 0; i<n; i++){
		tmp = data[2*i];
		if(i<n/2+1){
			k = 2.00L*PI*(double)i/L;
		} else {
			k = 2.00L*PI*((double)(i-n))/L;
		}
		if(i!=0 && i!=(n/2)){
			data[2*i] = -k*data[2*i+1];
			data[2*i+1] = k*tmp;
		} else { //Fix oddball and zero wavenumber
			data[2*i] = 0.00L;
			data[2*i+1] = 0.00L;
		}
	}
}

//2D FFT
void jdh_math_FFT2D(double** data, int n, int m, int dir){
	int i,j;
	double* ax = malloc(2*n*sizeof(double));
	double* ay = malloc(2*m*sizeof(double));
	//FFT in n
	for(j = 0; j<m; j++){
		for(i = 0; i<n; i++){
			ax[2*i] = data[i][2*j];
			ax[2*i+1] = data[i][2*j+1];
		}
		jdh_math_FFT(ax, n, dir);	
		for(i = 0; i<n; i++){
			data[i][2*j] = ax[2*i];
			data[i][2*j+1] = ax[2*i+1];
		}
	}
	//FFT in m	
	for(i = 0; i<n; i++){
		for(j = 0; j<m; j++){
			ay[2*j] = data[i][2*j];
			ay[2*j+1] = data[i][2*j+1];
		}
		jdh_math_FFT(ay, m, dir);	
		for(j = 0; j<m; j++){
			data[i][2*j] = ay[2*j];
			data[i][2*j+1] = ay[2*j+1];
		}
	}
	free(ax);
	free(ay);
}

//2D Fourier Derivative
void jdh_math_Fourier_Derivative_2D(double** data, int nx, int ny, int dir, double L){
	// derivative in first index
	if(dir == 0){
		double* temp1 = malloc(2*nx*sizeof(double));
		for(int j = 0; j<ny; j++){
			for(int i = 0; i<nx; i++){
				temp1[2*i] = data[i][2*j];
				temp1[2*i+1] = data[i][2*j+1];
			}
			jdh_math_Fourier_Derivative(temp1, nx, L);
			for(int i = 0; i<nx; i++){
				data[i][2*j] = temp1[2*i];
				data[i][2*j+1] = temp1[2*i+1];
			}
		}
	}
	//derivative in second index
	if(dir == 1){
		for(int i = 0; i<nx; i++){
			jdh_math_Fourier_Derivative(data[i],ny,L);
		}
	}
}

//Fast Sine Transform
//TODO: Fast Sine Transform

//Fast Cosine Transform
//TODO: Fast Cosine Transform

void jdh_math_Zero_Pad_FFT_2D(double** ndata, double** data, int nx, int ny, int nxn, int nyn){
	for(int i = 0; i<nxn; i++){
		if(i<(nx/2)){	
			jdh_math_Zero_Pad_FFT(ndata[i], data[i], ny, nyn);
		} else if (i<(nxn-nx/2)){
			for(int j = 0; j<nyn; j++){
				ndata[i][2*j] = 0.00L;
				ndata[i][2*j+1] = 0.00L;
			}
		} else {
			jdh_math_Zero_Pad_FFT(ndata[i], data[(nx-(nxn-i))], ny, nyn);
		}	
	}

}

	
void jdh_math_Zero_Pad_FFT(double* ndata, double* data, int n, int nn){
	for(int j = 0; j<nn; j++){
		if(j<(n/2)){	
			ndata[2*j] = data[2*j];
			ndata[2*j+1] = data[2*j+1];
		} else if (j<(nn-n/2)){
			ndata[2*j] = 0.00L;
			ndata[2*j+1] = 0.00L;
		} else {
			ndata[2*j] = data[2*(n-(nn-j))];
			ndata[2*j+1] = data[2*(n-(nn-j))+1];
		}
	}
}

void jdh_math_Shrink_FFT_2D(double** ndata, double** data, int nx, int ny, int nxn, int nyn){
	for(int i = 0; i<nxn; i++){
		if(i<(nx/2)){	
			jdh_math_Shrink_FFT(ndata[i], data[i], ny, nyn);
		
		} else if (i<(nxn-nx/2)){
		} else {
			jdh_math_Shrink_FFT(ndata[(nx-(nxn-i))], data[i], ny, nyn);
		}	
	}

}

void jdh_math_Shrink_FFT(double* ndata, double* data, int n, int nn){
	double L = 2.00L* (double) nn / (double) n;
	for(int i = 0; i<nn; i++){
		if(i<(n/2-1)){	
			ndata[2*i] = L*data[2*i];
			ndata[2*i+1] = L*data[2*i+1];
		} else if (i<(nn-n/2)){
		} else {
			ndata[2*(n-(nn-i))] = L*data[2*i];
			ndata[2*(n-(nn-i))+1] = L*data[2*i+1];
		}
	}
}




/************************************************
*************************************************
************** Array Operations *****************
*************************************************
************************************************/
//2D Complex Array Copy
void jdh_math_Array_Copy_2D_Complex(double** target, double** source, int nx, int ny){
	size_t data_size = sizeof(double)*2*ny;
	for(int i = 0; i<nx; i++){
		memcpy(target[i], source[i], data_size);
	}
} 

//2D Array Copy
void jdh_math_Array_Copy_2D(double** target, double** source, int nx, int ny){
	size_t data_size = sizeof(double)*ny;
	for(int i = 0; i<nx; i++){
		memcpy(target[i], source[i], data_size);
	}
}

//2D Array Sum
void jdh_math_Array_Sum_2D(double** ar1, double** ar2, double** ar3, int nx, int ny){
	for(int i = 0; i < nx; i++){
		for(int j = 0; j<ny; j++){
			ar3[i][j] = ar1[i][j]+ar2[i][j];	
		}
	}
}

//2D Array Sum Complex
void jdh_math_Array_Sum_2D_Complex(double** ar1, double** ar2, double** ar3, int nx, int ny){
	for(int i = 0; i < nx; i++){
		for(int j = 0; j<ny; j++){
			ar3[i][2*j] = ar1[i][2*j]+ar2[i][2*j];	
			ar3[i][2*j+1] = ar1[i][2*j+1]+ar2[i][2*j+1];	
		}
	}
}

//2D Array Product
void jdh_math_Array_Product_2D(double** ar1, double** ar2, double** ar3, int nx, int ny){	
	for(int i = 0; i < nx; i++){
		for(int j = 0; j<ny; j++){
			ar3[i][j] = ar1[i][j]*ar2[i][j];	
		}
	}
}

//2D Array Product Complex
void jdh_math_Array_Product_2D_Complex(double** ar1, double** ar2, double** ar3, int nx, int ny){
	double t1;
	double t2;
	for(int i = 0; i < nx; i++){
		for(int j = 0; j<ny; j++){
			t1 = ar1[i][2*j];
			t2 = ar2[i][2*j];
			ar3[i][2*j] = ar1[i][2*j]*ar2[i][2*j]-ar1[i][2*j+1]*ar2[i][2*j+1];
			ar3[i][2*j+1] = t1*ar2[i][2*j+1]+t2*ar1[i][2*j+1];
		}
	}
}

//2D Array Product No Complex
void jdh_math_Array_Product_2D_No_Complex(double** ar1, double** ar2, double** ar3, int nx, int ny){
	double t1;
	double t2;
	for(int i = 0; i < nx; i++){
		for(int j = 0; j<ny; j++){
			t1 = ar1[i][2*j];
			t2 = ar2[i][2*j];
			ar3[i][2*j] = ar1[i][2*j]*ar2[i][2*j];
			ar3[i][2*j+1] = 0.00L;
		}
	}
}

//2D Array Complex to Zero
void jdh_math_Array_Set_Complex_Zero_2D(double** ar, int nx, int ny){
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			ar[i][2*j+1] = 0.00L;
		}	
	}
}


/************************************************
*************************************************
**************** Chebyshev **********************
*************************************************
************************************************/


//Discrete Chebyshev Transform
void jdh_math_DCT(double* fx, double* am, int nx){
	double ci,cj;
	for(int i = 0; i<nx; i++){
		if(i==0 || i==(nx-1)){
			ci = 2.00;
		} else {
			ci = 1.00;
		}
		am[i] = 0.00;
		for(int j = 0; j<nx; j++){
			if(j == 0 || j == (nx-1)){
				cj = 2.00;
			} else {
				cj = 1.00;
			}
				am[i] += 2.00/((nx-1)*ci*cj)*fx[j]*cos((double)(i*j*PI/(nx-1)));
		}
	}
}

//Inverse Discrete Chebyshev Transform
void jdh_math_IDCT(double* fx, double* am, int nx){
	for(int i = 0; i<nx; i++){
		fx[i] = 0.00;
		for(int j = 0; j<nx; j++){
				fx[i] += am[j]*cos((double)(i*j*PI/(nx-1)));
		}
	}
}

//Derivative in Chebyshev Space
void jdh_math_Chebyshev_Derivative(double* am, double* bm, int nx){
	double ck;
	bm[nx] = 0.00;
	bm[nx-1] = 0.00;
	for(int k = nx-1; k>0; k--){
		if(k == 1){
			ck = 2.00;
		} else {
			ck = 1.00;
		}
		bm[k-1] = (2.00*(k)*am[k]+bm[k+1])/ck;
	}
} 

/************************************************
*************************************************
********* Domain Distributions ******************
*************************************************
************************************************/

//Domain Linear Distribution
double* jdh_math_LinDist(double xmin, double xmax, int nx){
	double* x = malloc(nx*sizeof(double));
	for(int i=0; i<nx; i++){
		x[i] = xmin+(xmax-xmin)/((double)(nx-1))*(double)i;	
	}
	return x;
}

//Domain Cosine Distribution
double* jdh_math_CosDist(double xmin, double xmax, int nx){
	double* x = malloc(nx*sizeof(double));
	for(int i = 0; i<nx; i++){
		x[i] = (cos(PI*(double)i/((double)(nx-1)))+1.00)/2.00*(xmin-xmax)+xmax;	
	}
	return x;
}
	
//Domain Chebyshev Distribution
double* jdh_math_ChbshvDist(double xmin, double xmax, int nx){
	return jdh_math_CosDist(xmax, xmin, nx);
}
/************************************************
*************************************************
************* Windowing *************************
*************************************************
************************************************/

//Apply hanning window to data
void jdh_math_Hanning(double* data, int n) {
	for(int i = 0; i < n; i++){
		data[i]= data[i]*0.5*(1-cos(2*PI*(double)i/(double)(n-1)));
	}
}

//Zero pad array of length n to new length nn
void jdh_math_Zero_Pad(double* data, double* ndat, int n, int nn){
	for(int i = 0; i<nn; i++){
		if (i < n){
			ndat[i] = data[i];
		} else {
			ndat[i] = 0.00;
		}
	}	
}

//Zero pad complex array of length 2*n to new length 2*nn
void jdh_math_Zero_Pad_Complex(double* data, double* ndat, int n, int nn){		
	for(int i = 0; i<nn; i++){
		if (i < n){
			ndat[2*i] = data[2*i];
			ndat[2*i+1] = data[2*i+1];
		} else {
			ndat[2*i] = 0.00;
			ndat[2*i+1] = 0.00;
		}
	}	
}

/************************************************
*************************************************
************* File Operations *******************
*************************************************
************************************************/


//Write 1D array to file
void jdh_math_fwrite(char* fname, double* x, int nx){
	FILE *fref = fopen(fname,"w");
	size_t count = fwrite(x,sizeof(double),nx,fref);
//	fclose(fref);	
}

//Write 2D array to file
void jdh_math_fwrite_2D(char* fname, double** x, int nx, int ny){
	FILE *fref = fopen(fname,"w");
	for(int i = 0; i<nx; i++){
		size_t count = fwrite(x[i],sizeof(double),ny,fref);
	}	
	fclose(fref);	
}

//Write 2D array to file
void jdh_math_fwrite_2D_Complex(char* fname, double** x, int nx, int ny){
	FILE *fref = fopen(fname,"w");
	for(int i = 0; i<nx; i++){
		size_t count = fwrite(x[i],2*sizeof(double),ny,fref);
	}	
	fclose(fref);	
}

