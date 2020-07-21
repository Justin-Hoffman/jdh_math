//Math Library - Justin Hoffman

//FFT Routine
extern void jdh_math_FFT(double* c, int nn, int isign);

//Derivative of Fourier Array
extern void jdh_math_Fourier_Derivative(double* data, int n, double L);

// 2D FFT Routine
extern void jdh_math_FFT2D(double** data, int n, int m, int dir);

//2D Fourier Derivative
extern void jdh_math_Fourier_Derivative_2D(double** data, int nx, int ny, int dir, double L);

//2D FFT Zero Pad
extern void jdh_math_Zero_Pad_FFT_2D(double** ndata, double** data, int nx, int ny, int nxn, int nyn);

//FFT Zero Pad
extern void jdh_math_Zero_Pad_FFT(double* ndata, double* data, int n, int nn);

//Shrink Zero Padded FFT
void jdh_math_Shrink_FFT_2D(double** ndata, double** data, int nx, int ny, int nxn, int nyn);

//Shrink Zero Padded FFT
void jdh_math_Shrink_FFT(double* ndata, double* data, int n, int nn);

//2D Array Copy
extern void jdh_math_Array_Copy_2D(double** target, double** source, int nx, int ny);

//2D Array Copy Complex
extern void jdh_math_Array_Copy_2D_Complex(double** target, double** source, int nx, int ny);

//2D Array Sum
extern void jdh_math_Array_Sum_2D(double** ar1, double** ar2, double** ar3, int nx, int ny);

//2D Array Sum Complex
extern void jdh_math_Array_Sum_2D_Complex(double** ar1, double** ar2, double** ar3, int nx, int ny);

//2D Array Product
void jdh_math_Array_Product_2D(double** ar1, double** ar2, double** ar3, int nx, int ny);

//2D Array Product Complex
void jdh_math_Array_Product_2D_Complex(double** ar1, double** ar2, double** ar3, int nx, int ny);

//2D Array Product Complex
void jdh_math_Array_Product_2D_No_Complex(double** ar1, double** ar2, double** ar3, int nx, int ny);

//2D Array Complex to Zero
void jdh_math_Array_Set_Complex_Zero_2D(double** ar, int nx, int ny);

//Discrete Chebyshev Transform
extern void jdh_math_DCT(double* fx, double* am, int nx);

//Inverse Discrete Chebyshev Transform
extern void jdh_math_IDCT(double* fx, double* am, int nx);

//Derivative of Chebyshev Array
extern void jdh_math_Chebyshev_Derivative(double* am, double* bm, int nx);

//Linear Distribution
extern double* jdh_math_LinDist(double xmin, double xmax, int nx);

//Cosine Distribution
extern double* jdh_math_CosDist(double xmin, double xmax, int nx);

//Chebyshev Distribution
extern double* jdh_math_ChbshvDist(double xmin, double xmax, int nx);

//Hanning Window Operator
extern void jdh_math_Hanning(double* data, int n);

//Zero Pad Array
extern void jdh_math_Zero_Pad(double* data, double* ndat, int n, int nn);

//Zero Pad Complex array
extern void jdh_math_Zero_Pad_Complex(double* data, double* ndat, int n, int nn);

//Write array to file
extern void jdh_math_fwrite(char* fname, double* x, int nx);

//Write array to file
extern void jdh_math_fwrite_2D(char* fname, double** x, int nx, int ny);
//Write array to file
extern void jdh_math_fwrite_2D_Complex(char* fname, double** x, int nx, int ny);
