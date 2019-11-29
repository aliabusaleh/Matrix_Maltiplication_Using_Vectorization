

#include <immintrin.h>  // portable to all x86 compilers
#include <stdio.h>
#include <time.h>

#define DATA float


#include <iostream>
using namespace std;
/* The size of 2D array or vector */
const int SIZE = 40;
/*
1- Array_      is the one will use to do matrix-multiplication and first one on matrix-matrix multiplication
2- Array_B     will use as the second Matrix for matrix-matrix multiplication
3- Array_out   will story the result for our array-matrix multiplication in normal way
4- Array_out_B will story the result for our array-matrix multiplication using SSE
4- Vector      will use to do the matrix-vector multiplication
5- Vector_out  will use to store the matrix-vector multiplication result
6- A           will use as the first Matrix for matrix-matrix multiplication (recursive)
6- B           will use as the second Matrix for matrix-matrix multiplication (recursive)
*/
DATA __attribute__((aligned(16))) Array_[SIZE][SIZE] ;
DATA __attribute__((aligned(16))) Array_B[SIZE][SIZE] ;
DATA __attribute__((aligned(16))) Array_out[SIZE][SIZE] ;
DATA __attribute__((aligned(16))) Array_out_B[SIZE][SIZE] ;
DATA __attribute__((aligned(16))) Vector_[SIZE] ;
DATA __attribute__((aligned(16))) Vector_out[SIZE] ;
int A[SIZE][SIZE] ;
int B[SIZE][SIZE] ;

// this function will print the matrix you give it.
void PrintMatrix(int *A, int MSIZE)
{
    printf("[ ");
    for (int i = 0; i < MSIZE; i++)
    {

        for (int j = 0; j < MSIZE; j++)
        {
            printf("%d ", A[i*MSIZE + j]);
        }

    }
    printf("]\n");
}
// first quarter for matrix multiplication recursively
void ProduceQuarter_00(int *A, int MSIZE, int *quarter)
{
    int n = MSIZE * 2;
    int m = 0;
    for (int i = 0; i < MSIZE; i++)
        for (int j = 0; j < MSIZE; j++)
            quarter[m++] = A[i*n + j];
}
// second quarter for matrix multiplication recursively
void ProduceQuarter_01(int *A, int MSIZE, int *quarter)
{
    int n = MSIZE * 2;
    int m = 0;
    for (int i = 0; i < MSIZE; i++)
        for (int j = 0; j < MSIZE; j++)
            quarter[m++] = A[i*n + (j + MSIZE)];
}
// third quarter for matrix multiplication recursively
void ProduceQuarter_10(int *A, int MSIZE, int *quarter)
{
    int n = MSIZE * 2;
    int m = 0;
    for (int i = 0; i < MSIZE; i++)
        for (int j = 0; j < MSIZE; j++)
            quarter[m++] = A[(i + MSIZE)*n + j];
}
// fourth quarter for matrix multiplication recursively
void ProduceQuarter_11(int *A, int MSIZE, int *quarter)
{
    int n = MSIZE * 2;
    int m = 0;
    for (int i = 0; i < MSIZE; i++)
        for (int j = 0; j < MSIZE; j++)
            quarter[m++] = A[(i + MSIZE)*n + (j + MSIZE)];
}
// this function will combine the Quarters
void CombineQuaters(int *quarter_00, int *quarter_01, int *quarter_10, int *quarter_11, int *C, int MSIZE)
{
    int n = MSIZE * 2;
    int m = 0;
    for (int i = 0; i < MSIZE; i++)
        for (int j = 0; j < MSIZE; j++)
            C[i*n + j] = quarter_00[m++];

    m = 0;
    for (int i = 0; i < MSIZE; i++)
        for (int j = 0; j < MSIZE; j++)
            C[i*n + (j + MSIZE)] = quarter_01[m++];

    m = 0;
    for (int i = 0; i < MSIZE; i++)
        for (int j = 0; j < MSIZE; j++)
            C[(i + MSIZE)*n + j] = quarter_10[m++];

    m = 0;
    for (int i = 0; i < MSIZE; i++)
        for (int j = 0; j < MSIZE; j++)
            C[(i + MSIZE)*n + (j + MSIZE)] = quarter_11[m++];

}
// add matricies
void AddMatrix(int *A, int *B, int *C, int MSIZE)
{
    for (int i = 0; i < MSIZE; i++)
        for (int j = 0; j < MSIZE; j++)
            C[i*MSIZE + j] = A[i*MSIZE + j] + B[i*MSIZE + j];

}
// for recursive mmultiplication
int *MatrixMulti(int *A, int *B, int size) {

    int *mat;
    mat = new int[size*size];
    int *Quad00_A, *Quad01_A, *Quad10_A, *Quad11_A;
    int *Quad00_B, *Quad01_B, *Quad10_B, *Quad11_B;
    int *mat1, *mat2, *mat3, *mat4;

    int msize=(size / 2);

    Quad00_A=new int[msize*msize];
    Quad01_A=new int[msize*msize];
    Quad10_A=new int[msize*msize];
    Quad11_A=new int[msize*msize];

    Quad00_B=new int[msize*msize];
    Quad01_B=new int[msize*msize];
    Quad10_B=new int[msize*msize];
    Quad11_B=new int[msize*msize];

    mat1=new int[msize*msize];
    mat2=new int[msize*msize];
    mat3=new int[msize*msize];
    mat4=new int[msize*msize];



    if (size == 1) {

        mat[0] = A[0] * B[0];


    }
    else{
    /*make quarter*/
    ProduceQuarter_00(A, size/2, Quad00_A);
    ProduceQuarter_00(B, size/2, Quad00_B);

    ProduceQuarter_01(A, size/2, Quad01_A);
    ProduceQuarter_01(B, size/2, Quad01_B);

    ProduceQuarter_10(A, size/2, Quad10_A);
    ProduceQuarter_10(B, size/2, Quad10_B);

    ProduceQuarter_11(A, size/2, Quad11_A);
    ProduceQuarter_11(B, size/2, Quad11_B);


    ProduceQuarter_10(A, size/2, mat1);
    ProduceQuarter_10(B, size/2, mat2);

    ProduceQuarter_11(A, size/2, mat3);
    ProduceQuarter_11(B, size/2, mat4);

    /*Recursive calling*/
    AddMatrix(MatrixMulti(Quad00_A, Quad00_B, size / 2), MatrixMulti(Quad01_A, Quad10_B, size / 2), mat1, size/2);
    AddMatrix(MatrixMulti(Quad00_A, Quad01_B, size / 2), MatrixMulti(Quad01_A, Quad11_B, size / 2), mat2, size/2);
    AddMatrix(MatrixMulti(Quad10_A, Quad00_B, size / 2), MatrixMulti(Quad11_A, Quad10_B, size / 2), mat3, size/2);
    AddMatrix(MatrixMulti(Quad10_A, Quad01_B, size / 2), MatrixMulti(Quad11_A, Quad11_B, size / 2), mat4, size/2);
    //Combine Quaters
    CombineQuaters(mat1, mat2, mat3, mat4, mat, size/2);

    }
    return mat;


}
// Return the current time in seconds
double seconds()
{

  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return now.tv_sec + now.tv_nsec / 1000000000.0;
}
// initilize the Arrays that will be used
void initialize_array(DATA a[][SIZE])
{
        for (int i = 0 ;  i < SIZE ; i++)
        for(int j=0; j<SIZE; j++)
        {
                a[i][j] = i+j;
            }

}
/* initilize the vector that will be used */
void initialize_vector(DATA V[SIZE])
{

        for(int j=0; j<SIZE; j++)
        {
                V[j] = j;
        }

}
/* matrix- vector multiplication normal way  */
void mat_vec_printf(DATA A[][SIZE],DATA V[SIZE], DATA Vout[SIZE])
{
        for (int i = 0 ; i < SIZE ; i++)
        {
                for (int j = 0 ; j < SIZE ; j++)
                {
                        Vout[i] += A[i][j]*V[j];
                }

        }

}
/* Matrix- Vector multiplication using SSE*/
void vec_sse(DATA m1[][SIZE], DATA m2[SIZE], DATA mout[SIZE])
{

    DATA prod = 0;

    __m128 X, Y, Z;



    for(int i=0; i< SIZE; i=i+1){
        Z[0] = Z[1] = Z[2] = Z[3] = 0;
        prod=0;
        for(int j=0; j<SIZE; j=j+4){
            X = _mm_load_ps(&m1[i][j]);
            Y = _mm_load_ps(&m2[j]);
            X = _mm_mul_ps(X, Y);
            Z = _mm_add_ps(X, Z);
        }

        for(int i_in=0; i_in<4; i_in++)
        {
        prod += Z[i_in];
        }
        mout[i] = prod;

    }

    return ;

}
/*Matrix-Matrix multiplication in a normal way  */
void mat_normal(DATA m1[][SIZE], DATA m2[][SIZE], DATA mout[][SIZE]) {
    for(int i=0; i<SIZE; i++)
        for(int k=0; k<SIZE; k++)
        for(int j=0; j<SIZE; j++){
        mout[i][j] += m1[i][k] * m2[k][j];
    }

}
/* matrix-matrix multiplication  using SSE (SIMD)*/
void mat_sse(DATA m1[][SIZE], DATA m2[][SIZE], DATA mout[][SIZE])
{

    for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j += 4) {
        __m128 sum = _mm_setzero_ps();
        for (int k = 0; k < SIZE; k++) {
            __m128 entry = _mm_set1_ps(m1[i][k]);
            __m128 row  = _mm_load_ps(&m2[k][j]);
            sum = _mm_add_ps(sum, _mm_mul_ps(entry, row));
        }
        _mm_store_ps(&mout[i][j], sum);
    }
}
    return ;

}


int main()
{


         /* This before and  after will be used to calculate time of execution for each method   */
    double before,after;
        /* Initilize the Arrays we will use  */
    initialize_array(Array_);
    initialize_array(Array_B);
    initialize_vector(Vector_);

    /**********************************************************************************************************************/
                                      /* Matrix - Vector Multiplication using normal method  */
    /**********************************************************************************************************************/
    printf("----------------------------------- matrix * Vector normal ------------------------------------\n");
        before = seconds();
    mat_vec_printf(Array_,Vector_,Vector_out);
        after = seconds();
        printf(" Time:%f for matrix-Vector multiplication\n",after-before);
        for(int j=0; j<SIZE; j++){
            printf("%f\n",Vector_out[j] );
        }

     /**********************************************************************************************************************/
                                      /* Matrix - Vector Multiplication using SSE method  */
    /**********************************************************************************************************************/
    printf("----------------------------------- matrix * Vector SSE ------------------------------------\n");
    before = seconds();
    vec_sse(Array_,Vector_,Vector_out);
        after = seconds();
        printf("Time:%f for Matrix-Vector SSE \n",after-before);
        //  for(int j=0; j<SIZE; j++){
        //     printf("%f\n",Vector_out[j] );
        // }

     /**********************************************************************************************************************/
                                      /* Matrix - Matrix Multiplication using normal method  */
    /**********************************************************************************************************************/
     printf("----------------------------------- matrix * Matix normal ------------------------------------\n");
        before = seconds();
     mat_normal(Array_, Array_B, Array_out);
    after = seconds();
    printf("Time:%f for Matrix-Matrix Normal \n",after-before);
      /*       Uncomment the folowing loop to print the output array    */

    // for(int j=0; j< SIZE; j++){
       //   for(int i=0; i< SIZE; i++)
    // {
      //   printf("%f ",Array_out[j][i]);
     //}
      //  printf("\n");
   // }

     /**********************************************************************************************************************/
                                      /* Matrix - Matrix Multiplication using SSE method  */
    /**********************************************************************************************************************/

    printf("----------------------------------- matrix * Matix SSE ------------------------------------\n");
        before = seconds();
    mat_sse(Array_, Array_B, Array_out_B);
    after = seconds();
    printf("Time:%f for Matrix-Matrix SSE \n",after-before);

      /*       Uncomment the folowing loop to print the output array    */
    // for(int j=0; j< SIZE; j++){
      //    for(int i=0; i< SIZE; i++)
    // {
      //   printf("%f ",Array_out_B[j][i]);
    // }
      //  printf("\n");
   // }

     /**********************************************************************************************************************/
                                      /* Matrix - Matrix Multiplication using SSE method  */
    /**********************************************************************************************************************/
    printf("----------------------------------- matrix * Matix Recursive ------------------------------------\n");
        before = seconds();
     int *Arr1 = new int[(SIZE*SIZE)];
    //Array number 2
    int *Arr2 = new int[(SIZE*SIZE)];
    // reading the matrix's 1  elements
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            Arr1[i*SIZE + j] = i+j;
        }
    }
    // reading the matrix's 2  elements
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
           Arr2[i*SIZE + j] =i+j;
        }
    }
    int *x;

     x =MatrixMulti(Arr1, Arr2, SIZE);
    after = seconds();
    printf("Time:%f for Matrix-Matrix Recursive \n",after-before);

   //

    return 0;
}
