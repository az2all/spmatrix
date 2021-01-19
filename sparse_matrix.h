/*******************************************************************************
 SMLC++ -- Sparse Matrix Library in C++ Using Old Yale Format.

 This header contains prototypes of the functions for storing and manipulating
 sparse matrix using Old Yale format.
  
 The toolbox is a free and open-source software, which is distributed under
 the BSD 2-Clause License.
 Copyright (c) <2021>

 Author: 	Anton Zaicenco, anton.az@gmail.com
 Date: 	    January 18, 2021
 Version:   v1.0

 
 [1]. Classes:
 sparse_Mtrx					sparse matrix class (compressed row format)
 full_Vector					full vector class
 coord_vectors					sparse matrix class (coordinate format)
 
 [2]. Functions for sparse matrix:
 sparse_assign					assign one value to sp. matrix (replace or add)
 matrix_coord_vectors_assign	assign values directly to JA, IA, A
 sparse_convert_to_full_print	convert sparse matrix to full and print
 
 [3]. Functions that operate on 1D arrays:
 sort__ascend					sort the elements of 3 arrays into ascending
 								order using int array A1 then A2
 sort__ascend_one				Sorts the elements of 3 arrays into ascending
 								order using int array A1
 int_array_capacity_up			boost 1D int array capacity
 
 [4]. Functions for "coord_vectors" class:
 coord_vectors_assign			assign values to 3 vectors (coordinates of non-zero)
 								of "coord_vectors" object
 
 [5]. Functions for matrix operations:
 sparse_matrix_multiply			multiply 2 sparse matrices
 sparse_dot_product				dot product of 2 sparse matrices
 sparse_matrix_vect_multiply	multiply sparse matrix and vector (1D array)
 sparse_A_times_b				sparse matrix * vector (full_Vector class)
 sparse_matrix_sum				sum of 2 sparse matrices
 sparse_matrix_scale_out		sparse matrix * number

*/


#ifndef _sparse_matrix_h_included_
#define _sparse_matrix_h_included_

#include <iostream>
#include <fstream>
#include <new>

using namespace std;

template <typename W>
class sparse_Mtrx
{
    public:
           int* IA;
           int* JA;
           W*    A;
           int cia;  // capacity of IA (nr of zero elements at the end)
           int cja;  // capacity of JA (nr of zero elements at the end)
           
           // constructor:
           sparse_Mtrx(int el_, int rs_, int cs_)
           {
              int i, i2;
              i2 = rs_+1;
              if (el_<1) {el_=1;};    // nr of elements cannot be <1
              
              IA = new (nothrow) int [el_];   // nr of elements
              JA = new  int [el_];
              A  = new (nothrow) W   [el_];
              eL = el_;
              rs = rs_;
              cs = cs_;
              for(i=0;i<el_;i++)    // length of JA and A = nr of elements
              {
                  JA[i]=0;
                   A[i]=0;
              }
              
              for (i=0;i<el_;i++)  // length of IA = nr of elements
              {
                  IA[i]=0;
              }
              
              cia = el_;         // IA capacity
              sia = cia;           // IA size
              cja = el_;           // JA capacity
              sja = el_;           // JA and A size
           }
           
           // destructor:
           virtual ~sparse_Mtrx()
           {
           // printf("  -- matrix destructor \n");
              delete [] IA;
              delete [] JA;
              delete []  A;
           }
           
           // functions:
           // int   get_nnz();
           void   print_Mtrx();             // print matrix
           int    get_sia(){return sia;};
           int    get_sja(){return sja;};
           int    get_rs(){return rs;};
           int    get_cs(){return cs;};
           int    get_matrix_nnz(){return eL;};    // return nr of non-zero elements
           void   set_matrix_nnz(int nv){eL=nv;};  // return nr of non-zero elements
           void   set_sia(int nv){sia=nv;};
           void   set_sja(int nv){sja=nv;};
           void   ia_capacity_up(int inc);   // boost IA capacity
           void   ia_capacity_down(int inc); // down  IA capacity
           void   ja_capacity_up(int inc);   // boost JA and A capacity
           // get data:
           W      get_entry(int r, int c);               // get an entry
           void   get_matrix_row(sparse_Mtrx *G, int r); // get a row
           void   get_matrix_col(sparse_Mtrx *G, int c); // get a column
           void   get_matrix_block(sparse_Mtrx *G, int r1, int r2, int c1, int c2);
           // compress row functions:
           void   compress_row();
           int    matrix_sort_vectors_ascend();    // sort using 1st and 2nd vectors
           int    matrix_Qsort_vectors_ascend();   // same sort but using quicksort
    protected:
           int rs;   // nr of rows
           int cs;   // nr of columns
           int eL;   // total nr. of elements
           int sia;  // size of IA
           int sja;  // size of JA = size of A
};
// -----------------------------------------------------------------------------



// #############################################################################
// Full vector class. Used with sparse matrix class
template <typename W>
class full_Vector
{
    public:
           W* V;
           // constructor:
           full_Vector(int elm)
           {
              int i;
              if (elm<1) {elm=1;};    // nr of elements cannot be <1
              V  = new (nothrow) W [elm];
              for(i=0;i<elm;i++)     // length of A = nr of elements
              {
                   V[i] = 0;
              }
           }
           // destructor:
           virtual ~full_Vector()
           {
              delete []  V;
           }
};
// -----------------------------------------------------------------------------


// #############################################################################
// input for sparse matrix class. Coordinate format
template <typename W>
class coord_vectors
{
    public:
           int*     iA;   // row index
           int*     jA;   // column index
           W*       aA;   // non-zero entry
           // constructor:
           coord_vectors(int el_)
           {
              int i;
              if (el_<1) {el_=1;};    // nr of elements cannot be <1
              iA  = new int    [el_];
              jA  = new int    [el_];
              aA  = new W      [el_];
              eL  = 0;
              sia = el_;
              for(i=0;i<el_;i++)    // zero fill
              {
                  iA[i]=0;
                  jA[i]=0;
                  aA[i]=0;
              }
           }
           // destructor:
           virtual ~coord_vectors()
           {
           // printf("  -- vector destructor \n");
              delete [] iA;
              delete [] jA;
              delete [] aA;
           }
           
           // functions:
           void  print_coord_vectors();            // print 3 vectors
           int   get_vectors_sia(){return sia;};   // return size of vectors
           int   get_vectors_nnz(){return eL;};    // return nr of non-zero elements
           void  vect_capacity_up(int inc);        // boost iA, jA and aA capacity
           void  set_vectors_nnz(int nv){eL=nv;};  // update nr of non-zero elements
           void  set_vectors_sia(int nv){sia=nv;}; // update length of vectors
           int   sort_vectors_ascend();            // sort using 1st and 2nd vectors
    protected:
           int eL;                                 // total nr. of non-zero elements
           int sia;                                // size of vectors
};

// -----------------------------------------------------------------------------



// #############################################################################
// basic matrix functions
// assign one value to sparse matrix (replace or add):
template <typename W>
void  sparse_assign(sparse_Mtrx<W> *G, W a, int j, int k, int SoA);

// assign values directly to 3 vectors of a sparse matrix (JA, IA, A):
template <typename W>
void  matrix_coord_vectors_assign(sparse_Mtrx<W> *G, W a, int j, int k, int SoA);

// convert sparse matrix to full and print:
template <typename W>
void  sparse_convert_to_full_print(sparse_Mtrx<W> *G);


// #############################################################################
// functions that operate on arrays:
// Sorts the elements of 3 arrays into ascending order using int array A1 then A2:
int   sort__ascend(int A1[], int A2[], double A3[], int length_vect);

// Sorts the elements of 3 arrays into ascending order using int array A1:
int   sort__ascend_one(int A1[], int A2[], double A3[], int length_vect);

// assign values to 3 vectors (coordinates of non-zero):
template <typename W>
void  coord_vectors_assign(coord_vectors<W> *G, W a, int j, int k, int SoA);

// boost 1D int array capacity:
void  int_array_capacity_up(int A[], int length_vect, int N);


// #############################################################################
// matrix operations
// matrix * matrix:
template <typename W>
void  sparse_matrix_multiply(sparse_Mtrx<W> *Rez, sparse_Mtrx<W> *B1, sparse_Mtrx<W> *B2);

// matrix: dot product (can be any of: row*column, row*row, etc.):
template <typename W>
W  sparse_dot_product(sparse_Mtrx<W> *M1, sparse_Mtrx<W> *M2);

// matrix * vector (using 1D array):
template <typename W>
void  sparse_matrix_vect_multiply(W *Rez, sparse_Mtrx<W> *B1, W *b);

// matrix * vector (using full_Vector class):
template <typename W>
void  sparse_A_times_b(full_Vector<W> *Rez, sparse_Mtrx<W> *B1, full_Vector<W> *b);

// matrix + matrix:
template <typename W>
void  sparse_matrix_sum(sparse_Mtrx<W> *Rez, sparse_Mtrx<W> *B1, sparse_Mtrx<W> *B2);

// matrix * number (result saved to another matrix -- Rez):
template <typename W>
void  sparse_matrix_scale_out(sparse_Mtrx<W> *Rez, sparse_Mtrx<W> *B1, W a);
// -----------------------------------------------------------------------------

#endif

