/*******************************************************************************
 SMLC++ -- Sparse Matrix Library in C++ Using Old Yale Format.
 
 Author: 	Anton Zaicenco, anton.az@gmail.com
 Date: 	    January 18, 2021
 Version:   v1.0

 
 This header contains prototypes of the functions for solving a system A*x=b
 where A is stored in Old Yale sparse matrix format.
 
 sparseJacobi		Jacobi
 sparseGaussSeidel	Gauss-Seidel
 sparseSOR			SOR
 sparseCG			conjugate gradient
 
 The toolbox is a free and open-source software, which is distributed under
 the BSD 2-Clause License.
 Copyright (c) <2021>
 
*/

#ifndef _sparse_linearsys_h_included_
#define _sparse_linearsys_h_included_

#include <iostream>
#include <fstream>
#include "sparse_matrix.h"

using namespace std;


// solve the system of linear equations Ax=b -----------------------------------
// stationary iterative methods:

/* Jacobi method. Careful! A must be diagonally dominant.
 Inputs:
 	Rez			result x
 	B1			coefficient matrix A (sparse, symm. pos. def.)
 	b			known vecor b
 	pr			precision to stop iterations: abs(sum(error))<pr
 	max_iter	max number of iterations
*/
template <typename W>
void  sparseJacobi(full_Vector<W> *Rez, sparse_Mtrx<W> *B1, full_Vector<W> *b, double pr, int max_iter);

/* Gauss-Seidel method. Careful! A must be diagonally dominant.
 Inputs:
 	Rez			result: vector x
 	B1			coefficient matrix A (sparse, symm. pos. def.)
 	b			known vecor b
 	pr			precision to stop iterations: abs(sum(error))<pr
 	max_iter	max number of iterations
*/
template <typename W>
void  sparseGaussSeidel(full_Vector<W> *Rez, sparse_Mtrx<W> *B1, full_Vector<W> *b, double pr, int max_iter);

/* SOR - succesive overrelaxation method. Careful! A must be diagonally dominant.
 Inputs:
 	Rez			result: vector x
 	B1			coefficient matrix A (sparse, symm. pos. def.)
 	b			known vecor b
 	pr			precision to stop iterations: abs(sum(error))<pr
 	max_iter	max number of iterations
*/
template <typename W>
void  sparseSOR(full_Vector<W> *Rez, sparse_Mtrx<W> *B1, full_Vector<W> *b, double pr, int max_iter);

// Krylov subspace methods:

/* CG - conjugate gradient method.
 Inputs:
 	Rez			result: vector x
 	B1			coefficient matrix A (sparse, symm. pos. def.)
 	b			known vecor b
 	pr			precision to stop iterations: abs(sum(error))<pr
 	max_iter	max number of iterations
*/
template <typename W>
int  sparseCG(full_Vector<W> *Rez, sparse_Mtrx<W> *B1, full_Vector<W> *b, double pr, int max_iter);

// -----------------------------------------------------------------------------


#endif

