/*******************************************************************************
   Functions for solving a system of linear equations with sparse matrix library

	Jacobi:        
    sparseJacobi(full_Vector<EL> *x2, sparse_Mtrx<EL> *B1, full_Vector<EL> *b, double pr)
	
	Gauss-Seidel:
    sparseGaussSeidel(full_Vector<EL> *x2, sparse_Mtrx<EL> *B1, full_Vector<EL> *b, double pr)
    
    SOR:
    sparseSOR(full_Vector<EL> *x2, sparse_Mtrx<EL> *B1, full_Vector<EL> *b, double pr)
    
    CG:
    sparseCG(full_Vector<EL> *x2, sparse_Mtrx<EL> *B1, full_Vector<EL> *b, double pr, int max_iter)


	author :	Anton Zaicenco, anton.az@gmail.com
	Date:		January 18, 2021
*/

#include <complex>
#include <math.h>
#include "sparse_linearsys.h"

using namespace std;

//------------------------------------------------------------------------------

template <typename EL>
EL absVal(EL a)
{
  EL c = 0;
  if(a>=c) {
   return a;
  }
  else {
   return -a;
  }
}

template <typename EL>
EL absVal(complex<EL> a)
{
  return abs(a);
  //return sqrt(a.imag()*a.imag() + a.real()*a.real());
}



// Several functions for CG method:
//------------------------------------------------------------------------------
// dot product of two 1D arrays. Used in iteration loop
template <typename EL>
EL dotArrays(EL *a, EL *b, int ln)
{
	int i;
	EL Rez = 0;
	for(i=0; i<ln; i++) {
		Rez += a[i]*b[i];
	}
	return Rez;
}
//------------------------------------------------------------------------------
// array1 = array2 + array3 * scalar. Used in iteration loop
template <typename EL>
void scaleSumArrays(EL *a1, EL *a2, EL *a3, EL sc, int ln)
{
	int i;
	for(i=0; i<ln; i++) {
		a1[i] = a2[i] + a3[i]*sc;
	}
}
//------------------------------------------------------------------------------
// full_vector = full_vector + array * scalar. Used in iteration loop
template <typename EL>
void scaleSumArray__Vect(full_Vector<EL> *a1, full_Vector<EL> *a2, EL *a3, EL sc, int ln)
{
	int i;
	for(i=0; i<ln; i++) {
		a1->V[i] = a2->V[i] + a3[i]*sc;
	}
}
//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
// solve A*x=b using Jacobi method
// matrix A must be symmetric positive-definite with non-zero diagonal terms
//
// since: A = D + L + U, we can write A*x=b in the form:
//        D*x + (L+U)*x = b,    and solve it in iteration:
//        x_{i+1} = -inv(D)*(L+U)*x_i + inv(D)*b
//
template <typename EL>
void  sparseJacobi(full_Vector<EL> *x2, sparse_Mtrx<EL> *B1, full_Vector<EL> *b, double pr, int max_iter)
{
	//int max_iter = 5000;            // max nr of iterations
	int k, j, i  = 0;               // iteration loop index
	int r1 = B1->get_rs();          // nr of rows in B1
	int nrnzB1;                     // nnz in row "j" of B1
	EL  nl = 0;                     // used for "if"
	EL  tmp_, tmp2;                 // temp variable to store B1(j,j)
	
	EL *xi    =   new EL[r1];       // it is faster than full_Vector class by ~5%
	EL *xitmp =   new EL[r1];
	
	for(k=0; k<r1; k++) {
    	x2->V[k] = nl;                    // clean output vector
	}
	
	while(i<max_iter) {
    for(j=0; j<r1; j++) {                 // go thru rows
        xi[j] = nl;
        //xi.V[j] = nl;
        nrnzB1 = B1->IA[j+1]-B1->IA[j];   // nnz in this row of B1
        if (nrnzB1!=0)                    // any non-zero entries for this row in B1?
        {
           for(k=0; k<nrnzB1; k++) {      // go thru columns
                if((B1->JA[B1->IA[j]+k])!=j) {
                   xi[j] = xi[j] + B1->A[B1->IA[j]+k] * x2->V[B1->JA[B1->IA[j]+k]];
                }
           }
           tmp_     = B1->get_entry(j,j); // get diagonal entry
           if ( tmp_!=nl ) {
              xi[j] =  (b->V[j] - xi[j]) / tmp_;
           }
        }
        x2->V[j] = xi[j];
	}
	
	sparse_matrix_vect_multiply(xitmp, B1, xi);
	
	tmp_ = nl;
	for(k=0; k<r1; k++) {
    	tmp2  = xitmp[k] - b->V[k];
    	tmp_ += tmp2*tmp2;
	}
	if ( sqrt((double) abs(tmp_) )<pr ) {
    	i = max_iter;
	}
	i++;
	}
	delete[] xi;
	delete[] xitmp;
}


//------------------------------------------------------------------------------
// solve A*x=b using Gauss-Seidel method
// matrix A must be symmetric positive-definite with non-zero diagonal terms
//
// since: A = D + L + U, we can write A*x=b in the form:
//        U*x + (L+D)*x = b,    and solve it in iteration:
//        x_{i+1} = -inv(L+D)*U*x_i + inv(L+D)*b
//
template <typename EL>
void  sparseGaussSeidel(full_Vector<EL> *x2, sparse_Mtrx<EL> *B1, full_Vector<EL> *b, double pr, int max_iter)
{
	//int max_iter = 5000;            // max nr of iterations
	int k, j, m, i  = 0;               // iteration loop index
	int r1 = B1->get_rs();          // nr of rows in B1
	int nrnzB1;                     // nnz in row "j" of B1
	EL  nl = 0;                     // used for "if"
	EL  tmp_, tmp2;                 // temp variable to store B1(j,j)
	
	EL *xi    =   new EL[r1];       // it is faster than full_Vector class by ~5%
	EL *xitmp =   new EL[r1];
	
	for(k=0; k<r1; k++) {
    	x2->V[k]    = nl;                    // clean output vector
	}
	
	while(i<max_iter) {
     for(j=0; j<r1; j++) {                 // go thru rows
         xi[j] = nl;
         //xi.V[j] = nl;
         nrnzB1 = B1->IA[j+1]-B1->IA[j];   // nnz in this row of B1
         if (nrnzB1!=0)                    // any non-zero entries for this row in B1?
         {
            for(k=0;  (B1->JA[B1->IA[j]+k])<j; k++) {      // go thru columns
                 xi[j] = xi[j] + B1->A[B1->IA[j]+k] * xi[B1->JA[B1->IA[j]+k]];
            }
            for(m=k+1;  m<nrnzB1; m++) {      // go thru columns
                 xi[j] = xi[j] + B1->A[B1->IA[j]+m] * x2->V[B1->JA[B1->IA[j]+m]];
            }
            tmp_     = B1->get_entry(j,j); // get diagonal entry
            if ( tmp_!=nl ) {
               xi[j] =  (b->V[j] - xi[j]) / tmp_;
            }
         }
         x2->V[j] = xi[j];
	 }
 	
	 sparse_matrix_vect_multiply(xitmp, B1, xi);
 	
	 tmp_ = nl;
	 for(k=0; k<r1; k++) {
    	 tmp2  = xitmp[k] - b->V[k];
    	 tmp_ += tmp2*tmp2;
     }

	 if ( sqrt((double) abs(tmp_) )<pr ) {
       i = max_iter;
	 }
	 i++;
	}
	delete[] xi;
	delete[] xitmp;
}


//------------------------------------------------------------------------------
// solve A*x=b using SOR method
// matrix A must be symmetric positive-definite with non-zero diagonal terms
//
// refines Gauss-Seidel method by introducing over-relaxation parameter "omega".
//
template <typename EL>
void  sparseSOR(full_Vector<EL> *x2, sparse_Mtrx<EL> *B1, full_Vector<EL> *b, double pr, int max_iter)
{
  //int max_iter = 11000;            // max nr of iterations
  int k, j, m, i  = 0;               // iteration loop index
  int r1 = B1->get_rs();          // nr of rows in B1
  int nrnzB1;                     // nnz in row "j" of B1
  EL  nl = 0;                     // used for "if"
  EL  tmp_, tmp2;                 // temp variable to store B1(j,j)
  EL  omega_ = (EL) 0.4;          // 0 .. 2
  
  EL *xi    =   new EL[r1];       // it is faster than full_Vector class by ~5%
  EL *xitmp =   new EL[r1];
  int btmp;
  EL  b2tmp;
  
  for(k=0; k<r1; k++) {
     x2->V[k]    = nl;                    // clean output vector
  }
  
  while(i<max_iter) {
    for(j=0; j<r1; j++) {                 // go thru rows
        xi[j]  = nl;
        btmp   = B1->IA[j];
        nrnzB1 = B1->IA[j+1]-btmp;        // nnz in this row of B1
        if (nrnzB1!=0)                    // any non-zero entries for this row in B1?
        {
           for(k=0;  (B1->JA[btmp+k])<j; k++) {      // go thru columns
                xi[j] = xi[j] + B1->A[btmp+k] * xi[B1->JA[btmp+k]];
           }
           for(m=k+1;  m<nrnzB1; m++) {      // go thru columns
                xi[j] = xi[j] + B1->A[btmp+m] * x2->V[B1->JA[btmp+m]];
           }
           tmp_     = B1->get_entry(j,j); // get diagonal entry
           if ( tmp_!=nl ) {
              xi[j] =  (b->V[j] - xi[j]) / tmp_;
           }
        }
        b2tmp = x2->V[j];
        x2->V[j] = b2tmp + omega_*(xi[j]-b2tmp);
   }
   sparse_matrix_vect_multiply(xitmp, B1, xi);
   tmp_ = nl;
   for(k=0; k<r1; k++) {
   		tmp2  = xitmp[k] - b->V[k];
    	tmp_ += tmp2*tmp2;
   }
   
   if ( sqrt((double) abs(tmp_) )<pr ) {
      i = max_iter;
   }
   i++;
  }
  
  delete[] xi;
  delete[] xitmp;
}


/*------------------------------------------------------------------------------
 solve A*x=b using CG method
 A does not to be diagonally dominant, but is symmetric and positive-definite.
 Inputs:
 	x2			solution x
 	B1			sparse matrix A
 	b			rhs full vector b
 	pr			precision (to stop iterations)
 	max_iter	max number of iterations
*/

template <typename EL>
int  sparseCG(full_Vector<EL> *x2, sparse_Mtrx<EL> *B1, full_Vector<EL> *b, double pr, int max_iter)
{
	int k, idx  = 0;                 // iteration loop index
	int r1		= B1->get_rs();      // nr of rows in B1
	full_Vector<EL> fv(r1);
	EL *p		= new EL[r1];        // search direction
	EL *q		= new EL[r1];
	EL *resid	= new EL[r1];        // residual vector
	
	EL nl = 0;                       // used for "if"
	EL rho_1=0, rho_2=0, beta_, alpha_, var1, tmp1;
	EL pr2 = (EL)(pr*pr);
	
	for(k=0; k<r1; k++) {
		x2->V[k] = nl;            // initial guess = 0
		resid[k] = b->V[k];       // residual = b-A*xi
	}
	
	while(idx<max_iter) {
		rho_2 = rho_1;
		rho_1 =  dotArrays(resid,resid,r1);
		if (idx==0) {
			for(k=0; k<r1; k++) {
				p[k] = resid[k];
			}
		}
		else {
			beta_ = rho_1/rho_2;
   			scaleSumArrays(p,resid,p,beta_,r1);   // p = resid + beta*p;
		}
		sparse_matrix_vect_multiply(q,B1,p);      // q = A*p;
   		alpha_ = rho_1/dotArrays(p,q,r1);         // alpha = rho_1/(p'*q);
   		scaleSumArray__Vect(x2,x2,p,alpha_,r1);   // x2 = x2 + alpha*p;
   		scaleSumArrays(resid,resid,q,-alpha_,r1); // res = res - alpha*q;
   		
   		sparse_A_times_b(&fv, B1, x2);            // A*x
   		var1 = 0;
		for(k=0; k<r1; k++) {
			tmp1 =  fv.V[k] - b->V[k];
			var1 += tmp1*tmp1;          // (A*xi-b)^2
		}
		
   		if ( var1<pr2 ) {
   			break; // terminate the main loop if norm<precision
		}
		idx++;
	}
	// cout << "iterations: " << idx << endl;
	
	delete[] p;
	delete[] q;
	delete[] resid;
	
	return idx;
}


//******************************************************************************

// The explicit instantiation part:
template void   sparseJacobi(full_Vector<int>     *Rez, sparse_Mtrx<int>     *A, full_Vector<int>     *b, double pr, int max_iter);
template void   sparseJacobi(full_Vector<float>   *Rez, sparse_Mtrx<float>   *A, full_Vector<float>   *b, double pr, int max_iter);
template void   sparseJacobi(full_Vector<double>  *Rez, sparse_Mtrx<double>  *A, full_Vector<double>  *b, double pr, int max_iter);

template void   sparseJacobi(full_Vector<complex<float> >  *Rez, sparse_Mtrx<complex<float> >  *A, full_Vector<complex<float> >  *b, double pr, int max_iter);
template void   sparseJacobi(full_Vector<complex<double> > *Rez, sparse_Mtrx<complex<double> > *A, full_Vector<complex<double> > *b, double pr, int max_iter);


template void   sparseGaussSeidel(full_Vector<int>     *Rez, sparse_Mtrx<int>     *A, full_Vector<int>     *b, double pr, int max_iter);
template void   sparseGaussSeidel(full_Vector<float>   *Rez, sparse_Mtrx<float>   *A, full_Vector<float>   *b, double pr, int max_iter);
template void   sparseGaussSeidel(full_Vector<double>  *Rez, sparse_Mtrx<double>  *A, full_Vector<double>  *b, double pr, int max_iter);

template void   sparseGaussSeidel(full_Vector<complex<float> >  *Rez, sparse_Mtrx<complex<float> >  *A, full_Vector<complex<float> >  *b, double pr, int max_iter);
template void   sparseGaussSeidel(full_Vector<complex<double> > *Rez, sparse_Mtrx<complex<double> > *A, full_Vector<complex<double> > *b, double pr, int max_iter);

template void   sparseSOR(full_Vector<int>    *x2, sparse_Mtrx<int>    *A, full_Vector<int>    *b, double pr, int max_iter);
template void   sparseSOR(full_Vector<float>  *x2, sparse_Mtrx<float>  *A, full_Vector<float>  *b, double pr, int max_iter);
template void   sparseSOR(full_Vector<double> *x2, sparse_Mtrx<double> *A, full_Vector<double> *b, double pr, int max_iter);

template void   sparseSOR(full_Vector<complex<float> >  *x2, sparse_Mtrx<complex<float> >  *A, full_Vector<complex<float> >  *b, double pr, int max_iter);
template void   sparseSOR(full_Vector<complex<double> > *x2, sparse_Mtrx<complex<double> > *A, full_Vector<complex<double> > *b, double pr, int max_iter);

template int   sparseCG(full_Vector<int>    *x2, sparse_Mtrx<int>    *A, full_Vector<int>    *b, double pr, int max_iter);
template int   sparseCG(full_Vector<float>  *x2, sparse_Mtrx<float>  *A, full_Vector<float>  *b, double pr, int max_iter);
template int   sparseCG(full_Vector<double> *x2, sparse_Mtrx<double> *A, full_Vector<double> *b, double pr, int max_iter);


//end
