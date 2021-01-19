/*******************************************************************************
	Example showing how to solve A*x=b with the sparse matrix library 

	author :	Anton Zaicenco, anton.az@gmail.com
	Date:		January 18, 2021
*/


#include <fstream>
#include <new>
#include "sparse_matrix.h"
#include "sparse_linearsys.h"
#include <time.h>




int main()
{
   clock_t tStart = clock();
   
   int i;
   int maxIter = 5000;
   
   sparse_Mtrx<double>   AN1(1,4,4);  // matrix A
   full_Vector<double>   b(4);        // vector b
   full_Vector<double>   x(4);        // vector x
   
   double pr = 0.00001;
   
   // entries of the sparse matrix A:
   int    rrow4[]  = {  0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3};
   int    ccol4[]  = {  0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3};
   double vval4[]  = {  3, 1,.1, 1, 2,-2,.1,-2,10, 7, 7,15};
   
   // entries of vector b:
   double vecb[] = { 4, 2, 5, 6};
   
   int arrSize_A = sizeof(rrow4) / sizeof(int);    // nnz of matrix A
   int arrSize_b = sizeof(vecb) / sizeof(double);  // nnz of vector b
   // 1. Assign vectors IA,JA,A in coordinate form:   
   for (i=0;i<arrSize_A;i++) {
      matrix_coord_vectors_assign<double>(&AN1, vval4[i], rrow4[i], ccol4[i], 0);
   }
   // assign vector b:
   for (i=0;i<arrSize_b;i++) {
      b.V[i] = vecb[i];
   }
   
   // 2. sort the vectors:
   AN1.matrix_sort_vectors_ascend();
   // 3. compress IA to create row-compressed format:
   AN1.compress_row();
   
   
   // AN1.print_Mtrx();
   
   cout << " -------- matrix A:" << endl;
   sparse_convert_to_full_print<double>(&AN1);
   cout << endl;
   
   cout << " -------- vector b:" << endl;
   for(i=0; i<4; i++) {
      printf("%5.3f \n", b.V[i]);
   }
   
   
   // cout << " -------- solution of A*x=b -------- " << endl;
     // sparseJacobi(&x, &AN1, &b, pr, maxIter);
     // sparseGaussSeidel(&x, &AN1, &b, pr, maxIter);
     // sparseSOR(&x, &AN1, &b, pr, maxIter);
     sparseCG(&x, &AN1, &b, pr, maxIter);
   //}
   
   
   cout << " -------- unknown vector x:" << endl;
   for(i=0;i<4;i++) {
            printf("%5.8f \n", x.V[i]);
   }
   
   printf("\n CPU time: %.15fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
   return 1;
}

