/*******************************************************************************
 SMLC++ -- Sparse Matrix Library in C++ Using Old Yale Format.
 
 Author: 	Anton Zaicenco, anton.az@gmail.com
 Date: 	    January 18, 2021
 Version:   v1.0

 
 Classes and functions for manipulating sparse matrices in Old Yale format.
 
 The toolbox is a free and open-source software, which is distributed under
 the BSD 2-Clause License.
 Copyright (c) <2021>
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <complex>
#include "sparse_matrix.h"


using namespace std;


/*******************************************************************************
// sparse matrix class functions
*******************************************************************************/


//------------------------------------------------------------------------------
// check the presence of an element in a vector : return FALSE if found
bool check_presence(int *JJ, int startP, int endP, int seekC)
{
    bool t=true;
    for (int i=startP; i<endP; i++)
    {
        if (JJ[i]==seekC) {
           t = false;
        }
    }
    return t;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// find place where to insert the element in a vector (index of preceding element)
int check_location(int *JJ, int startP, int endP, int seekC)
{
    int t=0;
    for (int i=startP; i<endP; i++)
    {
        if (JJ[i]<seekC) {
           t++;
        }
    }
    return t;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// down IA capacity. size of JA and A stays same
template <typename EL>
void  sparse_Mtrx<EL>::ia_capacity_down(int new_size)
{
    // Resize by creating a new array
    // cout << " ------ IA up \n" ;
    int *IA_save = new int[new_size];   // create new smaller array
    // Copy the data
    for ( int i = 0; i < new_size; i++ )
        IA_save[i] = IA[i];
    // Change the size to match
    set_sia(new_size);
    cia = 0;
    // Destroy the old array
    delete [] IA;
    // Reset to the new array
    IA = IA_save;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// boost IA capacity. size of JA and A stays same
template <typename EL>
void  sparse_Mtrx<EL>::ia_capacity_up(int inc)
{
    // Resize by creating a new array
    // cout << " ------ IA up \n" ;
    int N = get_sia();                  // get current size of IA
    int *IA_save = new int[N + inc];    // create new larger array
    // Copy the data
    for ( int i = 0; i < N; i++ )
        IA_save[i] = IA[i];
    // Change the size to match
    set_sia(N+inc);
    cia = cia + inc;
    
    
   // Destroy the old array
   delete [] IA;
   // Reset to the new array
   IA = IA_save;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// boost JA and A capacity. size of IA stays same
template <typename EL>
void  sparse_Mtrx<EL>::ja_capacity_up(int inc)
{
    // Resize by creating a new array
    // cout << "JA, A up \n" ;
    int i, N = get_sja();                  // get current size of JA
    int *JA_save = new int[N + inc];    // create new larger array
    EL  *A_save  = new EL [N + inc];    // create new larger array
    // Copy the data
    if (N!=0) // check if size of JA,A is >0
    {
       for ( i = 0; i < N; i++ )
       {
           JA_save[i] = JA[i];
            A_save[i] =  A[i];
           //printf("A :  %4.1f", A[i]);
       }
       for ( i = N; i < N+inc; i++ )
       {
           JA_save[i] = 0;
            A_save[i] = 0;
       }
    }
    // Change the size to match
    set_sja(N+inc);
    cja = cja + inc;
   // Destroy the old array
   delete [] JA;
   delete []  A;
   // Reset to the new array
   JA = JA_save;
    A =  A_save;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// assign values directly to 3 vectors of a sparse matrix (JA, IA, A).
// This turnes matrix into coordinate form. Row compression is needed afterwards.
// Inputs:
//    G         sparse matrix object
//    a         new value
//    j         new row
//    k         new column
//    SoA       add at the end / sum with existent [1 / !=1]
template <typename EL>
void matrix_coord_vectors_assign(sparse_Mtrx<EL> *G, EL a, int j, int k, int SoA)
{
	 int nrEl = 50000;
     int i, ii, stop_loop, siaN, sjaN;
     siaN = G->get_sia();
     sjaN = G->get_sja();
     // cout << " ------ assign \n" ;
     if (G->cja<2)                 // allocate larger space for JA and A
     {
        G->ja_capacity_up(nrEl);
        // mexPrintf("JA up 1\n");
     }
     if (G->cia<2)                 // allocate larger space if capacity is 0
     {
        G->ia_capacity_up(nrEl);
        // mexPrintf("JA up 2\n");
     }
     //ii = G->get_matrix_nnz();
     ii = G->get_sia() - G->cia;
     
     if (SoA==1)                  // option "add at the end"
     {
        G->IA[ii]  = j;
        G->JA[ii]  = k;
        G->A[ii]   = a;
        G->set_matrix_nnz(ii+1);
        G->cja--;
        G->cia--;
     }
     else                         // option "sum with existent"
     {
        i         = 0;
        stop_loop = 0;
        while ( stop_loop==0 && i<ii ) {
           if (G->IA[i]==j &&  G->JA[i]==k) {
                 G->A[i]  = G->A[i]+a;
                 stop_loop=1; // exit loop
              }
           i++;
        }
        if (stop_loop==0) {
           // nothing to sum -- add at the end
           G->IA[ii]  = j;
           G->JA[ii]  = k;
           G->A[ii]   = a;
           G->set_matrix_nnz(ii+1);
           G->cja--;
           G->cia--;
        }
     }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Sorts the elements of 3 vectors (IA, JA, A) in matrix into ascending order
// using 1st and 2nd vector.
// Used before converting coordinate form into row-compressed form.
// Use "selection sort" algorithm
template <typename EL>
int sparse_Mtrx<EL>::matrix_sort_vectors_ascend()
{
    int    length_vect = get_matrix_nnz();
    // cout << length_vect << endl;
    int    i, j, m=0, jj=0, idx2, idx=0, imin=0, amin=0, a2min = 0;
    EL     a3min;
    
    for (i=0; i<(length_vect-1); i++)
    {
        imin  = i;
        amin  = IA[imin];
        idx2 = i + 1;
        for (j = idx2; j < length_vect; j++) {
             if (IA[j] < amin) {
                 imin  = j;        // new index of minimum
                 amin  = IA[imin]; // new value of minimum
             } // if
        } // for j
        if (imin==i) { continue; } // element is already in place
        a2min = JA[imin];
        a3min =  A[imin];
        IA[imin] = IA[i]; // overwrite cell with minimum, with the value of boundary cell
        JA[imin] = JA[i];
        A[imin]  =  A[i];
        IA[i]    =  amin;  // replace boundary cell with newly found minimum
        JA[i]    = a2min;
        A[i]     = a3min;
    }  // for i
    
    // now sort using 2nd vector:
    
    while (jj<length_vect)
    {
        idx = IA[jj];
        m   = jj;
        while (idx==IA[m] && m<(eL-0)) {
            m++;
        }
        
        for (i=jj; i<(m-1); i++)
        {
            imin  = i;
            amin  = JA[imin];
            // a3min = A[imin];
            idx2 = i + 1;
            for (j = idx2; j < m; j++) {
                 if (JA[j] < amin) {
                     imin  = j;        // new index of minimum
                     amin  = JA[imin]; // new value of minimum
                 } // if
            } // for j
            if (imin==i) { continue; } // element is already in place
            a3min    = A[imin];
            JA[imin] = JA[i]; // overwrite sell with minimum, with the value of boundary cell
            A[imin]  = A[i];
            JA[i]    = amin;  // replace boundary cell with newly found minimum
            A[i]     = a3min;
        }  // for i
        jj++;
    }
    return 1;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// compress row (IA) of a sparse matrix (from coordinate form).
// This turnes matrix into row-compressed form.
template <typename EL>
void sparse_Mtrx<EL>::compress_row()
{
     int i, ii; // siaN, ciaN;
     int j=0, k, nr_nz=0;   // bg=0, en=0;
     int siaN = get_rs();   // get size of IA
     int *nnz_;            // temp array to store nnz for each row
     nnz_ = new int [siaN];
     
     ii = get_matrix_nnz();
     
     for(i=0; i<siaN; i++)   // get_rs()
     {
        for(k=j; IA[k]<i+1 && k<ii; k++) {
           if (IA[k]==i) {
              nr_nz++;
              j++;
           }
        }
        //cout << "      nr of non-zeros in this row: " << nr_nz << "  row =  "  <<  i <<  "  j= " << j << endl;
        nnz_[i]   = nr_nz;
        nr_nz    = 0;
     }
     
     IA[0] = 0;
     for(i=0; i<siaN; i++) {
         IA[i+1] = IA[i] + nnz_[i];
     }
     delete[] nnz_;
     ia_capacity_down(get_rs()+1);
}
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// assign one value to sparse matrix (replace or add)
// Inputs:
//    G         sparse matrix object
//    a         new value
//    j         new row
//    k         new column
//    RoA       replace / add to the existing entry if overlap [1/0]
template <typename EL>
void sparse_assign(sparse_Mtrx<EL> *G, EL a, int j, int k, int RoA)
{
     int ii, tmpa, NZ;
     //printf(" JA 1 2 3 ------  %4d %4d %4d %4d \n", j, G->IA[j] , G->get_sja(), G->JA[G->IA[j]] );
     // new entry comes after existing 1st non-zero in row j
     // ====================================================
     if ( (G->JA[G->IA[j]]<k) && ((G->IA[j+1]-G->IA[j])!=0) && check_presence(G->JA,G->IA[j],G->IA[j+1],k)  )
     // row's 1st non-zero is earlier && row is not empty && !! no value is in this place
     {
     //printf("after : row=%3d : column=%3d : %4.1f \n", j, k, a);
            if (G->cja<2) // allocate larger space if capacity is 0
            {
               G->ja_capacity_up(1000);
               //printf("JA capacity up \n");
            }
            // where to insert in current row
            int cl = check_location(G->JA,G->IA[j], G->IA[j+1],k);

            tmpa = G->get_sja() - G->cja; // nr of useful entries in JA
            // printf(" ******* tmpa = %+6d min = %4d \n", tmpa, G->IA[j]+cl+0);
            for (ii=tmpa+1; ii>G->IA[j]+cl+0; ii--)
            //for (ii=tmpa+1; ii<0; ii--)
            {
               G->JA[ii] = G->JA[ii-1];
               G->A[ii]  = G->A[ii-1];
            }
            G->JA[G->IA[j]+cl+0] = k;
            G->A[G->IA[j]+cl+0]  = a;
            G->cja = G->cja-1;
            NZ = G->get_matrix_nnz();
            G->set_matrix_nnz(NZ+1);
            // update IA:
            for (ii=j+1; ii<G->get_rs()+1; ii++)
            {
               G->IA[ii] = G->IA[ii] + 1;
            }
     }
     // new entry comes before existing 1st non-zero in row j
     // ====================================================
     else if ((G->JA[G->IA[j]]>k) && (G->IA[j+1]-G->IA[j]!=0))
     // row's 1st non-zero begins later && row is not empty
     {
        //printf("before : row=%4d : column%4d : %6.1f \n", j,k,a);
        if (G->cja<2) // allocate larger space if capacity is <2
        {
           G->ja_capacity_up(1000);
        }
        tmpa = G->get_sja() - G->cja;
        for (ii=tmpa; ii>=(G->IA[j]); ii--)
        {
           G->JA[ii] = G->JA[ii-1];
           G->A[ii]  = G->A[ii-1];
        }
        G->JA[G->IA[j]] = k;
        G->A[G->IA[j]]  = a;
        G->cja = G->cja-1;
        NZ = G->get_matrix_nnz();
        G->set_matrix_nnz(NZ+1);
        // update IA
        for (ii=j+1; ii<G->get_rs()+1; ii++)
        {
           G->IA[ii] = G->IA[ii] + 1;
        }
     }
     // new entry comes in place of existing 1st non-zero in row j
     // ====================================================
     else if ((G->JA[G->IA[j]]==k) && (G->IA[j+1]-G->IA[j]!=0))
     // row's 1st non-zero is in this place && row is not empty
     {
     //printf("in existing place  : row=%4d : column%4d : %6.1f \n", j,k,a);
     if (RoA==1) {
           G->A[G->IA[j]] = a;                  // replace
     }
     if (RoA==0) {
           G->A[G->IA[j]] = G->A[G->IA[j]] + a; // add
     }
        //G->A[G->IA[j]+k] = a;
     }
     // new entry comes in an empty row j
     // ====================================================
     else if (G->JA[G->IA[j+1]]-G->JA[G->IA[j]]==0 && (G->IA[j+1]-G->IA[j]==0))
     {
        //printf("in empty row  : row=%4d : column%4d : %6.1f \n", j,k,a);
        // printf("capacity JA : %3d \n", G->cja);
        if (G->cja<2) // allocate larger space if JA capacity is 0
        {
           G->ja_capacity_up(1000);
           //printf("JA capacity UP \n");
        }
        
        for (int i=j+1; i<G->get_rs()+1; i++)
        {
            G->IA[i] = G->IA[i] + 1;
        }
        
        tmpa = G->get_sja() - G->cja;
        // printf("tmpa G1 G2 | %3d  %3.2f %3.2f \n", tmpa, G->A[G->IA[j]],G->A[1]);
        
        for (int ii=tmpa; ii>G->IA[j+1]-1; ii--)
        {
           G->JA[ii] = G->JA[ii-1];
           G->A[ii]  = G->A[ii-1];
        }
        
        G->JA[G->IA[j]] = k;
        G->A[G->IA[j]]  = a;
        G->cja = G->cja-1;
        NZ = G->get_matrix_nnz();
        G->set_matrix_nnz(NZ+1);
     }
     // new entry goes into a clean matrix
     // ====================================================
     else if ((G->cja==G->get_sja()))
     {
         //printf("in clean matrix  : row=%4d : column%4d : %6.1f \n", j,k,a);
         for (int i=j+1; i<=G->get_rs(); i++)
         {
             G->IA[i] = 1;
         }
         G->JA[0] = k;
         G->A[0] = a;
         G->cja = G->cja-1;
         NZ = G->get_matrix_nnz();
         G->set_matrix_nnz(NZ+1);
         // printf("capacity JA : %3d", G->cja);
     }
     // new entry goes into existing entry
     // ====================================================
     else if ( (G->JA[G->IA[j]]<k) && (G->IA[j+1]-G->IA[j]!=0) && check_presence(G->JA,G->IA[j],G->IA[j+1],k)==false  )
     {
        //printf(" into existing element  : row=%4d : column%4d : %6.1f \n", j,k,a);
        // where to insert in current row
        int cl = check_location(G->JA,G->IA[j], G->IA[j+1],k);
        if (RoA==1) {
            G->A[G->IA[j]+cl+0]  = a;                         // replace
        }
        if (RoA==0) {
            G->A[G->IA[j]+cl+0]  = G->A[G->IA[j]+cl+0]+ a;    // add
        }
     }
}

//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
// convert sparse matrix to full and print
// Inputs:
//    G         sparse matrix object
//    a         new value
//    j         new row
//    k         new column
template <typename EL>
void sparse_convert_to_full_print (sparse_Mtrx<EL> *G)
{
    EL M [G->get_rs()][G->get_cs()];
    for (int i=0; i<G->get_rs(); i++)
    {
        for (int j=0; j<G->get_cs(); j++)
        {
            M[i][j]=0;
        }
    }
    
    // printf("rows = %5d  columns : %5d\n", G->get_rs(), G->get_cs());
    for (int i=0; i<G->get_rs(); i++)
    {
        for (int j=G->IA[i]; j<=G->IA[i+1]-1; j++)
        {
            M[i][G->JA[j]] = G->A[j];
            // printf("--- %3d  %3d %5.2f ", i, G->JA[j], G->A[j]);
        }
        // printf("\n");
    }
    
    for (int i=0; i<G->get_rs(); i++)
    {
        for (int j=0; j<G->get_cs(); j++)
        {
            // printf("(%2d,%2d) %5.2f ",i,j, M[i][j]);
            //printf(" %6.2f ", M[i][j]);
                cout << setprecision(3) <<  M[i][j] << "  " ;
        }
        printf("\n");
    }
    
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// print IA, JA, A vectors from a sparse matrix object
template <typename EL>
void sparse_Mtrx<EL>::print_Mtrx()
{
    int k;
    cout << " vector IA: \n";
    for(k=0;k<get_sia()-cia+0;k++)
    {
       // printf(" %4d  ", IA[k]);
       cout << IA[k] << "  ";
    }
    cout << "\n ========================================================= \n";
    cout << " vector JA: \n";
    for(k=0;k<get_matrix_nnz();k++)
    {
       // printf(" %4d  ", JA[k]);
       cout << JA[k] << "  ";
    }
    cout << "\n ========================================================= \n";
    cout << " vector A: \n";
    for(k=0;k<get_matrix_nnz();k++)
    {
       // printf("%+6.2f  ", A[k]);
       cout << A[k] << "  " ;
    }
    cout << "\n ========================================================= \n";
};
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// get entry (row,column) from the sparse matrix
template <typename EL>
EL sparse_Mtrx<EL>::get_entry(int r, int c)
{
    int j;
    EL  ab=0;
    if (IA[r]!=IA[r+1]) {  // any non-zero entries for this row?
         //for (j=IA[r]; j<=IA[r+1]-1; j++)
         j=IA[r];
         while (j<=IA[r+1]-1)
         {
             if (JA[j]==c) {
                 ab = A[j];
                 j=IA[r+1]+1; // exit while loop
             }
             j++;
         }
    }
    return ab;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// return a specified row from a sparse matrix (as a sparse matrix)
template <typename EL>
void sparse_Mtrx<EL>::get_matrix_row(sparse_Mtrx<EL> *G, int r)
{
    int j,k=0;
    if (IA[r]!=IA[r+1]) {  // any non-zero entries for this row?
         for (j=IA[r]; j<=IA[r+1]-1; j++)
         {
            // ab = A[j];
            sparse_assign(G, A[j], 0, JA[j], 1);
            k++;
         }
    }
    G->set_matrix_nnz(k);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// return a specified column from a sparse matrix (as a sparse matrix)
template <typename EL>
void sparse_Mtrx<EL>::get_matrix_col(sparse_Mtrx<EL> *G, int c)
{
    int i,j,k=0;
    for (i=0;i<get_rs();i++)
    {
       if (IA[i]!=IA[i+1])   // any non-zero entries for this row?
       {
         j = IA[i];
         while (j<=IA[i+1]-1)
         {
            if (JA[j]==c) {
               sparse_assign(G, A[j], i, 0, 1);
               k++;
               j = IA[i+1]+1; // exit while loop
            }
            j++;
         }
       }
    }
    G->set_matrix_nnz(k);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// return a block from a sparse matrix (as a sparse matrix)
template <typename EL>
void sparse_Mtrx<EL>::get_matrix_block(sparse_Mtrx<EL> *G, int r1, int r2, int c1, int c2)
{
    int i,j,ii=-1,k=0;
    for (i=0;i<get_rs();i++)
    {
      if (i>=r1 && i<=r2) {
        ii++;
        if (IA[i]!=IA[i+1])   // any non-zero entries for this row?
        {
           j = IA[i];
           while (j<=IA[i+1]-1 && JA[j]<=c2)
           {
               if (JA[j]>=c1 && JA[j]<=c2) {
                  sparse_assign(G, A[j], ii, JA[j]-c1, 1);
                  k++;
                  //printf("... %4.1f \n",A[j]);
               }
               j++;
           }
        }
       }
    }
    G->set_matrix_nnz(k);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// multiply 2 sparse matrices: B1*B2
// r1xc1 * r2xc2 = r1xc2
template <typename EL>
void  sparse_matrix_multiply(sparse_Mtrx<EL> *Rez, sparse_Mtrx<EL> *B1, sparse_Mtrx<EL> *B2)
{
   int i, j, k, r1, c1, r2, c2;
   EL rez_1, zero_=0;
   // coord_vectors<EL> Avec(5);  // create 3 coord. vectors, each with length 5
   r1 = B1->get_rs();       // nr of rows in B1
   c1 = B1->get_cs();       // nr of columns in B1
   r2 = B2->get_rs();       // nr of rows in B2
   c2 = B2->get_cs();       // nr of columns in B2

   if (Rez->get_rs()!=r1 || Rez->get_cs()!=c2) {
          return;           // return if rezulting matrix is not of correct size
   }
   
   sparse_Mtrx<EL>   B1_row(c1,1,c1); // row of B1
   sparse_Mtrx<EL>   B2_col(r2,r2,1); // column of B2
   
   // flash result matrix:
   for (i=0; i<Rez->get_sia(); i++) {
         Rez->IA[i]=0;
   }
   Rez->cia=Rez->get_rs()+1;
   Rez->IA[0]=0;
   Rez->IA[1]=0;
   Rez->JA[0]=0;
   Rez->A[0]=0;
   Rez->set_matrix_nnz(1);
   Rez->cja = Rez->get_sja(); // full capacity
   //cout << "NNZ = " << Rez->get_matrix_nnz() << "   SJA = " << Rez->get_sja() << endl;
   
   //for (i=0;i<Rez->get_rs();i++)
    if (B1->get_matrix_nnz()!=0 && B2->get_matrix_nnz()!=0)   // any non-zero entries in B1 or B2?
    {
       for (i=0; i<r1; i++) // go thru rows of B1
       {
          if (B1->IA[i]!=B1->IA[i+1])   // any non-zero entries for this row in B1?
          {
            B1->get_matrix_row(&B1_row,  i);
              for (j=0; j<c2; j++) // go thru columns of B2
              {
                  B2->get_matrix_col(&B2_col,  j);
                  if (B2_col.get_matrix_nnz()!=0)  // any non-zero entries for this column in B2?
                  {
                     rez_1 = sparse_dot_product(&B1_row, &B2_col);
                     // sparse_convert_to_full_print(&B1_row);
                     // sparse_convert_to_full_print(&B2_col);
                     if (rez_1!=zero_)
                     {
                        // coord_vectors_assign<EL>(&Avec, rez_1, i, j, 1);
                        matrix_coord_vectors_assign(Rez, rez_1, i, j, 0);
                        //printf("assigning: row %3d  column: %3d value: %4.1f \n",i,j,rez_1);
                     }
                  }
                  // flash column matrix:
                  for (k=0; k<=B2_col.get_rs(); k++)
                  {
                     B2_col.IA[k]=0;
                  }
                  B2_col.IA[1]=0;
                  B2_col.JA[0]=0;
                  B2_col.A[0] =0;
                  B2_col.set_matrix_nnz(1);
                  // cout << "B2_col.cja = " << B2_col.cja << "   B2_col.cia = " << B2_col.cia << endl;
                  B2_col.cja = r2;
                  B2_col.cia = r2+1;
              }
            }
         // flash row matrix:
         B1_row.IA[0]=0;
         B1_row.IA[1]=0;
         B1_row.JA[0]=0;
         B1_row.A[0] =0;
         B1_row.set_matrix_nnz(1);
         // cout << "B1_row.cja = " << B1_row.cja << "   B1_row.cia = " << B1_row.cia << endl;
         B1_row.cja = c1;
         B1_row.cia = 2;
       }
    }
   // sort the vectors:
   // Avec.sort_vectors_ascend();
   Rez->matrix_sort_vectors_ascend();
   Rez->compress_row();
   // cout << " === 1" << endl;
  /* for (i=0;i<Avec.get_vectors_nnz();i++) {
      sparse_assign<EL>(Rez, Avec.aA[i], Avec.iA[i], Avec.jA[i], 1);
   }
  */
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// multiply sparse matrix by a vector: B1*b
// 
template <typename EL>
void  sparse_matrix_vect_multiply(EL *Rez, sparse_Mtrx<EL> *B1, EL *b)
{
   int i, j, r1, c1;
   int nrnzB1;              // nnz in row of B1
   int tmp, tmpA;
   
   r1 = B1->get_rs();       // nr of rows in B1
   c1 = B1->get_cs();       // nr of columns in B1
   
   
   // flash result matrix:
   //for (i=0; i<r1; i++) {
   //      Rez[i]=0;
   //}
   
   
   if (B1->get_matrix_nnz()!=0)    // any non-zero entries in B1
   {
       for (i=0; i<r1; i++)        // go thru rows of B1
       {
	      Rez[i]=0;
	      tmpA = B1->IA[i];
          //nrnzB1 = B1->IA[i+1]-B1->IA[i];
          nrnzB1 = B1->IA[i+1]-tmpA;
          if (nrnzB1!=0)           // any non-zero entries for this row in B1?
          {
              for (j=0; j<nrnzB1; j++) // go thru columns of B1
              {
                 //Rez[i] = Rez[i] + b[B1->JA[B1->IA[i]+j]]*B1->A[B1->IA[i]+j];
                 tmp = tmpA+j;
                 Rez[i] += b[B1->JA[tmp]]*B1->A[tmp];
              }
            }
       }
   }
   
}
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// multiply sparse matrix by a vector: B1*b using vector class
// 
template <typename EL>
void  sparse_A_times_b(full_Vector<EL> *Rez, sparse_Mtrx<EL> *B1, full_Vector<EL> *b)
{
   int i, j, r1, c1, tmp, tmpA;
   int nrnzB1;              // nnz in row of B1
   
   r1 = B1->get_rs();       // nr of rows in B1
   c1 = B1->get_cs();       // nr of columns in B1
   
   
   if (B1->get_matrix_nnz()!=0)    // any non-zero entries in B1
   {
       for (i=0; i<r1; i++)        // go thru rows of B1
       {
	      Rez->V[i]=0;             // flash result matrix:
	      tmpA = B1->IA[i];
          nrnzB1 = B1->IA[i+1]-tmpA;
          if (nrnzB1!=0)           // any non-zero entries for this row in B1?
          {
              for (j=0; j<nrnzB1; j++) // go thru columns of B1
              {
	             tmp = tmpA+j;
                 Rez->V[i] += b->V[B1->JA[tmp]]*B1->A[tmp];
              }
            }
       }
   }
   
}
//------------------------------------------------------------------------------









//------------------------------------------------------------------------------
// dot product. (Can be row*column, row*row, etc.)
template <typename EL>
EL  sparse_dot_product (sparse_Mtrx<EL> *B1, sparse_Mtrx<EL> *B2)
{
  int  i, idx=0, r1, c1, r2, c2;
  EL   rez=0;
  r1 = B1->get_rs();       // nr of rows in B1
  c1 = B1->get_cs();       // nr of columns in B1
  r2 = B2->get_rs();       // nr of rows in B2
  c2 = B2->get_cs();       // nr of columns in B2
  // 1. column * row:
  if (c1==1 && r2==1 && r1==c2) {
    for (i=0;i<r1;i++) {
       if (B1->IA[i]!=B1->IA[i+1])   // any non-zero entries for this row in B1?
       {
          //cout << "1st::: " << B1->A[B1->IA[i]] << endl;
          if (i==B2->JA[idx]) {
             rez = rez + B1->A[B1->IA[i]]*B2->A[idx];
             idx++;
          }

       }
       else if (B2->JA[idx]==i) {
         idx++;
       }
    }      
  }
  
  // 2. row * row:
  if (r1==1 && r2==1 && c1==c2) {
    if (B1->IA[0]!=B1->IA[1] && B2->IA[0]!=B2->IA[1])   // any non-zero entries for these matrices?
    {
       for (i=0;i<B1->get_matrix_nnz();i++) {
          if (B1->JA[i]>B2->JA[idx]) {
             while (B1->JA[i]>B2->JA[idx] && idx<B2->get_matrix_nnz())
             {
                 idx++;
             }
          }
          if (B1->JA[i]==B2->JA[idx]) {
             //cout << "multiplying:  " << B1->A[i] << " * " << B2->A[idx] << endl;
             rez = rez + B1->A[i]*B2->A[idx];
             idx++;
          }
       }
    }
  }
  
  // 3. column * column:
  if (c1==1 && c2==1 && r1==r2) {
    for (i=0;i<r1;i++) {
       if (B1->IA[i]!=B1->IA[i+1])    // any non-zero entries for this row in B1?
       {
          if (B2->IA[i]!=B2->IA[i+1]) // any non-zero entries for this row in B2?
          {
             //cout << "row=" << i << "  1st: " << B1->A[B1->IA[i]] << "  2nd: " <<  B2->A[B2->IA[i]]  << endl;
             rez = rez + B1->A[B1->IA[i]]*B2->A[B2->IA[i]];
          }

       }
    }      
  }
  
  // 4. row * column:
  if (r1==1 && c2==1 && c1==r2) {
    if (B1->IA[0]!=B1->IA[1] && B2->get_matrix_nnz()!=0)   // any non-zero entries in B1?
    {
       for (i=0;i<B1->get_matrix_nnz();i++) {
          if (   -B2->IA[ B1->JA[i] ] + B2->IA[ B1->JA[i]+1 ]>0  )
          {
             //cout << "multiplying:  " << B1->A[i] << " * " << B2->A[ B2->IA[ B1->JA[i] ] ] << endl;
             rez = rez + B1->A[i]*B2->A[ B2->IA[ B1->JA[i] ] ];
          }
       }
    }
  }
  
  return rez;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Sum of 2 sparse matrices (must be equal size): Rez = B1+B2
template <typename EL>
void  sparse_matrix_sum(sparse_Mtrx<EL> *Rez, sparse_Mtrx<EL> *B1, sparse_Mtrx<EL> *B2)
{
   int i,   j, k=0, idx, nrnzB1, nrnzB2, mx;
   int r1, c1, r2, c2;
   EL  tmpv;
   r1 = B1->get_rs();       // nr of rows in B1
   c1 = B1->get_cs();       // nr of columns in B1
   r2 = B2->get_rs();       // nr of rows in B2
   c2 = B2->get_cs();       // nr of columns in B2
   
   if (r1!=r2 || c1!=c2) {
          return;           // return if matrices are not of equal size
   }

   // cout << "rows in A: " << r1 << "   Rows in B: " << r2 << endl;
   // flash result matrix:
   for (i=0; i<Rez->get_sia(); i++) {
         Rez->IA[i]=0;
   }
   Rez->cia=Rez->get_rs()+1;
   Rez->IA[0]=0;
   Rez->IA[1]=0;
   Rez->JA[0]=0;
   Rez->A[0] =0;
   Rez->set_matrix_nnz(1);
   Rez->cja = Rez->get_sja(); // full capacity

    if (B1->get_matrix_nnz()!=0 && B2->get_matrix_nnz()!=0)   // any non-zero entries in B1 or B2?
    {
       for (i=0; i<r1; i++) // go thru rows of B1 (and B2)
       {
              nrnzB1 = B1->IA[i+1]-B1->IA[i]; // nnz in this row of B1
              nrnzB2 = B2->IA[i+1]-B2->IA[i]; // nnz in this row of B2
              k = 0;
              if (nrnzB1>nrnzB2)
              {
                   // cout << "row " << i << ". Option 1 " << endl;
                   // nnz in row of B1 is longer than in B2
                   mx   = nrnzB1;
                   j    = 0;
                   while (j<mx)        // go thru non-zeros in this row of B1
                   {
                        if (nrnzB2!=0 && k<nrnzB2)    // B2 has non-zeros
                        {
                            idx = B2->JA[B2->IA[i]+k];
                        }
                        else
                        {
                            idx = c1-0;
                        }
                        if (B1->JA[B1->IA[i]+j] < idx) {
                            matrix_coord_vectors_assign(Rez, B1->A[B1->IA[i]+j], i, B1->JA[B1->IA[i]+j], 1);
                            j++;
                        }
                        if (B1->JA[B1->IA[i]+j] > idx) {
                            matrix_coord_vectors_assign(Rez, B2->A[B2->IA[i]+k], i, idx, 1);
                            k++;
                        }
                        if (B1->JA[B1->IA[i]+j] == idx) {
                            tmpv = B2->A[B2->IA[i]+k] + B1->A[B1->IA[i]+j];
                            matrix_coord_vectors_assign(Rez, tmpv, i, idx, 1);
                            j++;
                            k++;
                        }
                   }
              }
              
              if (nrnzB1<nrnzB2)
              {
                   // nnz in row of B2 is longer than in B1
                   mx   = nrnzB2;
                   j    = 0;
                   while (j<mx)        // go thru non-zeros in this row of B2
                   {
                        if (nrnzB1!=0 && k<nrnzB1)    // B1 has non-zeros
                        {
                            idx = B1->JA[B1->IA[i]+k];
                        }
                        else
                        {
                            idx = c1-0;
                        }
                        if (B2->JA[B2->IA[i]+j] < idx) {
                            matrix_coord_vectors_assign(Rez, B2->A[B2->IA[i]+j], i, B2->JA[B2->IA[i]+j], 1);
                            j++;
                        }
                        if (B2->JA[B2->IA[i]+j] > idx) {
                            matrix_coord_vectors_assign(Rez, B1->A[B1->IA[i]+k], i, idx, 1);
                            k++;
                        }
                        if (B2->JA[B2->IA[i]+j] == idx) {
                            tmpv = B1->A[B1->IA[i]+k] + B2->A[B2->IA[i]+j];
                            matrix_coord_vectors_assign(Rez, tmpv, i, idx, 1);
                            j++;
                            k++;
                        }
                   }

              }

              if (nrnzB1==nrnzB2)
              {
                   mx   = nrnzB2;
                   j    = 0;
                   while (j<mx)        // go thru non-zeros in this row of B2
                   {
                        if (nrnzB1!=0 && k<nrnzB1)    // B1 has non-zeros
                        {
                            idx = B1->JA[B1->IA[i]+k];
                        }
                        else
                        {
                            idx = c1-0;
                        }
                        if (B2->JA[B2->IA[i]+j] < idx) {
                            matrix_coord_vectors_assign(Rez, B2->A[B2->IA[i]+j], i, B2->JA[B2->IA[i]+j], 1);
                            j++;
                        }
                        if (B2->JA[B2->IA[i]+j] > idx) {
                            matrix_coord_vectors_assign(Rez, B1->A[B1->IA[i]+k], i, idx, 1);
                            k++;
                        }
                        if (B2->JA[B2->IA[i]+j] == idx) {
                            tmpv = B1->A[B1->IA[i]+k] + B2->A[B2->IA[i]+j];
                            matrix_coord_vectors_assign(Rez, tmpv, i, idx, 1);
                            j++;
                            k++;
                        }
                        if (j==mx && idx!=c1) {
                            matrix_coord_vectors_assign(Rez, B1->A[B1->IA[i]+k], i, B1->JA[B1->IA[i]+k], 1);
                        }
                   }

              }

       }
    }

    Rez->compress_row();
    
}


//------------------------------------------------------------------------------
// Matrix * scalar
// return another resulting matrix Rez
template <typename EL>
void  sparse_matrix_scale_out(sparse_Mtrx<EL> *Rez, sparse_Mtrx<EL> *B1, EL a)
{

   int i, nrnzB1, nrrB1, nrnzRez;
   int r1, c1, rz, cz, siaRez;
   r1     = B1->get_rs();       // nr of rows in B1
   c1     = B1->get_cs();       // nr of columns in B1
   rz     = Rez->get_rs();      // nr of rows in Rez
   cz     = Rez->get_cs();      // nr of columns in Rez
   siaRez = Rez->get_sia();     // size of IA in Rez
   
   if (r1!=rz || c1!=cz) {
          return;           // return if matrices are not of equal size
   }
   nrnzB1  = B1->get_matrix_nnz();
   nrnzRez = Rez->get_matrix_nnz();
   nrrB1   = B1->get_rs();
   
   if (nrnzRez+1<nrnzB1) {
      Rez->ja_capacity_up(nrnzB1-nrnzRez+1); // increase size of JA if necessary
   }
   if (siaRez+1<r1) {
      Rez->ia_capacity_up(r1-siaRez+1);      // increase size of IA if necessary
   }

   for(i=0;i<nrnzB1;i++) {
      Rez->JA[i] = B1->JA[i];
      Rez->A[i]  = B1->A[i]*a;
   }
   for(i=0;i<=nrrB1;i++) {
      Rez->IA[i] = B1->IA[i];
   }
   Rez->set_matrix_nnz(nrnzB1);

}




/*******************************************************************************
// "coordinate vector" class functions
*******************************************************************************/


//------------------------------------------------------------------------------
// boost iA, jA and aA capacity
template <typename EL>
void  coord_vectors<EL>::vect_capacity_up(int inc)
{
    // Resize by creating a new array
    int N = get_vectors_sia();              // get current length of vectors
    int *iA_save = new int[N + inc];        // create new larger array iA
    int *jA_save = new int[N + inc];        // create new larger array jA
    EL  *aA_save = new EL [N + inc];  // create new larger array aA
    // Copy the data
    for ( int i = 0; i < N; i++ ) {
        iA_save[i] = iA[i];
        jA_save[i] = jA[i];
        aA_save[i] = aA[i];
    }
    // Change the size to match
    set_vectors_sia(N+inc);
    // Destroy the old array
    delete [] iA;
    delete [] jA;
    delete [] aA;
    // Reset to the new array
    iA = iA_save;
    jA = jA_save;
    aA = aA_save;
}
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// assign values to 3 vectors (coordinates of non-zero)
// Inputs:
//    G         vector object
//    a         new value
//    j         new row
//    k         new column
//    SoA       add at the end / sum with existent [1 / !=1]
template <typename EL>
void coord_vectors_assign(coord_vectors<EL> *G, EL a, int j, int k, int SoA)
{
     int i, ii, tmpa, stop_loop;
     tmpa = G->get_vectors_sia() - G->get_vectors_nnz(); // nr of free entries
     if (tmpa<2)                   // allocate larger space if capacity is 0
     {
        G->vect_capacity_up(10);
     }
     ii = G->get_vectors_nnz();
     if (SoA==1)                  // option "add at the end"
     {
        G->iA[ii]  = j;
        G->jA[ii]  = k;
        G->aA[ii]  = a;
        G->set_vectors_nnz(ii+1);
     }
     else                         // option "sum with existent"
     {
        i         = 0;
        stop_loop = 0;
        while ( stop_loop==0 && i<G->get_vectors_nnz() ) {
           if (G->iA[i]==j &&  G->jA[i]==k) {
                 // printf(" sum with existing \n" );
                 G->aA[i]  = G->aA[i]+a;
                 stop_loop=1; // exit loop
              }
           i++;
        }
        if (stop_loop==0) {
           // nothing to sum -- add at the end
           G->iA[ii]  = j;
           G->jA[ii]  = k;
           G->aA[ii]  = a;
           G->set_vectors_nnz(ii+1);
        }
     }     
}
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Sorts the elements of 3 vectors into ascending order using 1st and 2nd vector
// Use "selection sort" algorithm
template <typename EL>
int coord_vectors<EL>::sort_vectors_ascend()
{
    int length_vect = eL;
    int    i, j, m=0, jj=0, idx=0, imin=0, amin=0, a2min = 0;
    EL     a3min;
    for (i=0; i<(length_vect-1); i++)
    {
        imin = i;
        amin = iA[imin];
        a2min= jA[imin];
        a3min= aA[imin];
        for (j = (i + 1); j < length_vect; j++) {
             if (iA[j] < amin) {
                 imin  = j;        // new index of minimum
                 amin  = iA[imin]; // new value of minimum
                 a2min = jA[imin];
                 a3min = aA[imin];
             } // endif
        } // end for j
        // swap operation: optimized as minimum is stored in variable amin
        if (imin==i) { continue; } // element is already in place
        iA[imin] = iA[i]; // overwrite sell with minimum, with the value of boundary cell
        jA[imin] = jA[i];
        aA[imin] = aA[i];
        iA[i]    = amin;  // replace boundary cell with newly found minimum
        jA[i]    = a2min;
        aA[i]    = a3min;
    }  // end for i

    // now sort using 2nd vector:
    
    while (jj<length_vect)
    {
        idx = iA[jj];
        m   = jj;
        while (idx==iA[m] && m<(eL-1)) {
            m++;
        }
        
        for (i=jj; i<(m-1); i++)
        {
            imin = i;
            amin = jA[imin];
            a3min= aA[imin];
            for (j = (i + 1); j < m; j++) {
                 if (jA[j] < amin) {
                     imin  = j;        // new index of minimum
                     amin  = jA[imin]; // new value of minimum
                     a3min = aA[imin];
                 } // if
            } // for j
            if (imin==i) { continue; } // element is already in place
            jA[imin] = jA[i]; // overwrite sell with minimum, with the value of boundary cell
            aA[imin] = aA[i];
            jA[i]    = amin;  // replace boundary cell with newly found minimum
            aA[i]    = a3min;
        }  // for i
        jj++;
    }
    return 1;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// print iA, jA, aA vectors
template <typename EL>
void coord_vectors<EL>::print_coord_vectors()
{
    int k;
    cout << "\n ========================================================= ";
    cout << "\n ==================== vector class ======================= \n";
    for(k=0;k<get_vectors_sia();k++)
    {
       cout << iA[k] << "  " << jA[k] << "  " << aA[k] << endl;
    }
    cout << "\n ========================================================= \n";
};
//------------------------------------------------------------------------------





/*******************************************************************************
// other functions that operate on arrays (vectors)
*******************************************************************************/


//------------------------------------------------------------------------------
// Sorts the elements of 3 arrays into ascending order using int array A1 then
// int array A2
// Use "selection sort" algorithm
int sort__ascend(int A1[], int A2[], double A3[], int length_vect)
{
    int    i, j, m, jj=0, idx, imin, amin, a2min;
    double a3min;
    for (i=0; i<(length_vect-1); i++)
    {
        imin = i;
        amin = A1[imin];
        a2min= A2[imin];
        a3min= A3[imin];
        for (j = (i + 1); j < length_vect; j++) {
             if (A1[j] < amin) {
                 imin  = j;        // new index of minimum
                 amin  = A1[imin]; // new value of minimum
                 a2min = A2[imin];
                 a3min = A3[imin];
             } // if
        } // for j
        // swap operation: optimized as minimum is stored in variable amin
        if (imin==i) { continue; } // element is already in place
        A1[imin] = A1[i]; // overwrite sell with minimum, with the value of boundary cell
        A2[imin] = A2[i];
        A3[imin] = A3[i];
        A1[i]    = amin;  // replace boundary cell with newly found minimum
        A2[i]    = a2min;
        A3[i]    = a3min;
    }  // for i

    // now sort using A2:

    while (jj<length_vect)
    {
        idx = A1[jj];
        m   = jj;
        while (idx==A1[m]) { m++; }
        
        for (i=jj; i<(m-1); i++)
        {
            imin = i;
            amin = A2[imin];
            a3min= A3[imin];
            for (j = (i + 1); j < m; j++) {
                 if (A2[j] < amin) {
                     imin  = j;        // new index of minimum
                     amin  = A2[imin]; // new value of minimum
                     a3min = A3[imin];
                 } // if
            } // for j
            if (imin==i) { continue; } // element is already in place
            A2[imin] = A2[i]; // overwrite sell with minimum, with the value of boundary cell
            A3[imin] = A3[i];
            A2[i]    = amin;  // replace boundary cell with newly found minimum
            A3[i]    = a3min;
        }  // for i
        jj++;
    }
    return 1;
}
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Sorts the elements of 3 arrays into ascending order using int array A1
// Use "selection sort" algorithm
int sort__ascend_one(int A1[], double A3[], int length_vect)
{
    int    i, j, imin, amin;
    double a3min;
    for (i=0; i<(length_vect-1); i++)
    {
        imin = i;
        amin = A1[imin];
        a3min= A3[imin];
        for (j = (i + 1); j < length_vect; j++) {
             if (A1[j] < amin) {
                 imin  = j;        // new index of minimum
                 amin  = A1[imin]; // new value of minimum
                 a3min = A3[imin];
             } // end if
        } // end for j
        // swap operation: optimized as minimum is stored in variable amin
        if (imin==i) { continue; } // element is already in place
        A1[imin] = A1[i]; // overwrite sell with minimum, with the value of boundary cell
        A3[imin] = A3[i];
        A1[i]    = amin;  // replace boundary cell with newly found minimum
        A3[i]    = a3min;
    }  // for i
    
    return 1;
}
//------------------------------------------------------------------------------





//------------------------------------------------------------------------------
// boost int array capacity
void  int_array_capacity_up(int *A, int length_vect, int N)
{
    // Resize by creating a new array
    int *A_save = new int[N + length_vect];    // create new larger array
    // Copy the data
    for ( int i = 0; i < length_vect; i++ ) {
        A_save[i] = A[i];
    }
   // Reset to the new array
   A = A_save;
}
//------------------------------------------------------------------------------





//------------------------------------------------------------------------------
// Quicksort algorithm functions:

/**
 * Swap the parameters.
 * @param a - The first parameter.
 * @param b - The second parameter.
*/
void swap(int& a, int& b)
{
    int temp = a;
    a = b;
    b = temp;
}

template <typename EL>
void swap2(EL& a, EL& b)
{
    EL temp = a;
    a = b;
    b = temp;
}

/**
 * Find and return the index of pivot element.
 * @param first - The start of the sequence.
 * @param last - The end of the sequence.
 * @return - the pivot element
*/
template <typename EL>
int pivot(int* IA, int* JA, EL* A, int first, int last) 
{
    int p = first;
    int pivotElement = IA[first];
 
    for(int i = first+1 ; i <= last ; i++)
    {
        /* If you want to sort the list in the other order, change "<=" to ">" */
        if(IA[i] <= pivotElement)
        {
            p++;
            swap(IA[i], IA[p]);
            swap(JA[i], JA[p]);
            swap2(A[i],  A[p]);
        }
    }
	
    swap(IA[p], IA[first]);
    swap(JA[p], JA[first]);
    swap2(A[p], A[first]);
 
    return p;
}


template <typename EL>
void quickSort( int* IA, int* JA, EL* A, int first, int last ) 
{
    int pivotElement;
    if(first < last)
    {
        pivotElement = pivot(IA, JA, A, first, last);
        quickSort(IA, JA, A, first, pivotElement-1);
        quickSort(IA, JA, A, pivotElement+1, last);
    }
}


template <typename EL>
int sparse_Mtrx<EL>::matrix_Qsort_vectors_ascend()
{
	int    length_vect = get_matrix_nnz();
	quickSort(IA, JA, A, 0, length_vect-1);
	return 0;
}



//******************************************************************************

// The explicit instantiation part:
template class sparse_Mtrx<int>; 
template class sparse_Mtrx<float>;
template class sparse_Mtrx<double>;
template class sparse_Mtrx<complex<int> >;
template class sparse_Mtrx<complex<float> >;
template class sparse_Mtrx<complex<double> >;

template class full_Vector<int>;
template class full_Vector<float>;
template class full_Vector<double>;
template class full_Vector<complex<int> >;
template class full_Vector<complex<float> >;
template class full_Vector<complex<double> >;

template class coord_vectors<int>;
template class coord_vectors<float>;
template class coord_vectors<double>;
template class coord_vectors<complex<int> >;
template class coord_vectors<complex<float> >;
template class coord_vectors<complex<double> >;


// force instantiation of sparse matrix class functions :

template void      sparse_assign(sparse_Mtrx<int> *G, int a, int j, int k, int SoA);
template void      sparse_assign(sparse_Mtrx<float> *G, float a, int j, int k, int SoA);
template void      sparse_assign(sparse_Mtrx<double> *G, double a, int j, int k, int SoA);
template void      sparse_assign(sparse_Mtrx<complex<int> > *G, complex<int> a, int j, int k, int SoA);
template void      sparse_assign(sparse_Mtrx<complex<float> > *G, complex<float> a, int j, int k, int SoA);
template void      sparse_assign(sparse_Mtrx<complex<double> > *G, complex<double> a, int j, int k, int SoA);



template void      sparse_convert_to_full_print(sparse_Mtrx<int> *G);
template void      sparse_convert_to_full_print(sparse_Mtrx<float> *G);
template void      sparse_convert_to_full_print(sparse_Mtrx<double> *G);



template void      sparse_matrix_multiply(sparse_Mtrx<int> *Rez, sparse_Mtrx<int> *B1, sparse_Mtrx<int> *B2);
template void      sparse_matrix_multiply(sparse_Mtrx<float> *Rez, sparse_Mtrx<float> *B1, sparse_Mtrx<float> *B2);
template void      sparse_matrix_multiply(sparse_Mtrx<double> *Rez, sparse_Mtrx<double> *B1, sparse_Mtrx<double> *B2);
template void      sparse_matrix_multiply(sparse_Mtrx<complex<int> > *Rez, sparse_Mtrx<complex<int> > *B1, sparse_Mtrx<complex<int> > *B2);
template void      sparse_matrix_multiply(sparse_Mtrx<complex<float> > *Rez, sparse_Mtrx<complex<float> > *B1, sparse_Mtrx<complex<float> > *B2);
template void      sparse_matrix_multiply(sparse_Mtrx<complex<double> > *Rez, sparse_Mtrx<complex<double> > *B1, sparse_Mtrx<complex<double> > *B2);

template void      sparse_matrix_sum(sparse_Mtrx<int> *Rez, sparse_Mtrx<int> *B1, sparse_Mtrx<int> *B2);
template void      sparse_matrix_sum(sparse_Mtrx<float> *Rez, sparse_Mtrx<float> *B1, sparse_Mtrx<float> *B2);
template void      sparse_matrix_sum(sparse_Mtrx<double> *Rez, sparse_Mtrx<double> *B1, sparse_Mtrx<double> *B2);

template int       sparse_dot_product(sparse_Mtrx<int> *A1, sparse_Mtrx<int> *A2);
template float     sparse_dot_product(sparse_Mtrx<float> *A1, sparse_Mtrx<float> *A2);
template double    sparse_dot_product(sparse_Mtrx<double> *A1, sparse_Mtrx<double> *A2);
template complex<int>       sparse_dot_product(sparse_Mtrx<complex<int> > *A1, sparse_Mtrx<complex<int> > *A2);
template complex<float>     sparse_dot_product(sparse_Mtrx<complex<float> > *A1, sparse_Mtrx<complex<float> > *A2);
template complex<double>    sparse_dot_product(sparse_Mtrx<complex<double> > *A1, sparse_Mtrx<complex<double> > *A2);

template void      coord_vectors_assign(coord_vectors<int>    *G, int    a, int j, int k, int SoA);
template void      coord_vectors_assign(coord_vectors<float>  *G, float  a, int j, int k, int SoA);
template void      coord_vectors_assign(coord_vectors<double> *G, double a, int j, int k, int SoA);
template void      coord_vectors_assign(coord_vectors<complex<int> >    *G, complex<int>    a, int j, int k, int SoA);
template void      coord_vectors_assign(coord_vectors<complex<float> >  *G, complex<float>  a, int j, int k, int SoA);
template void      coord_vectors_assign(coord_vectors<complex<double> > *G, complex<double> a, int j, int k, int SoA);

template void      matrix_coord_vectors_assign(sparse_Mtrx<int>    *G, int    a, int j, int k, int SoA);
template void      matrix_coord_vectors_assign(sparse_Mtrx<float>  *G, float  a, int j, int k, int SoA);
template void      matrix_coord_vectors_assign(sparse_Mtrx<double> *G, double a, int j, int k, int SoA);
template void      matrix_coord_vectors_assign(sparse_Mtrx<complex<int> >    *G, complex<int>    a, int j, int k, int SoA);
template void      matrix_coord_vectors_assign(sparse_Mtrx<complex<float> >  *G, complex<float>  a, int j, int k, int SoA);
template void      matrix_coord_vectors_assign(sparse_Mtrx<complex<double> > *G, complex<double> a, int j, int k, int SoA);

template void      sparse_matrix_scale_out(sparse_Mtrx<int>    *Rez, sparse_Mtrx<int>    *B1, int    a);
template void      sparse_matrix_scale_out(sparse_Mtrx<float>  *Rez, sparse_Mtrx<float>  *B1, float  a);
template void      sparse_matrix_scale_out(sparse_Mtrx<double> *Rez, sparse_Mtrx<double> *B1, double a);
template void      sparse_matrix_scale_out(sparse_Mtrx<complex<int> >    *Rez, sparse_Mtrx<complex<int> >    *B1, complex<int>    a);
template void      sparse_matrix_scale_out(sparse_Mtrx<complex<float> >  *Rez, sparse_Mtrx<complex<float> >  *B1, complex<float>  a);
template void      sparse_matrix_scale_out(sparse_Mtrx<complex<double> > *Rez, sparse_Mtrx<complex<double> > *B1, complex<double> a);

template void      sparse_matrix_vect_multiply(int    *Rez, sparse_Mtrx<int>    *B1, int    *b);
template void      sparse_matrix_vect_multiply(float  *Rez, sparse_Mtrx<float>  *B1, float  *b);
template void      sparse_matrix_vect_multiply(double *Rez, sparse_Mtrx<double> *B1, double *b);
template void      sparse_matrix_vect_multiply(complex<float>  *Rez, sparse_Mtrx<complex<float> >  *B1, complex<float>  *b);
template void      sparse_matrix_vect_multiply(complex<double>  *Rez, sparse_Mtrx<complex<double> >  *B1, complex<double>  *b);

template void      sparse_A_times_b(full_Vector<int>    *Rez, sparse_Mtrx<int>    *B1, full_Vector<int>    *b);
template void      sparse_A_times_b(full_Vector<float>  *Rez, sparse_Mtrx<float>  *B1, full_Vector<float>  *b);
template void      sparse_A_times_b(full_Vector<double> *Rez, sparse_Mtrx<double> *B1, full_Vector<double> *b);

template void      sparse_A_times_b(full_Vector<complex<float> >  *Rez, sparse_Mtrx<complex<float> >  *B1, full_Vector<complex<float> >  *b);
template void      sparse_A_times_b(full_Vector<complex<double> > *Rez, sparse_Mtrx<complex<double> > *B1, full_Vector<complex<double> > *b);

//******************************************************************************




