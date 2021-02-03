// Header inclusions, if any...

#include <mpi.h>
#include <cstring>
#include "lib/gemm.h"
#include "lib/common.h"
// You can directly use aligned_alloc
// with lab2::aligned_alloc(...)

#define block_size 64
// Using declarations, if any...

void GemmParallelBlocked(const float a[kI][kK], const float b[kK][kJ],
                         float c[kI][kJ])
{
    MPI_Status status;
    int rank;
    int Size;
    int rows;
    int i, j, k;
    MPI_Comm_size(MPI_COMM_WORLD, &Size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
     rows = kI / Size;
    //  int n = kI;
    //  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    float *a_local = (float *)aligned_alloc(1024, rows * kK * sizeof(float));
    float *c_local = (float *)aligned_alloc(1024, rows * kJ * sizeof(float));
    float *c_local2 = (float *)aligned_alloc(1024, rows * kJ * sizeof(float));
    float *b_local = (float *)lab2::aligned_alloc(1024, kK * kJ * sizeof(float));
   // float *b_local =  (float*)malloc(sizeof(float) * kK * kJ);
    //memcpy(b_local, b, kK * kJ* sizeof(float));
    
    

     MPI_Scatter(a, kK * rows, MPI_FLOAT, a_local, rows * kK, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //memcpy(a_local, a, rows * kK * sizeof(float));
    //memcpy(b_local, b, kJ * kK * sizeof(float));


   
      if(rank > 0) {
         //float *c_local = (float *)aligned_alloc(1024, rows * kJ * sizeof(float));
         //float *c_local2 = (float *)aligned_alloc(1024, rows * kJ * sizeof(float));
    
        float *b_local = (float*)lab2::aligned_alloc(1024,sizeof(float) * kK * kJ);
      //  memcpy(b_local, b, kK * kJ* sizeof(float)); 
    }
    else {
        // float *c_local = (float *)aligned_alloc(1024, rows * kJ * sizeof(float));
        //float *c_local2 = (float *)aligned_alloc(1024, rows * kJ * sizeof(float));
        memcpy(b_local, b, kK * kJ* sizeof(float));
    }
    MPI_Bcast(b_local, kK * kJ, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    
    
     for(i = 0; i < rows * kJ; i++) {
        c_local[i] = 0;
          c_local2[i] = 0;
    }

    int bi, bj, bk;
    int i_blocks = rows / block_size;
    int j_blocks = n / block_size;
    int k_blocks = n / block_size;

    for(i = 0; i < ( rows / block_size); i++) {
        for(k = 0; k < (kK / block_size); k++) {
            for(j = 0; j < (kJ / block_size); j++) {
                for(bi = 0; bi < block_size; bi++) {
                    for(bk = 0; bk < block_size; bk++) {
                        for(bj = 0; bj < block_size; bj++) {
                            c_local[(i * block_size + bi) * kJ + (j * block_size + bj)] += a_local[(i * block_size + bi) * kI + k * block_size + bk] * b_local[(k * block_size + bk) * kK + (j * block_size + bj)];
                        }
                    }
                }
            }
        }
    }
  
    //  for (i = 0; i < rows; i++)
    //  {
    //     for (k = 0; k < kK; k++)
    //     {
    //        for (j = 0; j < kJ; j++)
    //        {
    //            c_local[i * kI + j] += a_local[i * kI + k] * b_local[k * kK + j];
    //        }
    //    }
    //  }
 
  
     MPI_Gather(c_local, (rows * n), MPI_FLOAT, c, (rows * n), MPI_FLOAT, 0, MPI_COMM_WORLD);
   // memcpy(c,c_local2, kJ*kK*sizeof(float));
      //  for (int i = 0; i < kI; i++)
      //    for (int j = 0; j < kJ; j++)
      //    {
      //        c[i][j] = c_local2[i * kJ + j];
      //     }




    //float *a_local = (float *)aligned_alloc(8, kK * kJ * sizeof(float));
    //memcpy(a_local,&a, kK*kJ*4);

    // float *b_local = (float *)aligned_alloc(8, kK * kJ * sizeof(float));
    // memcpy(b_local, &b, kK * kJ * 4);

    // float *c_local = (float *)aligned_alloc(8, kI * kJ * sizeof(float));
    // memcpy(c_local, &c, kI * kJ * 4);
}

// Your code goes here...