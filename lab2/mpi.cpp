// Header inclusions, if any...

#include <mpi.h>
#include <cstring>
#include "lib/gemm.h"
#include "lib/common.h"
// You can directly use aligned_alloc
// with lab2::aligned_alloc(...)


// Using declarations, if any...

void GemmParallelBlocked(const float a[kI][kK], const float b[kK][kJ],
                         float c[kI][kJ])
{
    MPI_Status status;
    int rank;
    int Size;
    int rows;
    
    MPI_Comm_size(MPI_COMM_WORLD, &Size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
     rows = kI / Size;
    //  int n = kI;
    //  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    float *a_local = (float *)lab2::aligned_alloc(4096, rows * kK * sizeof(float));
    float *c_local = (float *)lab2::aligned_alloc(4096, rows * kJ * sizeof(float));
    float *c_local2 = (float *)lab2::aligned_alloc(4096, rows * kJ * sizeof(float));
    float *b_local = (float *)lab2::aligned_alloc(4096, kK * kJ * sizeof(float));
   // float *b_local =  (float*)malloc(sizeof(float) * kK * kJ);
    //memcpy(b_local, b, kK * kJ* sizeof(float));
    
    

     MPI_Scatter(a, kK * rows, MPI_FLOAT, a_local, rows * kK, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //memcpy(a_local, a, rows * kK * sizeof(float));
    //memcpy(b_local, b, kJ * kK * sizeof(float));


   
      if(rank > 0) {
         //float *c_local = (float *)aligned_alloc(1024, rows * kJ * sizeof(float));
         //float *c_local2 = (float *)aligned_alloc(1024, rows * kJ * sizeof(float));
    
        float *b_local = (float*)lab2::aligned_alloc(4096,sizeof(float) * kK * kJ);
      //  memcpy(b_local, b, kK * kJ* sizeof(float)); 
    }
    else {
        // float *c_local = (float *)aligned_alloc(1024, rows * kJ * sizeof(float));
        //float *c_local2 = (float *)aligned_alloc(1024, rows * kJ * sizeof(float));
        memcpy(b_local, b, kK * kJ* sizeof(float));
    }
    MPI_Bcast(b_local, kK * kJ, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    
    
    //  for(int i = 0; i < rows * kJ; i++) {
    //     c_local[i] = 0;
    //       c_local2[i] = 0;
    // }

    float temp;
     for (int b_i = 0; b_i < rows; b_i += 64)
  { 
    for (int b_j = 0; b_j < kJ; b_j += 64)
    {

    float temp_buff[64][64] = {0};
    
      for (int b_k = 0; b_k < kK; b_k += 64)
      {
        for (int i = b_i; i < ((b_i + 64)); i++)
        {
         // int indexI = i - b_i;
        //  int indexJ=0;
          for (int k = b_k; k < ((b_k + 64) ); k++)
          {
            temp = a_local[i*kI+k];
            for (int j = b_j; j < ((b_j + 64) ); j++)
            {
            //  indexJ = j - b_j;
              temp_buff[i - b_i][j - b_j] += temp * b_local[k*kK+j];
             // temp_buff2[indexI][indexJ] += temp * b[k][j]; 
            }

          }
        }
      }
      //memcpy ( &temp_buff2, &temp_buff, sizeof(temp_buff) );
     //&temp_buff2 = &temp_buff;
    // #pragma omp parallel for schedule(static, 8)
      for (int i = 0; i < 64; i++)
      { 
        //int indexI = b_i + i;
        memmove(&c_local[(b_i + i) *kI+b_j], &temp_buff[i][0], sizeof(float) * 64);}
    }
  }
    
    
    // for(int i = 0; i < ( rows / 64); i++) {
    //     for(int k = 0; k < (kK / 64); k++) {
    //         for(int j = 0; j < (kJ / 64); j++) {
    //             for(int bi = 0; bi < 64; bi++) {
    //                 for(int bk = 0; bk < 64; bk++) {
    //                     int aIndex =(i * 64 + bi) * kI + k * 64 + bk;
    //                     for(int bj = 0; bj < 64; bj++) {
    //                         int cIndex = (i * 64 + bi) * kJ + (j * 64 + bj);
                            
    //                         c_local[cIndex] += a_local[aIndex] * b_local[(k * 64 + bk) * kK + (j * 64 + bj)];
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
  
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
 
  
     MPI_Gather(c_local, (rows * kJ), MPI_FLOAT, c, (rows * kJ), MPI_FLOAT, 0, MPI_COMM_WORLD);
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