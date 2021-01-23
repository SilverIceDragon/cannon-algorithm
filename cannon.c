#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include "mpi.h"
#include <math.h>
#define N 20




int main()
{

    time_t t;
    srand ((unsigned) time (&t));
    double start, end, runtime
    int i, j, k, loop;
    int **MatrixA, **MatrixB;
    int MatrixC[N][N];


    MatrixA  = (int **)malloc(sizeof(int *) * N);
    MatrixA[0] = (int *)malloc(sizeof(int) * N * N);

    for(i = 0; i < N; i++) {
        MatrixA[i] = (*MatrixA + N * i);
    }

    MatrixB  = (int **)malloc(sizeof(int *) * N);
    MatrixB[0] = (int *)malloc(sizeof(int) * N * N);
 
    for(i = 0; i < N; i++) {
        MatrixB[i] = (*MatrixB + N * i);
    }

    for (i = 0; i < N; i++) {
        for (k = 0; k < N; k++) {
            MatrixA[i][k] = rand()%10;
        }
    }

    for (k = 0; k < N; k++) {
        for (j = 0; j < N; j++) {
            MatrixB[k][j] = rand()%10;
        }
    }


//Amount of rows and cols within blocks
    const int Blockrows = 2;
    const int Blockcols = 2;

//Amount of blocks
    const int Subrows=N/Blockrows;
    const int Subcols=N/Blockcols;



    int tempA[Blockrows][Blockcols];
        for (i=0; i<Blockrows; i++) {
            for (j=0; j<Blockcols; j++) {
                tempA[i][j] = 0;
        }
    }

    int tempB[Blockrows][Blockcols];
        for (i=0; i<Blockrows; i++) {
            for (j=0; j<Blockcols; j++) {
                tempB[i][j] = 0;
        }
    }

    int tempC[Blockrows][Blockcols];
        for (i=0; i<Blockrows; i++) {
            for (j=0; j<Blockcols; j++) {
                tempC[i][j] = 0;
        }
    }


//adjacent processes
    int leftNgbr = 0, rightNgbr = 0, upperNgbr = 0, lowerNgbr = 0;


//each own processor rank and number of processors
    int rank;
    int size;

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


//basic data for cartesian grid
    int dims[2];

//dim[0] is coloum dimension and [1] is row dimension
    dims[0] = sqrt(size);
    dims[1] = sqrt(size);

//when periods = 1 (true), the end of the matrix connects again to the beginning
    int periods[2];
    periods[0] = 1;
    periods[1] = 1;

//amount of coordinates x and y
    int coordinates[2];


//create cartesian grid
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&cart_comm );


     if (size != Subrows*Subcols) {
            if (rank == 0){
                 printf("MHas to be started with %d Prozesses.\n",Subrows*Subcols);
            }
        MPI_Finalize();
        exit(-1);
    }


    MPI_Datatype blocktype;
    MPI_Datatype blocktype2;

//creating a vector so that contiguous block matrices can be sent out
    MPI_Type_vector(Blockrows, Blockcols, N, MPI_INT, &blocktype2);
    MPI_Type_create_resized( blocktype2, 0, sizeof(int), &blocktype);
    MPI_Type_commit(&blocktype);


    int disps[Subrows*Subcols];
    int counts[Subrows*Subcols];


    for (i=0; i<Subrows; i++) {
        for (j=0; j<Subcols; j++) {
            disps[i*Subcols+j] = i*N*Blockrows+j*Blockcols;
            counts [i*Subcols+j] = 1;
        }
    }


    start = MPI_Wtime();

//distribute matrices to processes block by block
    MPI_Scatterv(*MatrixA, counts, disps, blocktype, &tempA, Blockrows*Blockcols, MPI_INT, 0, cart_comm);
    MPI_Scatterv(*MatrixB, counts, disps, blocktype, &tempB, Blockrows*Blockcols, MPI_INT, 0, cart_comm);



//initial alignment: move tempA horizontally
         
      MPI_Cart_shift(cart_comm, 1, rank/sqrt(size), &leftNgbr,&rightNgbr);
      MPI_Sendrecv_replace(&tempA,Blockrows*Blockcols,MPI_INT,leftNgbr,0,rightNgbr,0,cart_comm,MPI_STATUS_IGNORE);

//initial alignment: move tempB vertically
      MPI_Cart_shift(cart_comm, 0, rank%N, &upperNgbr,&lowerNgbr);
      MPI_Sendrecv_replace(&tempB,Blockrows*Blockcols,MPI_INT,upperNgbr,0,lowerNgbr,0,cart_comm,MPI_STATUS_IGNORE);

//multiply
     for (i=0; i<Blockrows; i++) {
        for (j=0; j<Blockcols; j++) {
            for (k=0; k<Blockcols; k++) {
                tempC[i][j] += tempA[i][k]*tempB[k][j];
                }
            }
        }

//loop for shifts with multiplications
       for(loop=0;loop<(sqrt(size)-1);loop++) {
         MPI_Cart_shift(cart_comm, 1, 1, &leftNgbr,&rightNgbr);
         MPI_Cart_shift(cart_comm, 0, 1, &upperNgbr,&lowerNgbr);
         MPI_Sendrecv_replace(&tempA,4,MPI_INT,leftNgbr,1,rightNgbr,1,cart_comm,MPI_STATUS_IGNORE);
         MPI_Sendrecv_replace(&tempB,4,MPI_INT,upperNgbr,1,lowerNgbr,1,cart_comm,MPI_STATUS_IGNORE);
              for (i=0; i<Blockrows; i++) {
                for (j=0; j<Blockcols; j++) {
                    for (k=0; k<Blockcols; k++) {
                        tempC[i][j] += tempA[i][k]*tempB[k][j];
                    }
                }
            }
        }


//collect result block matrices of the individual processes on process 0
       MPI_Gatherv(&tempC, Blockrows*Blockcols, MPI_INT, *MatrixC, counts, disps, blocktype, 0, cart_comm);


       end = MPI_Wtime();
       runtime =  end-start;

//local matrices
    for (loop=0; loop<size; loop++) {
        if (loop == rank) {
            printf("Rank = %d\n", rank);
            printf("Local matrix:\n");
            for (i=0; i<Blockrows; i++) {
                for (j=0; j<Blockcols; j++) {
                    printf("%3d ",(int)tempC[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }
        MPI_Barrier(cart_comm);
    }


//end-matrix
     if (rank == 0){
       printf ("Result:");
        for (i = 0; i < N; i++) {
         printf ("\n\t");
            for (j = 0; j < N; j++)
              printf (" %5d",MatrixC[i][j]);
            }
         printf ("\n");
         printf ("\nDuration = %f seconds\n", runtime);
         printf ("\n");
    }

 MPI_Type_free(&blocktype);
 MPI_Comm_free(&cart_comm);


 MPI_Finalize();



 return 0;

} 
