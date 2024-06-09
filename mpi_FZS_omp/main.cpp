#include <iostream>
#include <pthread.h>
#include<windows.h>
#include"mpi.h"
#include<omp.h>


using namespace std;

const int n=3000;

const int NUM_THREADS=14;


double run(int argc, char* argv[]){
    double head = 0;
    double tail = 0;


    static float A[n][n];
    static float B[n];
    static float X[n];
    static float factor[n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            A[i][j]=0;
    }
    for(int i=0;i<n;i++)
    {
        B[i]=1.0;
        A[i][i]=1.0;
        for(int j=i+1;j<n;j++)
            {
                A[i][j]= (rand() % 100)-50;
                B[i]+=A[i][j];
            }
    }

    for(int k=0;k<n;k++)
        for(int i=k+1;i<n;i++)
            {
                for(int j=0;j<n;j++)
                    A[i][j]+=A[k][j];
                B[i]+=B[k];
            }


    MPI_Init(&argc, &argv);
    int size=0;
    int rank=0;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int begin = n / size * rank;
    int end;
    if(rank == size - 1)
        end=n;
    else
        end = n / size * (rank + 1);

    if (rank == 0) {
       MPI_Request* request = new MPI_Request[n - end];
        for (int j = 1; j < size; j++) {
            int beg =n / size *j;
            int en;
            if(j == size - 1)
                en=n;
            else
                en=n / size *(j + 1);

            for (int i = beg; i < en; i++) {
                MPI_Isend(&A[i][0], n, MPI_FLOAT, j, 1, MPI_COMM_WORLD, &request[i - end]);
            }

        }
        MPI_Waitall(n - end, request, MPI_STATUS_IGNORE);

    }
    else {
        MPI_Request* request = new MPI_Request[end - begin];
        for (int i = begin; i < end; i++) {
            MPI_Irecv(&A[i][0], n, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &request[i - begin]);
        }
        MPI_Waitall(end - begin, request, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    head = MPI_Wtime();

#pragma omp parrallel if(parallel),num_threads(NUM_THREADS),private(i,j,k),schedule(dynamic)
    for (int k = 0; k < n; k++) {

        if ((begin <= k && k < end)) {
            #pragma omp for
            for (int j = k + 1; j < n; j++) {
                A[j][k]=A[j][k]/A[k][k];
            }
            MPI_Request* request = new MPI_Request[size - 1 - rank];
            for (int j = rank + 1; j < size; j++) {
                MPI_Isend(&A[k][0], n, MPI_FLOAT, j, 0, MPI_COMM_WORLD, &request[j - rank - 1]);
            }
            MPI_Waitall(size - 1 - rank, request, MPI_STATUS_IGNORE);
            if (k == end - 1)
                break;
        }
        else {
            int src;
            if (k < n / size * size)
                src = k / (n / size);
            else
                src = size - 1;
            MPI_Request request;
            MPI_Irecv(&A[k][0], n, MPI_FLOAT, src, 0, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }

        #pragma omp for
        for (int i = max(begin, k + 1); i < end; i++) {
            for (int j = k + 1; j < n; j++) {
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            B[i]=B[i]-A[i][k]*B[k];
            A[i][k]=0;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == size-1) {
        tail = MPI_Wtime();
        cout<<(tail - head)*1000<<endl;
    }
    MPI_Finalize();


    /*X[n-1]=B[n-1]/A[n-1][n-1];
    for(int i=n-2;i>-1;i--)
    {
        float sum=B[i];
        for(int j=i+1;j<n;j++)
            sum=sum-A[i][j]*X[j];
        X[i]=sum/A[i][i];


    }*/
    /*for(int i=0;i<n;i++)
        cout<<X[i]<<endl;*/

    return (tail - head) ;
}

int main(int argc, char* argv[])
{
    double time_sum;

    time_sum=run(argc,argv);

    //cout<<time_sum*1000<<endl;//  ms


    return 0;
}
