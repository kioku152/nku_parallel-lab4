#include <iostream>
#include <pthread.h>
#include<windows.h>
#include"mpi.h"


using namespace std;

const int n=3000;




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

    if (rank == 0) {
        for (int j = 1; j < size; j++) {
            for (int i = j; i < n; i += size) {
                MPI_Send(&A[i][0], n, MPI_FLOAT, j, 1, MPI_COMM_WORLD);
            }
        }

    }
    else {
        for (int i = rank; i < n; i+=size) {
            MPI_Recv(&A[i][0], n, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    head = MPI_Wtime();
    for (int k = 0; k < n; k++) {
        if (k%size == rank) {
            for (int j = k + 1; j < n; j++) {
                A[j][k]=A[j][k]/A[k][k];
            }
            for (int j = 0; j < size; j++) {
                if (j != rank)
                    MPI_Send(&A[k][0], n, MPI_FLOAT, j, 0, MPI_COMM_WORLD);
            }
        }
        else {
            int src = k % size;
            MPI_Recv(&A[k][0], n, MPI_FLOAT, src, 0, MPI_COMM_WORLD, &status);
        }
        int begin = k;
        while(begin%size != rank){begin++;}
        for (int i = begin; i < n; i+=size) {
            for (int j = k + 1; j < n; j++) {
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            B[i]=B[i]-A[i][k]*B[k];
            A[i][k]=0;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
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
