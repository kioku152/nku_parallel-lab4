#include <iostream>
#include <pthread.h>
#include<windows.h>
#include"mpi.h"


using namespace std;

const int n=500;




float run(){
    long long head, tail, freq;


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




    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

    QueryPerformanceCounter((LARGE_INTEGER*)&head);

    for(int k=0;k<n;k++)
    {
        for(int j=k+1;j<n;j++)
            A[j][k]=A[j][k]/A[k][k];

        for(int i=k+1;i<n;i++)
        {
            for(int j=k+1;j<n;j++)
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            B[i]=B[i]-A[i][k]*B[k];
            A[i][k]=0;
        }

    }

    QueryPerformanceCounter((LARGE_INTEGER*)&tail);


    X[n-1]=B[n-1]/A[n-1][n-1];
    for(int i=n-2;i>-1;i--)
    {
        float sum=B[i];
        for(int j=i+1;j<n;j++)
            sum=sum-A[i][j]*X[j];
        X[i]=sum/A[i][i];


    }





    float answer=(float)((tail - head) * 1000.0 / freq);

    /*for(int i=0;i<n;i++)
        cout<<X[i]<<endl;*/

    return answer ;
}

int main()
{
    float time_sum;
    int count1=0;

    while(count1<3)
    {
        count1+=1;
        time_sum+=run();
        cout<<count1<<endl;
    }

    cout<<time_sum/3.0<<endl;


    return 0;
}
