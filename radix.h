// We use the radix sort algorithm
// C++ implementation of Radix Sort
#include<iostream>
using namespace std;

// A utility function to get maximum value in arr[]
int getMax(int arr[], int size)
{
    int max = arr[0];
    for (int i = 1; i < size; i++)
        if (arr[i] > max)
            max = arr[i];
    return max;
}

void CountingSort(int arr[],int idx[], int size, int div)
{
    int output[size];
    int idxout[size];
    int count[10] = {0};

    for (int i = 0; i < size; i++)
        count[ (arr[i]/div)%10 ]++;

    for (int i = 1; i < 10; i++)
        count[i] += count[i - 1];

    for (int i = size - 1; i >= 0; i--){
        output[count[ (arr[i]/div)%10] - 1] = arr[i];
        idxout[count[ (arr[i]/div)%10] - 1] = idx[i];
        count[(arr[i]/div)%10]--;
    }

    for (int i = 0; i < size; i++){
        arr[i] = output[i];
        idx[i] = idxout[i];
    }
}


void RadixSort(int arr[],int idx[], int size)
{
    int m = getMax(arr, size);
    for (int div = 1; m/div > 0; div *= 10){
            CountingSort(arr,idx, size, div);
    }
}


void RadixSPH(int N,int key[], int keyS[],int idx[]){
    for(int i=0; i<N; i++){
        keyS[i]=key[i];
    }
	RadixSort(keyS, idx, N);
}

void Discretx(int N, double h_hash, double xmin, double R[], int Rd[]){//In this case we use the h_hash same like h_smoothing lenght
		for(int i=0; i<N;i++){
			Rd[i]=0;
			Rd[i]=int((R[i]-xmin)/h_hash);
//			cout << "Para la particula "<< i <<" la posicion continua es " << R[i] << ", la posicion discreta es " << Rd[i] << "\n";
		}
}



