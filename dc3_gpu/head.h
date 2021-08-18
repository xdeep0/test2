//
//  head.h
//  dc3 algorithm on GPU
//
//  Created by gangliao on 12/22/14.
//  Copyright (c) 2014 gangliao. All rights reserved.
//


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/scan.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>

__global__ void		Init_d_s12(int*, int);
__global__ void		keybits(int*, int*, int*, int, int);
__global__ void		InitScan(int*, int*, int*, int);
__global__ void		Set_suffix_rank(int*, int*, int*, int, int);
__global__ void		InitScan2(int*, int*, int, int);
__global__ void		Set_S0(int*, int*, int*, int, int);

void print_arr(int*, int, int);
void suffixArray(thrust::device_vector<int>&, thrust::device_vector<int>&, int, int);
//void radixPass(int* , int* , int* , int , int );
__global__ void merge_suffixes(int *, int *, int *, int *, int *, int, int, int, int);
//int leq(int a1, int a2, int b1, int b2);
//int leq2(int a1, int a2, int a3, int b1, int b2, int b3);
void read_genome(char *filename, char *buffer, int num);
int to_i(char c);
void print_suffix(char *cc, int i);
//int set_suffix_rank(int *orig_str, int *set_rank_arr, int *sorted_suff_arr, int n02, int n0);

#define MAX_ALPHA 26

#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)

struct mapping {
	__host__ __device__ int operator()(const int& x) const
	{
		return x + x / 2 + 1;
	}
};

void print_suffix(char *cc, int i)
{
	int j = 0;
	printf("%d: ", i);
	for (j = i; j < strlen(cc); j++)
		printf("%c", cc[j]);
	printf("\n");
}



int to_i(char c)
{
	return (int)c;
}


__device__ int leq(int a1, int a2, int b1, int b2) {
	return (a1 < b1 || (a1 == b1 && a2 <= b2));
}

__device__  int  leq2(int a1, int a2, int a3, int b1, int b2, int b3) {


	return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
}


void print_arr(int *arr, int no_arr, int dimension)
{
	int i = 0;
	for (i = 0; i<no_arr; i++){
		printf("%d ", arr[i]);
		if ((i + 1) % dimension == 0)
			printf(", ");
	}
	printf("\n");
}


__global__ void Init_d_s12(int* s12, int n)

{

	int index = blockIdx.x*blockDim.x + threadIdx.x;

	if (index >= n)

		return;

	s12[index] = index + index / 2 + 1;

	//printf("%d %d\n", index, s12[index]);

}


//d_SA12, d_s12, n02
__global__ void keybits(int* SA12, int* s12, int* s, int n, int i)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;

	if (index >= n)
		return;
	SA12[index] = s[s12[index] + i];

}

//s, d_SA12, d_scan, n02
__global__ void  InitScan(int* s, int* SA12, int* scan, int n)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index >= n)
		return;
	if ((s[SA12[index]] == s[SA12[index + 1]]) && (s[SA12[index] + 1] == s[SA12[index + 1] + 1]) && (s[SA12[index] + 2] == s[SA12[index + 1] + 2]))
	{
		scan[index] = 0;
	}
	else
		scan[index] = 1;
}



__global__ void Set_suffix_rank(int*  s12, int* SA12, int*  scan, int n02, int n0)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index >= n02)
		return;

	s12[SA12[index] / 3 + ((SA12[index] % 3) - 1) * n0] = scan[index] + 1;
}


__global__ void Store_unique_ranks(int*  s12, int*  SA12, int n)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index >= n)
		return;

	s12[SA12[index]] = index + 1;

}


__global__ void Compute_SA_From_UniqueRank(int*  s12, int* SA12, int n)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index >= n)
		return;
	SA12[s12[index] - 1] = index;
}



__global__ void  InitScan2(int* SA12, int* scan, int n0, int n02)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index >= n02)
		return;
	if (SA12[index] < n0)
		scan[index] = 1;
	else
		scan[index] = 0;
}

__global__ void  Set_S0(int* s0, int* SA12, int* scan, int n0, int n02)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index >= n02)
		return;
	if (SA12[index] < n0)
		s0[scan[index]] = 3 * SA12[index];

}



__global__ void merge_suffixes(int* SA0, int* SA12, int* SA, int* s, int* s12, int n0, int n02, int n)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	int left, right, mid;
	int flag = 0;
	if (index >= n0 + n02)
		return;

	if (n != n0 + n02)
	{
		flag = 1;
		if (index == n0)
			return;	
	}
	//if (n == n0 + n02 && index == n0 + n02 - 1)
		//return;

	int i, j;
	if (index < n0)
	{
		i = SA0[index];
		left = n0;
		right = n0 + n02;

		while (left < right)
		{
			mid = (left + right) / 2;
			if (SA12[mid - n0] < n0)
			{
				j = SA12[mid - n0] * 3 + 1;

				if (leq(s[j], s12[(j + 1) / 3 + ((j + 1) % 3 - 1)*n0], s[i], s12[i / 3]))
					left = mid + 1;
				else
					right = mid;
			}
			else
			{
				j = (SA12[mid - n0] - n0) * 3 + 2;

				if (leq2(s[j], s[j + 1], s12[(j + 2) / 3 + ((j + 2) % 3 - 1)*n0], s[i], s[i + 1], s12[i / 3 + n0]))
					left = mid + 1;
				else
					right = mid;
			}

		}
		//if(i == 0)
	//	printf("%%%%SA[%d] = %d index %d\n", index - n0 + left, i, index);
		SA[index + left - n0 - flag] = i;

	}
	else
	{
		if (SA12[index - n0] < n0)
		{
			i = SA12[index - n0] * 3 + 1;
			left = 0;
			right = n0;
			//leq(s[i], s12[SA12[t] + n0], s[j], s12[j/3])
			while (left < right)
			{
				mid = (left + right) / 2;
				//a1 < b1 || (a1 == b1 && a2 <= b2)

				if (leq(s[SA0[mid]], s12[SA0[mid] / 3], s[i], s12[(i + 1) / 3 + ((i + 1) % 3 - 1)*n0]))
					left = mid + 1;
				else
					right = mid;
			}
			SA[index - n0 + left - flag] = i;
			//if(i == 0)
			//printf("@@@SA[%d] = %d (%c %d) index: %d\n", index - n0 + left, i, s[i], s12[(i + 1) / 3 + ((i + 1) % 3 - 1)*n0], index);
		}
		else
		{
			i = (SA12[index - n0] - n0) * 3 + 2;
			left = 0;
			right = n0;
			while (left < right)
			{
				mid = (left + right) / 2;

				if (leq2(s[SA0[mid]], s[SA0[mid] + 1], s12[SA0[mid] / 3 + n0], s[i], s[i + 1], s12[(i + 2) / 3 + ((i + 2) % 3 - 1)*n0]))
					left = mid + 1;
				else
					right = mid;
			}
			SA[index - n0 + left - flag] = i;

			//printf("xxxxSA[%d] = %d (%c %c %d) index %d \n", index - n0 + left, i, s[i], s[i + 1], s12[(i + 2) / 3 + ((i + 2) % 3 - 1)*n0], index);
		}
	}
}