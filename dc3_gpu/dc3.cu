#include "head.h"

void read_data(char *filename, char *buffer, int num){
	FILE *fh;
	fh = fopen(filename, "r");
	fread(buffer, 1, num, fh);
	buffer[num] = '\0';

	fclose(fh);
}

// shin
__global__ void construct_bwt(char* T, char* BWT, int* SA, int n) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= n) return;
	BWT[idx] = SA[idx] == 0 ? 'A' : T[SA[idx] - 1];
}
// shin
void bwt(char *data, thrust::device_vector<int>& d_SA, int n) {
	thrust::device_vector<char> d_T(data, data + n + 1);
	thrust::device_vector<char> d_BWT(n + 1);
	char *pd_T = thrust::raw_pointer_cast(&d_T[0]);
	char *pd_BWT = thrust::raw_pointer_cast(&d_BWT[0]);
	int *pd_SA = thrust::raw_pointer_cast(&d_SA[0]);
    dim3 block(8, 1);
    dim3 grid((n + block.x - 1) / block.x, 1);
	construct_bwt<<<grid, block>>>(pd_T, pd_BWT, pd_SA, n);
	thrust::host_vector<char> h_BWT = d_BWT;
	printf("T: %s\n", data);
	printf("SA:\n");
	for (i = 0; i < n; i++) {
		printf("%d ", h_SA[i]);
	}
	printf("\nBWT:\n");
	for (i = 0; i < n; i++) {
		printf("%c", h_BWT[i]);
	}
	putchar('\n');
}

int main(int argc, char* argv[])
{
	
	float milliseconds = 0;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	char* filename = "genome.txt";				//load the local data set

	int n;										//input size

	char *data;									//data set pointer
	int i = 0;									//index
	int *SA;									//Suffix Array pointer

	printf("Please input the size of dataset you want to evaluate (10 - 1000000): \t");
	scanf("%d", &n);

	data = (char *)malloc((n + 1)*sizeof(char));

	read_data(filename, data, n);				//read data set from the local file

	data[n - 1] = 'A';	// shin

	thrust::host_vector<int> h_inp(n + 3);
	thrust::host_vector<int> h_SA(n + 3, 0);
	thrust::device_vector<int>d_inp;
	thrust::device_vector<int>d_SA;

	for (i = 0; i<n; i++)							//Ascii 'A' -> integer 0 by 'A' - 65
	{
		h_inp[i] = to_i(data[i]);
	}

	h_inp[i] = 0; h_inp[i + 1] = 0; h_inp[i + 2] = 0;				//prepare for triples
	d_inp = h_inp;
	d_SA = h_SA;


	cudaEventRecord(start);
	suffixArray(d_inp, d_SA, n, MAX_ALPHA);	        //dc3/skew algorithm

	h_SA = d_SA;
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&milliseconds, start, stop);

	printf("GPU construct Suffix Array\nNUM: %d \t Time: %f Sec\n", n, milliseconds / 1000);

	bwt(data, d_SA, n);  // shin

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	free(data);									//free allocated memory
	return 0;
}




void suffixArray(thrust::device_vector<int>& s, thrust::device_vector<int>& SA, int n, int K) {

	int n0 = (n + 2) / 3, n1 = (n + 1) / 3, n2 = n / 3, n02 = n0 + n2;

	thrust::device_vector<int>d_s12(n02 + 3, 0);
	int *pd_s12 = thrust::raw_pointer_cast(&d_s12[0]);

	thrust::device_vector<int>d_SA12(n02 + 3, 0);
	int *pd_SA12 = thrust::raw_pointer_cast(&d_SA12[0]);

	thrust::device_vector<int>d_s0(n0, 0);
	int *pd_s0 = thrust::raw_pointer_cast(&d_s0[0]);

	thrust::device_vector<int>d_SA0(n0, 0);
	int *pd_SA0 = thrust::raw_pointer_cast(&d_SA0[0]);

	thrust::device_vector<int>d_scan(n02 + 3);
	int *pd_scan = thrust::raw_pointer_cast(&d_scan[0]);
	
	
	int *pd_s = thrust::raw_pointer_cast(&s[0]);
	int *pd_SA = thrust::raw_pointer_cast(&SA[0]);

	dim3 numThreads(1024, 1, 1);
	dim3 numBlocks((n02 - 1) / 1024 + 1, 1, 1);

	
	// S12 initialization:
	//thrust::sequence(d_s12.begin(), d_s12.begin() + n02);
	//thrust::transform(d_s12.begin(), d_s12.begin() + n02, d_s12.begin(), mapping());

	Init_d_s12 <<<numBlocks, numThreads >>>(pd_s12, n02);

	//radix sort - using SA12 to store keys
	keybits <<<numBlocks, numThreads >>>(pd_SA12, pd_s12, pd_s, n02, 2);

	thrust::sort_by_key(d_SA12.begin(), d_SA12.begin() + n02, d_s12.begin());

	keybits <<<numBlocks, numThreads >>>(pd_SA12, pd_s12, pd_s, n02, 1);
	thrust::sort_by_key(d_SA12.begin(), d_SA12.begin() + n02, d_s12.begin());

	keybits <<<numBlocks, numThreads >>>(pd_SA12, pd_s12, pd_s, n02, 0);
	thrust::sort_by_key(d_SA12.begin(), d_SA12.begin() + n02, d_s12.begin());

	d_SA12 = d_s12;

	// stably sort the mod 0 suffixes from SA12 by their first character
	// find lexicographic names of triples
	
	InitScan <<<numBlocks, numThreads>>>(pd_s, pd_SA12, pd_scan, n02);

	thrust::exclusive_scan(d_scan.begin(), d_scan.begin() + n02 + 1, d_scan.begin());


	Set_suffix_rank <<<numBlocks, numThreads >>>(pd_s12, pd_SA12, pd_scan, n02, n0);


	int max_rank = d_scan[n02];

	// if max_rank is less than the size of s12, we have a repeat. repeat dc3.
	// else generate the suffix array of s12 directly

	if (max_rank < n02)
	{
		suffixArray(d_s12, d_SA12, n02, max_rank);
		Store_unique_ranks <<<numBlocks, numThreads >>>(pd_s12, pd_SA12, n02);
	}
	else{
		Compute_SA_From_UniqueRank <<<numBlocks, numThreads >>>(pd_s12, pd_SA12, n02);
	}

	InitScan2 <<<numBlocks, numThreads >>>(pd_SA12, pd_scan, n0, n02);
	thrust::exclusive_scan(d_scan.begin(), d_scan.begin() + n02, d_scan.begin()); 
	Set_S0 <<<numBlocks, numThreads >>>(pd_s0, pd_SA12, pd_scan, n0, n02);


	dim3 numBlocks3((n0 - 1) / 1024 + 1);
	keybits <<<numBlocks3, numThreads >>>(pd_SA0, pd_s0, pd_s, n0, 0);
	thrust::sort_by_key(d_SA0.begin(), d_SA0.begin() + n0, d_s0.begin());
	d_SA0 = d_s0;

	
	// merge sorted SA0 suffixes and sorted SA12 suffixes
	dim3 numBlocks2((n - 1) / 1024 + 1);
	merge_suffixes <<<numBlocks, numThreads >>>(pd_SA0, pd_SA12, pd_SA, pd_s, pd_s12, n0, n02, n);

}
