#include "head.h"

//char *cc;

void read_data(char *filename, char *buffer, int num){
	FILE *fh;
	fh = fopen(filename, "r");
	fread(buffer, 1, num, fh);
	buffer[num] = '\0';
	fclose(fh);
}


//	merge_suffixes(SA0, SA12, SA, s, s12, n0, n1, n02, n);

void merge_suffixes(int * SA0, int * SA12, int * SA, int * s, int * s12, int n0, int n1, int n02, int n){
	int p,t,k,i,j;
	for (p=0,  t=n0-n1,  k=0;  k < n;  k++) {
		int i = GetI(); // pos of current offset 12 suffix
		int j = SA0[p]; // pos of current offset 0  suffix
		if (SA12[t] < n0 ? leq(s[i], s12[SA12[t] + n0], s[j], s12[j/3]) : leq2(s[i],s[i+1],s12[SA12[t]-n0+1],s[j],s[j+1],s12[j/3+n0]))
		{ // suffix from SA12 is smaller
			SA[k] = i;  t++;
			if (t == n02) { // done --- only SA0 suffixes left
				for (k++;  p < n0;  p++, k++) SA[k] = SA0[p];
			}
		} else {
			SA[k] = j;  p++;
			if (p == n0)  { // done --- only SA12 suffixes left
				for (k++;  t < n02;  t++, k++) SA[k] = GetI();
			}
		}
	}
}

//set_suffix_rank(s,s12,SA12,n02,n0)

int set_suffix_rank(int *orig_str, int *set_rank_arr, int *sorted_suff_arr, int n02, int n0){

	int name = 0, c0 = -1, c1 = -1, c2 = -1,i;

	for (i = 0;  i < n02;  i++) {
		if (orig_str[sorted_suff_arr[i]] != c0 || orig_str[sorted_suff_arr[i]+1] != c1 || orig_str[sorted_suff_arr[i]+2] != c2) {
			name++;
			c0 = orig_str[sorted_suff_arr[i]];
			c1 = orig_str[sorted_suff_arr[i]+1];
			c2 = orig_str[sorted_suff_arr[i]+2];
		}
		if (sorted_suff_arr[i] % 3 == 1) {
			set_rank_arr[sorted_suff_arr[i]/3] = name;
		} // left half
		else{
			set_rank_arr[sorted_suff_arr[i]/3 + n0] = name;
		} // right half

	}
	return name;
}

//radixPass(s12 , SA12, s+2, n02, K);
//radixPass(SA12, s12 , s+1, n02, K);
//radixPass(s12 , SA12, s  , n02, K);

//radixPass(s0, SA0, s, n0, K);

void radixPass(int* to_be_sorted, int* sorted_suf_arr, int* orig_str, int n, int K)
{ // count occurrences
	int *count = (int *)malloc((K + 1)*sizeof(int)); // counter array
	int i=0,t=0,sum=0;
	for (i = 0;  i <= K;  i++) count[i] = 0;


	// reset counters
	for (i = 0;  i < n;  i++){
		count[orig_str[to_be_sorted[i]]]++;
	}

	// count occurrences
	for (i = 0, sum = 0;  i <= K;  i++) { // exclusive prefix sums
		t = count[i];  count[i] = sum;  sum += t;
	}
	for (i = 0;  i < n;  i++)
		sorted_suf_arr[count[orig_str[to_be_sorted[i]]]++] = to_be_sorted[i];      // sort
}




int main(int argc, char* argv[])
{
	//freopen("data","r",stdin);
	//freopen("output.txt","w",stdout);


	clock_t start, end;						    //record time
	double runTime;


	char *filename = "genome.txt";				//load the local data set


	int n;										//input size

	char *data;									//data set pointer
	int i = 0;									//index
	int *inp;									//transformed data pointer
	int *SA;									//Suffix Array pointer

	printf("Please input the size of dataset you want to evaluate (10 - 1000000): \t");
	scanf("%d", &n);

	data = (char *) malloc((n+1)*sizeof(char));

	read_data(filename, data, n);				//read data set from the local file

	start = clock();							//record the start time


	inp = (int *)malloc( (n+3)*sizeof(int) );	//dynamic allocate memory
	SA  = (int *)malloc( (n+3)*sizeof(int) );


	for(i=0;i<n;i++)							//Ascii 'A' -> integer 0 by 'A' - 65
	{
		inp[i] = to_i(data[i]);
		//inp[i] = data[i];
	}

	inp[i]=0;inp[i+1]=0;inp[i+2]=0;				//prepare for triples

    memset(SA,0,sizeof(int)*(n+3));      		//initialize the SA array

	suffixArray(inp,SA,n,MAX_ALPHA);	        //dc3/skew algorithm

	end = clock();								//record the end time
	runTime = (end - start) / (double) CLOCKS_PER_SEC ;   //run time


	//for(i = 0 ; i < n ; i++)					//print sorted suffixes from data set
	//{
	///	printf("No.%d Index.", i);
	//	print_suffix(data, SA[i]);
	//}

	printf("CPU linear construct Suffix Array\nNUM: %d \t Time: %f Sec\n", n, runTime);


	free(data);									//free allocated memory
	return 0;
}

void suffixArray(int* s, int* SA, int n, int K) {

	int n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2;
	//n0:# of i mode 3 = 0
	//n1:# of i mode 3 = 1
	//n2:# of i mode 3 = 2
	int i=0,j=0;
	int *s12 = (int *)malloc((n02 + 3)*sizeof(int));
	int *SA12 =(int *)malloc((n02 + 3)*sizeof(int));
	int *s0 = (int *)malloc(n0*sizeof(int));
	int *SA0 = (int *)malloc(n0*sizeof(int));

	/////////////////Parallel Part /////////////////////////////////
	for(i=0;i<n02+3;i++)
	{
		SA12[i]=0;
		s12[i]=0;
	}
	///////////////////////////////////////////////////////////////
	//synchronization!!!!!!!!!!!!!!!


	// generate positions of mod 1 and mod  2 suffixes
	// the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1

   /////////////////parallel part////////////////////////////////
	for (i=0, j=0;  i < n+(n0-n1);  i++)
		if (i%3 != 0)
			s12[j++] = i;
   /////////////////////////////////////////////////////////////
   ///synchronization!!!!!!!

	radixPass(s12 , SA12, s+2, n02, K);

	radixPass(SA12, s12 , s+1, n02, K);

	radixPass(s12 , SA12, s  , n02, K);

	///////////////////////////////

	// stably sort the mod 0 suffixes from SA12 by their first character

	// find lexicographic names of triples
	int max_rank = set_suffix_rank(s,s12,SA12,n02,n0);

	// if max_rank is less than the size of s12, we have a repeat. repeat dc3.
	// else generate the suffix array of s12 directly

	if(max_rank < n02)
	{
		suffixArray(s12,SA12,n02,max_rank);
		for(i = 0;  i < n02;  i++)
			s12[SA12[i]] = i + 1;
	}else{
		for(i = 0;  i < n02;  i++)
			SA12[s12[i] - 1] = i;
	}

	for (i=0, j=0;  i < n02;  i++)
		if (SA12[i] < n0)
			s0[j++] = 3*SA12[i];

	radixPass(s0, SA0, s, n0, K);

	// merge sorted SA0 suffixes and sorted SA12 suffixes
	merge_suffixes(SA0, SA12, SA, s, s12, n0, n1, n02, n);



	//printf("End of suffix array !!\n");
}
