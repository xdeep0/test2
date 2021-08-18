#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>

void print_arr(int*, int, int);
void suffixArray(int* , int*, int , int ) ;
void radixPass(int* , int* , int* , int , int );
void merge_suffixes(int *, int *, int *, int *, int *, int, int, int, int);
int CPU_left_boundary(int* , int , char* , char* );
int CPU_right_boundary(int* , int , char* , char* );
void trans_querySeq_to_digits(char*, int*);
inline int leq(int a1, int a2, int b1, int b2);
inline int leq2(int a1, int a2, int a3, int b1, int b2, int b3);
void read_genome(char *filename, char *buffer, int num);
int to_i(char c);
void print_suffix(char *cc, int i);
int set_suffix_rank(int *orig_str, int *set_rank_arr, int *sorted_suff_arr, int n02, int n0);

#define MAX_ALPHA 26

#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)


void print_suffix(char *cc, int i)
{
	int j = 0;
	printf("%d: ",i);
	for(j=i; j < strlen(cc);j++)
		printf("%c",cc[j]);
	printf("\n");
}


int CPU_left_boundary(int* SA, int n, char* query_seq, char* refer_seq)
{
	int m, temp;
	if (strncmp(query_seq,&refer_seq[SA[0]],strlen(query_seq)) == 0)
		return 0;

	if(strncmp(query_seq,&refer_seq[SA[0]],strlen(query_seq)) < 0)
		return -1;
	else if(strncmp(query_seq,&refer_seq[SA[n-1]],strlen(query_seq)) > 0)
		return -1;
	else
	{
		int l = 0;
		int r = n;
		while(r > l+1)
		{
			m = (l + r)/2;
			temp = strncmp(query_seq,&refer_seq[SA[m]],strlen(query_seq));

			if(temp <= 0)
				r = m;
			else
				l = m;
		}
		return r;
	}
}

int CPU_right_boundary(int* SA, int n, char* query_seq, char* refer_seq)
{
	if(strncmp(query_seq,&refer_seq[SA[n-1]],strlen(query_seq)) == 0)
		return n-1;
	int m, temp;

	int l = 0;
	int r = n;
	while(r > l + 1 )
	{
		m = (l + r)/2;
		temp = strncmp(query_seq,&refer_seq[SA[m]],strlen(query_seq));

		if(temp >= 0)
			l = m;
		else
			r = m;
	}
	return l;

}


int to_i(char c)
{
	return (int)c - 65;
}


inline int leq(int a1, int a2, int b1, int b2) {
	return(a1 < b1 || (a1 == b1 && a2 <= b2));
}

inline int leq2(int a1, int a2, int a3, int b1, int b2, int b3) {
	return(a1 < b1 || (a1 == b1 && leq(a2,a3, b2,b3)));
}


void print_arr(int *arr, int no_arr,int dimension)
{
	int i = 0;
	for(i=0; i<no_arr;i++){
		printf("%d ", arr[i]);
		if( (i+1)%dimension == 0)
			printf(", ");
	}
	printf("\n");
}



void trans_querySeq_to_digits(char*querySeq, int* queryInt) //# query sequence is 1024bp
{   // A:00 C:01 G:10 T:11
	int trans[100];
	trans['A'] = 0;  // 00
	trans['C'] = 1;  // 01
	trans['G'] = 2;  // 10
	trans['T'] = 3;  // 11
	for(int i = 0; i < 48; i++) //16bp characters -> 16*2(bits) = 32 bits = (integer type)4 bytes -> 1024bp / 16bp = 48
	{
		int t = 30,result = 0;
		for(int j = 0; j < 16; j++)
		{
			result += trans[querySeq[i*16+j]] << t; //transform 16bp to an integer
			t -= 2;
		}
		queryInt[i] = result;
	}
}
