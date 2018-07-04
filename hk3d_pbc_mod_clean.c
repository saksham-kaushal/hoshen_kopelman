/*
gcc -o hk3_pbc_mod_clean hk3d_pbc_mod_clean.c -Wall
./hk3_pbc_mod_clean
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define max(a,b) (a>b?a:b)

void print_matrix(int ***matrix, int l, int m, int n);
int hoshen_kopelman(int boundary, int ***matrix, int l, int m, int n);
void check_labelling(int boundary, int ***matrix, int l, int m, int n);
void uf_initialize(int max_labels);
int uf_find(int x);
int uf_union(int x, int y);
int uf_make_set(void);
void uf_initialize(int max_labels);
void uf_done(void);
int get_up(int i, int m);
int get_left(int j, int n);
int get_down(int i, int m);
int get_right(int j, int n);
int get_front(int k, int l);
int get_back(int k, int l);

//##########################################################

int *labels;
int n_labels = 0;     /* length of the labels array */

//##########################################################

int main(void)
{
	int l,m,n;
	int ***matrix;
	
	if (scanf("%d %d %d",&l,&m,&n) == 3)  // m = rows, n = columns
	{
		matrix = (int ***)calloc(m, sizeof(int**));
		for (int i=0; i<l; i++)
		{
			matrix[i] = (int **)calloc(m, sizeof(int*));
			for (int j=0; j<m; j++)
			{
				matrix[i][j] = (int *)calloc(n,sizeof(int));
				for (int k=0; k<n; k++)
				{
					scanf("%d",&(matrix[i][j][k]));
				}
			}
		}
		printf(" --input-- \n");
		print_matrix(matrix,l,m,n);
		
		int clusters = hoshen_kopelman(0,matrix,l,m,n);
		int new_clusters = hoshen_kopelman(1,matrix,l,m,n);
		check_labelling(1,matrix,l,m,n);
		uf_done();
		
		printf(" --output-- \n");
		print_matrix(matrix,l,m,n);
		printf("HK reports %d clusters in pbc\n", new_clusters);
		
		for (int i=0; i<l; i++)
		{
			for (int j=0; j<m; j++)
			{
				free(matrix[i][j]);
			}
			free(matrix[i]);
		}
		free(matrix);
		return 0;
	}
}

//##########################################################

void print_matrix(int ***matrix, int l, int m, int n)
{
	for (int i=0; i<l; i++)
	{
		for (int j=0; j<m; j++)
		{
			for (int k=0; k<n; k++)
			{
				printf("%3d ",matrix[i][j][k]);
			}
			printf("\n");
		}
		printf("\n\n");
	}
}

//----------------------------------------------------------

int hoshen_kopelman(int boundary, int ***matrix, int l, int m, int n)
{
	if (boundary == 0)
	{
		uf_initialize(l*m*n / 2);
		for (int i=0; i<l; i++)
		{
			for (int j=0; j<m; j++)
			{
				for (int k=0; k<n; k++)
				{
					if (matrix[i][j][k])
					{
						int front = (i==0 ? 0 : matrix[i-1][j][k]);
						int up = (j==0 ? 0 : matrix[i][j-1][k]);
						int left = (k==0 ? 0 : matrix[i][j][k-1]);
						
						switch (!!front + !!up + !!left)
						{
							case 0:
								matrix[i][j][k] = uf_make_set();      // a new cluster
								break;
							case 1:                              // part of an existing cluster
								matrix[i][j][k] = max(front,max(up,left));       // whichever is nonzero is labelled
								break;
							case 2:                              // this site binds two clusters
								matrix[i][j][k] = (!!front == 0 ? uf_union(left,up) : (!!up == 0 ? uf_union(front,left) : uf_union(front,up)));
								break;
							case 3:                              // this site binds three clusters
	  							matrix[i][j][k] = uf_union(left,up);
	  							matrix[i][j][k] = uf_union(front,left);
	  							matrix[i][j][k] = uf_union(front,up);
	  							break;
	  					}
					}
				}
			}
		}
		int total_clusters = labels[0];
		return total_clusters;
	}
	else if (boundary == 1)
	{
		for (int i=0; i<l-1; i++)
		{
			for (int j=0; j<m-1; j++)
			{
				for (int k=0; k<n-1; k++)
				{
					if (matrix[i][j][k])		// if occupied ...
					{
						int front = matrix[get_front(i,l)][j][k];
						int up = matrix[i][get_up(j,m)][k];    //  look up  
						int left = matrix[i][j][get_left(k,n)];  //  look left  
					
						switch (!!front + !!up + !!left)
						{
							case 0:
								break;
							case 1:
								if (front) {matrix[i][j][k] = uf_union(front,matrix[i][j][k]);}    // part of an existing cluster
								if (up) {matrix[i][j][k] = uf_union(up,matrix[i][j][k]);}       // whichever is nonzero is labelled
								if (left) {matrix[i][j][k] = uf_union(left,matrix[i][j][k]);}
								break;
							case 2:                              // this site binds two clusters
								matrix[i][j][k] = (!!front == 0 ? uf_union(left,up) : (!!up == 0 ? uf_union(front,left) : uf_union(front,up)));
								break;
							case 3:                              // this site binds three clusters
	  							matrix[i][j][k] = uf_union(left,up);
	  							matrix[i][j][k] = uf_union(front,left);
	  							matrix[i][j][k] = uf_union(front,up);
	  							break;
						}
					}
				}
			}
		}
		int *new_labels_pbc = calloc(sizeof(int), n_labels); // allocate array, initialized to zero
		for (int i=0; i<l; i++)
		{
			for (int j=0; j<m; j++)
			{
				for (int k=0; k<n; k++)
				{
					if (matrix[i][j][k])
					{
						int x = uf_find(matrix[i][j][k]);
						if (new_labels_pbc[x] == 0)
						{
							new_labels_pbc[0]++;
							new_labels_pbc[x] = new_labels_pbc[0];
						}
						matrix[i][j][k] = new_labels_pbc[x];
					}
				}
			}
		}
		
		int total_clusters = new_labels_pbc[0];
		free(new_labels_pbc);
		return total_clusters;
	}
	else
	{
		printf("possible boundary conditions :\n 0-open\n1-periodic\n");
		return 0;
	}
}

//------------------------------------------------------------

void uf_initialize(int max_labels) {
  n_labels = max_labels;
  labels = calloc(sizeof(int), n_labels);
  labels[0] = 0;
}

//------------------------------------------------------------

int uf_find(int x) {
  int y = x;
  while (labels[y] != y)
    y = labels[y];
  
  while (labels[x] != x) {
    int z = labels[x];
    labels[x] = y;
    x = z;
  }
  return y;
}

//-----------------------------------------------------------

int uf_union(int x, int y) {
  return labels[uf_find(x)] = uf_find(y);
}

//-----------------------------------------------------------

int uf_make_set(void) {
  labels[0] ++;
  assert(labels[0] < n_labels);
  labels[labels[0]] = labels[0];
  return labels[0];
}

//-----------------------------------------------------------

void uf_done(void) {
  n_labels = 0;
  free(labels);
  labels = 0;
}

//-----------------------------------------------------------

int get_up(int i, int m)
{
	if (i==0) {return m-1;}
	else {return i-1;}
}

int get_left(int j, int n)
{
	if (j==0) {return n-1;}
	else {return j-1;}
}

int get_front(int k, int l)
{
	if (k==0) {return l-1;}
	else {return k-1;}
}

int get_down(int i, int m)
{
	return i%m;
}

int get_right(int j, int n)
{
	return j%n;
}

int get_back(int k, int l)
{
	return k%l;
}

//----------------------------------------------------------

void check_labelling(int boundary, int ***matrix, int l, int m, int n)
{
	int N,S,E,W,F,B;
	if(boundary == 0)
	{
		for (int i=0; i<l; i++)
		{
			for (int j=0; j<m; j++)
			{
				for (int k=0; k<n; k++)
				{
					if (matrix[i][j][k])
					{
						F = ( i==0 ? 0 : matrix[i-1][j][k] );
						B = ( i==l-1 ? 0 : matrix[i+1][j][k] );
						N = ( j==0 ? 0 : matrix[i][j-1][k] );
						S = ( j==m-1 ? 0 : matrix[i][j+1][k] );
						E = ( k==n-1 ? 0 : matrix[i][j][k+1] );
						W = ( k==0 ? 0 : matrix[i][j][k-1] );
					
						assert( N==0 || matrix[i][j][k]==N );
						assert( S==0 || matrix[i][j][k]==S );
						assert( E==0 || matrix[i][j][k]==E );
						assert( W==0 || matrix[i][j][k]==W );
						assert( F==0 || matrix[i][j][k]==F );
						assert( B==0 || matrix[i][j][k]==B );
      				}
      			}
			}
		}
	}
	else if(boundary == 1)
	{
		for (int i=0; i<l; i++)
		{
			for (int j=0; j<m; j++)
			{
				for (int k=0; k<n; k++)
				{
					if (matrix[i][j][k])
					{
						F = ( matrix[get_front(i,l)][j][k] );
						B = ( matrix[get_back(i,l)][j][k] );
						N = ( matrix[i][get_up(j,m)][k] );
						S = ( matrix[i][get_down(j,m)][k] );
						E = ( matrix[i][j][get_right(k,n)] );
						W = ( matrix[i][j][get_left(k,n)] );
					
						assert( F==0 || matrix[i][j][k]==F );
						assert( B==0 || matrix[i][j][k]==B );
						assert( N==0 || matrix[i][j][k]==N );
						assert( S==0 || matrix[i][j][k]==S );
						assert( E==0 || matrix[i][j][k]==E );
						assert( W==0 || matrix[i][j][k]==W );
					}
				}
			}
		}
	}
	else
	{
		printf("possible boundary conditions :\n 0 - open\n1 - periodic\n");
	}
}

/* The sample program reads in a matrix from standard input, runs the HK algorithm on
   it, and prints out the results.  The form of the input is two integers giving the
   dimensions of the matrix, followed by the matrix elements (with data separated by
   whitespace).  

a sample input file is the following:

3 8 8
1 1 1 1 1 1 1 1
0 0 0 0 0 0 0 1
1 0 0 0 0 1 0 1
1 0 0 1 0 1 0 1
1 0 0 1 0 1 0 1
1 0 0 1 1 1 0 1
1 1 1 1 0 0 0 1
0 0 0 1 1 1 0 1 
1 1 1 1 1 1 1 1
0 0 0 0 0 0 0 1
1 0 0 0 0 1 0 1
1 0 0 1 0 1 0 1
1 0 0 1 0 1 0 1
1 0 0 1 1 1 0 1
1 1 1 1 0 0 0 1
0 0 0 1 1 1 0 1 
1 1 1 1 1 1 1 1
0 0 0 0 0 0 0 1
1 0 0 0 0 1 0 1
1 0 0 1 0 1 0 1
1 0 0 1 0 1 0 1
1 0 0 1 1 1 0 1
1 1 1 1 0 0 0 1
0 0 0 1 1 1 0 1 

this sample input gives the following output:

--input-- 
  1   1   1   1   1   1   1   1 
  0   0   0   0   0   0   0   1 
  1   0   0   0   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   1   1   0   1 
  1   1   1   1   0   0   0   1 
  0   0   0   1   1   1   0   1 


  1   1   1   1   1   1   1   1 
  0   0   0   0   0   0   0   1 
  1   0   0   0   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   1   1   0   1 
  1   1   1   1   0   0   0   1 
  0   0   0   1   1   1   0   1 


  1   1   1   1   1   1   1   1 
  0   0   0   0   0   0   0   1 
  1   0   0   0   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   1   1   0   1 
  1   1   1   1   0   0   0   1 
  0   0   0   1   1   1   0   1 


 --output-- 
  1   1   1   1   1   1   1   1 
  0   0   0   0   0   0   0   1 
  1   0   0   0   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   1   1   0   1 
  1   1   1   1   0   0   0   1 
  0   0   0   1   1   1   0   1 


  1   1   1   1   1   1   1   1 
  0   0   0   0   0   0   0   1 
  1   0   0   0   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   1   1   0   1 
  1   1   1   1   0   0   0   1 
  0   0   0   1   1   1   0   1 


  1   1   1   1   1   1   1   1 
  0   0   0   0   0   0   0   1 
  1   0   0   0   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   1   1   0   1 
  1   1   1   1   0   0   0   1 
  0   0   0   1   1   1   0   1 


HK reports 1 clusters in pbc

+++++++++++++++++++++++++++++++++++++++++++++++++

2 10 15
1 0 1 1 0 0 1 1 1 1 1 1 1 1 0 
0 0 1 0 0 1 0 0 1 1 1 0 0 0 1 
0 1 1 1 1 0 0 0 0 1 0 1 1 1 1 
1 0 1 1 0 1 0 0 0 0 1 1 0 0 1 
0 0 1 1 0 0 0 0 1 0 0 0 1 0 1 
1 0 1 0 1 0 0 0 1 0 0 1 0 1 1 
1 1 0 0 1 0 1 1 1 1 1 1 1 1 0 
1 0 0 1 1 1 0 0 1 1 1 0 0 1 1 
1 1 0 1 1 1 0 1 0 1 0 0 1 0 0 
0 0 1 0 0 1 1 0 1 1 0 1 1 0 0 
1 0 1 1 0 0 1 1 1 1 1 1 1 1 0 
0 0 1 0 0 1 0 0 1 1 1 0 0 0 1 
0 1 1 1 1 0 0 0 0 1 0 1 1 1 1 
1 0 1 1 0 1 0 0 0 0 1 1 0 0 1 
0 0 1 1 0 0 0 0 1 0 0 0 1 0 1 
1 0 1 0 1 0 0 0 1 0 0 1 0 1 1 
1 1 0 0 1 0 1 1 1 1 1 1 1 1 0 
1 0 0 1 1 1 0 0 1 1 1 0 0 1 1 
1 1 0 1 1 1 0 1 0 1 0 0 1 0 0 
0 0 1 0 0 1 1 0 1 1 0 1 1 0 0 

this sample input gives the following output:

 --input-- 
  1   0   1   1   0   0   1   1   1   1   1   1   1   1   0 
  0   0   1   0   0   1   0   0   1   1   1   0   0   0   1 
  0   1   1   1   1   0   0   0   0   1   0   1   1   1   1 
  1   0   1   1   0   1   0   0   0   0   1   1   0   0   1 
  0   0   1   1   0   0   0   0   1   0   0   0   1   0   1 
  1   0   1   0   1   0   0   0   1   0   0   1   0   1   1 
  1   1   0   0   1   0   1   1   1   1   1   1   1   1   0 
  1   0   0   1   1   1   0   0   1   1   1   0   0   1   1 
  1   1   0   1   1   1   0   1   0   1   0   0   1   0   0 
  0   0   1   0   0   1   1   0   1   1   0   1   1   0   0 


  1   0   1   1   0   0   1   1   1   1   1   1   1   1   0 
  0   0   1   0   0   1   0   0   1   1   1   0   0   0   1 
  0   1   1   1   1   0   0   0   0   1   0   1   1   1   1 
  1   0   1   1   0   1   0   0   0   0   1   1   0   0   1 
  0   0   1   1   0   0   0   0   1   0   0   0   1   0   1 
  1   0   1   0   1   0   0   0   1   0   0   1   0   1   1 
  1   1   0   0   1   0   1   1   1   1   1   1   1   1   0 
  1   0   0   1   1   1   0   0   1   1   1   0   0   1   1 
  1   1   0   1   1   1   0   1   0   1   0   0   1   0   0 
  0   0   1   0   0   1   1   0   1   1   0   1   1   0   0 


 --output-- 
  1   0   2   2   0   0   3   3   3   3   3   3   3   3   0 
  0   0   2   0   0   4   0   0   3   3   3   0   0   0   3 
  0   2   2   2   2   0   0   0   0   3   0   3   3   3   3 
  3   0   2   2   0   5   0   0   0   0   3   3   0   0   3 
  0   0   2   2   0   0   0   0   3   0   0   0   6   0   3 
  3   0   2   0   3   0   0   0   3   0   0   3   0   3   3 
  3   3   0   0   3   0   3   3   3   3   3   3   3   3   0 
  3   0   0   3   3   3   0   0   3   3   3   0   0   3   3 
  3   3   0   3   3   3   0   7   0   3   0   0   3   0   0 
  0   0   2   0   0   3   3   0   3   3   0   3   3   0   0 


  1   0   2   2   0   0   3   3   3   3   3   3   3   3   0 
  0   0   2   0   0   4   0   0   3   3   3   0   0   0   3 
  0   2   2   2   2   0   0   0   0   3   0   3   3   3   3 
  3   0   2   2   0   5   0   0   0   0   3   3   0   0   3 
  0   0   2   2   0   0   0   0   3   0   0   0   6   0   3 
  3   0   2   0   3   0   0   0   3   0   0   3   0   3   3 
  3   3   0   0   3   0   3   3   3   3   3   3   3   3   0 
  3   0   0   3   3   3   0   0   3   3   3   0   0   3   3 
  3   3   0   3   3   3   0   7   0   3   0   0   3   0   0 
  0   0   2   0   0   3   3   0   3   3   0   3   3   0   0 


HK reports 7 clusters in pbc
*/

