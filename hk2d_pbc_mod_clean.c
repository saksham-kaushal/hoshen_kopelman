#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define max(a,b) (a>b?a:b)

void print_matrix(int **matrix, int m, int n);
int hoshen_kopelman(int boundary, int **matrix, int m, int n);
void check_labelling(int boundary, int **matrix, int m, int n);
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

//##########################################################

int *labels;
int n_labels = 0;     /* length of the labels array */

//##########################################################

int main(void)
{
	int m,n;
	int **matrix;
	
	if (scanf("%d %d",&m,&n) == 2)  // m = rows, n = columns
	{
		matrix = (int **)calloc(m, sizeof(int*));
		for (int i=0; i<m; i++)
		{
			matrix[i] = (int *)calloc(n, sizeof(int));
			for (int j=0; j<n; j++)
			{
				scanf("%d",&(matrix[i][j]));
			}
		}
		printf(" --input-- \n");
		print_matrix(matrix,m,n);
		
		int clusters = hoshen_kopelman(0,matrix,m,n);
		int new_clusters = hoshen_kopelman(1,matrix,m,n);
		check_labelling(1,matrix,m,n);
		uf_done();
		
		printf(" --output-- \n");
		print_matrix(matrix,m,n);
		printf("HK reports %d clusters in pbc\n", new_clusters);
		
		for (int i=0; i<m; i++)
		{
			free(matrix[i]);
		}
		free(matrix);
		return 0;
	}
}

//##########################################################

void print_matrix(int **matrix, int m, int n)
{
	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
		{
			printf("%3d ",matrix[i][j]);
		}
	printf("\n");
	}
}

//----------------------------------------------------------

int hoshen_kopelman(int boundary, int **matrix, int m, int n)
{
	if (boundary == 0)
	{
		uf_initialize(m * n / 2);
		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
			{
				if (matrix[i][j])
				{
					int up = (i==0 ? 0 : matrix[i-1][j]);
					int left = (j==0 ? 0 : matrix[i][j-1]);
					
					switch (!!up + !!left)
					{
						case 0:
							matrix[i][j] = uf_make_set();      // a new cluster
							break;
						case 1:                              // part of an existing cluster
							matrix[i][j] = max(up,left);       // whichever is nonzero is labelled
							break;
						case 2:                              // this site binds two clusters
							matrix[i][j] = uf_union(left,up);
							break;
					}
				}
			}
		}
		int total_clusters = labels[0];
		return total_clusters;
	}
	else if (boundary == 1)
	{
		for (int i=0; i<m-1; i++)
		{
			for (int j=0; j<n-1; j++)
			{
				if (matrix[i][j])		// if occupied ...
				{
					int up = matrix[get_up(i,m)][j];    //  look up  
					int left = matrix[i][get_left(j,n)];  //  look left  
					
					switch (!!up + !!left)
					{
						case 0:
							
							break;
						case 1:                              // part of an existing cluster
							if (up) {matrix[i][j] = uf_union(up,matrix[i][j]);}       // whichever is nonzero is labelled
							if (left) {matrix[i][j] = uf_union(left,matrix[i][j]);}
							break;
						case 2:                              // this site binds two clusters
							matrix[i][j] = uf_union(left,up);
							break;
					}
				}
			}
		}
		int *new_labels_pbc = calloc(sizeof(int), n_labels); // allocate array, initialized to zero
		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
			{
				if (matrix[i][j])
				{
					int x = uf_find(matrix[i][j]);
					if (new_labels_pbc[x] == 0)
					{
						new_labels_pbc[0]++;
						new_labels_pbc[x] = new_labels_pbc[0];
					}
					matrix[i][j] = new_labels_pbc[x];
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

int get_down(int i, int m)
{
	return i%m;
}

int get_right(int j, int n)
{
	return j%n;
}

//----------------------------------------------------------

void check_labelling(int boundary, int **matrix, int m, int n)
{
	int N,S,E,W;
	if(boundary == 0)
	{
		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
			{
				if (matrix[i][j])
				{
					N = ( i==0 ? 0 : matrix[i-1][j] );
					S = ( i==m-1 ? 0 : matrix[i+1][j] );
					E = ( j==n-1 ? 0 : matrix[i][j+1] );
					W = ( j==0 ? 0 : matrix[i][j-1] );
					
					assert( N==0 || matrix[i][j]==N );
					assert( S==0 || matrix[i][j]==S );
					assert( E==0 || matrix[i][j]==E );
					assert( W==0 || matrix[i][j]==W );
      			}
			}
		}
	}
	else if(boundary == 1)
	{
		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
			{
				if (matrix[i][j])
				{
					N = ( matrix[get_up(i,m)][j] );
					S = ( matrix[get_down(i,m)][j] );
					E = ( matrix[i][get_right(j,n)] );
					W = ( matrix[i][get_left(j,n)] );
					
					assert( N==0 || matrix[i][j]==N );
					assert( S==0 || matrix[i][j]==S );
					assert( E==0 || matrix[i][j]==E );
					assert( W==0 || matrix[i][j]==W );
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

8 8
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
 --output-- 
  1   1   1   1   1   1   1   1 
  0   0   0   0   0   0   0   1 
  1   0   0   0   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   0   1   0   1 
  1   0   0   1   1   1   0   1 
  1   1   1   1   0   0   0   1 
  0   0   0   1   1   1   0   1 
HK reports 1 clusters in pbc

10 15
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
HK reports 7 clusters in pbc
*/

