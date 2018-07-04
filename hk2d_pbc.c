/* Tobin Fricke's implementation of the
   Hoshen-Kopelman algorithm for
   cluster labeling.

   Copyright (c) September 9, 2000, by Tobin Fricke <tobin@splorg.org>
   Distributed under the terms of the GNU Public License.

   Modified 2002-03-09 Tobin Fricke
   Modified substantially 2004-04-21 by Tobin Fricke

   This program is written in the 1999 standard of the C language (C99).  Older C
   compilers will refuse to compile it.   You can use a C++ compiler, a C99 compiler,
   or you can modify this code to comply with a previous version of the C standard.
   The GCC compiler supports C99 as of version 3.0.  Compile this program with:

   gcc-3.0 -Wall -std=c99 hk.c -o hk

   http://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
/* #include "hk.h" */

/* Implementation of Union-Find Algorithm */


/* The 'labels' array has the meaning that labels[x] is an alias for the label x; by
   following this chain until x == labels[x], you can find the canonical name of an
   equivalence class.  The labels start at one; labels[0] is a special value indicating
   the highest label already used. */

int *labels;
int  n_labels = 0;     /* length of the labels array */

/*  uf_find returns the canonical label for the equivalence class containing x */

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

/*  uf_union joins two equivalence classes and returns the canonical label of the resulting class. */

int uf_union(int x, int y) {
  return labels[uf_find(x)] = uf_find(y);
}

/*  uf_make_set creates a new equivalence class and returns its label */

int uf_make_set(void) {
  labels[0] ++;
  assert(labels[0] < n_labels);
  labels[labels[0]] = labels[0];
  return labels[0];
}

/*  uf_intitialize sets up the data structures needed by the union-find implementation. */

void uf_initialize(int max_labels) {
  n_labels = max_labels;
  labels = calloc(sizeof(int), n_labels);
  labels[0] = 0;
}

/*  uf_done frees the memory used by the union-find data structures */

void uf_done(void) {
  n_labels = 0;
  free(labels);
  labels = 0;
}

/* End Union-Find implementation */

#define max(a,b) (a>b?a:b)
#define min(a,b) (a>b?b:a)

/* print_matrix prints out a matrix that is set up in the "pointer to pointers" scheme
   (aka, an array of arrays); this is incompatible with C's usual representation of 2D
   arrays, but allows for 2D arrays with dimensions determined at run-time */

void print_matrix(int **matrix, int m, int n) {
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++)
      printf("%3d ",matrix[i][j]);
    printf("\n");
  }
}

void print_aug_matrix(int **matrix, int m, int n) {
  for (int i=1; i<m; i++) {
    for (int j=1; j<n; j++)
      printf("%3d ",matrix[i][j]);
    printf("\n");
  }
}

/* Label the clusters in "matrix".  Return the total number of clusters found. */

int hoshen_kopelman(int **matrix, int m, int n) {
  
  uf_initialize(m * n / 2);
  
  /* scan the matrix */
  
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      if (matrix[i][j]) {                        // if occupied ...

	int up = (i==0 ? 0 : matrix[i-1][j]);    //  look up  
	int left = (j==0 ? 0 : matrix[i][j-1]);  //  look left
	
	switch (!!up + !!left) {
	  
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
      
  
  /* apply the relabeling to the matrix */

  /* This is a little bit sneaky.. we create a mapping from the canonical labels
     determined by union/find into a new set of canonical labels, which are 
     guaranteed to be sequential. */
  
//      print_matrix(matrix,m,n);
//      printf("\n");
  int *new_labels = calloc(sizeof(int), n_labels); // allocate array, initialized to zero
  
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      if (matrix[i][j]) {
	int x = uf_find(matrix[i][j]);
	if (new_labels[x] == 0) {
	  new_labels[0]++;
	  new_labels[x] = new_labels[0];
	}
	matrix[i][j] = new_labels[x];
      }
 
  int total_clusters = new_labels[0];

  free(new_labels);
  uf_done();
  return total_clusters;
}

/* This procedure checks to see that any occupied neighbors of an occupied site
   have the same label. */

void check_labelling(int **matrix, int m, int n) {
  int N,S,E,W;
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      if (matrix[i][j]) {
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

//----------------------------------------------------------

int main(int argc, char **argv) {

  int m,n;
  int **matrix;

  /* Read in the matrix from standard input

     The whitespace-deliminated matrix input is preceeded
     by the number of rows and number of columns */

  if (scanf("%d %d",&m,&n) == 2) {  // m = rows, n = columns
    
    matrix = (int **)calloc(m, sizeof(int*));
    
    for (int i=0; i<m; i++) {
      matrix[i] = (int *)calloc(n, sizeof(int));
      for (int j=0; j<n; j++)
	scanf("%d",&(matrix[i][j]));
    }
    
    printf(" --input-- \n");
    
    print_matrix(matrix,m,n);
    
    
    /* Process the matrix */
    
    int clusters = hoshen_kopelman(matrix,m,n);
//    check_labelling(matrix,m,n);
//	pbc addition

    int **aug_matrix;
    aug_matrix = (int **)calloc(m, sizeof(int*));
    for (int i=0; i<m+1; i++) {
      aug_matrix[i] = (int *)calloc(n, sizeof(int));
    }
    aug_matrix[0][0]=0;
    for (int i=1;i<n+1;i++)
    {
    	aug_matrix[0][i] = matrix[m-1][i-1];
    }
    for (int j=1;j<m+1;j++)
    {
    	aug_matrix[j][0] = matrix[j-1][n-1];
    }
    for (int ii=1;ii<m+1;ii++)
    {
    	for (int jj=1;jj<n+1;jj++)
    	{
    		aug_matrix[ii][jj] = matrix[ii-1][jj-1];
    	}
    }   
    
    /* Output the result */
    
    printf(" --augmented matrix-- \n");
    print_matrix(aug_matrix,m+1,n+1);
    printf(" --output-- \n");
    int new_clusters = hoshen_kopelman(aug_matrix,m+1,n+1);
    print_matrix(aug_matrix,m+1,n+1);
    printf("\n\n");

/*    for(int i=0;i<m+1;i++)
    {
    	for(int j=0;j<n+1;j++)
    	{
    		printf("%d  ",aug_matrix[i][j]);
    	}
    	printf("\n");
    }
*/
    check_labelling(aug_matrix,m+1,n+1);
    
    printf("HK reports %d clusters found\n", new_clusters);

    for (int i=0; i<m; i++)
      free(matrix[i]);
    free(matrix);
    for (int i=0; i<m; i++)
      free(aug_matrix[i]);
    free(aug_matrix);
  }
  else {
  	printf("The form of the input is two integers giving the dimensions of the matrix, followed by the matrix elements (with data separated by whitespace).");
  }
  
  return 0;
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
  2   0   0   0   0   2   0   1 
  2   0   0   2   0   2   0   1 
  2   0   0   2   0   2   0   1 
  2   0   0   2   2   2   0   1 
  2   2   2   2   0   0   0   1 
  0   0   0   2   2   2   0   1 
HK reports 2 clusters found

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
  0   0   2   0   0   4   0   0   3   3   3   0   0   0   5 
  0   2   2   2   2   0   0   0   0   3   0   5   5   5   5 
  6   0   2   2   0   7   0   0   0   0   5   5   0   0   5 
  0   0   2   2   0   0   0   0   5   0   0   0   8   0   5 
  9   0   2   0  10   0   0   0   5   0   0   5   0   5   5 
  9   9   0   0  10   0   5   5   5   5   5   5   5   5   0 
  9   0   0  10  10  10   0   0   5   5   5   0   0   5   5 
  9   9   0  10  10  10   0  11   0   5   0   0  12   0   0 
  0   0  13   0   0  10  10   0   5   5   0  12  12   0   0 
HK reports 13 clusters found


*/


