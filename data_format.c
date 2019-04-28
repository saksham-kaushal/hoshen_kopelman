#include <stdio.h>

int main(void)
{
	FILE* inputfp;
	inputfp = fopen("./input","r");
	int l,m,n;
	fscanf(inputfp,"%d %d %d", &l, &m, &n);
	int arr[l][m][n];
	printf("%d %d %d\n",l,m,n);
	int i=0,j=0,k=0,x=0;
	for(i=0;i<l;i++)
	{
		for(j=0;j<m;j++)
		{
			for(k=0;k<n;k++)
			{
				fscanf(inputfp,"%d",&x);
				arr[i][j][k]=(2*x)-1;
				int neg=1,a=0;
				for(a=0;a<i+j+k;a++)
				{
					neg*=(-1);
				}
				printf("%d ",((arr[i][j][k]*neg+1)/2));
															
// Uncomment the following section to label the clusters of unoccupied sites. Also, comment the line above for the same.
/*				if (arr[i][j][k]*neg==1)
				{
					printf("0 ");
				}
				else if (arr[i][j][k]*neg==-1)
				{
					printf("1 ");
				}
*/
			}
			printf("\n");
		}
	}
	fclose(inputfp);
	return 0;
}
