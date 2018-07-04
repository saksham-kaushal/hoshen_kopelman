#include<stdio.h>
#include<stdlib.h>

int main(void)
{
	printf("%d\n",(-1%5));
	printf("%d\n",(1%5));
	int *arr1;
	arr1 = (int *)calloc(5, sizeof(int));
	for(int i=0;i<5;i++)
	{
		arr1[i] = i;
	}
	int *arr2;
	arr2 = (int *)calloc(5, sizeof(int));
	for(int i=0;i<5;i++)
	{
		arr2[i] = i*i;
	}
	for(int i=0;i<5;i++)
	{
		printf("%d",arr2[i]);
	}
	printf("\n");
	arr2=arr1;
	for(int i=0;i<5;i++)
	{
		printf("%d",arr2[i]);
	}
	printf("\n");
	free(arr1);
	free(arr2);
	return 0;
}
