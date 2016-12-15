//IDCT
//IDCT C code
/*
  20161215 Lab3 
  Y=At*X*A
  At->A matrix transpot
  X->input num
  A->A matrix
  Y->IDCT Result output matrix
*/
#include <stdio.h>
#include <stdlib.h>
#pragma warning (disable:4996)
//Y=At*X
#define S 8
/*IDCT Coefficients*/
#define a 125
#define b 118
#define c 106
#define d 90
#define e 71
#define f 48
#define g 24
/*Coefficients Matrix of 8x8 2D IDCT*/
int At[S][S] = { { d, a, b, c, d, e, f, g },
{ d, c, f, -g, -d, -a, -b, -e },
{ d, e, -f, -a, -d, g, b, c },
{ d, g, -b, -e, d, c, -f, -a },
{ d, -g, -b, e, d, -c, -f, a },
{ d, -e, -f, a, -d, -g, b, -c },
{ d, -c, f, g, -d, a, -b, e },
{ d, -a, b, -c, d, -e, f, -g } };
int X[S][S] = { { 224, -380, -403, -220, 0, 135, 72, 0 },
{ 100, -100, 55, 0, 0, 0, 0, 0 },
{ 71, -55, -32, 67, 0, 0, 0, 0 },
{ 55, -27, -32, 33, 0, -42, 0, 0 },
{ -27, 0, 101, -72, 0, 43, 0, 0 },
{ 32, -33, 0, 40, 0, 0, 0, 0 },
{ -97, 101, -36, 0, 47, 0, 0, 0 },
{ -67, 72, -43, 0, 0, 0, 0, 1 } };
int main()
{
	int M1[S][S], Y[S][S], Yt[S][S], A[S][S];
	int i, j, k, M = 0, acc = 0;
	FILE*fp;
	fp = fopen("martix.xls", "w+");
	/*----------M1=AT*X---------*/
	for (i = 0; i<S; i++)
	{
		for (j = 0; j<S; j++)
		{
			for (k = 0; k<S; k++)
			{
				M = At[i][k] * X[k][j];
				acc += M;
				if (k == S - 1)
				{
					M1[i][j] = acc >> 8;
					acc = 0;
				}
			}
		}
	}
	acc = 0;
	M = 0;
	/*----------outputfile-----M-----*/
	fprintf(fp, "Array M\n");
	for (i = 0; i<S; i++)
	{
		for (j = 0; j<S; j++)
		{
			fprintf(fp, "%d\t", M1[i][j]);
			if (j == S - 1)
				fprintf(fp, "\n");
		}
	}
	fprintf(fp, "\n");
	/*--------Transpot---------*/
	for (i = 0; i < S; i++)
	{
		for (j = 0; j < S; j++)
		{
			A[j][i] = At[i][j];
		}
	}
	/*----------outputfile-----A-----*/
	fprintf(fp, "Array A\n");
	for (i = 0; i<S; i++)
	{
		for (j = 0; j<S; j++)
		{
			fprintf(fp, "%d\t", A[i][j]);
			if (j == S - 1)
				fprintf(fp, "\n");
		}
	}
	fprintf(fp, "\n");
	/*----------Y=M1*A---------*/
	for (i = 0; i<S; i++)
	{
		for (j = 0; j<S; j++)
		{
			for (k = 0; k<S; k++)
			{
				M = M1[i][k] * A[k][j];
				acc += M;
				if (k == S - 1)
				{
					Yt[i][j] = acc >> 8;
					acc = 0;
				}
			}
		}
	}
	/*--------Transpot---------*/
	for (i = 0; i < S; i++)
	{
		for (j = 0; j < S; j++)
		{
			Y[j][i] = Yt[i][j];
		}
	}
	/*----------outputfile-----Y-----*/
	fprintf(fp, "Array Y\n");
	for (i = 0; i<S; i++)
	{
		for (j = 0; j<S; j++)
		{
			fprintf(fp, "%d\t", Y[i][j]);
			if (j == S - 1)
				fprintf(fp, "\n");
		}
	}
	fprintf(fp, "\n");
	fclose(fp);
	printf("output alreadly done...\n");
	system("pause");
	return 0;
}
