/*wo way implement the IDCT*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#pragma warning (disable:4996)

//Y=At*X
#define S 8
#define H 4
/*IDCT Coefficients*/
#define a 125
#define b 118
#define c 106
#define d 90
#define e 71
#define f 48
#define g 24

void acc1(void);
int Ya[S][S];
/*plus acc*/
int p[H][H] = { { d, b, d, f },
{ d, f, -d, -b },
{ d, -f, -d, b },
{ d, -b, d, -f } };
/*minus acc*/
int m[H][H] = { { a, c, e, g },
{ c, -g, -a, -e },
{ e, -a, g, c },
{ g, -e, c, -a } };
/*Coefficients Matrix of 8x8 2D IDCT*/
int At[S][S] = { { d, a, b, c, d, e, f, g },
{ d, c, f, -g, -d, -a, -b, -e },
{ d, e, -f, -a, -d, g, b, c },
{ d, g, -b, -e, d, c, -f, -a },
{ d, -g, -b, e, d, -c, -f, a },
{ d, -e, -f, a, -d, -g, b, -c },
{ d, -c, f, g, -d, a, -b, e },
{ d, -a, b, -c, d, -e, f, -g } };
/*input matrix*/
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
	acc1();
	clock_t start, end;
	double duringtime;
	int M1[S][S], Y[S][S], Yt[S][S], A[S][S];
	int i, j, k, M = 0, acc = 0;
	start = clock();
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
	//fprintf(fp,"NO Accelerate duringtiime:%f s\n",duringtime);
	fclose(fp);
	//Don't put the counting time in front of "fclose()" cause can't count using time
	end = clock();
	duringtime = (double)(end - start);
	printf("NO Accelerate IDCT\n\tusing:%f ms\n", duringtime);
	printf("output alreadly done...\n");
	system("pause");
	return 0;
}

void acc1(void)
{
	printf("Now excute IDCT Accelerate\n");
	/*----------define---------*/
	clock_t start, end;
	double waste;
	int i, j, k = 0, k1 = 1, M = 0, M1 = 0, acc = 0, acc1 = 0, E;
	int I[H], I1[H], Iu[H], Id[H], Iu1[H], Id1[H], Iu2[H], Id2[H], Iu3[H], Id3[H], Iu4[H], Id4[H], Iu5[H], Id5[H], Iu6[H], Id6[H], Iu7[H], Id7[H];
	int Macc[S][S],Macct[S][S],Ya[S][S];
	int MM=0, acc2 = 0;
	start = clock();
	FILE *fp;
	fp = fopen("acc.xls", "w+");
	/*------Accelerate caululate------*/
	for (E = 0; E < S; E++)
	{
		for (i = 0; i<H; i++)
		{
			for (j = 0; j<H; j++)
			{
				M = p[i][j] * X[j + k][E];
				M1 = m[i][j] * X[j + k1][E];
				acc += M;
				acc1 += M1;
				k++;
				k1++;
				if (j == H - 1)
				{
					switch (E)
					{
						case 0:
						{
							I[i] = acc;
							I1[i] = acc1;
							Iu[i] = I[i] + I1[i];
							Id[i] = I[i] - I1[i];
							break;
						}
						case 1:
						{
							I[i] = acc;
							I1[i] = acc1;
							Iu1[i] = I[i] + I1[i];
							Id1[i] = I[i] - I1[i];
							break;
						}
						case 2:
						{
							I[i] = acc;
							I1[i] = acc1;
							Iu2[i] = I[i] + I1[i];
							Id2[i] = I[i] - I1[i];
							break;
						}
						case 3:
						{
							I[i] = acc;
							I1[i] = acc1;
							Iu3[i] = I[i] + I1[i];
							Id3[i] = I[i] - I1[i];
							break;
						}
						case 4:
						{
							I[i] = acc;
							I1[i] = acc1;
							Iu4[i] = I[i] + I1[i];
							Id4[i] = I[i] - I1[i];
							break;
						}
						case 5:
						{
							I[i] = acc;
							I1[i] = acc1;
							Iu5[i] = I[i] + I1[i];
							Id5[i] = I[i] - I1[i];
							break;
						}
						case 6:
						{
							I[i] = acc;
							I1[i] = acc1;
							Iu6[i] = I[i] + I1[i];
							Id6[i] = I[i] - I1[i];
							break;
						}
						case 7:
						{
							I[i] = acc;
							I1[i] = acc1;
							Iu7[i] = I[i] + I1[i];
							Id7[i] = I[i] - I1[i];
							break;
						}
						default:
							break;
					}
					acc = 0;
					acc1 = 0;
					k = 0;
					k1 = 1;
				}
			}
		}
	}
		/*-------------sorting to matrix M-------------*/
		int temp=0;
		for (i = 0; i < S; i++)
		{
			k = 1;
			for (j = 0; j < S; j++)
			{
				switch (i)
				{
					case 0:
					{
						if (j < 4)
							Macc[j][i] = Iu[j];
						else if (j>=4)
						{
							Macc[j][i] = Id[j - k];
							k += 2;
						}
						break;
					}
					case 1:
					{
							  if (j < 4)
								  Macc[j][i] = Iu1[j];
							  else if (j >= 4)
							  {
								  Macc[j][i] = Id1[j - k];
								  k += 2;
							  }
							  break;
					}
					case 2:
					{
							  if (j < 4)
								  Macc[j][i] = Iu2[j];
							  else if (j >= 4)
							  {
								  Macc[j][i] = Id2[j - k];
								  k += 2;
							  }
							  break;
					}
					case 3:
					{
							  if (j < 4)
								  Macc[j][i] = Iu3[j];
							  else if (j >= 4)
							  {
								  Macc[j][i] = Id3[j - k];
								  k += 2;
							  }
							  break;
					}
					case 4:
					{
							  if (j < 4)
								  Macc[j][i] = Iu4[j];
							  else if (j >= 4)
							  {
								  Macc[j][i] = Id4[j - k];
								  k += 2;
							  }
							  break;
					}
					case 5:
					{
							  if (j < 4)
								  Macc[j][i] = Iu5[j];
							  else if (j >= 4)
							  {
								  Macc[j][i] = Id5[j - k];
								  k += 2;
							  }
							  break;
					}
					case 6:
					{
							  if (j < 4)
								  Macc[j][i] = Iu6[j];
							  else if (j >= 4)
							  {
								  Macc[j][i] = Id6[j - k];
								  k += 2;
							  }
							  break;
					}
					case 7:
					{
							  if (j < 4)
								  Macc[j][i] = Iu7[j];
							  else if (j >= 4)
							  {
								  Macc[j][i] = Id7[j - k];
								  k += 2;
							  }
							  break;
					}
					default:
						break;
				}
			}
		}
		/*printf result*/
		fprintf(fp, "Array  Macc\n");
		for (i = 0; i<S; i++)
		{
			for (j = 0; j<S; j++)
			{
				fprintf(fp, "%d\t", Macc[i][j]);
				if (j == S - 1)
					fprintf(fp, "\n");
			}
		}
		/*---------A transpot--------*/
		for (i = 0; i < S; i++)
		{
			for (j = 0; j < S; j++)
			{
				Macct[j][i] = Macc[i][j];
			}
		}
		/*--------Y=Macct*At-------*/
		for (i = 0; i<S; i++)
		{
			for (j = 0; j<S; j++)
			{
				for (k = 0; k<S; k++)
				{
					 MM= At[i][k] * Macct[k][j];
					acc2 += MM;
					if (k == S - 1)
					{
						Ya[i][j] = acc2 >> 16;
						acc2 = 0;
					}
				}
			}
		}
		/*--------output Ya-------*/
		fprintf(fp, "Array  Ya\n");
		for (i = 0; i<S; i++)
		{
			for (j = 0; j<S; j++)
			{
				fprintf(fp, "%d\t", Ya[i][j]);
				if (j == S - 1)
					fprintf(fp, "\n");
			}
		}
		fprintf(fp, "\n");
		fclose(fp);
		end = clock();
		waste = (double)(end - start);
		printf("\tAccelerate using time: %f ms\n",waste);
		puts("\tLeave Accelarate function");
		puts("*****************************");
	}
		
