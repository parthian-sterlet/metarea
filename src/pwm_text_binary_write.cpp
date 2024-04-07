#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define MATLEN 50 // max length of motives
#define OLIGNUM 4// di 16 mono 4
double pwm[MATLEN][OLIGNUM];
double pfm[MATLEN][OLIGNUM];

void DelChar(char* str, char c)
{
	int i, lens, size;

	size = 0;
	lens = strlen(str);
	for (i = 0; i < lens; i++)
	{
		if (str[i] != c)str[size++] = str[i];
	}
	str[size] = '\0';
}
int StrNStr(char* str, char c, int n)
{
	int i, len = (int)strlen(str);
	int k = 0;
	for (i = 0; i < len; i++)
	{
		if (str[i] == c)
		{
			k++;
			if (k == n)return i;
		}
	}
	return -1;
}
int UnderStol(char* str, int nstol, char* ret, size_t size, char sep)
{
	memset(ret, 0, size);
	int p1, p2, len;
	if (nstol == 0)
	{
		p2 = StrNStr(str, sep, 1);
		if (p2 == -1)p2 = strlen(str);
		strncpy(ret, str, p2);
		ret[p2] = '\0';
		return 1;
	}
	else
	{
		p1 = StrNStr(str, sep, nstol);
		p2 = StrNStr(str, sep, nstol + 1);
		if (p2 == -1)
		{
			p2 = strlen(str);
		}
		if (p1 == -1 || p2 == -1) return -1;
		len = p2 - p1 - 1;
		strncpy(ret, &str[p1 + 1], len);
		ret[len] = '\0';
		return 1;
	}
}
int main(int argc, char* argv[])
{
	char d[300], s[300], head[300], file_out_log[300], file_in_pfm[300], file_in_pwm[300], file_in_tab[300], file_out_distb[300];
	int i, j, k;
	FILE* out_log, * in_pfm, *in_pwm, *in_tab, * out_distb;
	if (argc != 6)
	{
		printf("%s 1file in_pfm 2file in_pwm 3file in tab 4file out_binary 5file out_log\n", argv[0]);		
		return -1;
	}
	strcpy(file_in_pfm, argv[1]);
	strcpy(file_in_pwm, argv[2]);
	strcpy(file_in_tab, argv[3]);
	strcpy(file_out_distb, argv[4]);
	strcpy(file_out_log, argv[5]);

	if ((in_pfm = fopen(file_in_pfm, "rt")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_in_pfm);
		return -1;
	}
	int alfabet = 0;
	fgets(head, sizeof(head), in_pfm);
	fgets(head, sizeof(head), in_pfm);
	//homer, cis-bp
	int nseq = 1000000000;
	int shift_col=0;		
	{
		DelChar(head, '\n');
		DelChar(head, '\r');
		int headlen = strlen(head);
		if (head[headlen - 1] == '\t')
		{
			head[headlen - 1] = '\0';
			headlen--;
		}
		int shift_col, counttab = 0;
		for (i = 0; i < headlen; i++)
		{
			if (head[i] == '\t')counttab++;
		}
		if (counttab == 4 || counttab == 16)
		{
			shift_col = 1;
			alfabet = counttab;
		}
		else
		{
			if (counttab == 3 || counttab == 15)
			{
				shift_col = 0;
				alfabet = counttab + 1;
			}
			else
			{
				printf("Reading errorin Cisbp/Homer matrix file %s !\n", file_in_pfm);
				return -1;
			}
		}
	}
	rewind(in_pfm);

	int len = 0;
	fgets(head, sizeof(head), in_pfm);
	for (i = 0; i < MATLEN; i++)
	{
		if (fgets(d, sizeof(d), in_pfm) != NULL)
		{
			char c = d[0];
			if (isdigit(c) || (strchr("-ATGC", c) != 0))len++;
		}
		else break;
	}
	rewind(in_pfm);
	char sep = '\t';
	fgets(head, sizeof(head), in_pfm);
	for (i = 0; i < len; i++)
	{
		if (fgets(d, sizeof(d), in_pfm) != NULL)
		{
			char c = d[0];
			if (isdigit(c) || (strchr("-ATGC", c) != 0))
			{
				//	 printf("%s",d[i]);						
				for (j = 0; j < alfabet; j++)
				{
					int test = UnderStol(d, j + shift_col, s, sizeof(s), sep);
					if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
						//double score=atof(s);
					//pfm[i][j]=nseq*score;
					pfm[i][j] = atof(s);
				}
			}
			else break;
		}
		else break;
	}
	fclose(in_pfm);
	if ((in_pwm = fopen(file_in_pwm, "rt")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_in_pwm);
		return -1;
	}
	char cm = '-';
	fgets(head, sizeof(head), in_pwm);
	for (i = 0; i < MATLEN; i++)
	{
		if (fgets(d, sizeof(d), in_pwm) != NULL)
		{
			char c = d[0];
			if (isdigit(c) || c ==cm)
			{
				//	 printf("%s",d[i]);						
				for (j = 0; j < alfabet; j++)
				{
					int test = UnderStol(d, j, s, sizeof(s), sep);
					if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
					//double score=atof(s);
					//pwm[i][j]=nseq*score;
					pwm[i][j] = atof(s);
				}
			}
			else break;
		}
		else break;
	}
	if ((in_tab = fopen(file_in_tab, "rt")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_in_tab);
		return -1;
	}
	int n_thresh = 0;
	//fgets(d, sizeof(d), in_tab);//header
	while (fgets(d, sizeof(d), in_tab) != NULL)
	{
		char c = d[0];		
		if (c == '-' || isdigit(c))
		{
			char s[30];
			int test = UnderStol(d, 1, s, sizeof(s), sep);
			if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
			n_thresh++;
			double fprx = atof(s);			
		}
	}
	rewind(in_tab);
	double* thr_dist, * fpr_dist;
	thr_dist = new double[n_thresh];
	if (thr_dist == NULL) { puts("Out of memory..."); return -1; }
	fpr_dist = new double[n_thresh];
	if (fpr_dist == NULL) { puts("Out of memory..."); return -1; }
	k = 0;
	while (fgets(d, sizeof(d), in_tab) != NULL)
	{
		char c = d[0];
		if (c == '-' || isdigit(c))
		{
			char s[30];
			char sep = '\t';
			int test = UnderStol(d, 1, s, sizeof(s), sep);
			if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
			thr_dist[k] = atof(d);
			fpr_dist[k] = atof(s);			
			k++;
		}
	}
	fclose(in_tab);

	if ((out_distb = fopen(file_out_distb, "wb")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_out_distb);
		return -1;
	}
	int len4 = 4 * len;
	fwrite(&len, sizeof(int), 1, out_distb);
	fwrite(pfm, sizeof(double), len4, out_distb);
	fwrite(pwm, sizeof(double), len4, out_distb);
	fwrite(&n_thresh, sizeof(int), 1, out_distb);
	fwrite(thr_dist, sizeof(double), n_thresh, out_distb);
	fwrite(fpr_dist, sizeof(double), n_thresh, out_distb);
	fclose(out_distb);	
	fclose(out_distb);
	if ((out_log = fopen(file_out_log, "at")) == NULL)
	{
		fprintf(out_log, "Input file %s can't be opened!\n", file_out_log);
		exit(1);
	}
	fprintf(out_log, "%d\t%d\t%.18f\t%.18f\n", len, n_thresh, thr_dist[0], fpr_dist[0]);
	fclose(out_log);
	return 0;
}
