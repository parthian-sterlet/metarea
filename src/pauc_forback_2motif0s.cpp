#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 12000
#define MATLEN 50 // max length of motives
#define OLIGNUM 4// di 16 mono 4
#define ARGLEN 300 //max argv length
#define MAXPAR 1500 //max no. of motifs
#define DIM 100

struct juxt
{
	double err;
	int sta;
	int ind;
	int mod;
};
int compare_juxt(const void* X1, const void* X2)//increase
{
	struct juxt* S1 = (struct juxt*)X1;
	struct juxt* S2 = (struct juxt*)X2;
	if (S1->err - S2->err > 0)return 1;
	if (S1->err - S2->err < 0)return -1;
	return 0;
}

struct acc {
	double auc;
	int rank;
	int nth;
	int num;
	long p;
};
int StrNStr(char* str, char c, int n)
{
	int i, len = (int)strlen(str);
	int k = 1;
	for (i = 0; i < len; i++)
	{
		if (str[i] == c)
		{
			if (k == n)return i;
			k++;
		}
	}
	return -1;
}
void DelHole(char* str)
{
	char* hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
}
double UnderStol(char* str, int nstol, char razd)
{
	if (nstol == 0)return atof(str);
	char ret[100];
	memset(ret, 0, sizeof(ret));
	int p1 = StrNStr(str, razd, nstol);
	int p2 = StrNStr(str, razd, nstol + 1);
	if (p2 == -1)
	{
		p2 = (int)strlen(str);
	}
	if (p1 == -1 || p2 == -1) return -1;
	int len = p2 - p1 - 1;
	strncpy(ret, &str[p1 + 1], len);
	ret[len] = '\0';
	return atof(ret);
}
int UnderStolStr(char* str, int nstol, char* ret, size_t size, char sep)
{
	memset(ret, '\0', size);
	int p1, p2, len;
	if (nstol == 0)
	{
		p2 = StrNStr(str, sep, 1);
		if (p2 == -1)p2 = (int)strlen(str);
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
			p2 = (int)strlen(str);
		}
		if (p1 == -1 || p2 == -1) return -1;
		len = p2 - p1 - 1;
		strncpy(ret, &str[p1 + 1], len);
		ret[len] = '\0';
		return 1;
	}
}
struct due {
	double buf;
	int sta;
	int end;
	int num;
	void get_copy(due* a);
	//	void print_all(void);
};
void due::get_copy(due* a)
{
	a->num = num;
	a->sta = sta;
	a->buf = buf;
	a->end = end;
};
/*
void due::print_all(void)
{
printf("[%d;%d]%s\t", sta, end, s[num].oli);
}*/
//set of dinucleotides
struct city {
	char site[300];
	int size;
	int len;
	double min;
	double raz;
	struct due tot[DIM];
	void get_copy(city* a);
	void sort_all(void);
	int get_file(char* file);
	//void city::fprint_tab(char *file);
}sta;
int city::get_file(char* file)
{
	FILE* in;
	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", file);
		return -1;
	}
	char d[300];
	fgets(d, sizeof(d), in);
	DelHole(d);
	strcpy(site, d);
	fgets(d, sizeof(d), in);
	size = atoi(d);
	fgets(d, sizeof(d), in);
	len = atoi(d);
	fgets(d, sizeof(d), in);
	min = atof(d);
	fgets(d, sizeof(d), in);
	raz = atof(d);
	char sep = '\t', s[30];
	int i, test;
	for (i = 0; i < size; i++)
	{
		fgets(d, sizeof(d), in);
		tot[i].sta = atoi(d);
		test = UnderStolStr(d, 1, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].end = atoi(s);
		test = UnderStolStr(d, 2, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].buf = atof(s);
		test = UnderStolStr(d, 3, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].num = atoi(s);
	}
	fclose(in);
	return 1;
}
void city::get_copy(city* a)
{
	strcpy(a->site, site);
	a->size = size;
	a->min = min;
	a->len = len;
	a->raz = raz;
	int i;
	for (i = 0; i < size; i++)
	{
		tot[i].get_copy(&a->tot[i]);
	}
}
int compare_due(const void* X1, const void* X2)
{
	struct due* S1 = (struct due*)X1;
	struct due* S2 = (struct due*)X2;
	if (S1->sta - S2->sta > 0)return 1;
	if (S1->sta - S2->sta < 0)return -1;
	if (S1->end - S2->end > 0)return 1;
	if (S1->end - S2->end < 0)return -1;
	if (S1->num - S2->num > 0)return 1;
	if (S1->num - S2->num < 0)return -1;
	return 0;
}
void city::sort_all(void)
{
	qsort((void*)tot, size, sizeof(tot[0]), compare_due);
}

struct qbs {
	double err;//ERR score
	int tpr;
	int fpr;
	int fps;
};
int compare_qq(const void* X1, const void* X2)//increase
{
	double X = (*(double*)X1 - *(double*)X2);
	if (X > 0)return -1;
	if (X < 0)return 1;
	return 0;
}
int compare_qbs(const void* X1, const void* X2)//decrease
{
	struct qbs* S1 = (struct qbs*)X1;
	struct qbs* S2 = (struct qbs*)X2;
	if (S1->err - S2->err > 0)return -1;
	if (S1->err - S2->err < 0)return 1;
	return 0;
}
int compare_auc(const void* X1, const void* X2)//decrease
{
	struct acc* S1 = (struct acc*)X1;
	struct acc* S2 = (struct acc*)X2;
	if (S1->auc - S2->auc > 0)return -1;
	if (S1->auc - S2->auc < 0)return 1;
	return 0;
}
int compare_num(const void* X1, const void* X2)//decrease
{
	struct acc* S1 = (struct acc*)X1;
	struct acc* S2 = (struct acc*)X2;
	if (S1->num - S2->num > 0)return 1;
	if (S1->num - S2->num < 0)return -1;
	return 0;
}
char* TransStr(char* d)
{
	int i, c, lens;
	lens = (int)strlen(d);
	for (i = 0; i < lens; i++)
	{
		c = int(d[i]);
		if (c < 97) d[i] = char(c + 32);
		//else break;
	}
	return(d);
}
int CheckStr(char* file, char* d, int n, int print)
{
	int i, len, ret;
	len = (int)strlen(d);
	ret = 1;
	for (i = 0; i < len; i++)
	{
		int di = (int)d[i];
		if (strchr("atgcATGC", di) != NULL)continue;
		if (strchr("nN", di) != NULL)
		{
			ret = 0; continue;
		}
		if (print == 1)printf("File %s; sequence %d position %d (%c) bad. Sequence deleted!\n", file, n, i + 1, d[i]);
		ret = -1;
		break;
	}
	return(ret);
}
int IdeLet(char c)
{
	int ret;
	switch (c) {
	case 'a': ret = 0; break;
	case 'c': ret = 1; break;
	case 'g': ret = 2; break;
	case 't': ret = 3; break;
	case 'n': ret = -1; break;
	default: ret = -2;
	}
	return(ret);
}
void EvalSeq(char* file, int& nseq, int olen, int len_peak_max)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int fl = 0;
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("EvalSeq! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	nseq = 0;
	int n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d);
			int check = CheckStr(file, d, n, 1);
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)nseq++;
			if (fl == -1)
			{
				fclose(in);
				break;
			}
			n++;
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelHole(l);
			strcpy(d, l);
			fl = 1; continue;
		}
		if (strlen(d) + strlen(l) > sizeof(d))
		{
			printf("Size is large...");
			printf("l:%s\nstrlen(l):%zu\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%zu\n", d, strlen(d));
			exit(1);
		}
		DelHole(l);
		strcat(d, l);
	}
}
void EvalLen(char* file, int* len, int olen, int len_peak_max)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int fl = 0;
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("EvalLen! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	int nn = 0;
	int n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d);
			int check = CheckStr(file, d, nn, 0);
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)len[n++] = lenx;
			nn++;
			if (fl == -1)
			{
				fclose(in);
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelHole(l);
			strcpy(d, l);
			fl = 1; continue;
		}
		if (strlen(d) + strlen(l) > sizeof(d))
		{
			printf("Size is large...");
			printf("l:%s\nstrlen(l):%zu\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%zu\n", d, strlen(d));
			exit(1);
		}
		DelHole(l);
		strcat(d, l);
	}
}
int ComplStr(char* d)
{
	int i, len;
	len = (int)strlen(d);
	char d1[SEQLEN];
	strcpy(d1, d);
	//	memset(d,0,sizeof(d));
	for (i = 0; i < len; i++)
	{
		switch (d1[len - i - 1])
		{
		case 'a': { d[i] = 't'; break; }
		case 't': { d[i] = 'a'; break; }
		case 'c': { d[i] = 'g'; break; }
		case 'g': { d[i] = 'c'; break; }
		default: d[i] = 'n';
		}
	}
	return 1;
}
void ReadSeq(char* file, int nseq, int* len, char*** seq_real, int olen, int len_peak_max)
{
	char l[SEQLEN], d[2][SEQLEN], head[400];
	int fl = 0, i, j;
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("ReadSeq! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	int size = 1;
	int nn = 0, n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d[0]);
			int check = CheckStr(file, d[0], nn, 0);
			nn++;
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)
			{
				TransStr(d[0]);
				d[0][len[n]] = '\0';
				strcpy(d[1], d[0]);
				ComplStr(d[1]);
				d[1][len[n]] = '\0';
				for (j = 0; j < 2; j++)
				{
					for (i = 0; i < lenx; i++)
					{
						seq_real[j][n][i] = d[j][i];
					}
					seq_real[j][n][lenx] = '\0';
				}
				n++;
			}
			else
			{
				if (lenx < olen)
				{
					printf("ReadSeq! Short peak %d (Len %d) ignored\n", nn + 1, lenx);
				}
				if (lenx > len_peak_max)
				{
					printf("ReadSeq! Long peak %d (Len %d) ignored\n", nn + 1, lenx);
				}
				if (check == -1)
				{
					printf("ReadSeq! Unusual symbol, peak %d ignored\n%s\n", nn + 1, d[0]);
				}
			}
			if (fl == -1)
			{
				fclose(in);
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d[0], 0, sizeof(d[0]));
			DelHole(l);
			strcpy(d[0], l);
			fl = 1; continue;
		}
		if (strlen(d[0]) + strlen(l) > sizeof(d[0]))
		{
			printf("Size is large...");
			printf("l:%s\nstrlen(l):%zu\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%zu\n", d[0], strlen(d[0]));
			exit(1);
		}
		DelHole(l);
		strcat(d[0], l);
	}
}
int SGA_rec_real_one(city* sta, int nthr_dist, double* thr_all, double* fpr_all, int* tpr, int nseq, char*** seq)
{
	int i, j, n, m;
	int compl1;
	char d[MATLEN];
	int nthr_dist1 = nthr_dist - 1;
	double thr_cr = thr_all[nthr_dist1];
	for (n = 0; n < nseq; n++)
	{
		//if ((n + 1) % 500 == 0)printf("\b\b\b\b\b\b\b%7d", n + 1);
		int len_pro1 = (int)strlen(seq[0][n]);
		int len21 = len_pro1 - sta->len;
		int index = nthr_dist;
		for (i = 0; i <= len21; i++)
		{
			double sco2 = 0;
			int gom = 0;
			for (compl1 = 0; compl1 < 2; compl1++)
			{
				int ista;
				if (compl1 == 0)ista = i;
				else ista = len21 - i;
				strncpy(d, &seq[compl1][n][ista], sta->len);
				d[sta->len] = '\0';
				if (strstr(d, "n") != NULL) { gom = 1; break; }
				double score = 0;
				for (j = 0; j < sta->size; j++)
				{
					int rlenj = (sta->tot[j].end - sta->tot[j].sta + 1);
					double fm = 0;
					for (m = sta->tot[j].sta; m <= sta->tot[j].end; m++)
					{
						int cod = 4 * IdeLet(d[m]) + IdeLet(d[m + 1]);
						if (sta->tot[j].num == cod) { fm++; }
					}
					if (fm != 0)
					{
						fm /= rlenj;
						score += sta->tot[j].buf * fm;
					}
				}
				score = (score - sta->min) / sta->raz;
				if (score > sco2)sco2 = score;
			}
			if (gom == 0)
			{
				if (sco2 >= thr_cr)
				{
					if (sco2 >= thr_all[0])
					{
						index = 0;
						break;
					}
					else
					{
						int index_here = nthr_dist;
						for (j = 1; j < nthr_dist; j++)
						{
							if (sco2 >= thr_all[j] && sco2 < thr_all[j - 1])
							{
								index_here = j;
								break;
							}
						}
						if (index_here < index)index = index_here;
					}
				}
			}
			if (index == 0)break;
		}
		for (j = index; j <= nthr_dist; j++)tpr[j]++;
	}
	return 1;
}
int SGA_rec_back_one(city* sta, int nthr_dist, double* thr_all, double* fpr_all, int* fpr, int* fp_nsites, int nseq, char*** seq, int& all_pos)
{
	int i, j, n, m;
	int compl1;
	char d[MATLEN];
	int nthr_dist1 = nthr_dist - 1;
	double thr_cr = thr_all[nthr_dist1];

	for (n = 0; n < nseq; n++)
	{
		//if ((n+1) % 500 == 0)printf("\b\b\b\b\b\b\b%7d", n+1);
		int len_pro1 = (int)strlen(seq[0][n]);
		int len21 = len_pro1 - sta->len;
		int index_best = nthr_dist;
		for (i = 0; i <= len21; i++)
		{
			double sco2 = 0;
			int gom = 0;
			for (compl1 = 0; compl1 < 2; compl1++)
			{
				int ista;
				if (compl1 == 0)ista = i;
				else ista = len21 - i;
				strncpy(d, &seq[compl1][n][ista], sta->len);
				d[sta->len] = '\0';
				if (strstr(d, "n") != NULL) { gom = 1; break; }
				double score = 0;
				for (j = 0; j < sta->size; j++)
				{
					int rlenj = (sta->tot[j].end - sta->tot[j].sta + 1);
					double fm = 0;
					for (m = sta->tot[j].sta; m <= sta->tot[j].end; m++)
					{
						int cod = 4 * IdeLet(d[m]) + IdeLet(d[m + 1]);
						if (sta->tot[j].num == cod) { fm++; }
					}
					if (fm != 0)
					{
						fm /= rlenj;
						score += sta->tot[j].buf * fm;
					}
				}
				score = (score - sta->min) / sta->raz;
				if (score > sco2)sco2 = score;
			}
			if (gom == 0)
			{
				all_pos++;
				int index = nthr_dist;
				if (sco2 >= thr_cr)
				{
					if (sco2 >= thr_all[0])
					{
						index = 0;
					}
					else
					{
						for (j = 1; j < nthr_dist; j++)
						{
							if (sco2 >= thr_all[j] && sco2 < thr_all[j - 1])
							{
								index = j;
								break;
							}
						}
					}
				}
				for (j = nthr_dist; j >= index; j--)fp_nsites[j]++;
				if (index < index_best)index_best = index;
			}
		}
		for (j = index_best; j <= nthr_dist; j++)fpr[j]++;
	}
	return 1;
}

int SGA_rec_real(city *sta, int nthr_dist[2], double** thr_all, double** fpr_all, int** tp_two, int nseq, char*** seq)
{
	int i, j, k, n, m;
	int compl1;
	char d[MATLEN];	
	int nthr_dist1[2], olen1[2];
	double thr_cr[2];
	for (n = 0; n < 2; n++)
	{
		nthr_dist1[n] = nthr_dist[n] - 1;
		olen1[n] = sta[n].len - 1;
		thr_cr[n] = thr_all[n][nthr_dist1[n]];
	}
	for (n = 0; n < nseq; n++)
	{
		//if ((n + 1) % 500 == 0)printf("\b\b\b\b\b\b\b%7d", n + 1);
		int best_inx[2];
		for (k = 0; k < 2; k++)
		{
			int index;
			int len_pro1 = (int)strlen(seq[0][n]);
			int len21 = len_pro1 - sta[k].len;
			index = best_inx[k] = nthr_dist[k];
			for (i = 0; i <= len21; i++)
			{
				double sco2 = 0;
				int gom = 0;
				for (compl1 = 0; compl1 < 2; compl1++)
				{
					int ista;
					if (compl1 == 0)ista = i;
					else ista = len21 - i;
					strncpy(d, &seq[compl1][n][ista], sta[k].len);
					d[sta[k].len] = '\0';
					if (strstr(d, "n") != NULL) { gom = 1; break; }
					double score = 0;
					for (j = 0; j < sta[k].size; j++)
					{
						int rlenj = (sta[k].tot[j].end - sta[k].tot[j].sta + 1);
						double fm = 0;
						for (m = sta[k].tot[j].sta; m <= sta[k].tot[j].end; m++)
						{
							int cod = 4 * IdeLet(d[m]) + IdeLet(d[m + 1]);
							if (sta[k].tot[j].num == cod) { fm++; }
						}
						if (fm != 0)
						{
							fm /= rlenj;
							score += sta->tot[j].buf * fm;
						}
					}
					score = (score - sta->min) / sta->raz;
					if (score > sco2)sco2 = score;
				}
				if (gom == 0)
				{
					if (sco2 >= thr_cr[k])
					{
						if (sco2 >= thr_all[k][0])
						{
							index = 0;
							break;
						}
						else
						{
							for (j = 1; j < nthr_dist[k]; j++)
							{
								if (sco2 >= thr_all[k][j] && sco2 < thr_all[k][j - 1])
								{
									index = j;
									break;
								}
							}
						}
					}
					if (index < best_inx[k])best_inx[k] = index;
					if (index == 0)break;
				}
			}
		}
		if (fpr_all[0][best_inx[0]] > fpr_all[1][best_inx[1]])tp_two[0][best_inx[0]]++;
		else tp_two[1][best_inx[1]]++;
	}
	/*
	int sum = 0;
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < nthr_dist[i]; j++)
		{
			sum += tp_two[i][j];
		}
	}*/
	return 1;
}
int SGA_rec_back(city *sta, int nthr_dist[2], double** thr_all, double** fpr_all, int** fp_two, int** fp_nsites, int nseq, char*** seq)
{
	int i, j, k, n, m;
	int compl1;
	char d[MATLEN];	
	int nthr_dist1[2], olen1[2];
	double thr_cr[2];
	for (n = 0; n < 2; n++)
	{
		nthr_dist1[n] = nthr_dist[n] - 1;
		olen1[n] = sta[n].len - 1;
		thr_cr[n] = thr_all[n][nthr_dist1[n]];
	}
	for (n = 0; n < nseq; n++)
	{
		//if ((n + 1) % 500 == 0)printf("\b\b\b\b\b\b\b%7d", n + 1);
		int best_inx[2];
		for (k = 0; k < 2; k++)
		{
			int len_pro1 = (int)strlen(seq[0][n]);
			int len21 = len_pro1 - sta[k].len;
			best_inx[k] = nthr_dist[k];
			for (i = 0; i <= len21; i++)
			{
				double sco2 = 0;
				int gom = 0;
				for (compl1 = 0; compl1 < 2; compl1++)
				{
					int ista;
					if (compl1 == 0)ista = i;
					else ista = len21 - i;
					strncpy(d, &seq[compl1][n][ista], sta[k].len);
					d[sta[k].len] = '\0';
					if (strstr(d, "n") != NULL) { gom = 1; break; }
					double score = 0;					
					for (j = 0; j < sta[k].size; j++)
					{
						int rlenj = (sta[k].tot[j].end - sta[k].tot[j].sta + 1);
						double fm = 0;
						for (m = sta[k].tot[j].sta; m <= sta[k].tot[j].end; m++)
						{
							int cod = 4 * IdeLet(d[m]) + IdeLet(d[m + 1]);
							if (sta[k].tot[j].num == cod) { fm++; }
						}
						if (fm != 0)
						{
							fm /= rlenj;
							score += sta->tot[j].buf * fm;
						}
					}
					score = (score - sta->min) / sta->raz;
					if (score > sco2)sco2 = score;
				}
				if (gom == 0)
				{
					int index = nthr_dist[k];
					if (sco2 >= thr_cr[k])
					{
						if (sco2 >= thr_all[k][0])
						{
							index = 0;
						}
						else
						{
							for (j = 1; j < nthr_dist[k]; j++)
							{
								if (sco2 >= thr_all[k][j] && sco2 < thr_all[k][j - 1])
								{
									index = j;
									break;
								}
							}
						}
					}
					for (j = nthr_dist[k]; j >= index; j--)fp_nsites[k][j]++;
					if (index < best_inx[k])best_inx[k] = index;
				}
				//	if (index == 0)break;				
			}
		}
		if (fpr_all[0][best_inx[0]] > fpr_all[1][best_inx[1]])fp_two[0][best_inx[0]]++;
		else fp_two[1][best_inx[1]]++;
	}
	/*
	int sum = 0;
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < nthr_dist[i]; j++)
		{
			sum += tp_two[i][j];
		}
	}*/
	return 1;
}
int main(int argc, char* argv[])
{
	int i, j, k, n, mot1, mot2;
	char*** seq_real, *** seq_back;
	char file_for[ARGLEN], file_back[ARGLEN], partner_db[ARGLEN];// , tfclass[80];path_fasta[ARGLEN], pfile_for[ARGLEN], pfile_back[ARGLEN],
	char file_mat[ARGLEN], file_auc[ARGLEN], file_prc[ARGLEN], file_log1[ARGLEN], file_log2[ARGLEN], file_prc1[ARGLEN];
	FILE* out_log1, * out_auc, * out_roc, * out_log2, * out_mat, * in_sga;

	if (argc != 11)
	{
		printf("%s 1,2input fasta foreground,background 3input binary files pwm_thresholds 4no. of motifs 5double ERRthresh 6,7,8,9,10output files pAUC_matrix, pAUC_list, log1, log2, PRcurves", argv[0]);
		return -1;
	}
	//strcpy(path_fasta, argv[1]);
	strcpy(file_for, argv[1]);
	strcpy(file_back, argv[2]);
	//strcpy(pfile_for, path_fasta);
	//strcpy(pfile_back, path_fasta);
	//strcat(pfile_for, file_for);
	//strcat(pfile_back, file_back);	
	strcpy(partner_db, argv[3]); //h12hs, h12mm
	int n_motifs = atoi(argv[4]);// no. of top-scoring matrices
	double fp2 = atof(argv[5]); //ERR threshold for pAUPRC		
	strcpy(file_mat, argv[6]);
	strcpy(file_auc, argv[7]);
	strcpy(file_log1, argv[8]);
	strcpy(file_log2, argv[9]);
	strcpy(file_prc, argv[10]);
	
	int* len_real, * len_back, nseq_real = 0, nseq_back = 0;
	int olen_min = 8;
	int len_peak_max = 3000;

	double fp2_lg = -log10(fp2);
	if ((out_auc = fopen(file_auc, "wt")) == NULL)
	{
		fprintf(out_auc, "Input file %s can't be opened!\n", file_auc);
		exit(1);
	}
	if ((out_log1 = fopen(file_log1, "wt")) == NULL)
	{
		fprintf(out_log1, "Input file %s can't be opened!\n", file_log1);
		exit(1);
	}
	fprintf(out_log1, "#Motif\tpAUC\n");
	if ((out_log2 = fopen(file_log2, "wt")) == NULL)
	{
		fprintf(out_log2, "Input file %s can't be opened!\n", file_log2);
		exit(1);
	}
	fprintf(out_log2, "#Motif 1\t#Motif 2\tRank 1\tRank 2\tSimilarity\tpAUC 1\tpAUC 2\tpAUC 1&2\tRatio\n");
	if ((out_mat = fopen(file_mat, "wt")) == NULL)
	{
		fprintf(out_mat, "Input file %s can't be opened!\n", file_mat);
		exit(1);
	}
	EvalSeq(file_for, nseq_real, olen_min, len_peak_max);
	EvalSeq(file_back, nseq_back, olen_min, len_peak_max);
	len_real = new int[nseq_real];
	if (len_real == NULL) { puts("Out of memory..."); exit(1); }
	len_back = new int[nseq_back];
	if (len_back == NULL) { puts("Out of memory..."); exit(1); }
	//printf("EvalLen\n");
	EvalLen(file_for, len_real, olen_min, len_peak_max);
	EvalLen(file_back, len_back, olen_min, len_peak_max);
	seq_real = new char** [2];
	if (seq_real == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		seq_real[i] = new char* [nseq_real];
		if (seq_real[i] == NULL) { puts("Out of memory..."); exit(1); }
		for (j = 0; j < nseq_real; j++)
		{
			seq_real[i][j] = new char[len_real[j] + 1];
			if (seq_real[i][j] == NULL) { puts("Out of memory..."); exit(1); }
		}
	}
	seq_back = new char** [2];
	if (seq_back == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		seq_back[i] = new char* [nseq_back];
		if (seq_back[i] == NULL) { puts("Out of memory..."); exit(1); }
		for (j = 0; j < nseq_back; j++)
		{
			seq_back[i][j] = new char[len_back[j] + 1];
			if (seq_back[i][j] == NULL) { puts("Out of memory..."); exit(1); }
		}
	}
	ReadSeq(file_for, nseq_real, len_real, seq_real, olen_min, len_peak_max);
	ReadSeq(file_back, nseq_back, len_back, seq_back, olen_min, len_peak_max);
	if ((in_sga = fopen(partner_db, "rb")) == NULL)
	{
		printf("Input file %s can't be opened!\n", partner_db);
		return -1;
	}
	double nseq_fb = (double)nseq_real / nseq_back;
	double prec_exp = 0.5;
	int nthr_dist[2];
	city* sta;
	sta = new city[2];	
	if (sta == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	int nthr_dist_max = 20000;
	double** thr_all;
	thr_all = new double* [2];
	if (thr_all == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		thr_all[i] = new double[nthr_dist_max];
		if (thr_all[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	double** fpr_all;
	fpr_all = new double* [2];
	if (fpr_all == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		fpr_all[i] = new double[nthr_dist_max];
		if (fpr_all[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	int n_motifs1 = n_motifs - 1;
	double** auctwei;
	auctwei = new double* [n_motifs1];
	if (auctwei == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (i = 0; i < n_motifs1; i++)
	{
		auctwei[i] = new double[n_motifs1 - i];
		if (auctwei[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	int* index_thr;
	index_thr = new int[n_motifs];
	if (index_thr == NULL) { puts("Out of memory..."); exit(1); }
	int* all_pos_thr;
	all_pos_thr = new int[n_motifs];
	if (all_pos_thr == NULL) { puts("Out of memory..."); exit(1); }
	int* all_pos;
	all_pos = new int[n_motifs];
	if (all_pos == NULL) { puts("Out of memory..."); exit(1); }
	double* fp_thr_rest;
	fp_thr_rest = new double[n_motifs];
	if (fp_thr_rest == NULL) { puts("Out of memory..."); exit(1); }
	//n_motifs = 50;
	int* tp_one, * fp_one;
	tp_one = new int[nthr_dist_max];
	if (tp_one == NULL) { puts("Out of memory..."); exit(1); }
	fp_one = new int[nthr_dist_max];
	if (fp_one == NULL) { puts("Out of memory..."); exit(1); }
	double* recall;
	recall = new double[nthr_dist_max];
	if (recall == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	double* prec;
	prec = new double[nthr_dist_max];
	if (prec == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	int** fp_nsites;// frequencies
	fp_nsites = new int* [2];
	if (fp_nsites == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		fp_nsites[i] = new int[nthr_dist_max];
		if (fp_nsites[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	int nmot_max = 1500;
	acc* motifp;
	motifp = new acc[n_motifs];
	if (motifp == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }	
	motifp[0].p = ftell(in_sga);
	for (mot1 = 0; mot1 < n_motifs; mot1++)all_pos[mot1] = all_pos_thr[mot1] = index_thr[mot1] = 0;
	for (mot1 = 0; mot1 < n_motifs; mot1++)fp_thr_rest[mot1] = 0;
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		//printf("%d\t%s\n", mot1 + 1, motif_class[mot1]);		
		size_t mu = fread(&sta[0], sizeof(sta[0]), 1, in_sga);
		if (mu != 1)
		{
			printf("Error reading file %s!\n", partner_db);
			return -1;
		}
		int nthr_dist_all;
		fread(&nthr_dist_all, sizeof(int), 1, in_sga);
		fread(thr_all[0], sizeof(double), nthr_dist_all, in_sga);
		fread(fpr_all[0], sizeof(double), nthr_dist_all, in_sga);
		if (mot1 != n_motifs1)motifp[mot1 + 1].p = ftell(in_sga);
		nthr_dist[0] = 0;
		for (i = 0; i < nthr_dist_all; i++)
		{
			nthr_dist[0]++;
			if (fpr_all[0][i] < fp2_lg)break;
		}
		motifp[mot1].nth = nthr_dist[0];
		/*if (strstr(motif_class[mot1], tfclass) == NULL)
		{
			motifp[mot1].auc = 0;
			continue;
		}		*/				
		for (j = 0; j <= nthr_dist[0]; j++)
		{
			tp_one[j] = fp_one[j] = fp_nsites[0][j] = 0;
		}
		//printf("Real %d\n",mot + 1);
		int sga_check = SGA_rec_real_one(&sta[0], nthr_dist[0], thr_all[0], fpr_all[0], tp_one, nseq_real, seq_real);
		if (sga_check == -1)
		{
			printf("One motif recognition error, real motif %d\n", mot1 + 1);
			exit(1);
		}
		//printf("Back %d\n", mot + 1);
		sga_check = SGA_rec_back_one(&sta[0], nthr_dist[0], thr_all[0], fpr_all[0], fp_one, fp_nsites[0], nseq_back, seq_back, all_pos[mot1]);
		if (sga_check == -1)
		{
			printf("One motif recognition error, back %d\n", mot1 + 1);
			exit(1);
		}
		/*for (i = 0; i <= nthr_dist[0]; i++)
		{
			printf("%f %f %d %d\t", thr_all[0][i],fpr_all[0][i],tp_one[i], fp_one[i]);
			if ((i + 1) % 4 == 0)printf("%d\n", i + 1);
		}*/
		//printf("ROC %d\n", mot + 1);
		all_pos_thr[mot1] = (int)(all_pos[mot1] * fp2);
		index_thr[mot1] = nthr_dist[0] - 1;
		int count_one = 0;
		for (i = 0; i < nthr_dist[0]; i++)
		{
			count_one = fp_nsites[0][i];
			//	printf("FPsites %d FPpeak %d TPpeak %d\n", fp_nsites[0][i], fp_one[i],tp_one[i]);
			if (count_one >= all_pos_thr[mot1])
			{
				index_thr[mot1] = i;
				fp_thr_rest[mot1] = 1 - (double)(count_one - all_pos_thr[mot1]) / fp_nsites[0][i];
				break;
			}
		}
		double auc_one = 0;
		int nthr_dist1 = nthr_dist[0] - 1;
		double prec_pred = 1;
		prec[0] = 1, recall[0] = 0;
		int nthr_here = 1;
		for (i = 0; i <= index_thr[mot1]; i++)
		{
			int dtp = tp_one[i];
			int i1 = i - 1;
			if (i > 0)dtp -= tp_one[i1];
			if (dtp > 0)
			{
				double dtpi = (double)tp_one[i];
				double dfpi = (double)fp_one[i];
				if (i > 0)
				{
					dfpi -= fp_one[i1];
					dtpi -= tp_one[i1];
				}
				if (i == index_thr[mot1])
				{
					dtpi *= fp_thr_rest[mot1];
					dfpi *= fp_thr_rest[mot1];
				}
				double fp_onec = dfpi;
				double tp_onec = dtpi;
				if (i > 0)
				{
					fp_onec += (double)fp_one[i1];
					tp_onec += (double)tp_one[i1];
				}
				double prec_cur = tp_onec / (tp_onec + nseq_fb * fp_onec);
				double prec_av = (prec_pred + prec_cur) / 2;
				double dauc = dtpi * (prec_av - prec_exp);
				recall[nthr_here] = tp_onec / nseq_real;
				prec[nthr_here] = prec_cur;
				auc_one += dauc;
				prec_pred = prec_cur;
				nthr_here++;
			}
		}
		auc_one *= 2;
		auc_one /= nseq_real;
		printf("%d\t%f\n", mot1 + 1, auc_one);
		fprintf(out_log1, "%d\t%f\n", mot1 + 1, auc_one);
		//fprintf(out_roc, "%s\t%s\t%d\t%s\t%g\n", file_for, file_back, mot1 + 1, motif_name[mot1], auc_one);
		motifp[mot1].auc = auc_one;
		//	motifp[mot1].rank = -1;		
		motifp[mot1].num = mot1;
		{
			strcpy(file_prc1, file_prc);
			strcat(file_prc1, "_");
			char buf[5];
			sprintf(buf, "%d", mot1 + 1);
			strcat(file_prc1, buf);
		}
		if ((out_roc = fopen(file_prc1, "wt")) == NULL)
		{
			fprintf(out_roc, "Input file %s can't be opened!\n", file_prc1);
			exit(1);
		}
		fprintf(out_roc, "%s\t%s\t%d\t%f\n", file_for, file_back, mot1 + 1, motifp[mot1].auc);
		for (i = 0; i < nthr_here; i++)
		{
			fprintf(out_roc, "%f\t%f\n", recall[i], prec[i]);
		}
		fclose(out_roc);
	}
	qsort(motifp, n_motifs, sizeof(motifp[0]), compare_auc);
	for (i = 0; i < n_motifs; i++)motifp[i].rank = i + 1;
	qsort(motifp, n_motifs, sizeof(motifp[0]), compare_num);
	int** tp_two;
	tp_two = new int* [2];
	if (tp_two == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		tp_two[i] = new int[nthr_dist_max];
		if (tp_two[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	int** fp_two;
	fp_two = new int* [2];
	if (fp_two == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (i = 0; i < 2; i++)
	{
		fp_two[i] = new int[nthr_dist_max];
		if (fp_two[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	qbs* tab;
	tab = new qbs[2 * nthr_dist_max];
	if (tab == NULL) { puts("Out of memory..."); exit(1); }
	fprintf(out_mat, "%s\t%s\npAUC\t", file_for, file_back);
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		fprintf(out_mat, "\t%f", motifp[mot1].auc);
	}
	fprintf(out_mat, "\n\tRank");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		fprintf(out_mat, "\t%d", motifp[mot1].rank);
	}
	fprintf(out_mat, "\n");
	fprintf(out_auc, "#Motif 1\t#Motif 2\tRank 1\tRank 2\tSimilarity\tpAUC 1\tpAUC 2\tpAUC 1&2\tRatio\n");
	int mot1c = 0;
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		fseek(in_sga, motifp[mot1].p, SEEK_SET);
		size_t mu0 = fread(&sta[0], sizeof(sta[0]), 1, in_sga);
		if (mu0 != 1)
		{
			printf("Error reading file %s!\n", partner_db);
			return -1;
		}
		int nthr_dist_all1 = 0;
		fread(&nthr_dist_all1, sizeof(int), 1, in_sga);
		fread(thr_all[0], sizeof(double), nthr_dist_all1, in_sga);
		fread(fpr_all[0], sizeof(double), nthr_dist_all1, in_sga);
		//if (strstr(motif_class[mot1], tfclass) == NULL)continue;	
		nthr_dist[0] = motifp[mot1].nth;
		{
			fprintf(out_mat, "%f\t%d\t", motifp[mot1].auc, motifp[mot1].rank);
		}
		for (i = 0; i < mot1c; i++)fprintf(out_mat, "\t");
		fprintf(out_auc, "\n");
		int mot2c = 0;
		for (mot2 = mot1 + 1; mot2 < n_motifs; mot2++)
		{
			fseek(in_sga, motifp[mot2].p, SEEK_SET);
			size_t mu1 = fread(&sta[1], sizeof(sta[1]), 1, in_sga);
			if (mu1 != 1)
			{
				printf("Error reading file %s!\n", partner_db);
				return -1;
			}
			int nthr_dist_all2;
			fread(&nthr_dist_all2, sizeof(int), 1, in_sga);
			fread(thr_all[1], sizeof(double), nthr_dist_all2, in_sga);
			fread(fpr_all[1], sizeof(double), nthr_dist_all2, in_sga);
			nthr_dist[1] = motifp[mot2].nth;
			for (n = 0; n < 2; n++)
			{
				for (j = 0; j <= nthr_dist[n]; j++)
				{
					fp_two[n][j] = 0;
					tp_two[n][j] = 0;
					fp_nsites[n][j] = 0;
				}
			}
			int all_pos_thr_two = (int)((all_pos[mot1] + all_pos[mot2]) * fp2 / 2);
			//printf("Real %d\n",mot + 1);
			int sga_check = SGA_rec_real(sta, nthr_dist, thr_all, fpr_all, tp_two, nseq_real, seq_real);
			if (sga_check == -1)
			{
				printf("Two motifs recognition error, real %d %d\n", mot1 + 1, mot2 + 1);
				exit(1);
			}
			//printf("Back %d\n", mot + 1);				
			sga_check = SGA_rec_back(sta, nthr_dist, thr_all, fpr_all, fp_two, fp_nsites, nseq_back, seq_back);
			if (sga_check == -1)
			{
				printf("Two motifs recognition error, back %d %d\n", mot1 + 1, mot2 + 1);
				exit(1);
			}
			int nthr_dist_two = nthr_dist[0] + nthr_dist[1];
			qbs* tab;
			tab = new qbs[nthr_dist_two];
			for (j = 0; j < nthr_dist_two; j++) { tab[j].tpr = tab[j].fpr = tab[j].fps = 0; }
			k = 0;
			for (j = 0; j < 2; j++)
			{
				for (i = 0; i < nthr_dist[j]; i++)
				{
					tab[k].tpr = tp_two[j][i];
					tab[k].fpr = fp_two[j][i];
					tab[k].fps = fp_nsites[j][i];
					tab[k].err = fpr_all[j][i];
					//	if (i < 20)printf("%d\t%d\t%d\t%d\t%f\t%d\n", j, i, tab[k].tpr, tab[k].fpr, tab[k].err, tab[k].fps);
					k++;
				}
			}
			qsort(tab, nthr_dist_two, sizeof(tab[0]), compare_qbs);
			int count_two = 0, count_pred = 0;
			int nthr_dist_two1 = nthr_dist_two - 1;
			int index_thr_two = 0;
			double fp_rest_two = 0;
			for (i = 0; i < nthr_dist_two; i++)
			{
				count_two = tab[i].fps;
				if (tab[i].fps > 0 && (i == nthr_dist_two1 || tab[i + 1].err != tab[i].err))
				{
					//printf("ERR %f Count %d FPsites %d FPpeak %d TPpeak %d\n", tab[i].err, count_two, tab[i].fps, tab[i].fpr, tab[i].tpr);
					if (count_two >= all_pos_thr_two || i == nthr_dist_two1)
					{
						fp_rest_two = (double)(all_pos_thr_two - count_pred) / (count_two - count_pred);
						index_thr_two = i;
						break;
					}
					count_pred = count_two;
				}
			}
			double prec_pred = 1;
			double auc_two = 0;
			int n_here = 0;
			double tpsum = 0, fpsum = 0;
			double dtp = 0, dfp = 0;
			for (i = 0; i <= index_thr_two; i++)
			{
				dtp += (double)tab[i].tpr;
				dfp += (double)tab[i].fpr;
				if (dtp > 0 && (i == index_thr_two || tab[i + 1].err != tab[i].err))
				{
					if (i == index_thr_two)
					{
						dtp *= fp_rest_two;
						dfp *= fp_rest_two;
					}
					tpsum += dtp;
					fpsum += dfp;
					double prec_cur = tpsum / (tpsum + nseq_fb * fpsum);
					double prec_av = (prec_pred + prec_cur) / 2;
					double dauc = dtp * (prec_av - prec_exp);
					recall[n_here] = tpsum / nseq_real;
					prec[n_here] = prec_cur;
					prec_pred = prec_cur;
					n_here++;
					auc_two += dauc;
					dtp = 0;
					dfp = 0;
				}
			}
			auc_two *= 2;
			auc_two /= nseq_real;
			{
				strcpy(file_prc1, file_prc);
				strcat(file_prc1, "_");
				char buf[5];
				sprintf(buf, "%d", mot1 + 1);
				strcat(file_prc1, buf);
				strcat(file_prc1, "&");
				sprintf(buf, "%d", mot2 + 1);
				strcat(file_prc1, buf);
			}
			if ((out_roc = fopen(file_prc1, "wt")) == NULL)
			{
				fprintf(out_roc, "Input file %s can't be opened!\n", file_prc1);
				exit(1);
			}
			fprintf(out_roc, "%s\t%s\t%d\t%d\t%f\t%f\t%f\n", file_for, file_back, mot1 + 1, mot2 + 1, motifp[mot1].auc, motifp[mot2].auc, auc_two);
			for (i = 0; i < n_here; i++)
			{
				fprintf(out_roc, "%f\t%f\n", recall[i], prec[i]);
			}
			fclose(out_roc);
			printf("%d\t%d\t%f\t%f\t%f", mot1 + 1, mot2 + 1, motifp[mot1].auc, motifp[mot2].auc, auc_two);
			//if (sims[mot1c][mot2c] != 1)printf("\t%f", sims[mot1c][mot2c]);
			printf("\n");
			fprintf(out_log2, "%d\t%d\t%d\t%d\t\t%f\t%f\t%f\n", mot1 + 1, mot2 + 1, motifp[mot1].rank, motifp[mot2].rank, motifp[mot1].auc, motifp[mot2].auc, auc_two);
			double auc_max = Max(motifp[mot1].auc, motifp[mot2].auc);
			auctwei[mot1c][mot2c] = auc_two / auc_max;
			//if (auc_two > motifp[mot1].auc && auc_two > motifp[mot2].auc)
			{
				fprintf(out_auc, "%d\t%d\t%d\t%d\t", mot1 + 1, mot2 + 1, motifp[mot1].rank, motifp[mot2].rank);
				fprintf(out_auc, "\t % f\t % f\t % f\t % f\n", motifp[mot1].auc, motifp[mot2].auc, auc_two, auctwei[mot1c][mot2c]);
			}
			mot2c++;
			//	else fprintf(out_auc, "\t");
			fprintf(out_mat, "\t%f", auc_two);
			/*if (auc_two > motifp[mot1].auc && auc_two > motifp[mot2].auc)
			{
				fprintf(out_roc, "\t%s\t%s\t%s\t%s\t%g\t%g\t%g\n", file_for, file_back, motif_name[motifp[mot1].num], motif_name[motifp[mot2].num], motifp[mot1].auc, motifp[mot2].auc, auc_two);
				for (i = 0; i < n_here; i++)
				{
					fprintf(out_roc, "%g\t%f\n", prec[i], recall[i]);
				}
			}*/
		}
		mot1c++;
		fprintf(out_mat, "\n");
	}
	fprintf(out_mat, "\n");
	fprintf(out_mat, "pAUC\t");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		fprintf(out_mat, "\t%f", motifp[mot1].auc);
	}
	fprintf(out_mat, "\n\tRank");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		fprintf(out_mat, "\t%d", motifp[mot1].rank);
	}
	fprintf(out_mat, "\n");
	k = 0;
	mot1c = 0;
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		fprintf(out_mat, "%f\t%d\t", motifp[mot1].auc, motifp[mot1].rank);
		for (i = 0; i < mot1c; i++)fprintf(out_mat, "\t");
		int mot2c = 0;
		for (mot2 = mot1 + 1; mot2 < n_motifs; mot2++)
		{
			fprintf(out_mat, "\t%f", auctwei[mot1c][mot2c++]);
		}
		mot1c++;
		fprintf(out_mat, "\n");
	}
	fprintf(out_mat, "pAUC\t");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		fprintf(out_mat, "\t%f", motifp[mot1].auc);
	}
	fprintf(out_mat, "\n\tRank");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		fprintf(out_mat, "\t%d", motifp[mot1].rank);
	}
	fprintf(out_mat, "\n");	
	fclose(out_mat);
	fclose(out_auc);
	fclose(out_log1);
	fclose(out_log2);
	//printf("All\t");
	for (k = 0; k < 2; k++)
	{
		delete[] thr_all[k];
	}
	delete[] thr_all;
	for (k = 0; k < 2; k++)
	{
		delete[] fpr_all[k];
	}
	delete[] fpr_all;
	delete[] len_real;
	delete[] len_back;
	delete[] tp_two;
	delete[] recall;
	delete[] prec;
	delete[] tp_one;
	delete[] index_thr;
	delete[] all_pos;
	delete[] all_pos_thr;
	delete[] fp_thr_rest;
	delete[] fp_one;
	for (k = 0; k < 2; k++)
	{
		delete[] fp_nsites[k];
	}
	delete[] fp_nsites;
	for (k = 0; k < n_motifs1; k++)
	{
		delete[] auctwei[k];
	}
	delete[] auctwei;
	delete[] motifp;
	delete[] tab;
	for (k = 0; k < 2; k++)
	{
		delete[] fp_two[k];
	}
	delete[] fp_two;
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < nseq_back; i++)
		{
			delete[] seq_back[k][i];
		}
		delete[] seq_back[k];
	}
	delete[] seq_back;
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < nseq_real; i++)
		{
			delete[] seq_real[k][i];
		}
		delete[] seq_real[k];
	}
	delete[] seq_real;
	return 0;
	delete[]sta;
}