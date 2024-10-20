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
#define NUM_LIBRARY 7// 4islo bibliotek
#define DIM 100

int StrNStr(char* str, char c, int n)
{
	int i, len = strlen(str);
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
double UnderStol(char* str, int nstol, char razd)
{
	if (nstol == 0)return atof(str);
	char ret[100];
	memset(ret, 0, sizeof(ret));
	int p1 = StrNStr(str, razd, nstol);
	int p2 = StrNStr(str, razd, nstol + 1);
	if (p2 == -1)
	{
		p2 = strlen(str);
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
void Mix(double* a, double* b)
{
	double buf = *a;
	*a = *b;
	*b = buf;
}
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
	double c;
	double std;
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
	DelChar(d, '\n');
	strcpy(site, d);
	fgets(d, sizeof(d), in);
	size = atoi(d);
	fgets(d, sizeof(d), in);
	len = atoi(d);
	fgets(d, sizeof(d), in);
	c = atof(d);
	std = 0.05;
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
	a->std = std;
	a->len = len;
	a->c = c;
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
	int nfo;
	double fpr;
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
void DelHole(char* str)
{
	char* hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
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
		case 'a': {d[i] = 't'; break; }
		case 't': {d[i] = 'a'; break; }
		case 'c': {d[i] = 'g'; break; }
		case 'g': {d[i] = 'c'; break; }
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
// ras4et 4asot oligonukleotidov po stroke (zdes' - nukleotidov)
void GetSostPro(char* d, int word, int* sost)
{
	int i, j, k, i_sost, let;
	char letter[] = "acgt";
	int ten[6] = { 1, 4, 16, 64, 256, 1024 };
	int lens = strlen(d);
	int size = 1;
	for (k = 0; k < word; k++)size *= 4;
	for (i = 0; i < size; i++)sost[i] = 0;
	for (i = 0; i < lens - word + 1; i++)
	{
		i_sost = 0;
		let = -1;
		for (j = word - 1; j >= 0; j--)
		{
			for (k = 0; k < 4; k++)
			{
				if (d[i + j] == letter[k]) { let = k; break; }
			}
			i_sost += ten[word - 1 - j] * let;
		}
		sost[i] = i_sost;
	}
}
void PWMScore(double(&pwm)[MATLEN][OLIGNUM], double& min, double& raz, int len1)
{
	int i, j;
	for (i = 0; i < len1; i++)
	{
		double pwmmin = 100;
		double pwmmax = -100;
		for (j = 0; j < OLIGNUM; j++)
		{
			if (pwm[i][j] < pwmmin)pwmmin = pwm[i][j];
			if (pwm[i][j] > pwmmax)pwmmax = pwm[i][j];
		}
		raz += pwmmax;
		min += pwmmin;
	}
	raz -= min;
}
int PWM_rec_real_one(double(&pwm)[MATLEN][OLIGNUM], double min, double raz, int nthr_dist, double* thr_all, double* fpr_all, int* tpr, int olen, int nseq, char*** seq)
{
	int i, j, n;
	int compl1;
	int cod[MATLEN];
	char d[MATLEN];
	int word = 1;
	int nthr_dist1 = nthr_dist - 1;
	double thr_cr = thr_all[nthr_dist1];
	for (n = 0; n < nseq; n++)
	{
		//if ((n + 1) % 500 == 0)printf("\b\b\b\b\b\b\b%7d", n + 1);
		int len_pro1 = strlen(seq[0][n]);
		int len21 = len_pro1 - olen;
		int index = nthr_dist;
		for (i = 0; i <= len21; i++)
		{
			for (compl1 = 0; compl1 < 2; compl1++)
			{
				int ista;
				if (compl1 == 0)ista = i;
				else ista = len21 - i;
				strncpy(d, &seq[compl1][n][ista], olen);
				d[olen] = '\0';
				if (strstr(d, "n") != NULL) { continue; }
				GetSostPro(d, word, cod);
				double score = 0;
				for (j = 0; j < olen; j++)
				{
					score += pwm[j][cod[j]];
				}
				score -= min;
				score /= raz;
				if (score >= thr_cr)
				{
					if (score >= thr_all[0])
					{
						index = 0;
						break;
					}
					else
					{
						int index_here = nthr_dist;
						for (j = 1; j < nthr_dist; j++)
						{
							if (score >= thr_all[j] && score < thr_all[j - 1])
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
int PWM_rec_back_one(double(&pwm)[MATLEN][OLIGNUM], double min, double raz, int nthr_dist, double* thr_all, double* fpr_all, int* fpr, int* fp_nsites, int olen, int nseq, char*** seq, int& all_pos)
{
	int i, j, n;
	int compl1;
	int cod[MATLEN];
	char d[MATLEN];
	int word = 1;
	int nthr_dist1 = nthr_dist - 1;
	double thr_cr = thr_all[nthr_dist1];

	for (n = 0; n < nseq; n++)
	{
		//if ((n+1) % 500 == 0)printf("\b\b\b\b\b\b\b%7d", n+1);
		int len_pro1 = strlen(seq[0][n]);
		int len21 = len_pro1 - olen;
		int index_best = nthr_dist;
		for (i = 0; i <= len21; i++)
		{
			for (compl1 = 0; compl1 < 2; compl1++)
			{
				int ista;
				if (compl1 == 0)ista = i;
				else ista = len21 - i;
				strncpy(d, &seq[compl1][n][ista], olen);
				d[olen] = '\0';
				if (strstr(d, "n") != NULL) { continue; }
				all_pos++;
				GetSostPro(d, word, cod);
				double score = 0;
				for (j = 0; j < olen; j++)
				{
					score += pwm[j][cod[j]];
				}
				score -= min;
				score /= raz;
				int index = nthr_dist;
				if (score >= thr_cr)
				{
					if (score >= thr_all[0])
					{
						index = 0;
					}
					else
					{
						for (j = 1; j < nthr_dist; j++)
						{
							if (score >= thr_all[j] && score < thr_all[j - 1])
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
int SGA_rec_real_one(city *sta, int nthr_dist, double* thr_all, double* fpr_all, int* tpr, int nseq, char*** seq)
{
	int i, j, n, m;
	int compl1;	
	char d[MATLEN];	
	int nthr_dist1 = nthr_dist - 1;
	double thr_cr = thr_all[nthr_dist1];
	for (n = 0; n < nseq; n++)
	{
		//if ((n + 1) % 500 == 0)printf("\b\b\b\b\b\b\b%7d", n + 1);
		int len_pro1 = strlen(seq[0][n]);
		int len21 = len_pro1 - sta->len;
		int index = nthr_dist;
		for (i = 0; i <= len21; i++)
		{
			for (compl1 = 0; compl1 < 2; compl1++)
			{
				int ista;
				if (compl1 == 0)ista = i;
				else ista = len21 - i;
				strncpy(d, &seq[compl1][n][ista], sta->len);
				d[sta->len] = '\0';
				if (strstr(d, "n") != NULL) { continue; }				
				double score = sta->c;
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
				score = 1 - fabs(score - 1);
				if (score >= thr_cr)
				{
					if (score >= thr_all[0])
					{
						index = 0;
						break;
					}
					else
					{
						int index_here = nthr_dist;
						for (j = 1; j < nthr_dist; j++)
						{
							if (score >= thr_all[j] && score < thr_all[j - 1])
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
int SGA_rec_back_one(city *sta, int nthr_dist, double* thr_all, double* fpr_all, int* fpr, int* fp_nsites, int nseq, char*** seq, int& all_pos)
{
	int i, j, n, m;
	int compl1;
	char d[MATLEN];
	int nthr_dist1 = nthr_dist - 1;
	double thr_cr = thr_all[nthr_dist1];

	for (n = 0; n < nseq; n++)
	{
		//if ((n+1) % 500 == 0)printf("\b\b\b\b\b\b\b%7d", n+1);
		int len_pro1 = strlen(seq[0][n]);
		int len21 = len_pro1 - sta->len;
		int index_best = nthr_dist;
		for (i = 0; i <= len21; i++)
		{
			for (compl1 = 0; compl1 < 2; compl1++)
			{
				int ista;
				if (compl1 == 0)ista = i;
				else ista = len21 - i;
				strncpy(d, &seq[compl1][n][ista], sta->len);
				d[sta->len] = '\0';
				if (strstr(d, "n") != NULL) { continue; }
				all_pos++;
				double score = sta->c;
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
				score = 1 - fabs(score - 1);
				int index = nthr_dist;
				if (score >= thr_cr)
				{
					if (score >= thr_all[0])
					{
						index = 0;
					}
					else
					{
						for (j = 1; j < nthr_dist; j++)
						{
							if (score >= thr_all[j] && score < thr_all[j - 1])
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
int PWM_SGA_rec_real(double(&pwm)[MATLEN][OLIGNUM], double min, double raz, city *sta, int nthr_dist[2], double** thr_all, double** fpr_all, int** tp_two, int olen[2], int nseq, char*** seq)
{
	int i, j, k, n, m;
	int compl1;
	int cod[MATLEN];
	char d[MATLEN];
	int word = 1;
	int nthr_dist1[2], olen1[2];
	double thr_cr[2];
	for (n = 0; n < 2; n++)
	{
		nthr_dist1[n] = nthr_dist[n] - 1;
		olen1[n] = olen[n] - 1;
		thr_cr[n] = thr_all[n][nthr_dist1[n]];
	}
	for (n = 0; n < nseq; n++)
	{
		//if ((n + 1) % 500 == 0)printf("\b\b\b\b\b\b\b%7d", n + 1);
		int best_inx[2];
		for (k = 0; k < 2; k++)
		{
			int index;
			int len_pro1 = strlen(seq[0][n]);
			int len21 = len_pro1 - olen[k];
			index = best_inx[k] = nthr_dist[k];
			for (i = 0; i <= len21; i++)
			{
				for (compl1 = 0; compl1 < 2; compl1++)
				{
					int ista;
					if (compl1 == 0)ista = i;
					else ista = len21 - i;
					strncpy(d, &seq[compl1][n][ista], olen[k]);
					d[olen[k]] = '\0';
					if (strstr(d, "n") != NULL) { continue; }
					double score = 0;
					if (k == 0)
					{
						GetSostPro(d, word, cod);						
						for (j = 0; j < olen[k]; j++)
						{
							score += pwm[j][cod[j]];
						}
						score -= min;
						score /= raz;
					}
					else
					{
						score += sta->c;
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
						score = 1 - fabs(score - 1);
					}
					if (score >= thr_cr[k])
					{
						if (score >= thr_all[k][0])
						{
							index = 0;
							break;
						}
						else
						{
							for (j = 1; j < nthr_dist[k]; j++)
							{
								if (score >= thr_all[k][j] && score < thr_all[k][j - 1])
								{
									index = j;
									break;
								}
							}
						}
					}
				}
				if (index < best_inx[k])best_inx[k] = index;
				if (index == 0)break;
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
int main(int argc, char* argv[])
{
	int i, j, k, n;
	char*** seq_real, *** seq_back;
	char file_for[ARGLEN], file_back[ARGLEN], partner_pwm[ARGLEN], partner_sga[ARGLEN];// path_fasta[ARGLEN], pfile_for[ARGLEN], pfile_back[ARGLEN]
	char file_prc[3][ARGLEN], file_auc[ARGLEN];
	FILE * out_roc[3], * out_auc, * in_pwm, * in_sga;

	if (argc != 10)
	{
		printf("%s 1,2input fasta foreground,background  3,4input pwm,sga binary files 5double ERRthresh 6output file pAUC 7,8,9output files PR curves", argv[0]);
		return -1;
	}
	//strcpy(path_fasta, argv[1]);
	strcpy(file_for, argv[1]);
	strcpy(file_back, argv[2]);
//	strcpy(pfile_for, path_fasta);
//	strcpy(pfile_back, path_fasta);
	//strcat(pfile_for, file_for);
//	strcat(pfile_back, file_back);
	strcpy(partner_pwm, argv[3]); //h12hs, h12mm
	strcpy(partner_sga, argv[4]); //h12hs, h12mm
	double fp2 = atof(argv[5]); //ERR threshold for pAUC-PR		
	strcpy(file_auc, argv[6]);
	strcpy(file_prc[0], argv[7]);
	strcpy(file_prc[1], argv[8]);
	strcpy(file_prc[2], argv[9]);
	int* len_real, * len_back, nseq_real = 0, nseq_back = 0;
	int olen_min = 8;
	int len_peak_max = 3000;
//	double fp2 = 0.005; //FPR threshold for pAUC	
	double fp2_lg = -log10(fp2);
	/*if ((outlog = fopen(file_log, "wt")) == NULL)
	{
		fprintf(outlog, "Input file %s can't be opened!\n", file_log);
		exit(1);
	}
	*/
	//	printf("EvalSeq\n");
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
	
	int index_thr[2], fp_thr[2];
	double fp_thr_rest[2];
	double nseq_fb = (double)nseq_real / nseq_back;
	double prec_exp = 0.5;
	int len_partner[2], nthr_dist[2];
	double pwm[MATLEN][OLIGNUM];
	double pfm[MATLEN][OLIGNUM];
	city sta;
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
	int* tp_one, * fp_one;
	tp_one = new int[nthr_dist_max];
	if (tp_one == NULL) { puts("Out of memory..."); exit(1); }
	fp_one = new int[nthr_dist_max];
	if (fp_one == NULL) { puts("Out of memory..."); exit(1); }
	double** recall1;
	recall1 = new double* [2];
	if (recall1 == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		recall1[i] = new double[nthr_dist_max];
		if (recall1[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	double** prec1;
	prec1 = new double* [2];
	if (prec1 == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		prec1[i] = new double[nthr_dist_max];
		if (prec1[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
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
	if (fp_two == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		fp_two[i] = new int[nthr_dist_max];
		if (fp_two[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	int** fp_nsites;// frequencies
	fp_nsites = new int* [2];
	if (fp_nsites == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		fp_nsites[i] = new int[nthr_dist_max];
		if (fp_nsites[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	if ((in_pwm = fopen(partner_pwm, "rb")) == NULL)
	{
		printf("Input file %s can't be opened!\n", partner_pwm);
		return -1;
	}
	if ((in_sga = fopen(partner_sga, "rb")) == NULL)
	{
		printf("Input file %s can't be opened!\n", partner_sga);
		return -1;
	}
	{
		int mu = 0;
		mu = fread(&len_partner[0], sizeof(int), 1, in_pwm);
		mu = fread(pfm, sizeof(double), 4 * len_partner[0], in_pwm);
		mu = fread(pwm, sizeof(double), 4 * len_partner[0], in_pwm);
		int nthr_dist_all[2];
		mu = fread(&nthr_dist_all[0], sizeof(int), 1, in_pwm);
		mu = fread(thr_all[0], sizeof(double), nthr_dist_all[0], in_pwm);
		mu = fread(fpr_all[0], sizeof(double), nthr_dist_all[0], in_pwm);
		fclose(in_pwm);
		mu = fread(&sta, sizeof(sta), 1, in_sga);
		mu = fread(&nthr_dist_all[1], sizeof(int), 1, in_sga);
		mu = fread(thr_all[1], sizeof(double), nthr_dist_all[1], in_sga);
		mu = fread(fpr_all[1], sizeof(double), nthr_dist_all[1], in_sga);
		fclose(in_sga);
		len_partner[1] = sta.len;
		for (n = 0; n < 2; n++)
		{
			nthr_dist[n] = 0;
			for (i = 0; i < nthr_dist_all[n]; i++)
			{
				nthr_dist[n]++;
			//	fpr_all[n][i] = -log10(fpr_all[n][i]);
				if (fpr_all[n][i] < fp2_lg)break;
			}
		}
	}
	double min = 0, raz = 0;
	PWMScore(pwm, min, raz, len_partner[0]);		
	double auc_one[2] = { 0,0 };
	int n_here1[2] = { 0,0 };	
	for(n=0;n<2;n++)
	{
		for (j = 0; j <= nthr_dist[n]; j++)tp_one[j] = fp_one[j] = fp_nsites[n][j] = 0;
		int all_pos = 0;
		if (n == 0)
		{
			int check = PWM_rec_real_one(pwm, min, raz, nthr_dist[0], thr_all[0], fpr_all[0], tp_one, len_partner[0], nseq_real, seq_real);
			if (check == -1)
			{
				printf("Motif pwm recognition error, real\n");
				exit(1);
			}
			check = PWM_rec_back_one(pwm, min, raz, nthr_dist[0], thr_all[0], fpr_all[0], fp_one, fp_nsites[0], len_partner[0], nseq_back, seq_back, all_pos);
			if (check == -1)
			{
				printf("Motif pwm recognition error, real\n");
				exit(1);
			}
		}
		else
		{
			int check = SGA_rec_real_one(&sta, nthr_dist[1], thr_all[1], fpr_all[1], tp_one, nseq_real, seq_real);
			if (check == -1)
			{
				printf("Motif sga recognition error, real\n");
				exit(1);
			}
			check = SGA_rec_back_one(&sta, nthr_dist[1], thr_all[1], fpr_all[1], fp_one, fp_nsites[1], nseq_back, seq_back, all_pos);
			if (check == -1)
			{
				printf("Motif sga recognition error, real\n");
				exit(1);
			}
		}
		int all_pos_thr = (int)(all_pos * fp2);
		index_thr[n] = nthr_dist[n] - 1;
		int count_one = 0;
		for (i = 0; i < nthr_dist[n]; i++)
		{
			count_one += fp_nsites[n][i];
			printf("FPsites %d FPpeak %d TPpeak %d\n", fp_nsites[n][i], fp_one[i],tp_one[i]);
			if (count_one >= all_pos_thr)
			{
				index_thr[n] = i;
				fp_thr[n] = fp_one[i];
				fp_thr_rest[n] = 1 - (double)(count_one - all_pos_thr) / fp_nsites[n][i];
				break;
			}
		}
		/*for (i = 0; i <= nthr_dist[n]; i++)
		{
			if ((i + 1) % 10 == 0)printf("%d\n", i + 1);
			printf("%5d %d\t", tp_one[i], fp_one[i]);
		}*/
		//printf("ROC %d\n", mot + 1);
		double prec_pred = 1;// (double)tp_one[0] / ((double)tp_one[0] + nseq_fb * fp_one[0]);
		recall1[n][0] = 0;//(double)tp_one[0]/nseq_real, 
		prec1[n][0] = prec_pred;
		n_here1[n] = 1;
		int index1 = index_thr[n] - 1;
		for (i = 0; i <= index_thr[n]; i++)
		{
			int dtp;
			if (i > 0)dtp = tp_one[i] - tp_one[i - 1];
			else dtp = tp_one[i];
			if (dtp > 0 || (i == index1 || i == index_thr[n]))
			{
				double dtpi = (double)tp_one[i];
				double prec_cur = dtpi / (dtpi + nseq_fb * fp_one[i]);
				double prec_av = (prec_pred + prec_cur) / 2;
				double dauc = dtp * (prec_av - prec_exp);
				recall1[n][n_here1[n]] = dtpi / nseq_real;
				prec1[n][n_here1[n]] = prec_cur;
				if (i == index_thr[n])dauc *= fp_thr_rest[n];
				auc_one[n] += dauc;
				n_here1[n]++;
			}
		}
		auc_one[n] *= 2;
		auc_one[n] /= nseq_real;
		if(n==0)printf("%s\t%s\t%d\t%s\t%g\n", file_for, file_back, n + 1, partner_pwm, auc_one[n]);
		else printf("%s\t%s\t%d\t%s\t%g\n", file_for, file_back, n + 1, partner_sga, auc_one[n]);
	}
	for (n = 0; n < 2; n++)for (j = 0; j <= nthr_dist[n]; j++)tp_two[n][j] = 0;
	for (n = 0; n < 2; n++)for (j = 0; j <= nthr_dist[n]; j++)fp_two[n][j] = 0;
	//printf("Real %d\n",mot + 1);
	int check = PWM_SGA_rec_real(pwm, min, raz, &sta, nthr_dist, thr_all, fpr_all, tp_two, len_partner, nseq_real, seq_real);
	if (check == -1)
	{
		printf("Motifs pwm&sga recognition error, real\n");
		exit(1);
	}
	//printf("Back %d\n", mot + 1);	
	int all_pos_two[2] = { 0,0 };
	check = PWM_SGA_rec_real(pwm, min, raz, &sta, nthr_dist, thr_all, fpr_all, fp_two, len_partner, nseq_back, seq_back);
	if (check == -1)
	{
		printf("Motifs pwm&sga recognition error, back\n");
		exit(1);
	}
	double fp_thr_two = (double)(fp_thr[0] + fp_thr[1]) / nseq_back / 2;
	double fp_thr_rest2 = (fp_thr_rest[0] + fp_thr_rest[1]) / 2;
	int all_pos;
	if (all_pos_two[0] > all_pos_two[1])all_pos = all_pos_two[0];
	else all_pos = all_pos_two[1];
	int nthr_dist_two = nthr_dist[0] + nthr_dist[1];
	qbs* tab;
	tab = new qbs[nthr_dist_two];
	for (j = 0; j < nthr_dist_two; j++) { tab[j].nfo = 0; tab[j].fpr = 0; }
	k = 0;
	for (j = 0; j < 2; j++)
	{
		for (i = 0; i < nthr_dist[j]; i++)
		{
			tab[k].nfo = tp_two[j][i];
			tab[k].fpr = fp_two[j][i];
			tab[k].err = fpr_all[j][i];
			k++;
		}
	}
	qsort(tab, nthr_dist_two, sizeof(tab[0]), compare_qbs);
	for (i = 1; i < nthr_dist_two; i++)
	{
		int i1 = i - 1;
		tab[i].nfo += tab[i1].nfo;
		tab[i].fpr += tab[i1].fpr;
	}
	double* recall;
	recall = new double[nthr_dist_two];
	if (recall == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	double* prec;
	prec = new double[nthr_dist_two];
	if (prec == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	double auc_two = 0;
	int n_here = 1;
	recall[0] = 0;
	prec[0] = 1;
	double prec_pred = 1;
	for (i = 0; i < nthr_dist_two; i++)
	{
		int dtp;
		if (i > 0)dtp = tab[i].nfo - tab[i - 1].nfo;
		else dtp = tab[i].nfo;
		if (dtp > 0)
		{
			double dtpi = (double)tab[i].nfo;
			double prec_cur = dtpi / (dtpi + nseq_fb * tab[i].fpr);
			double prec_av = (prec_pred + prec_cur) / 2;
			double dauc = dtp * (prec_av - prec_exp);
			recall[n_here] = dtpi / nseq_real;
			prec[n_here] = prec_cur;
			n_here++;
			double fproc_cur = (double)tab[i].fpr / nseq_back;
			if (fproc_cur >= fp_thr_two)
			{
				dauc *= fp_thr_rest2;
				auc_two += dauc;
				break;
			}
			else auc_two += dauc;
		}
	}
	auc_two *= 2;
	auc_two /= nseq_real;
	double auc_max = Max(auc_one[0], auc_one[1])
	double rat = auc_two / auc_max;
	printf("%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\n", file_for, file_back, partner_pwm, partner_sga, auc_one[0], auc_one[1], auc_two,rat);
	//fprintf(outlog, "%s\t%s\t%s\t%s\t%g\t%g\t%g\n", file_for, file_back, partner_pwm, partner_sga, auc_one[0], auc_one[1], auc_two);
	if ((out_auc = fopen(file_auc, "wt")) == NULL)
	{
		fprintf(out_auc, "Input file %s can't be opened!\n", file_auc);
		exit(1);
	}
	fprintf(out_auc, "%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\n", file_for, file_back, partner_pwm, partner_sga, auc_one[0], auc_one[1], auc_two,rat);
	fclose(out_auc);
	for (n = 0; n < 3; n++)
	{
		if ((out_roc[n] = fopen(file_prc[n], "wt")) == NULL)
		{
			fprintf(out_roc[n], "Input file %s can't be opened!\n", file_prc[n]);
			exit(1);
		}
	}
	for (n = 0; n < 2; n++)
	{
		if(n==0)fprintf(out_roc[n], "%s\t%s\t%s\t%g\n", file_for, file_back, partner_pwm, auc_one[n]);
		else fprintf(out_roc[n], "%s\t%s\t%s\t%g\n", file_for, file_back, partner_sga, auc_one[n]);
		for (i = 0; i < n_here1[n]; i++)
		{
			fprintf(out_roc[n], "%f\t%f\n", recall1[n][i], prec1[n][i]);
			if (prec1[n][i] == fp2)break;
		}
	}
	fprintf(out_roc[2], "%s\t%s\t%s\t%s\t%g\t%g\t%g\n", file_for, file_back, partner_pwm, partner_sga, auc_one[0], auc_one[1], auc_two);
	for (i = 0; i < n_here; i++)
	{
		fprintf(out_roc[2], "%f\t%f\n", recall[i], prec[i]);
		if (prec[i] == fp2)break;
	}
	for (n = 0; n < 3; n++)fclose(out_roc[n]);
	//fclose(outlog);
	fclose(in_pwm);
	fclose(in_sga);
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
	delete[] tab;
	delete[] tp_one;
	delete[] fp_one;
	delete[] recall;
	delete[] prec;
	for (k = 0; k < 2; k++)
	{
		delete[] fp_two[k];
	}
	delete[] fp_two;
	for (k = 0; k < 2; k++)
	{
		delete[] tp_two[k];
	}
	delete[] tp_two;
	for (k = 0; k < 2; k++)
	{
		delete[] recall1[k];
	}
	delete[] recall1;
	for (k = 0; k < 2; k++)
	{
		delete[] prec1[k];
	}
	delete[] prec1;
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
}