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

struct acc {
	double auc;
	int rank;
	int nth;
	int num;
	long p;
};

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
int PWM_rec_back_one(double(&pwm)[MATLEN][OLIGNUM], double min, double raz, int nthr_dist, double* thr_all, double* fpr_all, int* fpr, int olen, int nseq, char*** seq, int& all_pos)
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
				for (j = nthr_dist; j >= index; j--)fpr[j]++;
			}
		}
	}
	return 1;
}
int PWM_rec_real(double(&pwm)[2][MATLEN][OLIGNUM], double min[2], double raz[2], int nthr_dist[2], double** thr_all, double** fpr_all, int** tp_two, int olen[2], int nseq, char*** seq)
{
	int i, j, k, n;
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
					GetSostPro(d, word, cod);
					double score = 0;
					for (j = 0; j < olen[k]; j++)
					{
						score += pwm[k][j][cod[j]];
					}
					score -= min[k];
					score /= raz[k];
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
	int sum = 0;
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < nthr_dist[i]; j++)
		{
			sum += tp_two[i][j];
		}
	}
	/*for (i = 0; i < 2; i++)
	{
		printf("TP%d\t",i+1);
		for (j = 0; j < 10; j++)
		{
			printf("%d ", tp_two[i][j]);
		}
		printf("\n");
	}*/
	return 1;
}

int PWM_rec_back(double(&pwm)[2][MATLEN][OLIGNUM], double min[2], double raz[2], int nthr_dist[2], double** thr_all, double** fpr_all, double** fpr, int* olen, int nseq, char*** seq, int* all_pos_two)
{
	int i, j, k, m, n;
	int compl1, compl2;
	int cod[MATLEN];
	char d[MATLEN];
	int word = 1;
	double koef = (double)(olen[0] + olen[1]) / olen[0] / olen[1] / 2;
	int nthr_dist1[2], olen1[2];
	double thr_cr[2];
	for (n = 0; n < 2; n++)
	{
		nthr_dist1[n] = nthr_dist[n] - 1;
		olen1[n] = olen[n] - 1;
		thr_cr[n] = thr_all[n][nthr_dist1[n]];
	}
	int olenmin = Min(olen[0], olen[1]);
	int psco[2][2][SEQLEN];
	for (n = 0; n < nseq; n++)
	{
		//	int rec_pos[2] = { 0,0 };
		double over = 0;
		int len1 = strlen(seq[0][n]);
		int len21[2];
		for (i = 0; i < 2; i++)len21[i] = len1 - olen[i];
		for (i = 0; i < 2; i++)for (j = 0; j < 2; j++)for (k = 0; k <= len21[i]; k++)psco[i][j][k] = nthr_dist[i];
		for (i = 0; i < 2; i++)
		{
			//if ((n+1) % 500 == 0)printf("\b\b\b\b\b\b\b%7d", n+1);					
			for (k = 0; k <= len21[i]; k++)
			{
				for (compl1 = 0; compl1 < 2; compl1++)
				{
					int ista;
					if (compl1 == 0)ista = k;
					else ista = len21[i] - k;
					strncpy(d, &seq[compl1][n][ista], olen[i]);
					d[olen[i]] = '\0';
					if (strstr(d, "n") != NULL) { continue; }
					all_pos_two[i]++;
					GetSostPro(d, word, cod);
					double score = 0;
					for (j = 0; j < olen[i]; j++)
					{
						score += pwm[i][j][cod[j]];
					}
					score -= min[i];
					score /= raz[i];
					int index = nthr_dist[i];
					if (score >= thr_cr[i])
					{
						if (score >= thr_all[i][0])
						{
							index = 0;
						}
						else
						{
							for (j = 1; j < nthr_dist[i]; j++)
							{
								if (score >= thr_all[i][j] && score < thr_all[i][j - 1])
								{
									index = j;
									break;
								}
							}
						}
						//	rec_pos[i]++;
					}
					psco[i][compl1][k] = index;
					fpr[i][index]++;
				}
			}
		}
		//printf("%d\t%d\t%d\t\t", n + 1, rec_pos[0], rec_pos[1]);
		for (i = 0; i < 1; i++)
		{
			j = 1 - i;
			for (k = 0; k <= len21[i]; k++)
			{
				for (compl1 = 0; compl1 < 2; compl1++)
				{
					int psco1 = psco[i][compl1][k];
					if (psco1 != nthr_dist[i])
					{
						int pmin = Max(0, k - olen1[j]);
						int pmax = Min(k + olen1[i], len21[j]);
						for (m = pmin; m < pmax; m++)
						{
							for (compl2 = 0; compl2 < 2; compl2++)
							{
								int psco2 = psco[j][compl2][m];
								if (psco2 != nthr_dist[j])
								{
									int right = Min(k + olen1[i], m + olen1[j]);
									int left = Max(k, m);
									double dfpr = 0;
									int overlap = right - left + 1;
									if (overlap > 0)
									{
										if (overlap == olenmin)dfpr = 1;
										else
										{
											if (overlap > 0 && overlap < olenmin)dfpr = koef * overlap;
										}
										if (fpr_all[i][psco1] > fpr_all[j][psco2])fpr[j][psco2] -= dfpr;
										else fpr[i][psco1] -= dfpr;
										over += dfpr;
									}
								}
							}
						}
					}
				}
			}
		}
		//		printf("%f\n", over);
			/*	if ((n + 1) % 50 == 0)
				{
					int yy = 0;
				}*/
	}
	/*for (i = 0; i < 2; i++)
	{
		printf("FP%d\t", i + 1);
		for (j = 0; j < 10; j++)
		{
			printf("%.2f ", fpr[i][j]);
		}
		printf("\n");
	}*/
	return 1;
}
// position frequency mattrix (PFM), position weight matrix (PWM)
struct matrices {
	int len;
	double** fre;
	void get_copy(matrices* a);
	int mem_in(int len);
	void mem_out(int len);
	void norm(void);
};

int matrices::mem_in(int length)
{
	int i;
	len = length;
	fre = new double* [len];
	if (fre == NULL) return -1;
	for (i = 0; i < len; i++)
	{
		fre[i] = new double[OLIGNUM];
		if (fre[i] == NULL) return -1;
	}
	return 1;
}
void matrices::mem_out(int len)
{
	int i;
	for (i = 0; i < len; i++) delete[] fre[i];
	delete[] fre;
}
void matrices::get_copy(matrices* a)
{
	a->mem_in(len);
	a->len = len;
	int i, j;
	for (i = 0; i < len; i++)
	{
		for (j = 0; j < OLIGNUM; j++)
		{
			a->fre[i][j] = fre[i][j];
		}
	}
}
#include "pfm_similarity.h" //permutation test for anchor/partner motif comparison (separate task for all algorithm)

int main(int argc, char* argv[])
{
	int i, j, k, n, mot1, mot2;
	char*** seq_real, *** seq_back;
	char file_for[ARGLEN], file_back[ARGLEN], partner_db[ARGLEN];// , tfclass[80];path_fasta[ARGLEN], pfile_for[ARGLEN], pfile_back[ARGLEN],
	char file_mat[ARGLEN], file_auc[ARGLEN], file_roc[ARGLEN], file_log1[ARGLEN], file_log2[ARGLEN], file_roc1[ARGLEN];
	FILE* out_log1, * out_auc, * out_roc, * out_log2, * out_mat, * in_pwm;

	if (argc != 10)
	{
		printf("%s 1,2input fasta foreground,background 3input binary files pwm_thresholds 4no. of motifs 5,6,7,8,9output files pAUC_matrix, pAUC_list, log1, log2, ROC", argv[0]);
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
	strcpy(file_mat, argv[5]);
	strcpy(file_auc, argv[6]);
	strcpy(file_log1, argv[7]);
	strcpy(file_log2, argv[8]);	
	strcpy(file_roc, argv[9]);
	int s_overlap_min = 6, s_ncycle_small = 1000, s_ncycle_large = 10000;//for permutation(motif_comparison) min_length_of_alignment, no. of permutation (test & detailed)
	double s_granul = 0.001;//for permutation(motif_comparison) okruglenie 4astotnyh matric	
	double p_crit = 0.05;// threshold for similarity of matrices
	double pval_sim[4];
	matrices matrix[2];

	int* len_real, * len_back, nseq_real = 0, nseq_back = 0;
	int olen_min = 8;
	int len_peak_max = 3000;
	double fp2 = 0.001; //FPR threshold for pAUC		
	double fp2_lg = -log10(fp2);
	int kauc = 10000;
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
	if ((in_pwm = fopen(partner_db, "rb")) == NULL)
	{
		printf("Input file %s can't be opened!\n", partner_db);
		return -1;
	}
	int len_partner[2], nthr_dist[2];
	double pwm[2][MATLEN][OLIGNUM];
	double pfm[2][MATLEN][OLIGNUM];
	//double pwm1[MATLEN][OLIGNUM];
	//double pfm1[MATLEN][OLIGNUM];

	int nthr_dist_max = 5000;
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
	double** sims;
	sims = new double* [n_motifs1];
	if (sims == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (i = 0; i < n_motifs1; i++)
	{
		sims[i] = new double[n_motifs1 - i];
		if (sims[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	for (mot1 = 0; mot1 < n_motifs1; mot1++)for (mot2 = 0; mot2 < n_motifs1 - mot1; mot2++)sims[mot1][mot2] = 1;
	//n_motifs = 50;
	int* tp_one, * fp_one;
	tp_one = new int[nthr_dist_max];
	if (tp_one == NULL) { puts("Out of memory..."); exit(1); }
	fp_one = new int[nthr_dist_max];
	if (fp_one == NULL) { puts("Out of memory..."); exit(1); }
	double* tp_here;
	tp_here = new double[nthr_dist_max];
	if (tp_here == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	double* fp_here;
	fp_here = new double[nthr_dist_max];
	if (fp_here == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	int nmot_max = 1500;
	acc* motifp;
	motifp = new acc[n_motifs];
	if (motifp == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	double min[2], raz[2];
	motifp[0].p = ftell(in_pwm);
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		//printf("%d\t%s\n", mot1 + 1, motif_class[mot1]);
		fread(&len_partner[0], sizeof(int), 1, in_pwm);
		int len_partner4 = len_partner[0] * 4;
		fread(pfm[0], sizeof(double), len_partner4, in_pwm);
		fread(pwm[0], sizeof(double), len_partner4, in_pwm);
		int nthr_dist_all;
		fread(&nthr_dist_all, sizeof(int), 1, in_pwm);
		fread(thr_all[0], sizeof(double), nthr_dist_all, in_pwm);
		fread(fpr_all[0], sizeof(double), nthr_dist_all, in_pwm);
		if (mot1 != n_motifs1)motifp[mot1 + 1].p = ftell(in_pwm);
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
		min[0] = raz[0] = 0;
		PWMScore(pwm[0], min[0], raz[0], len_partner[0]);
		for (j = 0; j <= nthr_dist[0]; j++)
		{
			tp_one[j] = fp_one[j] = 0;
		}
		int all_pos = 0;
		//printf("Real %d\n",mot + 1);
		int pwm_check = PWM_rec_real_one(pwm[0], min[0], raz[0], nthr_dist[0], thr_all[0], fpr_all[0], tp_one, len_partner[0], nseq_real, seq_real);
		if (pwm_check == -1)
		{
			printf("One motif recognition error, real motif %d\n", mot1 +1);
			exit(1);
		}
		//printf("Back %d\n", mot + 1);
		pwm_check = PWM_rec_back_one(pwm[0], min[0], raz[0], nthr_dist[0], thr_all[0], fpr_all[0], fp_one, len_partner[0], nseq_back, seq_back, all_pos);
		if (pwm_check == -1)
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
		double auc_one = 0;
		int tproc_pred = 0;
		double fproc_pred = 0;
		int nthr_dist1 = nthr_dist[0] - 1;
		fp_here[0] = 0, tp_here[0] = 0;
		int nthr_here = 1;
		for (i = 0; i < nthr_dist[0]; i++)
		{
			if (fp_one[i] > fproc_pred && (i == nthr_dist1 || fp_one[i + 1] > fp_one[i]))
			{
				int tproc_cur = tp_one[i];
				double fproc_cur = (double)fp_one[i] / all_pos, fproc_cur_pauc = fproc_cur;
				if (fproc_cur >= fp2 || i == nthr_dist1)fproc_cur_pauc = fp2;
				double dauc = (tproc_cur + tproc_pred) * (fproc_cur_pauc - fproc_pred) / nseq_real / 2;
				tp_here[nthr_here] = (double)tp_one[i] / nseq_real;
				fp_here[nthr_here] = fproc_cur_pauc;
				nthr_here++;
				//fprintf(out,"%d\t%d\t%g\t%g\t%g\n", tproc_cur, tproc_pred, fproc_cur, fproc_pred, dauc);
				//fprintf(outq, "%g\t%f\n", fproc_cur, (double)tproc_cur / n_cnt_tot);
				if (fproc_cur <= fp2)auc_one += dauc;
				if (fproc_cur >= fp2)break;
				tproc_pred = tproc_cur;
				fproc_pred = fproc_cur;
			}
		}
		auc_one *= kauc;
		printf("%d\t%f\n", mot1+1, auc_one);
		fprintf(out_log1, "%d\t%f\n", mot1 + 1,  auc_one);
		//fprintf(out_roc, "%s\t%s\t%d\t%s\t%g\n", file_for, file_back, mot1 + 1, motif_name[mot1], auc_one);
		motifp[mot1].auc = auc_one;
		//	motifp[mot1].rank = -1;		
		motifp[mot1].num = mot1;
		{
			strcpy(file_roc1, file_roc);
			strcat(file_roc1, "_");
			char buf[5];
			sprintf(buf, "%d", mot1 + 1);
			strcat(file_roc1, buf);
		}
		if ((out_roc = fopen(file_roc1, "wt")) == NULL)
		{
			fprintf(out_roc, "Input file %s can't be opened!\n", file_roc1);
			exit(1);
		}
		for (i = 0; i < nthr_here; i++)
		{
			fprintf(out_roc, "%g\t%f\n", fp_here[i], tp_here[i]);
			if (fp_here[i] == fp2)break;
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
	double** fp_two;
	fp_two = new double* [2];
	if (fp_two == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (i = 0; i < 2; i++)
	{
		fp_two[i] = new double[nthr_dist_max];
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
		fseek(in_pwm, motifp[mot1].p, SEEK_SET);
		fread(&len_partner[0], sizeof(int), 1, in_pwm);
		int len_partner40 = len_partner[0] * 4;
		fread(pfm[0], sizeof(double), len_partner40, in_pwm);
		fread(pwm[0], sizeof(double), len_partner40, in_pwm);
		int nthr_dist_all1;
		fread(&nthr_dist_all1, sizeof(int), 1, in_pwm);
		fread(thr_all[0], sizeof(double), nthr_dist_all1, in_pwm);
		fread(fpr_all[0], sizeof(double), nthr_dist_all1, in_pwm);
		//if (strstr(motif_class[mot1], tfclass) == NULL)continue;	
		nthr_dist[0] = motifp[mot1].nth;
		{
			fprintf(out_mat, "%f\t%d\t", motifp[mot1].auc, motifp[mot1].rank);
		}
		for (i = 0; i < mot1c; i++)fprintf(out_mat, "\t");
		min[0] = raz[0] = 0;
		PWMScore(pwm[0], min[0], raz[0], len_partner[0]);
		matrix[0].mem_in(len_partner[0]);
		for (i = 0; i < len_partner[0]; i++)for (j = 0; j < OLIGNUM; j++)matrix[0].fre[i][j] = pfm[0][i][j];
		fprintf(out_auc, "\n");
		int mot2c = 0;
		for (mot2 = mot1 + 1; mot2 < n_motifs; mot2++)
		{
			fseek(in_pwm, motifp[mot2].p, SEEK_SET);
			fread(&len_partner[1], sizeof(int), 1, in_pwm);
			int len_partner41 = len_partner[1] * 4;
			fread(pfm[1], sizeof(double), len_partner41, in_pwm);
			fread(pwm[1], sizeof(double), len_partner41, in_pwm);
			int nthr_dist_all2;
			fread(&nthr_dist_all2, sizeof(int), 1, in_pwm);
			fread(thr_all[1], sizeof(double), nthr_dist_all2, in_pwm);
			fread(fpr_all[1], sizeof(double), nthr_dist_all2, in_pwm);
			nthr_dist[1] = motifp[mot2].nth;
			min[1] = raz[1] = 0;
			PWMScore(pwm[1], min[1], raz[1], len_partner[1]);
			matrix[1].mem_in(len_partner[1]);
			for (i = 0; i < len_partner[1]; i++)for (j = 0; j < OLIGNUM; j++)matrix[1].fre[i][j] = pfm[1][i][j];
			for (i = 0; i < 4; i++)pval_sim[i] = 1;
			double pvalue_similarity_tot = pfm_similarity(&matrix[0], &matrix[1], s_granul, s_overlap_min, s_ncycle_small, s_ncycle_large, pval_sim);
			if (pvalue_similarity_tot < p_crit)
			{
				pvalue_similarity_tot = -log10(pvalue_similarity_tot);
			}			
			else pvalue_similarity_tot = 1;
			sims[mot1c][mot2c] = pvalue_similarity_tot;
			for (n = 0; n < 2; n++)
			{
				for (j = 0; j <= nthr_dist[n]; j++)
				{
					fp_two[n][j] = 0;
					tp_two[n][j] = 0;
				}
			}
			//printf("Real %d\n",mot + 1);
			int pwm_check = PWM_rec_real(pwm, min, raz, nthr_dist, thr_all, fpr_all, tp_two, len_partner, nseq_real, seq_real);
			if (pwm_check == -1)
			{
				printf("Two motifs recognition error, real %d %d\n", mot1+1, mot2+1);
				exit(1);
			}
			//printf("Back %d\n", mot + 1);	
			int all_pos_two[2] = { 0,0 };
			pwm_check = PWM_rec_back(pwm, min, raz, nthr_dist, thr_all, fpr_all, fp_two, len_partner, nseq_back, seq_back, all_pos_two);
			if (pwm_check == -1)
			{
				printf("Two motifs recognition error, back %d %d\n", mot1 + 1, mot2 + 1);
				exit(1);
			}
			int all_pos;
			if (all_pos_two[0] > all_pos_two[1])all_pos = all_pos_two[0];
			else all_pos = all_pos_two[1];
			int nthr_dist_two = nthr_dist[0] + nthr_dist[1];
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
			/*for (i = 0; i < 30; i++)
			{
				printf("%d %g\t", tab[i].nfo,tab[i].fpr);
				if ((i + 1) % 10 == 0)printf("%d\n",i+1);
			}*/
			//printf("ROC %d\n", mot + 1);
			//qsort(tab, nseq_real, sizeof(tab[0]), compare_qbs);

			/*sum = 0;
			for (i = 0; i < nseq_razn; i++)
			{
				sum += tab[i].nfo;
				printf("%f\t%f\t%g\n", tab[i].err, (double)sum / nseq_real, tab[i].fpr / all_pos);
			}*/
			int tproc_pred = 0;
			double fproc_pred = 0;
			double auc_two = 0;
			int nthr_dist_two1 = nthr_dist_two - 1;
			int n_here = 0;
			for (i = 0; i < nthr_dist_two; i++)
			{
				if (tab[i].fpr > fproc_pred && (i == nthr_dist_two1 || tab[i + 1].fpr > tab[i].fpr))
				{
					int tproc_cur = tab[i].nfo;
					double fproc_cur = tab[i].fpr / all_pos, fproc_cur_pauc = fproc_cur;
					if (fproc_cur >= fp2 || i == nthr_dist_two1)fproc_cur_pauc = fp2;
					double dauc = (tproc_cur + tproc_pred) * (fproc_cur_pauc - fproc_pred) / nseq_real / 2;
					tp_here[n_here] = (double)tproc_cur / nseq_real;
					fp_here[n_here] = fproc_cur_pauc;
					n_here++;
					//fprintf(out,"%d\t%d\t%g\t%g\t%g\n", tproc_cur, tproc_pred, fproc_cur, fproc_pred, dauc);
					//fprintf(outq, "%g\t%f\n", fproc_cur, (double)tproc_cur / n_cnt_tot);
					if (fproc_cur <= fp2)auc_two += dauc;
					if (fproc_cur >= fp2)break;
					tproc_pred = tproc_cur;
					fproc_pred = fproc_cur;
				}
			}
			auc_two *= kauc;
			{
				strcpy(file_roc1, file_roc);
				strcat(file_roc1, "_");
				char buf[5];
				sprintf(buf, "%d", mot1 + 1);
				strcat(file_roc1, buf);
				strcat(file_roc1, "&");
				sprintf(buf, "%d", mot2 + 1);
				strcat(file_roc1, buf);
			}
			if ((out_roc = fopen(file_roc1, "wt")) == NULL)
			{
				fprintf(out_roc, "Input file %s can't be opened!\n", file_roc1);
				exit(1);
			}
			for (i = 0; i < n_here; i++)
			{
				fprintf(out_roc, "%g\t%f\n", fp_here[i], tp_here[i]);
				if (fp_here[i] == fp2)break;
			}
			fclose(out_roc);
			printf("%d\t%d\t%f\t%f\t%f", mot1+1, mot2+1, motifp[mot1].auc, motifp[mot2].auc, auc_two);
			if (sims[mot1c][mot2c] != 1)printf("\t%f", sims[mot1c][mot2c]);
			printf("\n");
			fprintf(out_log2, "%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n", mot1 + 1, mot2 + 1, motifp[mot1].rank, motifp[mot2].rank,sims[mot1c][mot2c], motifp[mot1].auc, motifp[mot2].auc, auc_two);
			double auc_max = Max(motifp[mot1].auc, motifp[mot2].auc);
			auctwei[mot1c][mot2c] = auc_two / auc_max;
			//if (auc_two > motifp[mot1].auc && auc_two > motifp[mot2].auc)
			{
				fprintf(out_auc, "%d\t%d\t%d\t%d\t", mot1+1, mot2+1, motifp[mot1].rank, motifp[mot2].rank);
				fprintf(out_auc, "%f\t%f\t%f\t%f\t%f\n", sims[mot1c][mot2c], motifp[mot1].auc, motifp[mot2].auc, auc_two, auctwei[mot1c][mot2c]);				
			}
			mot2c++;
			//	else fprintf(out_auc, "\t");
			fprintf(out_mat, "\t%f", auc_two);
			/*if (auc_two > motifp[mot1].auc && auc_two > motifp[mot2].auc)
			{
				fprintf(out_roc, "\t%s\t%s\t%s\t%s\t%g\t%g\t%g\n", file_for, file_back, motif_name[motifp[mot1].num], motif_name[motifp[mot2].num], motifp[mot1].auc, motifp[mot2].auc, auc_two);
				for (i = 0; i < n_here; i++)
				{
					fprintf(out_roc, "%g\t%f\n", fp_here[i], tp_here[i]);
					if (fp_here[i] == fp2)break;
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
	mot1c = 0;
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		{
			fprintf(out_mat, "%f\t%d\t", motifp[mot1].auc, motifp[mot1].rank);
			for (i = 0; i < mot1c; i++)fprintf(out_mat, "\t");
			int mot2c = 0;
			for (mot2 = mot1 + 1; mot2 < n_motifs; mot2++)
			{
				{
					fprintf(out_mat, "\t");
					if (sims[mot1c][mot2c] != 1)fprintf(out_mat, "%f", sims[mot1c][mot2c]);
					mot2c++;
				}
			}
			mot1c++;
			fprintf(out_mat, "\n");
		}
	}
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
	delete[] tp_here;
	delete[] fp_here;
	delete[] tp_one;
	delete[] fp_one;
	for (k = 0; k < n_motifs1; k++)
	{
		delete[] auctwei[k];
	}
	delete[] auctwei;
	for (k = 0; k < n_motifs1; k++)
	{
		delete[] sims[k];
	}
	delete[] sims;
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
}