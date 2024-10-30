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

struct acc {
	double auc;
	int rank;
	int nth;
	int num;
	long p;
};
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
int PWM_rec_back(double(&pwm)[2][MATLEN][OLIGNUM], double min[2], double raz[2], int nthr_dist[2], double** thr_all, double** fpr_all, int** fp_two, int** fp_nsites, int olen[2], int nseq, char*** seq)
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
			int len_pro1 = strlen(seq[0][n]);
			int len21 = len_pro1 - olen[k];
			best_inx[k] = nthr_dist[k];
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
					int index = nthr_dist[k];
					if (score >= thr_cr[k])
					{
						if (score >= thr_all[k][0])
						{
							index = 0;							
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
#include "pfm_families.h"
#include "pfm_classes.h"
#include "pfm_list.h" //motifs
#include "pfm_tfs.h" //TFs

int main(int argc, char* argv[])
{
	int i, j, k, n, mot1, mot2;
	char*** seq_real, *** seq_back;
	char file_for[ARGLEN], file_back[ARGLEN], partner_db[ARGLEN], partner_sim[ARGLEN];// , tfclass[80];path_fasta[ARGLEN], pfile_for[ARGLEN], pfile_back[ARGLEN],
	char file_mat[ARGLEN], file_auc[ARGLEN], file_log1[ARGLEN], file_log2[ARGLEN];
	FILE* out_log1, * out_auc, * out_log2, * out_mat, * in_pwm;

	if (argc != 11)
	{
		printf("%s 1,2input fasta foreground,background 3,4input binary files pwm_thresholds,pwm_similarity 5ntop matrices 6double ERRthresh 7,8,9,10output files pAUC_matrix, pAUC_list, log1, log2", argv[0]);
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
	strcpy(partner_sim, argv[4]); //h12hs, h12mm
	int ntop = atoi(argv[5]);// no. of top-scoring matrices
	double fp2 = atof(argv[6]); //ERR threshold for pAUC-PR		
	strcpy(file_mat, argv[7]);
	strcpy(file_auc, argv[8]);
	strcpy(file_log1, argv[9]);
	strcpy(file_log2, argv[10]);
	int* len_real, * len_back, nseq_real=0, nseq_back = 0;
	int olen_min = 8;
	int len_peak_max = 3000;
	//double fp2 = 0.005; //FPR threshold for pAUC		
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
	fprintf(out_log1, "#Motif\tTF\tMotif\tpAUC\tClass\tFamily\n");
	if ((out_log2 = fopen(file_log2, "wt")) == NULL)
	{
		fprintf(out_log2, "Input file %s can't be opened!\n", file_log2);
		exit(1);
	}
	fprintf(out_log2, "#Motif 1\t#Motif 2\tTF 1\tTF 2\tMotif 1\tMotif 2\tSimilarity\tpAUC 1\tpAUC 2\tpAUC 1&2\t\tClass 1\tFamily 1\t\tClass 2\tFamily 2\n");
	if ((out_mat = fopen(file_mat, "wt")) == NULL)
	{
		fprintf(out_mat, "Input file %s can't be opened!\n", file_mat);
		exit(1);
	}
	FILE* out_sim;
	if ((out_sim = fopen(partner_sim, "rb")) == NULL)
	{
		printf("Out file %s can't be opened!\n", partner_sim);
		return -1;
	}
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
	//motif library
	char  motif_name[MAXPAR][40], motif_class[MAXPAR][80], motif_family[MAXPAR][80], motif_class_inx[MAXPAR][8], motif_family_inx[MAXPAR][10], motif_tf[MAXPAR][30];
	int n_motifs, motif_library = -1;
	{
		char library_tag[NUM_LIBRARY][30] = { "h12core_hg38" , "h12core_mm10" , "h11core_hg38" , "h11core_mm10" , "dapseq", "jaspar24_at10", "jaspar24_dm6" };
		int motif_count_library[NUM_LIBRARY] = { 1420,1142,391,346,510,556,151 };
		for (i = 0; i < NUM_LIBRARY; i++)
		{
			if (strstr(partner_db, library_tag[i]) != NULL)
			{
				motif_library = i;
				n_motifs = motif_count_library[i];
				break;
			}
		}
		if (motif_library == -1)
		{
			printf("Wrong motif library %s\tAllowed motif library options:\n", partner_db);
			int i1 = NUM_LIBRARY - 1;
			for (i = 0; i < NUM_LIBRARY; i++)
			{
				printf("%s", library_tag[i]);
				if (i == i1)printf("\n");
				else printf("\t");
			}
			exit(1);
		}
		switch (motif_library) {
		case 0: {
			for (i = 0; i < n_motifs; i++)
			{			
				strcpy(motif_name[i], hs_core_12_names[i]);
				strcpy(motif_class[i], h12hs_classes[i]);
				strcpy(motif_family[i], h12hs_families[i]);
				strcpy(motif_tf[i], h12hs_tfs[i]);
			}
			break;
		}
		case 1: {
			for (i = 0; i < n_motifs; i++)
			{				
				strcpy(motif_name[i], mm_core_12_names[i]);
				strcpy(motif_class[i], h12mm_classes[i]);
				strcpy(motif_family[i], h12mm_families[i]);
				strcpy(motif_tf[i], h12mm_tfs[i]);
			}
			break;
		}
		case 5: {
			for (i = 0; i < n_motifs; i++)
			{
				strcpy(motif_name[i], jaspar24_plants_names[i]);
				strcpy(motif_class[i], jaspar24_plants_classes[i]);
				strcpy(motif_family[i], jaspar24_plants_families[i]);
				strcpy(motif_tf[i], jaspar24_plants_tfs[i]);				
			}
			break;
		}
		case 6: {
			for (i = 0; i < n_motifs; i++)
			{
				strcpy(motif_name[i], jaspar24_insects_names[i]);
				strcpy(motif_class[i], jaspar24_insects_classes[i]);
				strcpy(motif_family[i], jaspar24_insects_families[i]);
				strcpy(motif_tf[i], jaspar24_insects_tfs[i]);				
			}
			break;
		}
		default:
			break;
		}
	}
	for (i = 0; i < n_motifs; i++)
	{
		int sta, clen;
		clen = strlen(motif_class[i]);
		sta = 0;
		char st = '{', en = '}';
		for (j = 0; j < clen; j++)
		{
			if (motif_class[i][j] == st)
			{
				sta = j + 1;
				break;
			}
		}
		memset(motif_class_inx[i], '\0', sizeof(motif_class_inx[i]));
		k = 0;
		for (j = sta; motif_class[i][j] != en; j++)motif_class_inx[i][k++] = motif_class[i][j];
		motif_family_inx[i][k] = '\0';
		clen = strlen(motif_family[i]);
		sta = 0;
		for (j = 0; j < clen; j++)
		{
			if (motif_family[i][j] == st)
			{
				sta = j + 1;
				break;
			}
		}
		if (motif_library <= 1 || motif_library == 6)
		{
			memset(motif_family_inx[i], '\0', sizeof(motif_family_inx[i]));
			k = 0;
			for (j = sta; motif_family[i][j] != en; j++)motif_family_inx[i][k++] = motif_family[i][j];
			motif_family_inx[i][k] = '\0';
		}
		else
		{
			memset(motif_family_inx[i], '\0', sizeof(motif_family_inx[i]));
		}
	}
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
	int ntop1 = ntop - 1;
	double** auctwei;
	auctwei = new double* [ntop1];
	if (auctwei == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (i = 0; i < ntop1; i++)
	{
		auctwei[i] = new double[ntop1 - i];
		if (auctwei[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	double** sims;
	sims = new double*[ntop1];
	if (sims == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (i = 0; i < ntop1; i++)
	{
		sims[i] = new double[ntop1-i];
		if (sims[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	for (mot1 = 0; mot1 < ntop1; mot1++)for (mot2 = 0; mot2 < ntop1-mot1; mot2++)sims[mot1][mot2] = 1;
	//n_motifs = 50;
	int* tp_one, * fp_one;
	tp_one = new int[nthr_dist_max];
	if (tp_one == NULL) { puts("Out of memory..."); exit(1); }
	fp_one = new int[nthr_dist_max];
	if (fp_one == NULL) { puts("Out of memory..."); exit(1); }	
	double* recall;
	recall = new double[2*nthr_dist_max];
	if (recall == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	double* prec;
	prec = new double[2*nthr_dist_max];
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
	double min[2], raz[2];
	motifp[0].p = ftell(in_pwm);
	int n_motifs1 = n_motifs - 1;
	double nseq_fb = (double)nseq_real / nseq_back;
	double prec_exp = 0.5;
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
	for (mot1 = 0; mot1 < n_motifs; mot1++)all_pos[mot1] = all_pos_thr[mot1] = index_thr[mot1] = 0;
	for (mot1 = 0; mot1 < n_motifs; mot1++)fp_thr_rest[mot1] = 0;
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		//printf("%d\t%s\n", mot1 + 1, motif_class[mot1]);
		fread(&len_partner[0], sizeof(int), 1, in_pwm);
		int len_partner4 = len_partner[0] * 4;
		fread(pfm[0], sizeof(double), len_partner4, in_pwm);
		fread(pwm[0], sizeof(double), len_partner4, in_pwm);
		int nthr_dist_all=0;
		fread(&nthr_dist_all, sizeof(int), 1, in_pwm);
		fread(thr_all[0], sizeof(double), nthr_dist_all, in_pwm);
		fread(fpr_all[0], sizeof(double), nthr_dist_all, in_pwm);
		if(mot1!=n_motifs1)motifp[mot1+1].p = ftell(in_pwm);
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
			tp_one[j] = fp_one[j] = fp_nsites[0][j] = 0;
		}		
		//printf("Real %d\n",mot + 1);
		int pwm_check = PWM_rec_real_one(pwm[0], min[0], raz[0], nthr_dist[0], thr_all[0], fpr_all[0], tp_one, len_partner[0], nseq_real, seq_real);
		if (pwm_check == -1)
		{
			printf("One motif recognition error, real %s\n", motif_name[mot1]);
			exit(1);
		}
		//printf("Back %d\n", mot + 1);
		pwm_check = PWM_rec_back_one(pwm[0], min[0], raz[0], nthr_dist[0], thr_all[0], fpr_all[0], fp_one, fp_nsites[0], len_partner[0], nseq_back, seq_back, all_pos[mot1]);
		if (pwm_check == -1)
		{
			printf("One motif recognition error, back %s\n", motif_name[mot1]);
			exit(1);
		}
		all_pos_thr[mot1] = (int)(all_pos[mot1] * fp2);
		index_thr[mot1] = nthr_dist[0] - 1;
		int count_one = 0;
		for (i = 0; i < nthr_dist[0]; i++)
		{
			count_one += fp_nsites[0][i];
			//printf("FPsites %d FPpeak %d TPpeak %d\n", fp_nsites[0][i], fp_one[i],tp_one[i]);
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
		printf("%d\t%s\t%s\t%f\n", mot1 + 1, motif_tf[mot1], motif_name[mot1], auc_one);
		fprintf(out_log1, "%d\t%s\t%s\t%f\t%s\t%s\n", mot1 +1, motif_tf[mot1], motif_name[mot1], auc_one,motif_class[mot1],motif_family[mot1]);
		//fprintf(out_roc, "%s\t%s\t%d\t%s\t%g\n", file_for, file_back, mot1 + 1, motif_name[mot1], auc_one);
		motifp[mot1].auc =  auc_one;
	//	motifp[mot1].rank = -1;		
		motifp[mot1].num = mot1;
	/*	for (i = 0; i < nthr_here; i++)
		{
			fprintf(out_roc, "%g\t%f\n", prec[i], recall[i]);
		}*/
	}	
	qsort(motifp, n_motifs, sizeof(motifp[0]), compare_auc);
	for (i = 0; i < n_motifs; i++)motifp[i].rank = i + 1;
	qsort(motifp, n_motifs, sizeof(motifp[0]), compare_num);
	{
		int* simn;
		simn = new int[n_motifs1];
		if (simn == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
		double* simv;
		simv = new double[n_motifs1];
		if (simv == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }		
		{
			int sim_here, mot1k=0;
			for (mot1 = 0; mot1 < n_motifs; mot1++)
			{								
				fread(&sim_here, sizeof(int), 1, out_sim);
				//printf("%d Rank %d\t%s %s\t%d similar\n", mot1 + 1, motifp[mot1].rank, motif_name[mot1], motif_tf[mot1], sim_here);
				if (sim_here > 0)
				{
					fread(simn, sizeof(double), sim_here, out_sim);
					fread(simv, sizeof(double), sim_here, out_sim);				
					if (motifp[mot1].rank <= ntop)
					{
						int mot2k = 0, is = 0, isend=sim_here-1;
						for (mot2 = mot1 + 1; mot2 < n_motifs; mot2++)
						{
							if (simn[is] == mot2)
							{
								if (motifp[mot2].rank <= ntop)
								{
									sims[mot1k][mot2k] = simv[is];
								}
								//printf("%d(%d,%d) Rank2nd %d\t%s %s\t%s %s\t%f\n", is+1,mot1k+1,mot2k+1, motifp[mot2].rank,motif_name[mot1], motif_name[mot2], motif_tf[mot1], motif_tf[mot2], simv[i]);
								if (is == isend)break;
								is++;
							}
							if (motifp[mot2].rank <= ntop)mot2k++;
						}
						//int yy = 1;
					}
				}				
				if (motifp[mot1].rank <= ntop)mot1k++;
			}
		}
		delete[] simv;
		delete[] simn;
	}
	fclose(out_sim);
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
	tab = new qbs[2*nthr_dist_max];
	if (tab == NULL) { puts("Out of memory..."); exit(1); }	
	fprintf(out_mat, "%s\t%s\t\t\t\tpAUC",file_for,file_back);
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)fprintf(out_mat, "\t%f", motifp[mot1].auc);
	}
	fprintf(out_mat, "\n\t\t\t\t\tRank");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)fprintf(out_mat, "\t%d", motifp[mot1].rank);
	}
	fprintf(out_mat, "\n\t\t\t\t\tClass index");			
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{			
			fprintf(out_mat, "\t%s", motif_class_inx[mot1]);
		}
	}
	fprintf(out_mat, "\n\t\t\t\t\tFamily index");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "\t%s", motif_family_inx[mot1]);
		}
	}
	fprintf(out_mat, "\n\t\t\t\t\tTF");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "\t%s", motif_tf[mot1]);
		}
	}
	fprintf(out_mat, "\npAUC\tRank\tClass index\tFamily index\tTF\tMotif");				
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "\t%s", motif_name[mot1]);
		}
	}
	fprintf(out_mat, "\n");
	fprintf(out_auc, "Rank 1\tRank 2\tTF 1\tTF 2\tMotif 1\tMotif 2\tSimilarity\tpAUC 1\tpAUC 2\tpAUC 1&2\tRatio\t\tTF class 1\tTF class 2\t\tTF family 1\tTF family 2\n");	
	int mot1c = 0;
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank > ntop)continue;
		fseek(in_pwm, motifp[mot1].p, SEEK_SET);
		fread(&len_partner[0], sizeof(int), 1, in_pwm);
		int len_partner40 = len_partner[0] * 4;
		fread(pfm[0], sizeof(double), len_partner40, in_pwm);
		fread(pwm[0], sizeof(double), len_partner40, in_pwm);
		int nthr_dist_all1=0;
		fread(&nthr_dist_all1, sizeof(int), 1, in_pwm);		
		fread(thr_all[0], sizeof(double), nthr_dist_all1, in_pwm);
		fread(fpr_all[0], sizeof(double), nthr_dist_all1, in_pwm);
		//if (strstr(motif_class[mot1], tfclass) == NULL)continue;	
		nthr_dist[0] = motifp[mot1].nth;
		{			
			fprintf(out_mat, "%f\t%d\t%s\t%s\t%s\t%s\t", motifp[mot1].auc, motifp[mot1].rank, motif_class_inx[mot1], motif_family_inx[mot1], motif_tf[mot1], motif_name[mot1]);
		}
		for(i=0;i<mot1c;i++)fprintf(out_mat, "\t");
		min[0] = raz[0] = 0;
		PWMScore(pwm[0], min[0], raz[0], len_partner[0]);
		int mot2c = 0;
		for (mot2 = mot1+1; mot2 < n_motifs; mot2++)
		{
			if (motifp[mot2].rank > ntop)continue;
			fseek(in_pwm, motifp[mot2].p, SEEK_SET);
			fread(&len_partner[1], sizeof(int), 1, in_pwm);
			int len_partner41 = len_partner[1] * 4;
			fread(pfm[1], sizeof(double), len_partner41, in_pwm);
			fread(pwm[1], sizeof(double), len_partner41, in_pwm);
			/*for (i = 0; i < len_partner[1]; i++)
			{
				for (j = 0; j < OLIGNUM; j++)
				{
					pfm[1][i][j] = pfm1[i][j];
					pwm[1][i][j] = pwm1[i][j];
				}
			}*/
			int nthr_dist_all2=0;
			fread(&nthr_dist_all2, sizeof(int), 1, in_pwm);			
			fread(thr_all[1], sizeof(double), nthr_dist_all2, in_pwm);
			fread(fpr_all[1], sizeof(double), nthr_dist_all2, in_pwm);
			//if (strstr(motif_class[mot2], tfclass) == NULL)continue;			
			nthr_dist[1] = motifp[mot2].nth;					
			min[1] = raz[1] = 0;
			PWMScore(pwm[1], min[1], raz[1], len_partner[1]);
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
			int pwm_check = PWM_rec_real(pwm, min,raz,nthr_dist, thr_all, fpr_all, tp_two, len_partner, nseq_real, seq_real);
			if (pwm_check == -1)
			{
				printf("Two motifs recognition error, real %s %s\n", motif_name[mot1], motif_name[mot2]);
				exit(1);
			}
			//printf("Back %d\n", mot + 1);				
			pwm_check = PWM_rec_back(pwm, min, raz, nthr_dist, thr_all, fpr_all, fp_two, fp_nsites, len_partner, nseq_back, seq_back);
			if (pwm_check == -1)
			{
				printf("Two motifs recognition error, back %s %s\n", motif_name[mot1], motif_name[mot2]);
				exit(1);
			}
			int nthr_dist_two = nthr_dist[0] + nthr_dist[1];
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
					//	if(i<20)printf("%d %d\t%d\t%d\n", j,i,tab[k].tpr, tab[k].fpr);
					k++;
				}
			}
			qsort(tab, nthr_dist_two, sizeof(tab[0]), compare_qbs);
			//printf("ROC %d\n", mot + 1);
			int count_two = 0, count_pred = 0;
			int nthr_dist_two1 = nthr_dist_two - 1;
			int index_thr_two = 0;
			double fp_rest_two = 0;
			for (i = 0; i < nthr_dist_two; i++)
			{
				count_two += tab[i].fps;
				if (tab[i].fps > 0 && (i == nthr_dist_two1 || tab[i + 1].err != tab[i].err))
				{
				//	printf("ERR %f Count %d FPsites %d FPpeak %d TPpeak %d\n", tab[i].err, count_two, tab[i].fps, tab[i].fpr, tab[i].tpr);
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
			printf("%s\t%s\t%s\t%s\t%f\t%f\t%f", motif_tf[mot1], motif_tf[mot2], motif_name[mot1], motif_name[mot2], motifp[mot1].auc, motifp[mot2].auc, auc_two);
			if(sims[mot1c][mot2c]!=1)printf("\t%f", sims[mot1c][mot2c]);
			printf("\n");			
			fprintf(out_log2, "%d\t%d\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t\t", mot1 +1, mot2 +1, motif_tf[mot1], motif_tf[mot2], motif_name[mot1], motif_name[mot2], sims[mot1c][mot2c],motifp[mot1].auc, motifp[mot2].auc, auc_two);
			fprintf(out_log2, "%s\t%s\t\t%s\t%s\n", motif_class[mot1], motif_family[mot1], motif_class[mot2], motif_family[mot2]);
			double auc_max = Max(motifp[mot1].auc, motifp[mot2].auc);
			auctwei[mot1c][mot2c] = auc_two / auc_max;
			if (auc_two > motifp[mot1].auc && auc_two > motifp[mot2].auc)
			{
				int mot1x, mot2x;
				if (motifp[mot1].auc >= motifp[mot2].auc) 
				{ 
					mot1x = mot1; mot2x = mot2; 					
				}
				else 
				{ 
					mot1x = mot2; mot2x = mot1; 
				}
				fprintf(out_auc, "%d\t%d\t%s\t%s\t%s\t%s\t", motifp[mot1x].rank, motifp[mot2x].rank, motif_tf[mot1x], motif_tf[mot2x], motif_name[mot1x], motif_name[mot2x]);
				fprintf(out_auc, "%f\t%f\t%f\t%f\t%f\t", sims[mot1c][mot2c], motifp[mot1x].auc, motifp[mot2x].auc, auc_two, auctwei[mot1c][mot2c]);
				fprintf(out_auc, "\t%s\t%s\t\t%s\t%s\n", motif_class[mot1x], motif_class[mot2x], motif_family[mot1x], motif_family[mot2x]);
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
	fprintf(out_mat, "%s\t%s\t\t\t\tpAUC", file_for, file_back);
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)fprintf(out_mat, "\t%f", motifp[mot1].auc);
	}
	fprintf(out_mat, "\n\t\t\t\t\tRank");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)fprintf(out_mat, "\t%d", motifp[mot1].rank);
	}
	fprintf(out_mat, "\n\t\t\t\t\tClass index");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{			
			fprintf(out_mat, "\t%s", motif_class_inx[mot1]);
		}
	}
	fprintf(out_mat, "\n\t\t\t\t\tFamily index");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "\t%s", motif_family_inx[mot1]);
		}
	}
	fprintf(out_mat, "\n\t\t\t\t\tTF");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "\t%s", motif_tf[mot1]);
		}
	}
	fprintf(out_mat, "\npAUC\tRank\tClass index\tFamily index\tTF\tMotif");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "\t%s", motif_name[mot1]);
		}
	}
	fprintf(out_mat, "\n");
	k = 0;
	mot1c = 0;
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank > ntop)continue;
		{
			fprintf(out_mat, "%f\t%d\t%s\t%s\t%s\t%s\t", motifp[mot1].auc, motifp[mot1].rank, motif_class_inx[mot1], motif_family_inx[mot1], motif_tf[mot1], motif_name[mot1]);
		}
		for (i = 0; i < mot1c; i++)fprintf(out_mat, "\t");	
		int mot2c = 0;
		for (mot2 = mot1 + 1; mot2 < n_motifs; mot2++)
		{
			if (motifp[mot2].rank <= ntop)
			{
				fprintf(out_mat, "\t%f", auctwei[mot1c][mot2c++]);
			}
		}
		mot1c++;
		fprintf(out_mat, "\n");
	}
	fprintf(out_mat, "%s\t%s\t\t\t\tpAUC", file_for, file_back);
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)fprintf(out_mat, "\t%f", motifp[mot1].auc);
	}
	fprintf(out_mat, "\n\t\t\t\t\tRank");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)fprintf(out_mat, "\t%d", motifp[mot1].rank);
	}
	fprintf(out_mat, "\n\t\t\t\t\tClass index");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "\t%s", motif_class_inx[mot1]);
		}
	}
	fprintf(out_mat, "\n\t\t\t\t\tFamily index");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "\t%s", motif_family_inx[mot1]);
		}
	}
	fprintf(out_mat, "\n\t\t\t\t\tTF");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "\t%s", motif_tf[mot1]);
		}
	}
	fprintf(out_mat, "\npAUC\tRank\tClass index\tFamily index\tTF\tMotif");
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "\t%s", motif_name[mot1]);
		}
	}
	fprintf(out_mat, "\n");
	mot1c = 0;
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		if (motifp[mot1].rank <= ntop)
		{
			fprintf(out_mat, "%f\t%d\t%s\t%s\t%s\t%s\t", motifp[mot1].auc, motifp[mot1].rank, motif_class_inx[mot1], motif_family_inx[mot1], motif_tf[mot1], motif_name[mot1]);
			for (i = 0; i < mot1c; i++)fprintf(out_mat, "\t");
			int mot2c = 0;
			for (mot2 = mot1 + 1; mot2 < n_motifs; mot2++)
			{
				if (motifp[mot2].rank <= ntop)
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
	for (k = 0; k < 2; k++)
	{
		delete[] fp_nsites[k];
	}
	delete[] fp_nsites;
	delete[] fpr_all;
	delete[] len_real;
	delete[] len_back;	
	delete[] tp_two;
	delete[] recall;
	delete[] prec;
	delete[] tp_one;
	delete[] fp_one;	
	delete[] index_thr;
	delete[] all_pos;
	delete[] all_pos_thr;
	delete[] fp_thr_rest;
	for (k = 0; k < ntop1; k++)
	{
		delete[] auctwei[k];
	}
	delete[] auctwei;
	for (k = 0; k < ntop1; k++)
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