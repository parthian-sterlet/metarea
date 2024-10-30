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
	int i, j, k, n;
	char*** seq_real, *** seq_back;
	char file_for[ARGLEN], file_back[ARGLEN], partner_db[2][ARGLEN];//path_fasta[ARGLEN], pfile_for[ARGLEN], pfile_back[ARGLEN],
	char file_roc[3][ARGLEN], file_auc[ARGLEN];
	FILE* out_roc[3], * out_auc, * in_pwm[2];

	if (argc != 10)
	{
		printf("%s 1,2input fasta foreground,background fasta 3,4input binary files library1,2 5double ERRthresh 6,7,8,9output files pAUC, PR curves PWM1, PWM2, PWM1&PWM2", argv[0]);
		return -1;
	}
	//strcpy(path_fasta, argv[1]);
	strcpy(file_for, argv[1]);
	strcpy(file_back, argv[2]);
	//strcpy(pfile_for, path_fasta);
	//strcpy(pfile_back, path_fasta);
	//strcat(pfile_for, file_for);
	//strcat(pfile_back, file_back);
	strcpy(partner_db[0], argv[3]); //h12hs, h12mm
	strcpy(partner_db[1], argv[4]); //h12hs, h12mm	
	double fp2 = atof(argv[5]); //ERR threshold for pAUC-PR	
	strcpy(file_auc, argv[6]);
	strcpy(file_roc[0], argv[7]);
	strcpy(file_roc[1], argv[8]);
	strcpy(file_roc[2], argv[9]);
	int* len_real, * len_back, nseq_real = 0, nseq_back = 0;
	int olen_min = 8;
	int len_peak_max = 3000;
	//double fp2 = 0.01; //FPR threshold for pAUC	
	double fp2_lg = -log10(fp2);
	int s_overlap_min = 6, s_ncycle_small = 1000, s_ncycle_large = 10000;//for permutation(motif_comparison) min_length_of_alignment, no. of permutation (test & detailed)
	double s_granul = 0.001;//for permutation(motif_comparison) okruglenie 4astotnyh matric	
	double p_crit = 0.05;// threshold for similarity of matrices
	double pval_sim[4];
	matrices matrix[2];

	/*if ((outlog = fopen(file_log, "at")) == NULL)
	{
		fprintf(outlog, "Input file %s can't be opened!\n", file_log);
		exit(1);
	}*/
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
	for (i = 0; i < 2; i++)
	{
		if ((in_pwm[i] = fopen(partner_db[i], "rb")) == NULL)
		{
			printf("Input file %s can't be opened!\n", partner_db[i]);
			return -1;
		}
	}
	int len_partner[2], nthr_dist[2];
	double pwm[2][MATLEN][OLIGNUM];
	double pfm[2][MATLEN][OLIGNUM];
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
	double** recall_1;
	recall_1 = new double* [2];
	if (recall_1 == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		recall_1[i] = new double[nthr_dist_max];
		if (recall_1[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	double** prec_1;
	prec_1 = new double* [2];
	if (prec_1 == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		prec_1[i] = new double[nthr_dist_max];
		if (prec_1[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
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

	for (n = 0; n < 2; n++)
	{
		fread(&len_partner[n], sizeof(int), 1, in_pwm[n]);
		fread(pfm[n], sizeof(double), 4 * len_partner[n], in_pwm[n]);
		fread(pwm[n], sizeof(double), 4 * len_partner[n], in_pwm[n]);
		int nthr_dist_all;
		fread(&nthr_dist_all, sizeof(int), 1, in_pwm[n]);
		fread(thr_all[n], sizeof(double), nthr_dist_all, in_pwm[n]);
		fread(fpr_all[n], sizeof(double), nthr_dist_all, in_pwm[n]);
		nthr_dist[n] = 0;
		for (i = 0; i < nthr_dist_all; i++)
		{
			nthr_dist[n]++;
			if (fpr_all[n][i] < fp2_lg)break;
		}
		fclose(in_pwm[n]);
		matrix[n].mem_in(len_partner[n]);
		for (i = 0; i < len_partner[n]; i++)for (j = 0; j < OLIGNUM; j++)matrix[n].fre[i][j] = pfm[n][i][j];
	}
	for (i = 0; i < 4; i++)pval_sim[i] = 1;
	double pvalue_similarity_tot = pfm_similarity(&matrix[0], &matrix[1], s_granul, s_overlap_min, s_ncycle_small, s_ncycle_large, pval_sim);
	if (pvalue_similarity_tot < p_crit)
	{
		pvalue_similarity_tot = -log10(pvalue_similarity_tot);
	}
	else pvalue_similarity_tot = 1;
	double min[2] = { 0, 0 }, raz[2] = { 0,0 };
	for (n = 0; n < 2; n++)
	{
		PWMScore(pwm[n], min[n], raz[n], len_partner[n]);
	}
	double auc_one[2] = { 0,0 }, prauc_one[2] = { 0,0 };
	int n_here1[2] = { 0,0 };
	int index_thr[2] = { 0,0 }, index_thr_two=0;
	double fp_thr_rest[2], fp_rest_two =0;
	double nseq_fb = (double)nseq_real / nseq_back;
	double prec_exp = 0.5;
	int all_pos[2] = { 0,0 };
	int all_pos_thr[2];
	for (n = 0; n < 2; n++)
	{
		for (j = 0; j <= nthr_dist[n]; j++)tp_one[j] = fp_one[j] = fp_nsites[0][j] = fp_nsites[1][j] = 0;		
		//printf("Real %d\n",mot + 1);
		int pwm_check = PWM_rec_real_one(pwm[n], min[n], raz[n], nthr_dist[n], thr_all[n], fpr_all[n], tp_one, len_partner[n], nseq_real, seq_real);
		if (pwm_check == -1)
		{
			printf("Motif %d recognition error, real\n", n + 1);
			exit(1);
		}
		//printf("Back %d\n", mot + 1);
		pwm_check = PWM_rec_back_one(pwm[n], min[n], raz[n], nthr_dist[n], thr_all[n], fpr_all[n], fp_one, fp_nsites[n], len_partner[n], nseq_back, seq_back, all_pos[n]);
		if (pwm_check == -1)
		{
			printf("Motif %d recognition error, back\n", n);
			exit(1);
		}
		all_pos_thr[n] = (int)(all_pos[n] * fp2);
		index_thr[n] = nthr_dist[n] - 1;
		int count_one = 0;		
		for (i = 0; i < nthr_dist[n]; i++)
		{
			count_one += fp_nsites[n][i];
			//printf("FPsites %d FPpeak %d TPpeak %d\n", fp_nsites[n][i], fp_one[i],tp_one[i]);
			if (count_one >= all_pos_thr[n]) //|| tp_one[i] >= nseq_real_thr
			{
				index_thr[n] = i;
				//fp_thr[n] = fp_one[i];				
				fp_thr_rest[n] = 1 - (double)(count_one - all_pos_thr[n]) / fp_nsites[n][i];
				break;
			}
		}
		/*for (i = 0; i <= index_thr[n]; i++)
		{
			if ((i + 1) % 10 == 0)printf("%d\n", i + 1);
			printf("%5d %d\t", tp_one[i], fp_one[i]);
		}*/
		//printf("ROC %d\n", mot + 1);
		int nthr_dist1 = nthr_dist[n] - 1;
		int index1 = index_thr[n] - 1;
		double prec_pred = 1;// (double)tp_one[0] / ((double)tp_one[0] + nseq_fb * fp_one[0]);
		recall_1[n][0] = 0;//(double)tp_one[0]/nseq_real, 
		prec_1[n][0] = prec_pred;
		n_here1[n] = 1;
		for (i = 0; i <= index_thr[n]; i++)
		{
			int dtp = tp_one[i];
			int i1 = i - 1;
			if (i > 0)dtp -= tp_one[i1];			
			if (dtp > 0)
			{				
				double dfpi = (double)fp_one[i];
				double dtpi = (double)tp_one[i];
				if (i > 0)
				{
					dfpi -= fp_one[i1];
					dtpi -= tp_one[i1];
				}
				if (i == index_thr[n])
				{
					dtpi *= fp_thr_rest[n];
					dfpi *= fp_thr_rest[n];
				}
				double fp_onec = dfpi;
				double tp_onec = dtpi;
				if (i > 0)
				{
					fp_onec += (double)fp_one[i1];
					tp_onec += (double)tp_one[i1];
				}
				double prec_cur = tp_onec / (tp_onec + nseq_fb * fp_onec);
				double prec_av = (prec_pred + prec_cur)/2;
				double dauc = dtpi * (prec_av - prec_exp);
				recall_1[n][n_here1[n]] = tp_onec / nseq_real;
				prec_1[n][n_here1[n]] = prec_cur;				
				prec_pred = prec_cur;
				prauc_one[n] += dauc;
				n_here1[n]++;
			}
		}
		prauc_one[n] *= 2;
		prauc_one[n] /= nseq_real;
		printf("%s\t%s\t%d\t%s\t%g\t%g\n", file_for, file_back, n + 1, partner_db[n], auc_one[n], prauc_one[n]);
	}
	int all_pos_thr_two = (int)((all_pos[0]+ all_pos[1]) * fp2 / 2);
	//printf("FPR thr - 1 %d 2 %d 12 %f\n", fp_thr[0],fp_thr[1], fp_thr_two);
	for (n = 0; n < 2; n++)for (j = 0; j <= nthr_dist[n]; j++)tp_two[n][j] = 0;
	for (n = 0; n < 2; n++)for (j = 0; j <= nthr_dist[n]; j++)fp_two[n][j] = 0;
	//printf("Real %d\n",mot + 1);
	int pwm_check = PWM_rec_real(pwm, min, raz, nthr_dist, thr_all, fpr_all, tp_two, len_partner, nseq_real, seq_real);
	if (pwm_check == -1)
	{
		printf("Motifs recognition error, real\n");
		exit(1);
	}
	//printf("Back %d\n", mot + 1);		
	pwm_check = PWM_rec_real(pwm, min, raz, nthr_dist, thr_all, fpr_all, fp_two, len_partner, nseq_back, seq_back);
	if (pwm_check == -1)
	{
		printf("Motifs recognition error, back\n");
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
			//	if(i<20)printf("%d %d\t%d\t%d\n", j,i,tab[k].tpr, tab[k].fpr);
			k++;
		}
	}
	qsort(tab, nthr_dist_two, sizeof(tab[0]), compare_qbs);
	//	printf(" %d\t%d\n", tab[0].tpr, tab[0].fpr);
	int count_two = 0, count_pred = 0;
	int nthr_dist_two1 = nthr_dist_two - 1;
	for (i = 0; i < nthr_dist_two; i++)
	{		
		count_two += tab[i].fps;
		if (tab[i].fps > 0 && (i == nthr_dist_two1 || tab[i+1].err != tab[i].err))
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
	double* recall;
	recall = new double[nthr_dist_two];
	if (recall == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	double* prec;
	prec = new double[nthr_dist_two];
	if (prec == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	double prauc_two = 0;
	int n_herepr = 1;
	recall[0] = 0;
	prec[0] = 1;
	double prec_pred = 1;
	double tpsum = 0, fpsum = 0;
	double dtp = 0, dfp = 0;	
	for (i = 0; i <= index_thr_two; i++)
	{		
		dtp += (double)tab[i].tpr;
		dfp += (double)tab[i].fpr;
		if (dtp > 0 && (i == index_thr_two || tab[i+1].err!= tab[i].err))
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
			recall[n_herepr] = tpsum / nseq_real;
			prec[n_herepr] = prec_cur;
			prec_pred = prec_cur;
			n_herepr++;
			prauc_two += dauc;
			dtp = 0;
			dfp = 0;			
		}
	}
	prauc_two *= 2;
	prauc_two /= nseq_real;
	{
		//double auc_max = Max(auc_one[0], auc_one[1]);
		double prauc_max = Max(prauc_one[0], prauc_one[1]);
		printf("%s\t%s\t%s\t%s\t%g\t%g\t%f\n", file_for, file_back, partner_db[0], partner_db[1], prauc_two, prauc_two / prauc_max, pvalue_similarity_tot);
	}
	//	fprintf(outlog, "%s\t%s\t%s\t%s\t%g\t%g\t%g\n", file_for, file_back, partner_db[0], partner_db[1], auc_one[0], auc_one[1], auc_two);
	if ((out_auc = fopen(file_auc, "wt")) == NULL)
	{
		fprintf(out_auc, "Input file %s can't be opened!\n", file_auc);
		exit(1);
	}
//	double auc_max = Max(auc_one[0], auc_one[1]);
	double prauc_max = Max(prauc_one[0], prauc_one[1]);
	int m0, m1;
	m0 = 0; m1 = 1;
	//if (auc_one[0] > auc_one[1]) { m0 = 0; m1 = 1; }
	//else { m0 = 1; m1 = 0; }
	fprintf(out_auc, "%s\t%s\t%s\t%s\t", file_for, file_back, partner_db[m0], partner_db[m1]);
	//fprintf(out_auc, "%g\t%g\t%g\t%f\t\t", auc_one[m0], auc_one[m1], auc_two, auc_two / auc_max);
	fprintf(out_auc, "%g\t%g\t%g\t%f\t", prauc_one[m0], prauc_one[m1], prauc_two, prauc_two / prauc_max);
	fprintf(out_auc, "%f\n", pvalue_similarity_tot);
	fclose(out_auc);
	for (n = 0; n < 3; n++)
	{
		if ((out_roc[n] = fopen(file_roc[n], "wt")) == NULL)
		{
			fprintf(out_roc[n], "Input file %s can't be opened!\n", file_roc[n]);
			exit(1);
		}
	}
	/*for (n = 0; n < 2; n++)
	{
		fprintf(out_roc[n], "%s\t%s\t%s\t%g\n", file_for, file_back, partner_db[n], auc_one[n]);
		for (i = 0; i < n_here1[n]; i++)
		{
			fprintf(out_roc[n], "%g\t%f\n", fp_here1[n][i], tp_here1[n][i]);
		}
	}
	*/
	for (n = 0; n < 2; n++)
	{
		fprintf(out_roc[n], "%s\t%s\t%s\t%g\n", file_for, file_back, partner_db[n], prauc_one[n]);
		for (i = 0; i < n_here1[n]; i++)
		{
			fprintf(out_roc[n], "%g\t%f\n", recall_1[n][i], prec_1[n][i]);
		}
	}
	/*fprintf(out_roc[2], "%s\t%s\t%s\t%s\t%g\t%g\t%g\n", file_for, file_back,partner_db[0],partner_db[1],auc_one[0],auc_one[1], auc_two);
	for (i = 0; i < n_here; i++)
	{
		fprintf(out_roc[2], "%g\t%f\n", fp_here[i], tp_here[i]);
	}*/
	fprintf(out_roc[2], "%s\t%s\t%s\t%s\t%g\t%g\t%g\n", file_for, file_back, partner_db[0], partner_db[1], prauc_one[0], prauc_one[1], prauc_two);
	for (i = 0; i < n_herepr; i++)
	{
		fprintf(out_roc[2], "%g\t%f\n", recall[i], prec[i]);
	}
	for (n = 0; n < 3; n++)fclose(out_roc[n]);
	//fclose(outlog);
	for (k = 0; k < 2; k++)fclose(in_pwm[k]);
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
	delete[] prec;
	delete[] recall;
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
		delete[] fp_nsites[k];
	}
	delete[] fp_nsites;
	for (k = 0; k < 2; k++)
	{
		delete[] recall_1[k];
	}
	delete[] recall_1;
	for (k = 0; k < 2; k++)
	{
		delete[] prec_1[k];
	}
	delete[] prec_1;
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