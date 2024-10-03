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

int main(int argc, char* argv[])
{
	int i, j, k, n;
	char**** seq_real, **** seq_back;
	char file_for[2][ARGLEN], file_back[2][ARGLEN], partner_db[2][ARGLEN];
	char file_out[ARGLEN];
	FILE * out_auc, * in_pwm[2];//* out_roc[4],

	if (argc != 8)
	{
		printf("%s 1,2,3input 1st motif foreground fasta,background fasta, binary file, 4,5,6input 1st motif foreground fasta,background fasta, binary file, 7output txt file", argv[0]);
		return -1;
	}
	//strcpy(path_fasta, argv[1]);
	strcpy(file_for[0], argv[1]);
	strcpy(file_back[0], argv[2]);
	strcpy(partner_db[0], argv[3]); //h12hs, h12mm
	strcpy(file_for[1], argv[4]);
	strcpy(file_back[1], argv[5]);
	strcpy(partner_db[1], argv[6]); //h12hs, h12mm	
	strcpy(file_out, argv[7]);
	//strcpy(file_roc[0], argv[6]);
	//strcpy(file_roc[1], argv[7]);
	//strcpy(file_roc[2], argv[8]);
	int** len_real, **len_back, nseq_real[2] = {0, 0}, nseq_back[2] = { 0, 0 };
	int olen_min = 8;
	int len_peak_max = 3000;
	double fp2 = 0.001; //FPR threshold for pAUC	
	double fp2_lg = -log10(fp2);
	matrices matrix[2];

	/*if ((outlog = fopen(file_log, "at")) == NULL)
	{
		fprintf(outlog, "Input file %s can't be opened!\n", file_log);
		exit(1);
	}*/
	len_real = new int*[2];
	if (len_real == NULL) { puts("Out of memory..."); exit(1); }
	len_back = new int* [2];
	if (len_back == NULL) { puts("Out of memory..."); exit(1); }
	//	printf("EvalSeq\n");
	for (i = 0; i < 2; i++)
	{
		EvalSeq(file_for[i], nseq_real[i], olen_min, len_peak_max);
		EvalSeq(file_back[i], nseq_back[i], olen_min, len_peak_max);
		len_real[i] = new int[nseq_real[i]];
		if (len_real[i] == NULL) { puts("Out of memory..."); exit(1); }
		len_back[i] = new int[nseq_back[i]];
		if (len_back[i] == NULL) { puts("Out of memory..."); exit(1); }
		//printf("EvalLen\n");
		EvalLen(file_for[i], len_real[i], olen_min, len_peak_max);
		EvalLen(file_back[i], len_back[i], olen_min, len_peak_max);
	}
	seq_real = new char*** [2];
	if (seq_real == NULL) { puts("Out of memory..."); exit(1); }
	for (n = 0; n < 2; n++)
	{
		seq_real[n] = new char** [2];
		if (seq_real[n] == NULL) { puts("Out of memory..."); exit(1); }
		for (k = 0; k < 2; k++)
		{
			seq_real[n][k] = new char* [nseq_real[n]];
			if (seq_real[n][k] == NULL) { puts("Out of memory..."); exit(1); }
			for (i = 0; i < nseq_real[n]; i++)
			{
				seq_real[n][k][i] = new char[len_real[n][i] + 1];
				if (seq_real[n][k][i] == NULL) { puts("Out of memory..."); exit(1); }				
			}
		}
	}
	seq_back = new char*** [2];
	if (seq_back == NULL) { puts("Out of memory..."); exit(1); }
	for (n = 0; n < 2; n++)
	{
		seq_back[n] = new char** [2];
		if (seq_back[n] == NULL) { puts("Out of memory..."); exit(1); }
		for (k = 0; k < 2; k++)
		{
			seq_back[n][k] = new char* [nseq_back[n]];
			if (seq_back[n][k] == NULL) { puts("Out of memory..."); exit(1); }
			for (i = 0; i < nseq_back[n]; i++)
			{
				seq_back[n][k][i] = new char[len_back[n][i] + 1];
				if (seq_back[n][k][i] == NULL) { puts("Out of memory..."); exit(1); }
			}
		}
	}
	for (n = 0; n < 2; n++)
	{
		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < nseq_real[n]; i++)
			{
				memset(seq_real[n][k][i], '\0', len_real[n][i] + 1);
			}
			for (i = 0; i < nseq_back[n]; i++)
			{
				memset(seq_back[n][k][i], '\0', len_back[n][i] + 1);
			}
		}
	}
	for (i = 0; i < 2; i++)
	{		
		ReadSeq(file_for[i], nseq_real[i], len_real[i], seq_real[i], olen_min, len_peak_max);
		ReadSeq(file_back[i], nseq_back[i], len_back[i], seq_back[i], olen_min, len_peak_max);
	}
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
	double** tp_here1;
	tp_here1 = new double* [4];
	if (tp_here1 == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 4; i++)
	{
		tp_here1[i] = new double[nthr_dist_max];
		if (tp_here1[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	double** fp_here1;
	fp_here1 = new double* [4];
	if (fp_here1 == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 4; i++)
	{
		fp_here1[i] = new double[nthr_dist_max];
		if (fp_here1[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	int** tp_one;
	tp_one = new int* [4];
	if (tp_one == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 4; i++)
	{
		tp_one[i] = new int[nthr_dist_max];
		if (tp_one[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	int** fp_one;
	fp_one = new int* [4];
	if (fp_one == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 4; i++)
	{
		fp_one[i] = new int[nthr_dist_max];
		if (fp_one[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
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
	for (k = 0; k < 2; k++)fclose(in_pwm[k]);
	double min[2] = { 0, 0 }, raz[2] = { 0,0 };
	for (n = 0; n < 2; n++)
	{
		PWMScore(pwm[n], min[n], raz[n], len_partner[n]);
	}
	double auc_one[4] = { 0,0,0,0 };	
	int nmodel[4] = { 0,1,0,1 };
	int nsequ[4] = { 0,1,1,0 };
	int n_here1[4] = { 0,0,0,0 };	
	for (n = 0; n < 4; n++)
	{
		int all_pos = 0;
		int nm = nmodel[n], ns = nsequ[n];
		for (j = 0; j <= nthr_dist[nm]; j++)tp_one[n][j] = fp_one[n][j] = 0;
		//printf("Real %d\n",mot + 1);
		int pwm_check = PWM_rec_real_one(pwm[nm], min[nm], raz[nm], nthr_dist[nm], thr_all[nm], fpr_all[nm], tp_one[n], len_partner[nm], nseq_real[ns], seq_real[ns]);
		if (pwm_check == -1)
		{
			printf("Motif %d recognition error, real\n", n + 1);
			exit(1);
		}
		//printf("Back %d\n", mot + 1);
		pwm_check = PWM_rec_back_one(pwm[nm], min[nm], raz[nm], nthr_dist[nm], thr_all[nm], fpr_all[nm], fp_one[n], len_partner[nm], nseq_back[ns], seq_back[ns], all_pos);
		if (pwm_check == -1)
		{
			printf("Motif %d recognition error, back\n", n);
			exit(1);
		}
		/*for (i = 0; i <= nthr_dist[n]; i++)
		{
			if ((i + 1) % 10 == 0)printf("%d\n", i + 1);
			printf("%5d %d\t", tp_one[i], fp_one[i]);
		}*/
		//printf("ROC %d\n", mot + 1);
		int tproc_pred = 0;
		double fproc_pred = 0;
		int nthr_dist1 = nthr_dist[nm] - 1;
		fp_here1[n][0] = 0, tp_here1[n][0] = 0;
		n_here1[n] = 1;
		for (i = 0; i < nthr_dist[nm]; i++)
		{
			double fproc_cur = (double)fp_one[n][i] / all_pos;
			if (fproc_cur > fproc_pred && (i == nthr_dist1 || fp_one[n][i + 1] > fp_one[n][i]))
			{
				int tproc_cur = tp_one[n][i];
				double fproc_cur_pauc = fproc_cur;
				if (fproc_cur >= fp2 || i == nthr_dist1)fproc_cur_pauc = fp2;
				double dauc = (tproc_cur + tproc_pred) * (fproc_cur_pauc - fproc_pred) / nseq_real[ns] / 2;
				tp_here1[n][n_here1[n]] = (double)tp_one[n][i] / nseq_real[ns];
				fp_here1[n][n_here1[n]] = fproc_cur_pauc;
				n_here1[n]++;
				//fprintf(out,"%d\t%d\t%g\t%g\t%g\n", tproc_cur, tproc_pred, fproc_cur, fproc_pred, dauc);
				//fprintf(outq, "%g\t%f\n", fproc_cur, (double)tproc_cur / n_cnt_tot);
				if (fproc_cur <= fp2)auc_one[n] += dauc;
				if (fproc_cur >= fp2)break;
				tproc_pred = tproc_cur;
				fproc_pred = fproc_cur;
			}
		}
		printf("%s\t%s\t%d\t%s\t%g\n", file_for[ns], file_back[ns], n + 1, partner_db[nm], auc_one[n]);
	}
	if ((out_auc = fopen(file_out, "wt")) == NULL)
	{
		fprintf(out_auc, "Input file %s can't be opened!\n", file_out);
		exit(1);
	}
	fprintf(out_auc, "%s\t%s\t%s\t%s\t", file_for[0], file_back[0], file_for[1], file_back[1]);
	for (n = 0; n < 4; n++)fprintf(out_auc, "\t%f", auc_one[n]);
	fprintf(out_auc, "\t\t%f\n", (auc_one[2] + auc_one[3])/(auc_one[0] + auc_one[1]));
	fclose(out_auc);
	//for (n = 0; n < 4; n++)fclose(out_roc[n]);
	//fclose(outlog);
	//printf("All\t");
	for (k = 0; k < 2; k++)
	{
		delete[] fpr_all[k];
		delete[] thr_all[k];
		delete[] len_real[k];
		delete[] len_back[k];
	}
	for (k = 0; k < 4; k++)
	{
		delete[] tp_one[k];
		delete[] fp_one[k];
	}
	delete[] len_back;
	delete[] len_real;
	delete[] fpr_all;
	delete[] thr_all;
	delete[] tp_one;
	delete[] fp_one;
	for (k = 0; k < 2; k++)
	{
		delete[] tp_here1[k];
	}
	delete[] tp_here1;
	for (k = 0; k < 2; k++)
	{
		delete[] fp_here1[k];
	}
	delete[] fp_here1;
	for (n = 0; n < 2; n++)
	{
		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < nseq_back[n]; i++)
			{
				delete[] seq_back[n][k][i];
			}
			delete[] seq_back[n][k];
		}
		delete[] seq_back[n];
	}
	delete[] seq_back;
	for (n = 0; n < 2; n++)
	{
		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < nseq_real[n]; i++)
			{
				delete[] seq_real[n][k][i];
			}
			delete[] seq_real[n][k];
		}
		delete[] seq_real[n];
	}
	delete[] seq_real;
	for (n = 0; n < 2; n++)matrix[n].mem_out(len_partner[n]);
	return 0;
}