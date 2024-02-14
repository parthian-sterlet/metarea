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
//#define MAXPAR 1500 //max no. of motifs
#define NUM_LIBRARY 7// 4islo bibliotek

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
	int i, j, k, mot1, mot2;
	char partner_db[ARGLEN];// , tfclass[80];path_fasta[ARGLEN], pfile_for[ARGLEN], pfile_back[ARGLEN],
	char file_nmat[ARGLEN], file_pv[ARGLEN], file_bin[ARGLEN];
	FILE* out_pv, * out_mat, * in_pwm, *out_bin;

	if (argc != 5)
	{
		printf("%s 1input binary file 2,3output text files nmat,pvalues 4 output binary file", argv[0]);
		return -1;
	}	
	strcpy(partner_db, argv[1]); //h12hs, h12mm	
	strcpy(file_nmat, argv[2]);
	strcpy(file_pv, argv[3]);
	strcpy(file_bin, argv[4]);
	int s_overlap_min = 6, s_ncycle_small = 1000, s_ncycle_large = 10000;//for permutation(motif_comparison) min_length_of_alignment, no. of permutation (test & detailed)
	double s_granul = 0.001;//for permutation(motif_comparison) okruglenie 4astotnyh matric	
	double p_crit = 0.05;// threshold for similarity of matrices
	double pval_sim[4];

	if ((out_mat = fopen(file_nmat, "wt")) == NULL)
	{
		fprintf(stderr, "Input file %s can't be opened!\n", file_nmat);
		exit(1);
	}
	if ((out_pv = fopen(file_pv, "wt")) == NULL)
	{
		fprintf(stderr, "Input file %s can't be opened!\n", file_pv);
		exit(1);
	}
	if ((out_bin = fopen(file_bin, "wb")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_bin);
		return -1;
	}
	//motif library
	//char  motif_name[MAXPAR][40], motif_class[MAXPAR][80], motif_family[MAXPAR][80], motif_class_inx[MAXPAR][8], motif_family_inx[MAXPAR][10], motif_tf[MAXPAR][30];
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
	}
	if ((in_pwm = fopen(partner_db, "rb")) == NULL)
	{
		printf("Input file %s can't be opened!\n", partner_db);
		return -1;
	}
	int len_partner[2];
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
	long *motifp;
	motifp = new long[n_motifs];
	if (motifp == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	int n_motifs1 = n_motifs - 1;
	int* simn;
	simn = new int[n_motifs1];
	if (simn == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	double* simv;
	simv = new double[n_motifs1];
	if (simv == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	matrices matrix[2];
	int sim_here;
	int sim_tot = 0;

	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		//printf("%d\t%s\n", mot1 + 1, motif_class[mot1]);
		motifp[mot1] = ftell(in_pwm);
		fread(&len_partner[0], sizeof(int), 1, in_pwm);
		int len_partner4 = len_partner[0] * 4;
		fread(pfm[0], sizeof(double), len_partner4, in_pwm);
		fread(pwm[0], sizeof(double), len_partner4, in_pwm);
		int nthr_dist_all;
		fread(&nthr_dist_all, sizeof(int), 1, in_pwm);
		fread(thr_all[0], sizeof(double), nthr_dist_all, in_pwm);
		fread(fpr_all[0], sizeof(double), nthr_dist_all, in_pwm);				
	}	
	for (mot1 = 0; mot1 < n_motifs; mot1++)
	{
		fseek(in_pwm, motifp[mot1], SEEK_SET);
		fread(&len_partner[0], sizeof(int), 1, in_pwm);
		int len_partner40 = len_partner[0] * 4;
		fread(pfm[0], sizeof(double), len_partner40, in_pwm);
		fread(pwm[0], sizeof(double), len_partner40, in_pwm);
		int nthr_dist_all1;
		fread(&nthr_dist_all1, sizeof(int), 1, in_pwm);
		fread(thr_all[0], sizeof(double), nthr_dist_all1, in_pwm);
		fread(fpr_all[0], sizeof(double), nthr_dist_all1, in_pwm);
		matrix[0].mem_in(len_partner[0]);
		for (i = 0; i < len_partner[0]; i++)for (j = 0; j < OLIGNUM; j++)matrix[0].fre[i][j] = pfm[0][i][j];		
		sim_here = 0;
		for (mot2 = mot1 + 1; mot2 < n_motifs; mot2++)
		{
			//fseek(in_pwm, motifp[mot2], SEEK_SET);
			fread(&len_partner[1], sizeof(int), 1, in_pwm);
			int len_partner41 = len_partner[1] * 4;
			fread(pfm[1], sizeof(double), len_partner41, in_pwm);
			fread(pwm[1], sizeof(double), len_partner41, in_pwm);
			int nthr_dist_all2;
			fread(&nthr_dist_all2, sizeof(int), 1, in_pwm);
			fread(thr_all[1], sizeof(double), nthr_dist_all2, in_pwm);
			fread(fpr_all[1], sizeof(double), nthr_dist_all2, in_pwm);
			matrix[1].mem_in(len_partner[1]);
			for (i = 0; i < len_partner[1]; i++)for (j = 0; j < OLIGNUM; j++)matrix[1].fre[i][j] = pfm[1][i][j];		
			for (i = 0; i < 4; i++)pval_sim[i] = 1;
			double pvalue_similarity_tot = pfm_similarity(&matrix[0], &matrix[1], s_granul, s_overlap_min, s_ncycle_small, s_ncycle_large, pval_sim);						
			if (pvalue_similarity_tot < p_crit)
			{
				simn[sim_here] = mot2;
				simv[sim_here] = -log10(pvalue_similarity_tot);
				sim_here++;
			}
			matrix[1].mem_out(len_partner[1]);
		}	
		matrix[0].mem_out(len_partner[0]);
		fwrite(&sim_here, sizeof(int), 1, out_bin);
		if (sim_here > 0)
		{
			sim_tot += sim_here;			
			fwrite(simn, sizeof(double), sim_here, out_bin);
			fwrite(simv, sizeof(double), sim_here, out_bin);			
			fprintf(out_mat, "%d", mot1 + 1);
			for (i = 0; i < sim_here; i++)fprintf(out_mat, "\t%d", simn[i]);
			fprintf(out_pv, "%d", mot1 + 1);
			for (i = 0; i < sim_here; i++)fprintf(out_pv, "\t%f", simv[i]);
		}		
		fprintf(out_pv,"\n");
		fprintf(out_mat, "\n");
		printf("%d\t%d\n", mot1 + 1, sim_tot);
	}
	fclose(out_pv);
	fclose(in_pwm);
	fclose(out_mat);
	fclose(out_bin);
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
	delete[] simv;
	delete[] simn;
	delete[] motifp;
	return 0;
}