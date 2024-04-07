#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define DIM 100

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
		test = UnderStol(d, 1, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].end = atoi(s);
		test = UnderStol(d, 2, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].buf = atof(s);
		test = UnderStol(d, 3, s, sizeof(s), sep);
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

int main(int argc, char* argv[])
{
	char d[300], file_out_log[300], file_in_sga[300], file_in_tab[300], file_out_distb[300];
	char s[100];
	int k;
	FILE* out_log, * in_tab, * out_distb;
	city sta;

	if (argc != 5)
	{
		printf("%s 1file in_pfm 2file in_pwm 3file in tab 4file out_binary 5file out_log\n", argv[0]);		//6file out_cpp_arr 7file out_cpp_struct_many 8char name
		return -1;
	}
	strcpy(file_in_sga, argv[1]);
	strcpy(file_in_tab, argv[2]);
	strcpy(file_out_distb, argv[3]);
	strcpy(file_out_log, argv[4]);

	if (sta.get_file(file_in_sga) == -1)
	{
		printf("Site %s function not found!", file_in_sga);
		exit(1);
	}
	int len1 = sta.len;
	if ((in_tab = fopen(file_in_tab, "rt")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_in_tab);
		return -1;
	}
	int n_thresh = 0;
	//fgets(d, sizeof(d), in_tab);//header
	char sep = '\t';
	while (fgets(d, sizeof(d), in_tab) != NULL)
	{
		char c = d[0];
		if (c == '-' || isdigit(c))
		{			
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
	fwrite(&sta, sizeof(sta), 1, out_distb);
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
	fprintf(out_log, "%d\t%d\t%.18f\t%.18f\n", sta.len, n_thresh, thr_dist[0], fpr_dist[0]);
	fclose(out_log);
	return 0;
}
