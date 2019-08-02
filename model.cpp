// /*************************************************************************
// * 
// * Baylor College of Medicine CONFIDENTIAL
// * __________________
// * 
// * [2016] - [20**] Baylor College of Medicine & Yongtao Guan
// * All Rights Reserved.
// * 
// * NOTICE: All information contained herein is, and remains
// * the property of Baylor College of Medicine and the creator. 
// * The intellectual and technical concepts contained
// * herein are proprietary to Baylor College of Medicine and the creator
// * and may be covered by U.S. and Foreign Patents, patents in process, 
// * and are protected by trade secret or copyright law.
// * Dissemination of this information or reproduction of this material
// * is strictly forbidden unless prior written permission is obtained
// * from Baylor College of Medicine.
//                  
// version 0.4
// This version uses a new partioning method proposed by Dr. Xie. 
// It implicitly imposes a limit on the number of bins.
#include <fstream> 
#include <iostream> 
#include <array>
#include <unordered_map>
// #include <algorithm>
// #include <iterator>
#include <string.h>
#include <zlib.h>
#include <string>
#include <cmath> 
#include <complex>
#include <vector> 
#include <map>
#include <set> 
#include <stdio.h>
#include <sys/vfs.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_permute.h"
//#include "gsl/gsl_sf_gamma.h"
//#include "gsl/gsl_multifit.h"
     
using namespace std; 
                          
#define dd 2 
void operation(double a[][dd], char); 
void operation(double a[][dd], int, double); 
void operation_rep1(double prod[2][2], char);
void operation_rep2(double prod[3][3], char);
void operation_rep3(double prod[4][4], char);
void operation_rep4(double prod[4][4], char);
void operation_rep5(std::complex<long double> prod[2][2], char nuc);
void operation_rep6(double prod[4][4], char);
void multiply2(double A[2][2], double B[2][2], double C[2][2]);
void complexmultiply2(std::complex<long double> A[2][2], std::complex<long double> B[2][2], std::complex<long double> C[2][2]);
void multiply3(double A[3][3], double B[3][3], double C[3][3]);
void multiply4(double A[4][4], double B[4][4], double C[4][4]);
void simreads(void); 
void compute_stats(vector<string> &, string, long, string, unsigned long, int withgc, int ncol);
void compute_stats_unequal_datasize(vector<string> &vfn, string fout, string method, int withgc, int ncol, int n_bin, string bin_file); 
void compute_stats2(vector<string> &, vector<string> &, string, long, string, unsigned long, int ncol); 
void compute_stats2_unequal_datasize(vector<string> &vfn1, vector<string> &vfn2, string fout, string method, int ncol, int n_bin, string bin_file);
void adio(string, double*, int *); 
int compare_doubles(const void * a, const void *b); 
int compare_doubles_2D(const void *pa, const void *pb);
double pairwiseL2(double * xx, double * yy, int nx, int ny); 
double * select_quantiles_old(double * xx, double * yy, int nx, int ny, int ng);
double * select_quantiles(double * xx, long nx, int ng);
double ** select_quantiles_all(double ** xx, long nx, int ng, int ncol);
int quickselect(double * brk, int ng, double num);
double hellinger(double * xx, double * yy, int nx, int ny);
double hellinger_2D(double ** xx, double ** yy, int nx, int ny); 
double hellinger_4D(double ** xx, double ** yy, int nx, int ny);
double hellinger_1D_quantiles_old(double * xx, double * yy, int nx, int ny);
double hellinger_1D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng);
double hellinger_2D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng);
double hellinger_3D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng);
double hellinger_4D_quantiles(unordered_map<string, double> wx, unordered_map<string, double> wy, double pseudo_freq);
double hellinger_5D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk);
double hamming_4D_quantiles(unordered_map<string, double> wx, unordered_map<string, double> wy);
double l1_4D_quantiles(unordered_map<string, double> wx, unordered_map<string, double> wy);
double adtest(double * xx, double * yy, int nx, int ny); 
double cosine(double * xx, double * yy, int nx, int ny);
double cosine_4D(double ** xx, double ** yy, int nx, int ny);
double cosine_1D_quantiles_old(double * xx, double * yy, int nx, int ny);
double cosine_1D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng);
double cosine_2D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng);
double cosine_3D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng);
double cosine_4D_quantiles(unordered_map<string, double> wx, unordered_map<string, double> wy);
double cosine_5D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk);
double canberra_4D(double ** xx, double ** yy, int nx, int ny);
double adtest1(double * xx, double * yy, int nx, int ny); 
double kstest(double * xx, double * yy, int nx, int ny); 
void twoalleles(void); 
void quantify_reads(string fin, string fout, int dim, int w, int step, int nm, int, int, double alpha); 
void regress_gc(string fin, string fout); 
void qc_reads(string fin, string fout, int dim, int nmiss, int tph); 
void rmdup_quad(string fin, string fout); 
double simpson(unordered_map<string, double> wx);
double shannon(unordered_map<string, double> wx);
void compute_alpha_unequal_datasize(vector<string> &vfn, string fout, string method, int withgc, int ncol, int n_bin, string bin_file);

//global gsl variable in use; 
const gsl_rng_type * gslType;
gsl_rng * gsl_r; 

//global log file; 
FILE * flog; 

//global length;
int seqLen; 

double now()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.;
}

int main(int argc, char ** argv)
{
	string fin; 
	string fin2; 
	string fout; 
	string fnlog("\0"); 
	string method("h2_4d"); 
	vector<string> vfn; 
	vector<string> vfn2; 
	int dim = 75; 
	int w = dim; //sliding window size
	int step = w; //sliding window step
	int nm = 1; 
	double alpha=1.0;
	long s = 0;  //sample this many reads; 
	int adaptor_size = 0; 
	int cad = 0; 
	int cad2 = 0; 
	int two = 0; 
	int qtf = 0; 
	int rgc = 0; 
	int within_div = 0; 
	int qc = 0; 
	int rmdup = 0; 
	int tphred = 20; 
	int tmiss = 2; 
	int tlen = 76; 
	int withgc = 0; //toggle for whether to compute ad distance for gc contents.
	int ncol=4; 
	int n_bin=15; 
	string bin_file = "";
	unsigned long seed = time(NULL); 

	for(int i = 1; i < argc; i++) 
	{   
		string str;  
		string opt; 
		
		if(i>1 && argv[i][0] != '-') 
			continue;
		if(i==1 && argv[i][0] == '-') 
		{
			printf("./nubeam cad or ./nubeam qtf \n"); 
			exit(0);
		}
		str.assign(argv[i]);
		opt.assign(str, 0, str.length());

		if(str.compare("cad") == 0) {			
			cad = 1; 
		}
		else if(str.compare("cad2") == 0) {			
			cad2 = 1; 
		}
		else if(str.compare("two") == 0) {			
			two = 1; 
		}
		else if(str.compare("quad") == 0 || str.compare("qtf") == 0) {			
			qtf = 1; 
		}
		else if(str.compare("nogc") == 0 || str.compare("rgc") == 0) {			
			rgc = 1; 
		}
		else if(str.compare("within_div") == 0) {			
			within_div = 1; 
		}
		else if(str.compare("qc") == 0 || str.compare("rgc") == 0) {			
			qc = 1; 
		}
		else if(str.compare("rmdup") == 0 ) {			
			rmdup = 1; 
		}
		else if(str.compare("-h") == 0) {			
			if(qtf == 1)
			{
				printf("./nubeam quad [-idonwh]\n"); 
				printf("compute quadriples for reads in (gzipped) fastq format.\n");
				printf("produces prefix.quad.gz (gc content is within) and prefix.quad.log.\n");  
				printf("-a : adaptor size (default 8)\n"); 
				printf("-i : input filename\n"); 
				printf("-o : output prefix\n"); 
				printf("-d : dimension, or length of the reads (default d=75).\n"); 
				printf("-n : number of missing nucleotide allowed.\n"); 
				printf("-w : sliding window size (default w=d).\n"); 
				printf("-S : sliding window step (default S=w).\n"); 
				printf("-f : value, plus 33 is the PHRED quality value of fastq reads.\n");
				printf("-A : alpha value (0,1].\n"); 
				printf("-h : print this help\n"); 
			}
			else if(cad == 1)
			{
				printf("./nubeam cad [-iomsh]\n"); 
				printf("compute pariwise distances of a set; the inputs are nubeam qtf outputs.\n"); 
				printf("produces prefix.cad.log.\n");  
				printf("-i : specifies input file which is output of nubeam.\n"); 
				printf("-o : output prefix (prefix.log contains pairwise distance matrix)\n"); 
				printf("-m : choice of methods: h2, h2_1d, h2_2d, h2_3d, h2_4d, h2_5d, ham_4d, l1_4d, l1_1d, cos, cos_1d, cos_2d, cos_3d, cos_4d, cos_5d, can_4d, ad, ks, ad1.\n"); 
				printf("-s : sampling s many entries to perform test\n");
				printf("-c : designating the number of columns of scores\n"); 
				printf("-b : designating the number of bins per column of scores\n");
//				printf("-a : allow negative entres or not.\n"); 
				printf("-h : print this help\n"); 
			}
			else if(cad2 == 1)
			{
				printf("./nubeam cad2 [-ijomsh]\n"); 
				printf("compute pariwise cross distances between two sets; the inputs are nubeam qtf outputs.\n"); 
				printf("produces prefix.cad2.log.\n");  
				printf("-i : specifies input file of first set, which is output of nubeam.\n"); 
				printf("-j : specifies input file of second set, which is output of nubeam.\n"); 
				printf("-o : output prefix (prefix.log contains pairwise distance matrix)\n"); 
				printf("-m : choice of methods: h2, ad, ks, ad1.\n"); 
				printf("-s : sampling s many entries to perform test\n");
				printf("-c : designating the number of columns of scores\n"); 
				printf("-b : designating the number of bins per column of scores\n");
				printf("-bf : the file describing how to partition the bins\n");
//				printf("-a : allow negative entres or not.\n"); 
				printf("-h : print this help\n"); 
			}
			else if(rgc == 1) 
			{
                printf("./nubeam nogc [-ioh]\n"); 
				printf("regress out gc contents from read quantification and output residuals.\n");  
				printf("produces prefix.nogc.gz and prefix.nogc.log.\n");  
				printf("-i : input file name.\n"); 
				printf("-o : output prefix.\n"); 
				printf("-h : print this help\n"); 
			}
			else if(within_div == 1) 
			{
                printf("./nubeam within_div [-iohN]\n"); 
				printf("calculate within-sample diversity.\n");  
				printf("produces prefix.alpha.gz\n");  
				printf("-i : input file name.\n"); 
				printf("-o : output prefix.\n"); 
				printf("-N : total number of reads.\n");
				printf("-m : choice of methods: simpson_4d, shannon_4d.\n"); 
				printf("-c : designating the number of columns of scores\n"); 
				printf("-b : designating the number of bins per column of scores\n");
				printf("-bf : the file describing how to partition the bins\n");
				printf("-h : print this help\n"); 
			}
			else if(qc == 1) 
			{
                printf("./nubeam qc [-ioh]\n"); 
				printf("regress out gc contents from read quantification and output residuals.\n");  
				printf("produces prefix.qc.gz and prefix.qc.log.\n");  
				printf("-i : input file name.\n"); 
				printf("-o : output prefix.\n"); 
				printf("-d : threshold of read length (default 76).\n"); 
				printf("-n : threshold of missing counts (default 2).\n"); 
				printf("-w : threshold of phred score (default 25).\n"); 
				printf("-h : print this help\n"); 
			}
//			else if(rmdup == 1) 
//			{
//                printf("./nubeam rmdup [-ioh]\n"); 
//				printf("remove duplicates from quad.gz file.\n");  
//				printf("produces prefix.rmdup.gz and prefix.rmdup.log.\n");  
//				printf("-i : input file name.\n"); 
//				printf("-o : output prefix.\n"); 
//				printf("-h : print this help\n"); 
//			}
//			else if(two == 1) 
//			{
//                printf("./nubeam two [-dsh]\n"); 
//				printf("simulate two similar haplotypes and reads and compuate ad distance.\n");  
//				printf("-s : sequence length.\n"); 
//				printf("-w : read length.\n"); 
//				printf("-d : sequence depth.\n"); 
//				printf("-n : mutations before two sequences.\n"); 
//				printf("-h : print this help\n"); 
//			}
			else {
                printf("./nubeam [qc,quad,nogc,cad,cad2]\n"); 
				printf("For example, use ./nubeam quad -h for more options.\n"); 
			}

			exit(0); 
		}
		else if (str.compare("-a") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                adaptor_size = atoi(argv[i+1]);
		}	
		else if (str.compare("-A") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                alpha = atof(argv[i+1]);
		}
		else if (str.compare("-in") == 0 || str.compare("-i") == 0) {
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fin.clear();
				fin.assign(argv[i+1]);
				vfn.push_back(fin); 
		}
		else if (str.compare("-j") == 0 ) {
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fin2.clear();
				fin2.assign(argv[i+1]);
				vfn2.push_back(fin2); 
		}
		else if (str.compare("-out") == 0 || str.compare("-o") == 0) {
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fout.clear();
				fout.assign(argv[i+1]);
				fnlog.assign(fout); 
		}
		else if (str.compare("-d") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                dim = w = step = atoi(argv[i+1]);
                tlen = atoi(argv[i+1]);
		}	
		else if (str.compare("-method") == 0 || str.compare("-m") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				method.clear(); 
                method.assign(argv[i+1]);
				if(method.compare("ad") != 0 && method.compare("ad1") != 0 && method.compare("ks") != 0 && method.compare("h2") !=0 && method.compare("h2_1d") !=0 && method.compare("h2_2d") !=0 && method.compare("h2_3d") !=0 && method.compare("h2_4d") !=0 && method.compare("h2_5d") !=0 && method.compare("ham_4d") !=0 && method.compare("l1_4d") !=0 && method.compare("l1_1d") !=0 && method.compare("cos") !=0 && method.compare("cos_1d") !=0 && method.compare("cos_2d") !=0 && method.compare("cos_3d") !=0 && method.compare("cos_4d") !=0 && method.compare("cos_5d") !=0 && method.compare("can_4d") !=0 && method.compare("simpson_4d") !=0 && method.compare("shannon_4d") !=0)
				{
                    printf("wrong augument after option.\n");
					exit(0); 
				}
		}	
		else if (str.compare("-nm") == 0 || str.compare("-n") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                nm = atoi(argv[i+1]);
                tmiss = atoi(argv[i+1]);
		}	
		else if (str.compare("-w") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                w = step = atoi(argv[i+1]);
		}
		else if (str.compare("-S") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                step = atoi(argv[i+1]);
		}	
		else if (str.compare("-f") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                tphred = atoi(argv[i+1]);
		}	
		else if (str.compare("-s") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                s = atol(argv[i+1]);
		}	
		else if (str.compare("-R") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                seed = atol(argv[i+1]);
		}	
		else if (str.compare("-l") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                tlen = atoi(argv[i+1]);
		}	
		else if (str.compare("-withgc") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                withgc = atoi(argv[i+1]);
		}	
		else if (str.compare("-c") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                ncol = atoi(argv[i+1]);
                }
        else if (str.compare("-b") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                n_bin = atoi(argv[i+1]);
                }
        else if (str.compare("-bf") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue; 
                bin_file.assign(argv[i+1]);
		}
		else 
		{
				fprintf(stderr,"Bad option %s\n", argv[i]);
				exit(0);
		}              
	}                          
	
	gsl_rng_env_setup();
	gslType = gsl_rng_default;
	gsl_r = gsl_rng_alloc (gslType);
	gsl_rng_set(gsl_r, seed);
	//set up random number generator; 

	if(fnlog.size()==0)
	{
		char temp[1000]; 
		sprintf(temp, "%ld", seed); 
		fnlog.assign(temp); 
	}

	seqLen = dim; 
	
	if(two) {
		twoalleles(); 
		return 1; 
	}
	if(qtf) {
		fnlog.append(".quad.log");
		printf("Log file name: %s \n", fnlog.c_str()); 
	   flog = fopen(fnlog.c_str(), "w");
	   if(flog == NULL) 
		   printf("can't open log file %s to write, proceed anyway.", fnlog.c_str());
       quantify_reads(fin, fout, dim, w, step, nm, adaptor_size, tphred, alpha); 
	   fclose(flog); 
	   return 1; 
	}
	if(rgc) {
		fnlog.append(".nogc.log");
		printf("Log file name: %s \n", fnlog.c_str()); 
	   flog = fopen(fnlog.c_str(), "w");
	   if(flog == NULL) 
		   printf("can't open log file %s to write, proceed anyway.", fnlog.c_str());
		regress_gc(fin, fout); 
		fclose(flog); 
		return 1; 
	}
	// if(within_div) {
	// 	alpha_diversity(fin, fout, N);  
	// 	return 1; 
	// }
	if(within_div) {
		fnlog.append(".alpha.log");
		printf("Log file name: %s \n", fnlog.c_str()); 
		flog = fopen(fnlog.c_str(), "w");
		if(flog == NULL) 
			printf("can't open log file %s to write, proceed anyway.", fnlog.c_str());
		compute_alpha_unequal_datasize(vfn, fout, method, withgc, ncol, n_bin, bin_file); 
		fclose(flog); 
		return 1; 
	}
	if(qc) {
		fnlog.append(".qc.log");
		printf("Log file name: %s \n", fnlog.c_str()); 
	   flog = fopen(fnlog.c_str(), "w");
	   if(flog == NULL) 
		   printf("can't open log file %s to write, proceed anyway.", fnlog.c_str());
		qc_reads(fin, fout, tlen, tmiss, tphred); 
		fclose(flog); 
		return 1; 
	}
	if(cad) {        
		fnlog.append(".cad.log");
		printf("Log file name: %s \n", fnlog.c_str()); 
	   flog = fopen(fnlog.c_str(), "w");
	   if(flog == NULL) 
		   printf("can't open log file %s to write, proceed anyway.", fnlog.c_str());
		//compute_stats(vfn, fout, s, method, seed, withgc, ncol);
		compute_stats_unequal_datasize(vfn, fout, method, withgc, ncol, n_bin, bin_file); 
		fclose(flog); 
		return 1;  
	}
	if(cad2) {        
		fnlog.append(".cad2.log");
		printf("Log file name: %s \n", fnlog.c_str()); 
	   flog = fopen(fnlog.c_str(), "w");
	   if(flog == NULL) 
		   printf("can't open log file %s to write, proceed anyway.", fnlog.c_str());
		//compute_stats2(vfn, vfn2, fout, s, method, seed, ncol); 
		compute_stats2_unequal_datasize(vfn, vfn2, fout, method, ncol, n_bin, bin_file); 
		fclose(flog); 
		return 1;  
	}
	if(rmdup) {        
		fnlog.append(".rmdup.log");
		printf("Log file name: %s \n", fnlog.c_str()); 
	   flog = fopen(fnlog.c_str(), "w");
	   if(flog == NULL) 
		   printf("can't open log file %s to write, proceed anyway.", fnlog.c_str());
		rmdup_quad(fin, fout); 
		fclose(flog); 
		return 1;  
	}
	return 0; 
}

void qc_reads(string fin, string fout, int dim, int nmiss, int tph) 
{
	fprintf(flog, "### qc -d %d -n %d -w %d -o %s -i %s \n",  dim, nmiss, tph, fout.c_str(), fin.c_str());  
	fout.append(".qc.gz"); 
	gzFile fp1 = gzopen(fout.c_str(), "wb"); 
	if(fp1 == NULL) 
	{
		printf("can't open %s file to write\n", fout.c_str()); 
		exit(0); 
    }
	fprintf(flog, "### -i %s -o %s -d %d -n %d -w %d \n", fin.c_str(), fout.c_str(), dim, nmiss, tph);  

	gzFile fgz = gzopen(fin.c_str(), "rb"); 
	if(fgz == NULL) {
		printf("can't open %s file to read. \n", fin.c_str()); 
		exit(0); 
	}
    
    char * tline = new char[1024]; 
    char * tscore = new char[1024]; 
    char * first = new char[1024]; 
    char * third = new char[1024]; 

	long c1, c2, c3; 
	c1=c2=c3=0; 
	long ca, cc, ct, cg, cn, cphred; 
	ca=cc=ct=cg=cn=cphred=0; 
	long total=0;  
	gzgets(fgz, first, 1024); 
	gzgets(fgz, tline, 1024); 
	gzgets(fgz, third, 1024); 
	gzgets(fgz, tscore, 1024); 
	while(1) 
	{
		total += 1; 
		char * line = tline; 
		char * score = tscore; 

		line[strcspn(line, "\r\n")] = 0; 
        int skip = 0; 
		if(strlen(line) < dim) 
		{
			skip = 1; 
			c1++; 
		}
		int nn = 0; 
		for(int k = 0; k < strlen(line); k++)
		{
			char ch = toupper(line[k]);
			if( ch == 'A') ca ++; 
			else if( ch == 'T') ct ++; 
			else if( ch == 'C') cc ++; 
			else if( ch == 'G') cg ++; 
			else { cn++; nn++;}
		}
		if(nn > nmiss ) 
		{
			skip = 2; 
			c2++; 
		}

		int lowscore = 0; 
		for(int k = 0; k < dim; k++)
		{
			int phred = (int) score[k] - 33;  
			if(phred < tph) 
				lowscore ++; 
		}
		cphred += lowscore; 
		if(lowscore > 1) 
		{
			skip = 3;
			c3++; 
		}

        if(skip == 0) 
		{
			gzprintf(fp1, "%s", first); 
			gzprintf(fp1, "%s\n", line); 
			gzprintf(fp1, "%s", third); 
			gzprintf(fp1, "%s", score); 
		}

		char * ret = Z_NULL; 
		ret = gzgets(fgz, first, 1024); 
		if(ret == Z_NULL) break; 
		ret = gzgets(fgz, tline, 1024); 
		if(ret == Z_NULL) break; 
		ret = gzgets(fgz, third, 1024); 
		if(ret == Z_NULL) break; 
		ret = gzgets(fgz, tscore, 1024); 
		if(ret == Z_NULL) break; 
	}
	gzclose(fp1); 
	gzclose(fgz); 
	delete[] first; 
	delete[] third; 
	delete[] tline; 
	delete[] tscore; 
	fprintf(flog, "total_reads %ld \n", total); 
	fprintf(flog, "fail_length_test %ld \n", c1); 
	fprintf(flog, "fail_missing_test %ld \n", c2); 
	fprintf(flog, "fail_phred_test %ld \n", c3); 
	fprintf(flog, "%ld many bp has phred schore less than %d \n", cphred, tph); 
	fprintf(flog, "A %ld \n", ca); 
	fprintf(flog, "T %ld \n", ct); 
	fprintf(flog, "C %ld \n", cc); 
	fprintf(flog, "G %ld \n", cg); 
	fprintf(flog, "N %ld \n", cn); 
}

void quantify_reads(string fin, string fout, int dim, int w, int step, int nm, int adaptor, int tphred, double alpha) 
{
	fprintf(flog, "### quad -d %d -w %d -s %d -n %d -a %d -A %f -o %s -i %s \n",  dim, w, step, nm, adaptor, alpha, fout.c_str(), fin.c_str());  
	map<int, char> i2c; 
	i2c[0] = 'A'; 
	i2c[1] = 'T'; 
	i2c[2] = 'C'; 
	i2c[3] = 'G'; 
	map<char, int> c2i; 
	c2i['A'] = 0; 
	c2i['T'] = 1; 
	c2i['C'] = 2; 
	c2i['G'] = 3; 
	c2i['N'] = -1; 
	map<char, char> c2c; 
	c2c['A'] = 'T'; 
	c2c['T'] = 'A'; 
	c2c['C'] = 'G'; 
	c2c['G'] = 'C'; 
	c2c['N'] = 'N'; 
	map<char, char> c2cp; 
	c2cp['A'] = 'C'; 
	c2cp['T'] = 'G'; 
	c2cp['C'] = 'A'; 
	c2cp['G'] = 'T'; 
	c2cp['N'] = 'N'; 
	// string spacer1 = ""; // nothing
	// string spacer1 = "GATC"; //4
	// string spacer1 = "GGAATTCC"; //8
	// string spacer1 = "GGGAAATTTCCC"; //12
	// string spacer1 = "GATCGATCGATC";
	// string spacer1 = "GGGGAAAATTTTCCCC"; //16
	// string spacer1 = "GGGGGAAAAATTTTTCCCCC"; //20
	// string spacer1 = "GGGGGGAAAAAATTTTTTCCCCCC"; //24
	// string spacer1 = "GGGGGGGGGGAAAAAAAAAATTTTTTTTTTCCCCCCCCCC"; //40
	// string spacer1 = "GGGGGGGGGGGGGGGAAAAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCCCCCCC"; //60
	// string spacer1 = "GGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCC"; //100
	// string spacer1 = "GATC"; 
	// string spacer2 = "CTAG"; //4
//	string spacer2; 
//	for (unsigned k = 0; k < spacer1.size(); k++)
//		spacer2[k] = c2c[spacer1[spacer1.size()-k]]; 
//	//spacer2 is the reverse compliment (self-defined) of spacer1; 

	string f1(fout); 
	f1.append(".quad.gz"); 
	gzFile fp1 = gzopen(f1.c_str(), "wb"); 
	if(fp1 == NULL) 
	{
		printf("can't open %s file to write\n", fout.c_str()); 
		exit(0); 
    }

//	string f2(fout); 
//	f2.append(".atcg.gz"); 
//	gzFile fp2 = gzopen(f2.c_str(), "wb"); 
//	if(fp2 == NULL) 
//	{
//		printf("can't open %s file to write\n", fout.c_str()); 
//		exit(0); 
//    }

	gzFile fgz = NULL; 
	fgz = gzopen(fin.c_str(), "rb"); 
	if(fgz == NULL) {
		printf("can't open %s file to read. \n", fin.c_str()); 
		exit(0); 
	}
    
	char * first = new char[1024];  
    char * tline = new char[1024]; 
    char * tscore = new char[1024]; 
	gzgets(fgz, first, 1024); 
	gzgets(fgz, tline, 1024); 
	gzgets(fgz, tscore, 1024); 
	gzgets(fgz, tscore, 1024); 

	int c1, c2, c3; 
	c1=c2=c3=0; 
	while (tline != NULL) 
	{
		char * line = tline + adaptor; 
		char * score = tscore + adaptor; 

        int skip = 0; 
		if(strlen(line) < dim+adaptor) 
			skip = 1; 
		if(skip == 0) {
            int temp = 0; 
			for(int k = 0; k < dim; k++)
				temp += (int) (toupper(line[k]) == 'N');
			if(temp > nm) 
				skip = 2; 
		}
		if(skip == 0) {
			int lowscore = 0; 
			for(int k = 0; k < dim; k++)
			{
            	int phred = (int) score[k] - 33;  
				if(phred < tphred) 
					lowscore ++; 
				if(lowscore > 1) break; 
			}
			if(lowscore > 1) 
				skip = 3; 
		}
        if(skip == 0) {
			// for(int s = 0; s <= (dim-w); s+=w)
			for(int s = 0; s <= (dim-w); s+=step)
			{
				string zz(""); 
				zz.append(line, s, w); 
				// write gc contents.
				int nc[4]; 
				int na = 0; 
				memset(nc, 0, 4*sizeof(int)); 
				// start the for loop
				for (unsigned k = 0; k < w; k++) {
					// all capital letters
					zz[k] = toupper(zz[k]);
					// get GC contents
					int t1=c2i[zz[k]];
					if(t1>=0) nc[t1]++;
					else na++;
				}
				// // reverse complement of zz
				// string zz_rc(w, 'N');
				// for (unsigned k = 0; k < w; k++) {
				// 	// get reverse complement
				// 	zz_rc[k] = c2c[zz[w - 1 - k]];
				// }
				// // use the alphabetically smaller read
				// if (zz.compare(zz_rc) > 0) zz = zz_rc;
				
				// string zz_complement(w, 'N');
				// for (unsigned k = 0; k < w; k++)
				// 	zz_complement[k] = c2c[zz[k]]; 
				// zz.append(zz_complement);
				// zz.append(spacer1); 
				// zz.insert(0, spacer1);
//				string t2(zz.size(), 'N'); 
//				for(int k = 0; k < w; k++)
//					t2[k]=c2c[toupper(zz[zz.size()-1-k])];  
//				zz.append(t2); 
//				string t3(zz.size(), 'N'); 
//				for(unsigned k = 0; k < zz.size(); k++)
//					t3[k]=c2cp[zz[zz.size()-1-k]];  
//				zz.append(t3);

 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 //						operation(prod, (int) (zz.at(k) == 'A' || zz.at(k) == 'C')); 
 						operation(prod, (int) (zz.at(k) == 'A'), alpha); 
 					}
 					// gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]));
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(7) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(10) * prod[1][1]));
 					// gzprintf(fp1, "%lf %lf %lf ",  log(0.5 * (prod[0][1] + prod[1][0]) + 0.5 * prod[0][0] + 0.5 * prod[1][1]), 
 					// 	log((sqrt(3) / 4) * (prod[0][1] + prod[1][0]) + 0.25 * prod[0][0] + 0.75 * prod[1][1]),
 					// 	log((sqrt(3) / 4) * (prod[0][1] + prod[1][0]) + 0.75 * prod[0][0] + 0.25 * prod[1][1]));
 					// gzprintf(fp1, "%lf ",  log(0.5 * (prod[0][1] + prod[1][0]) + 0.5 * prod[0][0] + 0.5 * prod[1][1]));
 				}
 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 //						operation(prod, (int) (zz.at(k) == 'A' || zz.at(k) == 'T')); 
 						operation(prod, (int) (zz.at(k) == 'T'), alpha); 
 					}
 					// gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]));
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(7) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(10) * prod[1][1]));
 					// gzprintf(fp1, "%lf %lf %lf ",  log(0.5 * (prod[0][1] + prod[1][0]) + 0.5 * prod[0][0] + 0.5 * prod[1][1]), 
 					// 	log((sqrt(3) / 4) * (prod[0][1] + prod[1][0]) + 0.25 * prod[0][0] + 0.75 * prod[1][1]),
 					// 	log((sqrt(3) / 4) * (prod[0][1] + prod[1][0]) + 0.75 * prod[0][0] + 0.25 * prod[1][1]));
 					// gzprintf(fp1, "%lf ",  log(0.5 * (prod[0][1] + prod[1][0]) + 0.5 * prod[0][0] + 0.5 * prod[1][1]));
 				}
 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 //						operation(prod, (int) (zz.at(k) == 'A' || zz.at(k) == 'C')); 
 						operation(prod, (int) (zz.at(k) == 'C'), alpha); 
 					}
 					// gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]));
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(7) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(10) * prod[1][1]));
 					// gzprintf(fp1, "%lf %lf %lf ",  log(0.5 * (prod[0][1] + prod[1][0]) + 0.5 * prod[0][0] + 0.5 * prod[1][1]), 
 					// 	log((sqrt(3) / 4) * (prod[0][1] + prod[1][0]) + 0.25 * prod[0][0] + 0.75 * prod[1][1]),
 					// 	log((sqrt(3) / 4) * (prod[0][1] + prod[1][0]) + 0.75 * prod[0][0] + 0.25 * prod[1][1]));
 					// gzprintf(fp1, "%lf ",  log(0.5 * (prod[0][1] + prod[1][0]) + 0.5 * prod[0][0] + 0.5 * prod[1][1]));
 				}
 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 //						operation(prod, (int) (zz.at(k) == 'A' || zz.at(k) == 'G')); 
 						operation(prod, (int) (zz.at(k) == 'G'), alpha); 
 					}
 					// gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]));
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(7) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(10) * prod[1][1]));
 					// gzprintf(fp1, "%lf %lf %lf ",  log(0.5 * (prod[0][1] + prod[1][0]) + 0.5 * prod[0][0] + 0.5 * prod[1][1]), 
 					// 	log((sqrt(3) / 4) * (prod[0][1] + prod[1][0]) + 0.25 * prod[0][0] + 0.75 * prod[1][1]),
 					// 	log((sqrt(3) / 4) * (prod[0][1] + prod[1][0]) + 0.75 * prod[0][0] + 0.25 * prod[1][1]));
 					// gzprintf(fp1, "%lf ",  log(0.5 * (prod[0][1] + prod[1][0]) + 0.5 * prod[0][0] + 0.5 * prod[1][1]));
 				}
 				// gzprintf(fp1, "%d %d\n", nc[0]+nc[1], nc[2]+nc[3]);
 				// {
 				// 	double prod[dd][dd] = {{1,0},{0,1}}; 
 				// 	for(unsigned k = 0; k < zz.size(); k++) {
 				// 		operation(prod, (int) (zz_rc.at(k) == 'A'), alpha); 
 				// 	}
 				// 	gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(7) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(10) * prod[1][1]));
 				// }
 				// {
 				// 	double prod[dd][dd] = {{1,0},{0,1}}; 
 				// 	for(unsigned k = 0; k < zz.size(); k++) {
 				// 		operation(prod, (int) (zz_rc.at(k) == 'T'), alpha); 
 				// 	}
 				// 	gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(7) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(10) * prod[1][1]));
 				// }
 				// {
 				// 	double prod[dd][dd] = {{1,0},{0,1}}; 
 				// 	for(unsigned k = 0; k < zz.size(); k++) {
 				// 		operation(prod, (int) (zz_rc.at(k) == 'C'), alpha); 
 				// 	}
 				// 	gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(7) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(10) * prod[1][1]));
 				// }
 				// {
 				// 	double prod[dd][dd] = {{1,0},{0,1}}; 
 				// 	for(unsigned k = 0; k < zz.size(); k++) { 
 				// 		operation(prod, (int) (zz_rc.at(k) == 'G'), alpha); 
 				// 	}
 				// 	gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(7) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(10) * prod[1][1]));
 				// }
 				// {
 				// 	double prod[dd][dd] = {{1,0},{0,1}}; 
 				// 	for(unsigned k = 0; k < zz.size(); k++)
 				// 		operation(prod, (int) (zz.at(k) == 'A'), alpha); 
 				// 	gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]));
 				// }
 				// {
 				// 	double prod[dd][dd] = {{1,0},{0,1}}; 
 				// 	for(unsigned k = 0; k < zz.size(); k++)
 				// 		operation(prod, (int) (zz.at(k) == 'C'), alpha); 
 				// 	gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]));
 				// 	gzprintf(fp1, "%d %d\n", nc[0]+nc[1], nc[2]+nc[3]);
 				// }
 				// {
 				// 	double prod[dd][dd] = {{1,0},{0,1}}; 
 				// 	for(unsigned k = 0; k < zz.size(); k++)
 				// 		operation(prod, (int) (zz.at(k) == 'T'), alpha); 
 				// 	gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]));
 				// }
 				// {
 				// 	double prod[dd][dd] = {{1,0},{0,1}}; 
 				// 	for(unsigned k = 0; k < zz.size(); k++)
 				// 		operation(prod, (int) (zz.at(k) == 'G'), alpha); 
 				// 	gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]));
 				// 	gzprintf(fp1, "%d %d\n", nc[0]+nc[1], nc[2]+nc[3]);
 				// }
				// {
				// 	double prod[2][2] = {{1, 0},{0, 1}};
				// 	for(unsigned k = 0; k < zz.size(); k++)
				// 		operation_rep1(prod, zz.at(k)); 
				// 	gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]));
				// }
				// {
				// 	double prod[3][3] = {{1, 0, 0},{0, 1, 0},{0, 0, 1}};
				// 	for(unsigned k = 0; k < zz.size(); k++)
				// 		operation_rep2(prod, zz.at(k)); 
				// 	gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]+prod[2][2]));
				// }
				// {
				// 	double prod[4][4] = {{1, 0, 0, 0},{0, 1, 0, 0},{0, 0, 1, 0},{0, 0, 0, 1}};
				// 	for(unsigned k = 0; k < zz.size(); k++)
				// 		operation_rep3(prod, zz.at(k)); 
				// 	gzprintf(fp1, "%lf ",  log(prod[0][0]+prod[1][1]+prod[2][2]+prod[3][3]));
				// }
				// {
				// 	double prod[4][4] = {{1, 0, 0, 0},{0, 1, 0, 0},{0, 0, 1, 0},{0, 0, 0, 1}};
				// 	for(unsigned k = 0; k < zz.size(); k++)
				// 		operation_rep4(prod, zz.at(k)); 
				// 	gzprintf(fp1, "%lf ",  (prod[0][0]+prod[1][1]+prod[2][2]+prod[3][3]));
				// }
				// {
				// 	std::complex<long double> prod[2][2] = {
    //     				{std::complex<long double> (1.0), std::complex<long double> (0.0)},
    //     				{std::complex<long double> (0.0), std::complex<long double> (1.0)}
    // 				};
				// 	for(unsigned k = 0; k < zz.size(); k++)
				// 		operation_rep5(prod, zz.at(k)); 
				// 	// gzprintf(fp1, "%Lf ",  (std::real(prod[0][0]+prod[1][1])));
				// 	gzprintf(fp1, "%Lf\n",  (std::real(prod[0][0]+prod[1][1])));
				// }
				// {
				// 	double prod[4][4] = {{1, 0, 0, 0},{0, 1, 0, 0},{0, 0, 1, 0},{0, 0, 0, 1}};
				// 	for(unsigned k = 0; k < zz.size(); k++)
				// 		operation_rep6(prod, zz.at(k)); 
				// 	gzprintf(fp1, "%lf ",  (prod[0][0]+prod[1][1]+prod[2][2]+prod[3][3]));
				// }

				// write gc contents. 
				// int nc[4]; 
				// int na = 0; 
				// memset(nc, 0, 4*sizeof(int)); 
				// for(unsigned k = 0; k < w; k++)
				// {
				// 	int t1=c2i[toupper(line[s+k])]; 
				// 	if(t1>=0) nc[t1]++;
				// 	else na++; 
				// }
				gzprintf(fp1, "%d %d\n", nc[0]+nc[1], nc[2]+nc[3]);
			}
		}
		else if (skip == 1) 
			c1++; 
		else if (skip == 2) 
			c2++; 
		else if (skip == 3) 
			c3++; 
		char * ret = Z_NULL; 
		ret = gzgets(fgz, first, 1024); 
		if(ret == Z_NULL) break; 
		ret = gzgets(fgz, tline, 1024); 
		if(ret == Z_NULL) break; 
		ret = gzgets(fgz, tscore, 1024); 
		if(ret == Z_NULL) break; 
		ret = gzgets(fgz, tscore, 1024); 
		if(ret == Z_NULL) break; 
	}
	fprintf(flog, " %d was skipped due to insufficient length \n", c1); 
	fprintf(flog, " %d was skipped due to missing rate \n", c2); 
	fprintf(flog, " %d was skipped due to low phred \n", c3); 
	gzclose(fp1); 
//	gzclose(fp2); 
	gzclose(fgz); 
	delete[] tline; 
	delete[] first; 
	delete[] tscore; 
}

void operation(double prod[][dd], int yes, double aa) 
{
//    static int ii = 1;
//	ii = ii % 7 + 1; 
//	double aa = ii * 0.01 + 0.1; 
//    double aa = 0.21; 
	if(yes == 1) {
		prod[0][1] += aa*prod[0][0]; 
		prod[1][1] += aa*prod[1][0]; 
	}
	else {
		prod[0][0] += aa*prod[0][1]; 
		prod[1][0] += aa*prod[1][1]; 
	}
}

void operation(double prod[][dd], char nuc) {
//	double temp1, temp2; 
//	capital A append on right; 
//	small a remove from left; 

	double mI[dd][dd] ={{1,0},{0,1}}; 
	if(nuc == 'A') {mI[0][1]=1;mI[1][0]=2;mI[1][1]=3;}
	if(nuc == 'T') {mI[0][1]=2;mI[1][0]=1;mI[1][1]=3;}
	if(nuc == 'C') { mI[0][0]=3; mI[0][1]=2; mI[1][0]=mI[1][1]=1;}
	if(nuc == 'G') { mI[0][0]=3; mI[0][1]=1; mI[1][0]=2; mI[1][1]=1; }

	double res[dd][dd]; 
	for (int i = 0; i < dd; i++) 
		for (int j = 0; j < dd; j++)
		{
			res[i][j] = 0; 
			for (int k = 0; k < dd; k++)
				res[i][j] += prod[i][k] * mI[k][j]; 
		}
	for (int i = 0; i < dd; i++) 
		for (int j = 0; j < dd; j++)
			prod[i][j] = res[i][j]; 
}

// This function multiplies A[][] and B[][], and stores the result in C[][]
void multiply2(double A[2][2], double B[2][2], double C[2][2])
{
    int i, j, k;
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            C[i][j] = 0;
            for (k = 0; k < 2; k++)
                C[i][j] += A[i][k]*B[k][j];
        }
    }
}

// This function multiplies A[][] and B[][], and stores the result in C[][]
void complexmultiply2(std::complex<long double> A[2][2], std::complex<long double> B[2][2], std::complex<long double> C[2][2])
{
    int i, j, k;
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            C[i][j] = 0;
            for (k = 0; k < 2; k++)
                C[i][j] += A[i][k]*B[k][j];
        }
    }
}

// This function multiplies A[][] and B[][], and stores the result in C[][]
void multiply3(double A[3][3], double B[3][3], double C[3][3])
{
    int i, j, k;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            C[i][j] = 0;
            for (k = 0; k < 3; k++)
                C[i][j] += A[i][k]*B[k][j];
        }
    }
}

// This function multiplies A[][] and B[][], and stores the result in C[][]
void multiply4(double A[4][4], double B[4][4], double C[4][4])
{
    int i, j, k;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            C[i][j] = 0;
            for (k = 0; k < 4; k++)
                C[i][j] += A[i][k]*B[k][j];
        }
    }
}

void operation_rep1(double prod[2][2], char nuc) {
	double mI[2][2] = {
    	{1, 0},
    	{0, 1}
    };

	if(nuc == 'A') {
		mI[0][1]=2;
	}
	if(nuc == 'T') {
		mI[1][0]=2;
	}
	if(nuc == 'C') { 
		// mI[0][0]=2; mI[0][1]=1; mI[1][0]=1;
		mI[0][0]=2; mI[1][0]=1;
	}
	if(nuc == 'G') {
		// mI[0][1]=1; mI[1][0]=1; mI[1][1]=2;
		mI[0][0]=2; mI[0][1]=1;
	}

	double res[2][2]; // store result
	multiply2(prod, mI, res);


	for (int i = 0; i < 2; i++) 
		for (int j = 0; j < 2; j++)
			prod[i][j] = res[i][j]; 
}

void operation_rep2(double prod[3][3], char nuc) {
	double mI[3][3] = {
    	{1, 0, 0},
    	{0, 1, 0},
    	{0, 0, 1}
    };
	if(nuc == 'A') {
		mI[0][0]=0.16; mI[1][0]=1.96; mI[2][0]=2.68; mI[2][1]=0.52; mI[2][2]=2.56;
	}
	if(nuc == 'T') {
		mI[0][0]=0.16; mI[0][1]=2.52; mI[0][2]=2.68; mI[1][2]=1.56; mI[2][2]=2.56;
	}
	if(nuc == 'C') {
		mI[0][0]=0.16; mI[0][1]=0.72; mI[0][2]=0.32; 
		mI[1][0]=0.56; mI[1][1]=3.52; mI[1][2]=1.72;
		mI[2][0]=0.48; mI[2][1]=2.36; mI[2][2]=3.64;
	}
	if(nuc == 'G') {
		mI[0][0]=3.64; mI[0][1]=1.96; mI[0][2]=1.28; 
		mI[1][0]=2.12; mI[1][1]=1.12; mI[1][2]=0.96;
		mI[2][0]=1.92; mI[2][1]=0.32; mI[2][2]=2.56;
	}

	double res[3][3]; // store result
	multiply3(prod, mI, res);


	for (int i = 0; i < 3; i++) 
		for (int j = 0; j < 3; j++)
			prod[i][j] = res[i][j]; 
}

void operation_rep3(double prod[4][4], char nuc) {
	double w1 = 0.5;
	double w2 = 0.3;
	double w3 = 1.0;
	double mI[4][4] = {
    	{1, 0, 0, 0},
    	{0, 1, 0, 0},
    	{0, 0, 1, 0},
    	{0, 0, 0, 1}
    };

	if(nuc == 'A') {
		mI[0][1]=w3; mI[0][2]=w1; mI[0][3]=w2;
	}
	if(nuc == 'T') {
		mI[1][0]=w3; mI[1][2]=w2; mI[1][3]=w1;
	}
	if(nuc == 'C') { 
		mI[2][0]=w1; mI[2][1]=w2; mI[2][3]=w3;
	}
	if(nuc == 'G') {
		mI[3][0]=w2; mI[3][1]=w1; mI[3][2]=w3;
	}

	double res[4][4]; // store result
	multiply4(prod, mI, res);


	for (int i = 0; i < 4; i++) 
		for (int j = 0; j < 4; j++)
			prod[i][j] = res[i][j]; 
}

void operation_rep4(double prod[4][4], char nuc) {
	double da=M_PI/13;
	double dc=M_PI/17;
	double mI[4][4] = {
    	{1, 0, 0, 0},
    	{0, 1, 0, 0},
    	{0, 0, 1, 0},
    	{0, 0, 0, 1}
    };

	if(nuc == 'A') {
		mI[0][0]=cos(da); mI[1][0]=sin(da); mI[0][1]=-sin(da); mI[1][1]=cos(da);
	}
	if(nuc == 'T') {
		mI[1][1]=cos(dc); mI[2][1]=sin(dc); mI[1][2]=-sin(dc); mI[2][2]=cos(dc);
	}
	if(nuc == 'C') { 
		mI[2][2]=cos(da); mI[3][2]=sin(da); mI[2][3]=-sin(da); mI[3][3]=cos(da);
	}
	if(nuc == 'G') {
		mI[0][0]=cos(dc); mI[3][0]=sin(dc); mI[0][3]=-sin(dc); mI[3][3]=cos(dc);
	}

	double res[4][4]; // store result
	multiply4(prod, mI, res);


	for (int i = 0; i < 4; i++) 
		for (int j = 0; j < 4; j++)
			prod[i][j] = res[i][j]; 
}

void operation_rep5(std::complex<long double> prod[2][2], char nuc) {
    std::complex<long double> mI[2][2] = {
        {std::complex<long double> (0.0), std::complex<long double> (0.0)},
        {std::complex<long double> (0.0), std::complex<long double> (0.0)}
    };
    if(nuc == 'A') {
            mI[0][0] = std::complex<long double> (1.0, 0.5); mI[0][1] = std::complex<long double> (-1.0/3, 0.2);
            mI[1][0] = std::complex<long double> (1.0/3, 0.2); mI[1][1] = std::complex<long double> (1.0, -0.5);
    }
    if(nuc == 'T') {
            // mI[0][0] = std::complex<long double> (0.5, 1.0/3); mI[0][1] = std::complex<long double> (-0.2, 1.0);
            // mI[1][0] = std::complex<long double> (0.2, 1.0); mI[1][1] = std::complex<long double> (0.5, -1.0/3);
    		mI[0][0] = std::complex<long double> (1.0, 0.5); mI[0][1] = std::complex<long double> (1.0/3, 0.2);
            mI[1][0] = std::complex<long double> (-1.0/3, 0.2); mI[1][1] = std::complex<long double> (1.0, -0.5);
    }
    if(nuc == 'C') { 
            // mI[0][0] = std::complex<long double> (1.0/3, 0.2); mI[0][1] = std::complex<long double> (-1.0, 0.5);
            // mI[1][0] = std::complex<long double> (1.0, 0.5); mI[1][1] = std::complex<long double> (1.0/3, -0.2);
    		mI[0][0] = std::complex<long double> (0.2, 1.0); mI[0][1] = std::complex<long double> (0.5, 1.0/3);
            mI[1][0] = std::complex<long double> (-0.5, 1.0/3); mI[1][1] = std::complex<long double> (0.2, -1.0);
    }
    if(nuc == 'G') {
            mI[0][0] = std::complex<long double> (0.2, 1.0); mI[0][1] = std::complex<long double> (-0.5, 1.0/3);
            mI[1][0] = std::complex<long double> (0.5, 1.0/3); mI[1][1] = std::complex<long double> (0.2, -1.0);
    }

    std::complex<long double> res[2][2]; // store result
    complexmultiply2(prod, mI, res);


    for (int i = 0; i < 2; i++) 
        for (int j = 0; j < 2; j++)
            prod[i][j] = res[i][j]; 
}

void operation_rep6(double prod[4][4], char nuc) {
	double da=M_PI/13;
	double dc=M_PI/17;
	double mI[4][4] = {
    	{1, 0, 0, 0},
    	{0, 1, 0, 0},
    	{0, 0, 1, 0},
    	{0, 0, 0, 1}
    };

	if(nuc == 'A') {
		mI[0][0]=std::pow(cos(da), 2); mI[1][0]=cos(da)*sin(da); mI[2][0]=-sin(da);
		mI[0][1]=-cos(da)*sin(da)+std::pow(sin(da), 2)*cos(da); mI[1][1]=std::pow(cos(da), 2)+std::pow(sin(da), 3); mI[2][1]=cos(da)*sin(da); 
		mI[0][2]=std::pow(sin(da), 2)+std::pow(cos(da), 2)*sin(da); mI[1][2]=-cos(da)*sin(da)+std::pow(sin(da), 2)*cos(da); mI[2][2]=std::pow(cos(da), 2);
	}
	if(nuc == 'T') {
		mI[1][1]=std::pow(cos(dc), 2); mI[2][1]=cos(dc)*sin(dc); mI[3][1]=-sin(dc);
		mI[1][2]=-cos(dc)*sin(dc)+std::pow(sin(dc), 2)*cos(dc); mI[2][2]=std::pow(cos(dc), 2)+std::pow(sin(dc), 3); mI[3][2]=cos(dc)*sin(dc);
		mI[1][3]=std::pow(sin(dc), 2)+std::pow(cos(dc), 2)*sin(dc); mI[2][3]=-cos(dc)*sin(dc)+std::pow(sin(dc), 2)*cos(dc); mI[3][3]=std::pow(cos(dc), 2);
	}
	if(nuc == 'C') { 
		mI[0][0]=std::pow(cos(da), 2); mI[2][0]=cos(da)*sin(da); mI[3][0]=-sin(da); 
		mI[0][2]=-cos(da)*sin(da)+std::pow(sin(da), 2)*cos(da); mI[2][2]=std::pow(cos(da), 2)+std::pow(sin(da), 3); mI[3][2]=cos(da)*sin(da);
		mI[0][3]=std::pow(sin(da), 2)+std::pow(cos(da), 2)*sin(da); mI[2][3]=-cos(da)*sin(da)+std::pow(sin(da), 2)*cos(da); mI[3][3]=std::pow(cos(da), 2);
	}
	if(nuc == 'G') {
		mI[0][0]=std::pow(cos(dc), 2); mI[1][0]=cos(dc)*sin(dc); mI[3][0]=-sin(dc);
		mI[0][1]=-cos(dc)*sin(dc)+std::pow(sin(dc), 2)*cos(dc); mI[1][1]=std::pow(cos(dc), 2)+std::pow(sin(dc), 3); mI[3][1]=cos(dc)*sin(dc); 
		mI[0][3]=std::pow(sin(dc), 2)+std::pow(cos(dc), 2)*sin(dc); mI[1][3]=-cos(dc)*sin(dc)+std::pow(sin(dc), 2)*cos(dc); mI[3][3]=std::pow(cos(dc), 2);
	}

	double res[4][4]; // store result
	multiply4(prod, mI, res);


	for (int i = 0; i < 4; i++) 
		for (int j = 0; j < 4; j++)
			prod[i][j] = res[i][j]; 
}

void simreads(void) {
	map<int, char> i2c; 
	i2c[0] = 'A'; 
	i2c[1] = 'T'; 
	i2c[2] = 'C'; 
	i2c[3] = 'G'; 
    int len = 8; 
	long int * div = new long int[len]; 
	
	div[len-1] = 1; 
	for (int i = 2; i <= len; i++)
		div[len-i] = div[len-i+1] * 4; 
	long int comb = div[0] * 4; 

    for (long int i = 0; i < comb; i++)
	{
		long int res = i; 
//		char * str = new char[len]; 
        for (int j = 0; j < len; j++)
		{
			int ii = (int) (res/div[j]); 
			res -= ii * div[j]; 
//			str[j] = (i2c[ii]); 
		}
	}
}

int compare_doubles(const void * a, const void * b)
{
	return ((int) (*(double *) a > *(double *) b));  
}

/*for how to write compare functions, see http://www.cplusplus.com/reference/cstdlib/qsort/
https://stackoverflow.com/questions/14092642/sorting-a-2-dimensional-array-in-c
https://stackoverflow.com/questions/17202178/c-qsort-with-dynamic-n-by-2-multi-dimensional-array*/
int compare_doubles_2D(const void *pa, const void *pb) {
    const double *a = *(const double **)pa;
    const double *b = *(const double **)pb;
    if(a[4] > b[4])
        return 1;
    else if(a[4] < b[4])
        return -1;
    else
        return 0;
}

double pairwiseL2(double * xx, double * yy, int nx, int ny)
{
	qsort((double *) xx, nx, sizeof(double), compare_doubles); 
	qsort((double *) yy, ny, sizeof(double), compare_doubles); 
	double * hh = new double[nx]; 
	double * ww = new double[nx]; 
	for (int i = 0; i < nx; i++)
	{
		ww[i] = (xx[i]+yy[i])/2; 
		hh[i] = abs(xx[i]-yy[i]); 
	}
	for (int i = 0; i < nx-1; i++)
		ww[i] = ww[i+1]-ww[i]; 
	ww[nx-1] = 0; 
	double sum = 0; 
	for (int i = 0; i < nx; i++)
		sum +=  ww[i] * hh[i]; 
	delete[] hh; 
	delete[] ww; 
	return(sum); 
}

double * select_quantiles_old(double * xx, double * yy, int nx, int ny, int ng)
{
	// combine the two pre-sorted vectors and sort it into one vector
	double * xy = new double[nx+ny];
	memset(xy, 1, sizeof(double)*(nx+ny));
	int i=0; int j=0; int k=0;
	while (i<nx && j<ny) {
		if (xx[i]<=yy[j]) {
			xy[k++]=xx[i++];
		}	else {
			xy[k++]=yy[j++];
		}
	}
	if (i==nx) {
		while (j<ny) {
			xy[k++]=yy[j++];
		}
	}	else {
		while (i<nx) {
			xy[k++]=xx[i++];
		}
	}

	// get the quantiles of the vector of numbers
	double * quantiles = new double[ng+1];
	for (int i=0; i<ng; i++) {
		quantiles[i]=xy[(i*(nx+ny))/ng];
	}
	quantiles[ng]=xy[nx+ny-1]+0.1;

	delete[] xy;
	return(quantiles);
}

double * select_quantiles(double * xx, long nx, int ng)
{
	double * temp = new double[nx];
	std::copy(xx, xx + nx, temp);
	double * brk = new double[ng+1]; // breaks
	qsort((double *) temp, nx, sizeof(double), compare_doubles);
	for (int j = 0; j < ng; j++) {
		brk[j] = temp[(nx / ng) * j];
	}
	// brk[ng] = temp[nx - 1] + 0.1;
	brk[ng] = INFINITY;
	brk[0] = -1 * INFINITY;
	delete[] temp;
	return(brk);
}

double ** select_quantiles_all(double ** xx, long nx, int ng, int ncol)
{
	double ** brk = new double*[ncol]; //for each of columns
	for (int i = 0; i < ncol; i++) {
		brk[i] = new double[ng+1]; // breaks
		qsort((double *) xx[i], nx, sizeof(double), compare_doubles);
		for (int j = 0; j < ng; j++) {
			brk[i][j] = xx[i][(nx / ng) * j];
		}
		brk[i][ng] = xx[i][nx - 1] + 0.1;
	}
	return(brk);
}

int quickselect(double * brk, int ng, double num)
{
	// ng number of bins, num the number to be placed in a certain bin
	// brk, the break points of bin for a single dimension
	int start = 0;
	int end = ng - 1;
	while (start < end) {
		int middle = start + ceil((end - start) / 2.0);
		if (num < brk[middle]) {
			end = middle - 1;
		}
		else {
			start = middle;
		}
	}
	return(start);
}

double hellinger(double * xx, double * yy, int nx, int ny)
{
//	double a1,a2,a3; 
//	a1=a2=a3=1; 
//	for (int i = 1; i < seqLen; i++) 
//	{
//		a3=a1+a2; 
//		a1=a2; 
//		a2=a3; 
////		printf("%lf %lf %lf \n", a1, a2, a3); 
//	}
//    printf("%lf %d \n",log(a3), seqLen);  
	int ng = 100; 
	double * brk = new double[ng+1]; 
	double * wx = new double[ng]; 
	double * wy = new double[ng]; 
	int initial_value = 1;
	memset(wx, initial_value, sizeof(double) * ng); 
	memset(wy, initial_value, sizeof(double) * ng); 
	qsort((double *) xx, nx, sizeof(double), compare_doubles); 
	qsort((double *) yy, ny, sizeof(double), compare_doubles); 

//	int xp1=(int) (nx*0.0001); 
//	int yp1=(int) (ny*0.0001); 
//	int xp9=(int) (nx*0.9999); 
//	int yp9=(int) (ny*0.9999); 

	brk[0] = min(xx[0], yy[0]); 
	brk[ng] = max(xx[nx-1], yy[ny-1]); 
	for (int i = 1; i < ng; i++)
		brk[i] = brk[0] + i * (brk[ng]-brk[0])/ng; 
	brk[ng] = max(xx[nx-1], yy[ny-1]) + 0.1; 

//    brk[0] = 28; 
//	brk[ng-1] = 45; 
//	brk[ng] = log(a3); 

//	for (int i = 0; i <= ng; i++)
//		brk[i] = i; 
//	brk[ng] += 0.1; 

	int bx = 0; 
	for (int i = 0; i < nx; i++)
	{
		if(xx[i] < brk[bx+1]) {
			wx[bx] += 1; 
		}
		else {
			wx[++bx] += 1; 
			//if(bx > ng) break; unnec
		}
	}
	int by = 0; 
	for (int i = 0; i < ny; i++)
	{
		if(yy[i] < brk[by+1]) {
			wy[by] += 1; 
		}
		else {
			wy[++by] += 1; 
			//if(by > ng) break; unnec
		}
	}

	for (int i = 0; i < ng; i++)
	{
		wx[i] /= (nx + ng * initial_value); 
		wy[i] /= (ny + ng * initial_value); 
	}

	double res = 0; 
	for (int i = 0; i < ng; i++)
	{
		double temp = sqrt(wx[i]) - sqrt(wy[i]); 
		res += temp * temp; 
	}
	res = sqrt(res); 
	res *= 1/sqrt(2); 
		
	delete[] brk; 
	delete[] wx; 
	delete[] wy; 

	return(res); 
}

double hellinger_1D_quantiles_old(double * xx, double * yy, int nx, int ny)
{
	int ng = 100; 
	int initial_value = 1;
	double * brk = new double[ng+1]; 
	double * wx = new double[ng]; 
	double * wy = new double[ng]; 
	memset(wx, initial_value, sizeof(double) * ng); 
	memset(wy, initial_value, sizeof(double) * ng); 
	qsort((double *) xx, nx, sizeof(double), compare_doubles); 
	qsort((double *) yy, ny, sizeof(double), compare_doubles); 

	brk = select_quantiles_old(xx, yy, nx, ny, ng);

	int bx = 0; 
	for (int i = 0; i < nx; i++)
	{
		if(xx[i] < brk[bx+1]) {
			wx[bx] += 1; 
		}
		else {
			wx[bx++] += 1; 
			if(bx > ng) break; 
		}
	}
	int by = 0; 
	for (int i = 0; i < ny; i++)
	{
		if(yy[i] < brk[by+1]) {
			wy[by] += 1; 
		}
		else {
			wy[by++] += 1; 
			if(by > ng) break; 
		}
	}

	for (int i = 0; i < ng; i++)
	{
		wx[i] /= (nx + ng * initial_value); 
		wy[i] /= (ny + ng * initial_value); 
	}

	double res = 0; 
	for (int i = 0; i < ng; i++)
	{
		double temp = sqrt(wx[i]) - sqrt(wy[i]); 
		res += temp * temp; 
	}
	res = sqrt(res); 
	res *= 1/sqrt(2); 
		
	delete[] brk; 
	delete[] wx; 
	delete[] wy; 

	return(res); 
}

double hellinger_1D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng) // xx and yy are arrays with 1 columns
{
	//int ng = 15; // number of bins, you can adjust it
	std::vector<double> wx; //1D vector
	std::vector<double> wy; //1D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	// Set up initial values
	int initial_value = 0;
	for (int i = 0; i < ng; i++) {
		wx[i] = initial_value;
		wy[i] = initial_value;
	}

	/* get the density for each 'volumn' */  
	for (int i = 0; i < nx; i++) { // for each line
		int bin_num[1]; //bin_num track which bin should the data point go in for 'one of the 1 columns'
		for ( int j = 0; j < 1; j++ ) { // for each column
			for ( int bx = 0; bx < ng; bx++ ) {
				if ( xx[i][j] < brk[j][bx + 1] ) {
					bin_num[j] = bx;
					break; // if an interval found, exit
				}
			}
		}
		wx[bin_num[0]] += 1;
	}
	for (int i = 0; i < ny; i++) { // for each line
		int bin_num[1]; //bin_num track which bin should the data point go in for 'one of the 1 columns'
		for ( int j = 0; j < 1; j++ ) { // for each column
			for ( int by = 0; by < ng; by++ ) {
				if ( yy[i][j] < brk[j][by + 1] ) {
					bin_num[j] = by;
					break; // if an interval found, exit
				}
			}
		}
		wy[bin_num[0]] += 1;
	}

	for (int i = 0; i < ng; i++) {
		wx[i] /= (nx + initial_value * std::pow(ng, 1));
		wy[i] /= (ny + initial_value * std::pow(ng, 1));
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	for (int i = 0; i < ng; i++) {
		double temp = sqrt(wx[i]) - sqrt(wy[i]);
		res += temp * temp;
	}
	res = sqrt(res);
	res *= 1/sqrt(2);
	return(res);
}

double hellinger_2D(double ** xx, double ** yy, int nx, int ny) // xx and yy are arrays with 5 columns
{
	int ng = 100; int nh = 100; // you can adjust it
	int volumns = ng * nh;
	double brk[ng+1]; 
	double wx[volumns][4]; 
	double wy[volumns][4]; 
	memset(wx, 0, sizeof(double) * volumns * 4); 
	memset(wy, 0, sizeof(double) * volumns * 4); 
	qsort((double **) xx, nx, sizeof(double), compare_doubles_2D); 
	qsort((double **) yy, ny, sizeof(double), compare_doubles_2D); 
	/* find the breaking points for the last column */
	brk[0] = min(xx[0][4], yy[0][4]); 
	brk[ng] = max(xx[nx-1][4], yy[ny-1][4]); 
	for (int i = 1; i < ng; i++)
		brk[i] = brk[0] + i * (brk[ng]-brk[0])/ng; 
	brk[ng] += 0.1; 
	/* find the min and max for first 4 columns */
	double max_xx[4];
	double min_xx[4];
	double max_yy[4];
	double min_yy[4];
	for ( int i = 0; i < 4 ; i++ ) { //for each column
		min_xx[i] = xx[0][i]; //Assume first element in i th Column is Minimum
		max_xx[i] = xx[0][i]; //Assume first element in i th Column is Maximum
		min_yy[i] = yy[0][i]; //Assume first element in i th Column is Minimum
		max_yy[i] = yy[0][i]; //Assume first element in i th Column is Maximum
		for ( int j = 1; j < nx; j++ ){ // for each row 
			if (max_xx[i] < xx[j][i])
				max_xx[i] = xx[j][i]; // update max		
			else if (min_xx[i] > xx[j][i])
				min_xx[i] = xx[j][i]; // update min
		}
		for ( int j = 1; j < ny; j++ ){ // for each row 
			if (max_yy[i] < yy[j][i])
				max_yy[i] = yy[j][i]; // update max		
			else if (min_yy[i] > yy[j][i])
				min_yy[i] = yy[j][i]; // update min
		}
	}
	/* get the bin size for each pair of first 4 columns */
	double bin_size[4];
	for ( int i = 0; i < 4 ; i++ ) { //for each pair of column
		bin_size[i]=((max(max_xx[i], max_yy[i]) - min(min_xx[i], min_yy[i]))/nh)+0.1;
	}	
	/* get the density for each volumn, for each of A, T, C, G */ 
	for (int atcg = 0; atcg < 4; atcg++) { // for each of A, T, C, G
		int bx = 0; int bin_num; // bx tracks which bin should the data point in according to the last column; bin_num tracks one of the first 4 columns
		for (int i = 0; i < nx; i++)
		{
			if(xx[i][4] < brk[bx+1]) {
				bin_num = int(xx[i][atcg]/bin_size[atcg]);
				wx[bx * nh + bin_num][atcg] += 1; 
			}
			else {
				bin_num = int(xx[i][atcg]/bin_size[atcg]);
				wx[(++bx) * nh + bin_num][atcg] += 1; 
			}
		}
	}
	for (int atcg = 0; atcg < 4; atcg++) { // for each of A, T, C, G
		int by = 0; int bin_num; // by tracks which bin should the data point in according to the last column; bin_num tracks one of the first 4 columns
		for (int i = 0; i < ny; i++)
		{
			if(yy[i][4] < brk[by+1]) {
				bin_num = int(yy[i][atcg]/bin_size[atcg]);
				wy[by * nh + bin_num][atcg] += 1; 
			}
			else {
				bin_num = int(yy[i][atcg]/bin_size[atcg]);
				wy[(++by) * nh + bin_num][atcg] += 1; 
			}
		}
	}
	for (int atcg = 0; atcg < 4; atcg++) { // for each of A, T, C, G
		for (int i = 0; i < volumns; i++)
		{
			wx[i][atcg] /= (nx); 
			wy[i][atcg] /= (ny); 
		}
	}
	/* calculate the pair-wise distance for each pair of columns */
	double res_final = 0;
	for (int atcg = 0; atcg < 4; atcg++) { // for each of A, T, C, G
		double res = 0; 
		for (int i = 0; i < volumns; i++)
		{
			double temp = sqrt(wx[i][atcg]) - sqrt(wy[i][atcg]); 
			res += temp * temp; 
		}
		res = sqrt(res); 
		res *= 1/sqrt(2); 
		res_final += res;
	}
	res_final /= 4;
	return(res_final); 
}

double hellinger_4D(double ** xx, double ** yy, int nx, int ny) // xx and yy are arrays with 4 columns
{
	int ng = 15; // number of bins, you can adjust it
	std::vector<std::vector<std::vector<std::vector<double> > > > wx; //4D vector
	std::vector<std::vector<std::vector<std::vector<double> > > > wy; //4D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	for (int i = 0; i < ng; i++) {
		wx[i].resize(ng);
		wy[i].resize(ng);
		for (int j = 0; j < ng; j++) {
			wx[i][j].resize(ng);
			wy[i][j].resize(ng);
			for (int k = 0; k < ng; k++) {
				wx[i][j][k].resize(ng);
				wy[i][j][k].resize(ng);
			}
		}
	}
	// Set up initial values
	int initial_value = 10;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					wx[i][j][k][l] = initial_value;
					wy[i][j][k][l] = initial_value;
				}
			}
		}
	}

	/* find the min and max for 4 columns of Nubeam scores*/
	double max_xx[4];
	double min_xx[4];
	double max_yy[4];
	double min_yy[4];
	for ( int i = 0; i < 4 ; i++ ) { //for each column, A, T, C, G
		min_xx[i] = xx[0][i]; //Assume first element in i th Column is Minimum
		max_xx[i] = xx[0][i]; //Assume first element in i th Column is Maximum
		min_yy[i] = yy[0][i]; //Assume first element in i th Column is Minimum
		max_yy[i] = yy[0][i]; //Assume first element in i th Column is Maximum
		for ( int j = 1; j < nx; j++ ){ // for each row 
			if (max_xx[i] < xx[j][i])
				max_xx[i] = xx[j][i]; // update max		
			else if (min_xx[i] > xx[j][i])
				min_xx[i] = xx[j][i]; // update min
		}
		for ( int j = 1; j < ny; j++ ){ // for each row 
			if (max_yy[i] < yy[j][i])
				max_yy[i] = yy[j][i]; // update max		
			else if (min_yy[i] > yy[j][i])
				min_yy[i] = yy[j][i]; // update min
		}
	}

	/* get the bin size for each pair of the 4 columns */
	double bin_size[4];
	double mins[4];
	for ( int i = 0; i < 4 ; i++ ) { //for each pair of column
		mins[i] = min(min_xx[i], min_yy[i]);
		bin_size[i] = ((max(max_xx[i], max_yy[i]) - mins[i])/ng) + 0.000001;
	}

	/* get the density for each 'volumn' */ 
	int bin_num1; int bin_num2; int bin_num3; int bin_num4; //bin_num track which bin should the data point go in for one of the 4 columns
	for (int i = 0; i < nx; i++)
	{
		bin_num1 = int((xx[i][0] - mins[0])/bin_size[0]);
		bin_num2 = int((xx[i][1] - mins[1])/bin_size[1]);
		bin_num3 = int((xx[i][2] - mins[2])/bin_size[2]);
		bin_num4 = int((xx[i][3] - mins[3])/bin_size[3]);
		wx[bin_num1][bin_num2][bin_num3][bin_num4] += 1;
	}
	for (int i = 0; i < ny; i++)
	{
		bin_num1 = int((yy[i][0] - mins[0])/bin_size[0]);
		bin_num2 = int((yy[i][1] - mins[1])/bin_size[1]);
		bin_num3 = int((yy[i][2] - mins[2])/bin_size[2]);
		bin_num4 = int((yy[i][3] - mins[3])/bin_size[3]);
		wy[bin_num1][bin_num2][bin_num3][bin_num4] += 1;
	}
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					wx[i][j][k][l] /= (nx + initial_value * std::pow(ng, 4));
					wy[i][j][k][l] /= (ny + initial_value * std::pow(ng, 4));
				}
			}
		}
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					double temp = sqrt(wx[i][j][k][l]) - sqrt(wy[i][j][k][l]);
					res += temp * temp;
				}
			}
		}
	}
	res = sqrt(res);
	res *= 1/sqrt(2);
	return(res);
}

double hellinger_2D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng) // xx and yy are arrays with 2 columns
{
	// int ng = 50; // number of bins, you can adjust it
	std::vector<std::vector<double>  > wx; //2D vector
	std::vector<std::vector<double>  > wy; //2D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	for (int i = 0; i < ng; i++) {
		wx[i].resize(ng);
		wy[i].resize(ng);
	}
	// Set up initial values
	int initial_value = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			wx[i][j] = initial_value;
			wy[i][j] = initial_value;
		}
	}

	// /* get the breaks for each of 4 columns according to the quantiles */
	// double ** brk = new double*[4];
	// for ( int j = 0; j < 4; j++ ) { // for each column of A, T, C and G
	// 	double * xx_col = new double[nx]; // replicate data for a column
	// 	for ( int i = 0; i < nx; i++ )
	// 		xx_col[i] = xx[i][j];
	// 	qsort((double *) xx_col, nx, sizeof(double), compare_doubles);
	// 	double * yy_col = new double[ny];
	// 	for ( int i = 0; i < ny; i++ )
	// 		yy_col[i] = yy[i][j];
	// 	qsort((double *) yy_col, ny, sizeof(double), compare_doubles); 
	// 	brk[j] = new double[ng+1]; // breaks
	// 	brk[j] = select_quantiles(xx_col, yy_col, nx, ny, ng);
	// 	delete[] xx_col;
	// 	delete[] yy_col;
	// }

	/* get the density for each 'volumn' */  
	for (int i = 0; i < nx; i++) { // for each line
		int bin_num[2]; //bin_num track which bin should the data point go in for one of the 2 columns
		for ( int j = 0; j < 2; j++ ) { // for each of the 2 columns
			for ( int bx = 0; bx < ng; bx++ ) {
				if ( xx[i][j] < brk[j][bx + 1] ) {
					bin_num[j] = bx;
					break; // if an interval found, exit
				}
			}
		}
		wx[bin_num[0]][bin_num[1]] += 1;
	}
	for (int i = 0; i < ny; i++) { // for each line
		int bin_num[2]; //bin_num track which bin should the data point go in for one of the 2 columns
		for ( int j = 0; j < 2; j++ ) { // for each of the 2 columns
			for ( int by = 0; by < ng; by++ ) {
				if ( yy[i][j] < brk[j][by + 1] ) {
					bin_num[j] = by;
					break; // if an interval found, exit
				}
			}
		}
		wy[bin_num[0]][bin_num[1]] += 1;
	}

	// for ( int j = 0; j < 4; j++ )
	// 	delete[] brk[j];
	// delete[] brk;

	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			wx[i][j] /= (nx + initial_value * std::pow(ng, 2));
			wy[i][j] /= (ny + initial_value * std::pow(ng, 2));
		}
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			double temp = sqrt(wx[i][j]) - sqrt(wy[i][j]);
			res += temp * temp;
		}
	}
	res = sqrt(res);
	res *= 1/sqrt(2);
	return(res);
}

double hellinger_3D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng) // xx and yy are arrays with 3 columns
{
	//int ng = 15; // number of bins, you can adjust it
	std::vector<std::vector<std::vector<double> > > wx; //3D vector
	std::vector<std::vector<std::vector<double> > > wy; //3D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	for (int i = 0; i < ng; i++) {
		wx[i].resize(ng);
		wy[i].resize(ng);
		for (int j = 0; j < ng; j++) {
			wx[i][j].resize(ng);
			wy[i][j].resize(ng);
		}
	}
	// Set up initial values
	int initial_value = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				wx[i][j][k] = initial_value;
				wy[i][j][k] = initial_value;
			}
		}
	}

	/* get the density for each 'volumn' */  
	for (int i = 0; i < nx; i++) { // for each line
		int bin_num[3]; //bin_num track which bin should the data point go in for one of the 3 columns
		for ( int j = 0; j < 3; j++ ) { // for each column
			for ( int bx = 0; bx < ng; bx++ ) {
				if ( xx[i][j] < brk[j][bx + 1] ) {
					bin_num[j] = bx;
					break; // if an interval found, exit
				}
			}
		}
		wx[bin_num[0]][bin_num[1]][bin_num[2]] += 1;
	}
	for (int i = 0; i < ny; i++) { // for each line
		int bin_num[3]; //bin_num track which bin should the data point go in for one of the 3 columns
		for ( int j = 0; j < 3; j++ ) { // for each column
			for ( int by = 0; by < ng; by++ ) {
				if ( yy[i][j] < brk[j][by + 1] ) {
					bin_num[j] = by;
					break; // if an interval found, exit
				}
			}
		}
		wy[bin_num[0]][bin_num[1]][bin_num[2]] += 1;
	}

	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				wx[i][j][k] /= (nx + initial_value * std::pow(ng, 3));
				wy[i][j][k] /= (ny + initial_value * std::pow(ng, 3));
			}
		}
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				double temp = sqrt(wx[i][j][k]) - sqrt(wy[i][j][k]);
				res += temp * temp;
			}
		}
	}
	res = sqrt(res);
	res *= 1/sqrt(2);
	return(res);
}


double l1_4D_quantiles(unordered_map<string, double> wx, unordered_map<string, double> wy) // wx and wy are unordered hash maps
{
	/* calculate the pair-wise distance for two samples */
	double res = 0;
	unordered_map<string, double>:: iterator itr; 
	for (itr = wx.begin(); itr != wx.end(); itr++)
	{
		// itr works as a pointer to pair<string, double> 
		// type itr->first stores the key part  and 
		// itr->second stroes the value part
		// string key = itr->first;
		if (wy.find(itr->first) == wy.end()) {
			res += itr->second;
		}
		else {
			res += abs(itr->second - wy[itr->first]);
		}
	}

	for (itr = wy.begin(); itr != wy.end(); itr++) 
	{
		// string key = itr->first;
		if (wx.find(itr->first) == wx.end()) {
			res += itr->second;
		}
	}

	return(res);
}

double hellinger_4D_quantiles(unordered_map<string, double> wx, unordered_map<string, double> wy, double pseudo_freq) // wx and wy are unordered hash maps
{
	/* calculate the pair-wise distance for two samples */
	double res = 0;
	double temp = 0;
	// int total_bin_num = 0;
	// int total_different_bin_num = 0;
	unordered_map<string, double>:: iterator itr; 
	for (itr = wx.begin(); itr != wx.end(); itr++)
	{
		// itr works as a pointer to pair<string, double> 
		// type itr->first stores the key part  and 
		// itr->second stroes the value part
		// string key = itr->first;
		if (wy.find(itr->first) == wy.end()) {
			// res += itr->second;
			temp = sqrt(itr->second) - sqrt(pseudo_freq);
			res += temp * temp;
			// total_different_bin_num++;
		}
		else {
			temp = sqrt(itr->second) - sqrt(wy[itr->first]);
			res += temp * temp;
			// if (temp != 0) {
			// 	total_different_bin_num++;
			// }
		}
		// total_bin_num++;
	}

	for (itr = wy.begin(); itr != wy.end(); itr++) 
	{
		// string key = itr->first;
		if (wx.find(itr->first) == wx.end()) {
			// res += itr->second;
			temp = sqrt(itr->second) - sqrt(pseudo_freq);
			res += temp * temp;
			// total_different_bin_num++;
		}
		// total_bin_num++;
	}

	res = sqrt(res);
	res *= 1/sqrt(2);
	// printf("%d\t%d\n", total_different_bin_num, total_bin_num);
	// cout << total_different_bin_num << "\n";
	// cout << total_bin_num << "\n";
	return(res);
}

double hamming_4D_quantiles(unordered_map<string, double> wx, unordered_map<string, double> wy) // wx and wy are unordered hash maps
{
	/* calculate the pair-wise distance for two samples */
	double res = 0;
	unordered_map<string, double>:: iterator itr; 
	for (itr = wx.begin(); itr != wx.end(); itr++)
	{
		// itr works as a pointer to pair<string, double> 
		// type itr->first stores the key part  and 
		// itr->second stroes the value part
		// string key = itr->first;
		if (wy.find(itr->first) == wy.end()) {
			res += 1;
		}
	}

	for (itr = wy.begin(); itr != wy.end(); itr++) 
	{
		// string key = itr->first;
		if (wx.find(itr->first) == wx.end()) {
			res += 1;
		}
	}

	res /= (wx.size() + wy.size());
	return(res);
}

double hellinger_5D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk) // xx and yy are arrays with 5 columns
{
	int ng = 5; // number of bins, you can adjust it
	std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > wx; //5D vector
	std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > wy; //5D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	for (int i = 0; i < ng; i++) {
		wx[i].resize(ng);
		wy[i].resize(ng);
		for (int j = 0; j < ng; j++) {
			wx[i][j].resize(ng);
			wy[i][j].resize(ng);
			for (int k = 0; k < ng; k++) {
				wx[i][j][k].resize(ng);
				wy[i][j][k].resize(ng);
				for (int l = 0; l < ng; l++) {
					wx[i][j][k][l].resize(ng);
					wy[i][j][k][l].resize(ng);
				}
			}
		}
	}
	// Set up initial values
	int initial_value = 10;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					for (int m = 0; m < ng; m++) {
						wx[i][j][k][l][m] = initial_value;
						wy[i][j][k][l][m] = initial_value;
					}
				}
			}
		}
	}

	/* get the density for each 'volumn' */
	for (int i = 0; i < nx; i++) { // for each line
		int bin_num[5]; //bin_num track which bin should the data point go in for one of the 5 columns
		for ( int j = 0; j < 5; j++ ) { // for each column
			for ( int bx = 0; bx < ng; bx++ ) {
				if ( xx[i][j] < brk[j][bx + 1] ) {
					bin_num[j] = bx;
					break; // if an interval found, exit
				}
			}
		}
		wx[bin_num[0]][bin_num[1]][bin_num[2]][bin_num[3]][bin_num[4]] += 1;
	}
	for (int i = 0; i < ny; i++) { // for each line
		int bin_num[5]; //bin_num track which bin should the data point go in for one of the 5 columns
		for ( int j = 0; j < 5; j++ ) { // for each column
			for ( int by = 0; by < ng; by++ ) {
				if ( yy[i][j] < brk[j][by + 1] ) {
					bin_num[j] = by;
					break; // if an interval found, exit
				}
			}
		}
		wy[bin_num[0]][bin_num[1]][bin_num[2]][bin_num[3]][bin_num[4]] += 1;
	}

	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					for (int m = 0; m < ng; m++) {
						wx[i][j][k][l][m] /= (nx + initial_value * std::pow(ng, 5));
						wy[i][j][k][l][m] /= (ny + initial_value * std::pow(ng, 5));
					}
				}
			}
		}
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					for (int m = 0; m < ng; m++) {
						double temp = sqrt(wx[i][j][k][l][m]) - sqrt(wy[i][j][k][l][m]);
						res += temp * temp;
					}
				}
			}
		}
	}
	res = sqrt(res);
	res *= 1/sqrt(2);
	return(res);
}

double cosine(double * xx, double * yy, int nx, int ny)
{
//	double a1,a2,a3; 
//	a1=a2=a3=1; 
//	for (int i = 1; i < seqLen; i++) 
//	{
//		a3=a1+a2; 
//		a1=a2; 
//		a2=a3; 
////		printf("%lf %lf %lf \n", a1, a2, a3); 
//	}
//    printf("%lf %d \n",log(a3), seqLen);  
	int ng = 100; 
	int initial_value = 1;
	double * brk = new double[ng+1]; 
	double * wx = new double[ng]; 
	double * wy = new double[ng]; 
	memset(wx, initial_value, sizeof(double) * ng); 
	memset(wy, initial_value, sizeof(double) * ng); 
	qsort((double *) xx, nx, sizeof(double), compare_doubles); 
	qsort((double *) yy, ny, sizeof(double), compare_doubles); 

//	int xp1=(int) (nx*0.0001); 
//	int yp1=(int) (ny*0.0001); 
//	int xp9=(int) (nx*0.9999); 
//	int yp9=(int) (ny*0.9999); 

	brk[0] = min(xx[0], yy[0]); 
	brk[ng] = max(xx[nx-1], yy[ny-1]); 
	for (int i = 1; i < ng; i++)
		brk[i] = brk[0] + i * (brk[ng]-brk[0])/ng; 
	brk[ng] = max(xx[nx-1], yy[ny-1]) + 0.1; 

//    brk[0] = 28; 
//	brk[ng-1] = 45; 
//	brk[ng] = log(a3); 

//	for (int i = 0; i <= ng; i++)
//		brk[i] = i; 
//	brk[ng] += 0.1; 

	int bx = 0; 
	for (int i = 0; i < nx; i++)
	{
		if(xx[i] < brk[bx+1]) {
			wx[bx] += 1; 
		}
		else {
			wx[++bx] += 1; 
			//if(bx > ng) break; unnec
		}
	}
	int by = 0; 
	for (int i = 0; i < ny; i++)
	{
		if(yy[i] < brk[by+1]) {
			wy[by] += 1; 
		}
		else {
			wy[++by] += 1; 
			//if(by > ng) break; unnec
		}
	}

	for (int i = 0; i < ng; i++)
	{
		wx[i] /= (nx + ng * initial_value); 
		wy[i] /= (ny + ng * initial_value); 
	}

	double res = 0; 
	double temp1 = 0;
	double temp2 = 0;
	double temp3 = 0;
	for (int i = 0; i < ng; i++)
	{
		temp1 += wx[i] * wy[i];
		temp2 += wx[i] * wx[i];
		temp3 += wy[i] * wy[i];
	}
	res = temp1 / (sqrt(temp2) * sqrt(temp3));
	res = 1 - res; 
	res *= 0.5; 
		
	delete[] brk; 
	delete[] wx; 
	delete[] wy; 

	return(res); 
}

double cosine_1D_quantiles_old(double * xx, double * yy, int nx, int ny)
{
	int ng = 100; 
	int initial_value = 1;
	double * brk = new double[ng+1]; 
	double * wx = new double[ng]; 
	double * wy = new double[ng]; 
	memset(wx, initial_value, sizeof(double) * ng); 
	memset(wy, initial_value, sizeof(double) * ng); 
	qsort((double *) xx, nx, sizeof(double), compare_doubles); 
	qsort((double *) yy, ny, sizeof(double), compare_doubles); 

	brk = select_quantiles_old(xx, yy, nx, ny, ng);

	int bx = 0; 
	for (int i = 0; i < nx; i++)
	{
		if(xx[i] < brk[bx+1]) {
			wx[bx] += 1; 
		}
		else {
			wx[++bx] += 1;
		}
	}
	int by = 0; 
	for (int i = 0; i < ny; i++)
	{
		if(yy[i] < brk[by+1]) {
			wy[by] += 1; 
		}
		else {
			wy[++by] += 1;
		}
	}

	for (int i = 0; i < ng; i++)
	{
		wx[i] /= (nx + ng * initial_value); 
		wy[i] /= (ny + ng * initial_value); 
	}

	double res = 0; 
	double temp1 = 0;
	double temp2 = 0;
	double temp3 = 0;
	for (int i = 0; i < ng; i++)
	{
		temp1 += wx[i] * wy[i];
		temp2 += wx[i] * wx[i];
		temp3 += wy[i] * wy[i];
	}
	res = temp1 / (sqrt(temp2) * sqrt(temp3));
	res = 1 - res; 
	res *= 0.5; 
		
	delete[] brk; 
	delete[] wx; 
	delete[] wy; 

	return(res); 
}

double cosine_1D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng) // xx and yy are arrays with 1 column
{
	std::vector<double> wx; //1D vector
	std::vector<double> wy; //1D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	// Set up initial values
	int initial_value = 0;
	for (int i = 0; i < ng; i++) {
		wx[i] = initial_value;
		wy[i] = initial_value;
	}

	/* get the density for each 'volumn' */  
	for (int i = 0; i < nx; i++) { // for each line
		int bin_num[1]; //bin_num track which bin should the data point go in for 'one of the 1 columns'
		for ( int j = 0; j < 1; j++ ) { // for each column
			for ( int bx = 0; bx < ng; bx++ ) {
				if ( xx[i][j] < brk[j][bx + 1] ) {
					bin_num[j] = bx;
					break; // if an interval found, exit
				}
			}
		}
		wx[bin_num[0]] += 1;
	}
	for (int i = 0; i < ny; i++) { // for each line
		int bin_num[1]; //bin_num track which bin should the data point go in for 'one of the 1 columns'
		for ( int j = 0; j < 1; j++ ) { // for each column
			for ( int by = 0; by < ng; by++ ) {
				if ( yy[i][j] < brk[j][by + 1] ) {
					bin_num[j] = by;
					break; // if an interval found, exit
				}
			}
		}
		wy[bin_num[0]] += 1;
	}

	for (int i = 0; i < ng; i++) {
		//printf("%f\t", wx[i]);
		wx[i] /= (nx + initial_value * std::pow(ng, 1));
		//printf("%f\n", wy[i]);
		wy[i] /= (ny + initial_value * std::pow(ng, 1));
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	double temp1 = 0;
	double temp2 = 0;
	double temp3 = 0;
	for (int i = 0; i < ng; i++) {
		temp1 += wx[i] * wy[i];
		temp2 += wx[i] * wx[i];
		temp3 += wy[i] * wy[i];
	}
	res = temp1 / (sqrt(temp2) * sqrt(temp3));
	res = 1 - res;
	res *= 0.5;
	return(res);
}


double cosine_4D(double ** xx, double ** yy, int nx, int ny) // xx and yy are arrays with 4 columns
{
	int ng = 15; // number of bins, you can adjust it
	std::vector<std::vector<std::vector<std::vector<double> > > > wx; //4D vector
	std::vector<std::vector<std::vector<std::vector<double> > > > wy; //4D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	for (int i = 0; i < ng; i++) {
		wx[i].resize(ng);
		wy[i].resize(ng);
		for (int j = 0; j < ng; j++) {
			wx[i][j].resize(ng);
			wy[i][j].resize(ng);
			for (int k = 0; k < ng; k++) {
				wx[i][j][k].resize(ng);
				wy[i][j][k].resize(ng);
			}
		}
	}
	// Set up initial values
	int initial_value = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					wx[i][j][k][l] = initial_value;
					wy[i][j][k][l] = initial_value;
				}
			}
		}
	}

	/* find the min and max for 4 columns of Nubeam scores*/
	double max_xx[4];
	double min_xx[4];
	double max_yy[4];
	double min_yy[4];
	for ( int i = 0; i < 4 ; i++ ) { //for each column, A, T, C, G
		min_xx[i] = xx[0][i]; //Assume first element in i th Column is Minimum
		max_xx[i] = xx[0][i]; //Assume first element in i th Column is Maximum
		min_yy[i] = yy[0][i]; //Assume first element in i th Column is Minimum
		max_yy[i] = yy[0][i]; //Assume first element in i th Column is Maximum
		for ( int j = 1; j < nx; j++ ){ // for each row 
			if (max_xx[i] < xx[j][i])
				max_xx[i] = xx[j][i]; // update max		
			else if (min_xx[i] > xx[j][i])
				min_xx[i] = xx[j][i]; // update min
		}
		for ( int j = 1; j < ny; j++ ){ // for each row 
			if (max_yy[i] < yy[j][i])
				max_yy[i] = yy[j][i]; // update max		
			else if (min_yy[i] > yy[j][i])
				min_yy[i] = yy[j][i]; // update min
		}
	}

	/* get the bin size for each pair of the 4 columns */
	double bin_size[4];
	double mins[4];
	for ( int i = 0; i < 4 ; i++ ) { //for each pair of column
		mins[i] = min(min_xx[i], min_yy[i]);
		bin_size[i] = ((max(max_xx[i], max_yy[i]) - mins[i])/ng) + 0.000001;
	}

	/* get the density for each 'volumn' */ 
	int bin_num1; int bin_num2; int bin_num3; int bin_num4; //bin_num track which bin should the data point in for one of the 4 columns
	for (int i = 0; i < nx; i++)
	{
		bin_num1 = int((xx[i][0] - mins[0])/bin_size[0]);
		bin_num2 = int((xx[i][1] - mins[1])/bin_size[1]);
		bin_num3 = int((xx[i][2] - mins[2])/bin_size[2]);
		bin_num4 = int((xx[i][3] - mins[3])/bin_size[3]);
		wx[bin_num1][bin_num2][bin_num3][bin_num4] += 1;
	}
	for (int i = 0; i < ny; i++)
	{
		bin_num1 = int((yy[i][0] - mins[0])/bin_size[0]);
		bin_num2 = int((yy[i][1] - mins[1])/bin_size[1]);
		bin_num3 = int((yy[i][2] - mins[2])/bin_size[2]);
		bin_num4 = int((yy[i][3] - mins[3])/bin_size[3]);
		wy[bin_num1][bin_num2][bin_num3][bin_num4] += 1;
	}
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					wx[i][j][k][l] /= (nx + initial_value * std::pow(ng, 4));
					wy[i][j][k][l] /= (ny + initial_value * std::pow(ng, 4));
				}
			}
		}
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	double temp1 = 0;
	double temp2 = 0;
	double temp3 = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					temp1 += wx[i][j][k][l] * wy[i][j][k][l];
					temp2 += wx[i][j][k][l] * wx[i][j][k][l];
					temp3 += wy[i][j][k][l] * wy[i][j][k][l];
				}
			}
		}
	}
	res = temp1 / (sqrt(temp2) * sqrt(temp3));
	res = 1 - res;
	res *= 0.5;
	return(res);
}

double cosine_2D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng) // xx and yy are arrays with 2 columns
{
	//int ng = 15; // number of bins, you can adjust it
	std::vector<std::vector<double> > wx; //2D vector
	std::vector<std::vector<double> > wy; //2D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	for (int i = 0; i < ng; i++) {
		wx[i].resize(ng);
		wy[i].resize(ng);
	}
	// Set up initial values
	int initial_value = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			wx[i][j] = initial_value;
			wy[i][j] = initial_value;
		}
	}

	/* get the density for each 'volumn' */  
	for (int i = 0; i < nx; i++) { // for each line
		int bin_num[2]; //bin_num track which bin should the data point go in for one of the 2 columns
		for ( int j = 0; j < 2; j++ ) { // for each of the 2 columns
			for ( int bx = 0; bx < ng; bx++ ) {
				if ( xx[i][j] < brk[j][bx + 1] ) {
					bin_num[j] = bx;
					break; // if an interval found, exit
				}
			}
		}
		wx[bin_num[0]][bin_num[1]] += 1;
	}
	for (int i = 0; i < ny; i++) { // for each line
		int bin_num[2]; //bin_num track which bin should the data point go in for one of the 2 columns
		for ( int j = 0; j < 2; j++ ) { // for each of the 2 columns
			for ( int by = 0; by < ng; by++ ) {
				if ( yy[i][j] < brk[j][by + 1] ) {
					bin_num[j] = by;
					break; // if an interval found, exit
				}
			}
		}
		wy[bin_num[0]][bin_num[1]] += 1;
	}

	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			wx[i][j] /= (nx + initial_value * std::pow(ng, 2));
			wy[i][j] /= (ny + initial_value * std::pow(ng, 2));
		}
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	double temp1 = 0;
	double temp2 = 0;
	double temp3 = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			temp1 += wx[i][j] * wy[i][j];
			temp2 += wx[i][j] * wx[i][j];
			temp3 += wy[i][j] * wy[i][j];
		}
	}
	res = temp1 / (sqrt(temp2) * sqrt(temp3));
	res = 1 - res;
	res *= 0.5;
	return(res);
}


double cosine_3D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk, int ng) // xx and yy are arrays with 3 columns
{
	//int ng = 15; // number of bins, you can adjust it
	std::vector<std::vector<std::vector<double> > > wx; //3D vector
	std::vector<std::vector<std::vector<double> > > wy; //3D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	for (int i = 0; i < ng; i++) {
		wx[i].resize(ng);
		wy[i].resize(ng);
		for (int j = 0; j < ng; j++) {
			wx[i][j].resize(ng);
			wy[i][j].resize(ng);
		}
	}
	// Set up initial values
	int initial_value = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				wx[i][j][k] = initial_value;
				wy[i][j][k] = initial_value;
			}
		}
	}

	/* get the density for each 'volumn' */  
	for (int i = 0; i < nx; i++) { // for each line
		int bin_num[3]; //bin_num track which bin should the data point go in for one of the 3 columns
		for ( int j = 0; j < 3; j++ ) { // for each column
			for ( int bx = 0; bx < ng; bx++ ) {
				if ( xx[i][j] < brk[j][bx + 1] ) {
					bin_num[j] = bx;
					break; // if an interval found, exit
				}
			}
		}
		wx[bin_num[0]][bin_num[1]][bin_num[2]] += 1;
	}
	for (int i = 0; i < ny; i++) { // for each line
		int bin_num[3]; //bin_num track which bin should the data point go in for one of the 3 columns
		for ( int j = 0; j < 3; j++ ) { // for each column
			for ( int by = 0; by < ng; by++ ) {
				if ( yy[i][j] < brk[j][by + 1] ) {
					bin_num[j] = by;
					break; // if an interval found, exit
				}
			}
		}
		wy[bin_num[0]][bin_num[1]][bin_num[2]] += 1;
	}

	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				//printf("%f\t", wx[i][j][k]);
				wx[i][j][k] /= (nx + initial_value * std::pow(ng, 3));
				//printf("%f\n", wy[i][j][k]);
				wy[i][j][k] /= (ny + initial_value * std::pow(ng, 3));
			}
		}
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	double temp1 = 0;
	double temp2 = 0;
	double temp3 = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				temp1 += wx[i][j][k] * wy[i][j][k];
				temp2 += wx[i][j][k] * wx[i][j][k];
				temp3 += wy[i][j][k] * wy[i][j][k];
			}
		}
	}
	res = temp1 / (sqrt(temp2) * sqrt(temp3));
	res = 1 - res;
	res *= 0.5;
	return(res);
}

double cosine_4D_quantiles(unordered_map<string, double> wx, unordered_map<string, double> wy) // wx and wy are unordered hash maps
{
	/* calculate the pair-wise distance for two samples */
	double res = 0;
	double temp1 = 0;
	double temp2 = 0;
	double temp3 = 0;
	unordered_map<string, double>:: iterator itr; 
	for (itr = wx.begin(); itr != wx.end(); itr++)
	{
		// itr works as a pointer to pair<string, double> 
		// type itr->first stores the key part  and 
		// itr->second stroes the value part
		// string key = itr->first;
		if (wy.find(itr->first) == wy.end()) {
			double x = itr->second;
			temp2 += x * x;
		}
		else {
			double x = itr->second;
			double y = wy[itr->first];
			temp1 += x * y;
			temp2 += x * x;
			temp3 += y * y;
		}
	}

	for (itr = wy.begin(); itr != wy.end(); itr++) 
	{
		// string key = itr->first;
		if (wx.find(itr->first) == wx.end()) {
			double y = itr->second;
			temp3 += y * y;
		}
	}

	res = temp1 / (sqrt(temp2) * sqrt(temp3));
	res = 1 - res;
	res *= 0.5;
	return(res);
}

double cosine_5D_quantiles(double ** xx, double ** yy, int nx, int ny, double ** brk) // xx and yy are arrays with 5 columns
{
	int ng = 5; // number of bins, you can adjust it
	std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > wx; //5D vector
	std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > wy; //5D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	for (int i = 0; i < ng; i++) {
		wx[i].resize(ng);
		wy[i].resize(ng);
		for (int j = 0; j < ng; j++) {
			wx[i][j].resize(ng);
			wy[i][j].resize(ng);
			for (int k = 0; k < ng; k++) {
				wx[i][j][k].resize(ng);
				wy[i][j][k].resize(ng);
				for (int l = 0; l < ng; l++) {
					wx[i][j][k][l].resize(ng);
					wy[i][j][k][l].resize(ng);
				}
			}
		}
	}
	// Set up initial values
	int initial_value = 10;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					for (int m = 0; m < ng; m++) {
						wx[i][j][k][l][m] = initial_value;
						wy[i][j][k][l][m] = initial_value;
					}
				}
			}
		}
	}

	/* get the density for each 'volumn' */
	for (int i = 0; i < nx; i++) { // for each line
		int bin_num[5]; //bin_num track which bin should the data point go in for one of the 5 columns
		for ( int j = 0; j < 5; j++ ) { // for each column
			for ( int bx = 0; bx < ng; bx++ ) {
				if ( xx[i][j] < brk[j][bx + 1] ) {
					bin_num[j] = bx;
					break; // if an interval found, exit
				}
			}
		}
		wx[bin_num[0]][bin_num[1]][bin_num[2]][bin_num[3]][bin_num[4]] += 1;
	}
	for (int i = 0; i < ny; i++) { // for each line
		int bin_num[5]; //bin_num track which bin should the data point go in for one of the 5 columns
		for ( int j = 0; j < 5; j++ ) { // for each column
			for ( int by = 0; by < ng; by++ ) {
				if ( yy[i][j] < brk[j][by + 1] ) {
					bin_num[j] = by;
					break; // if an interval found, exit
				}
			}
		}
		wy[bin_num[0]][bin_num[1]][bin_num[2]][bin_num[3]][bin_num[4]] += 1;
	}

	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					for (int m = 0; m < ng; m++) {
						wx[i][j][k][l][m] /= (nx + initial_value * std::pow(ng, 5));
						wy[i][j][k][l][m] /= (ny + initial_value * std::pow(ng, 5));
					}
				}
			}
		}
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	double temp1 = 0;
	double temp2 = 0;
	double temp3 = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					for (int m = 0; m < ng; m++) {
						temp1 += wx[i][j][k][l][m] * wy[i][j][k][l][m];
						temp2 += wx[i][j][k][l][m] * wx[i][j][k][l][m];
						temp3 += wy[i][j][k][l][m] * wy[i][j][k][l][m];
					}
				}
			}
		}
	}
	res = temp1 / (sqrt(temp2) * sqrt(temp3));
	res = 1 - res;
	res *= 0.5;
	return(res);
}

double canberra_4D(double ** xx, double ** yy, int nx, int ny) // xx and yy are arrays with 4 columns
{
	int ng = 15; // number of bins, you can adjust it
	std::vector<std::vector<std::vector<std::vector<double> > > > wx; //4D vector
	std::vector<std::vector<std::vector<std::vector<double> > > > wy; //4D vector
	// Set up sizes
	wx.resize(ng);
	wy.resize(ng);
	for (int i = 0; i < ng; i++) {
		wx[i].resize(ng);
		wy[i].resize(ng);
		for (int j = 0; j < ng; j++) {
			wx[i][j].resize(ng);
			wy[i][j].resize(ng);
			for (int k = 0; k < ng; k++) {
				wx[i][j][k].resize(ng);
				wy[i][j][k].resize(ng);
			}
		}
	}
	// Set up initial values
	int initial_value = 10;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					wx[i][j][k][l] = initial_value;
					wy[i][j][k][l] = initial_value;
				}
			}
		}
	}

	/* find the min and max for 4 columns of Nubeam scores*/
	double max_xx[4];
	double min_xx[4];
	double max_yy[4];
	double min_yy[4];
	for ( int i = 0; i < 4 ; i++ ) { //for each column, A, T, C, G
		min_xx[i] = xx[0][i]; //Assume first element in i th Column is Minimum
		max_xx[i] = xx[0][i]; //Assume first element in i th Column is Maximum
		min_yy[i] = yy[0][i]; //Assume first element in i th Column is Minimum
		max_yy[i] = yy[0][i]; //Assume first element in i th Column is Maximum
		for ( int j = 1; j < nx; j++ ){ // for each row 
			if (max_xx[i] < xx[j][i])
				max_xx[i] = xx[j][i]; // update max		
			else if (min_xx[i] > xx[j][i])
				min_xx[i] = xx[j][i]; // update min
		}
		for ( int j = 1; j < ny; j++ ){ // for each row 
			if (max_yy[i] < yy[j][i])
				max_yy[i] = yy[j][i]; // update max		
			else if (min_yy[i] > yy[j][i])
				min_yy[i] = yy[j][i]; // update min
		}
	}

	/* get the bin size for each pair of the 4 columns */
	double bin_size[4];
	double mins[4];
	for ( int i = 0; i < 4 ; i++ ) { //for each pair of column
		mins[i] = min(min_xx[i], min_yy[i]);
		bin_size[i] = ((max(max_xx[i], max_yy[i]) - mins[i])/ng) + 0.000001;
	}

	/* get the density for each 'volumn' */ 
	int bin_num1; int bin_num2; int bin_num3; int bin_num4; //bin_num track which bin should the data point in for one of the 4 columns
	for (int i = 0; i < nx; i++)
	{
		bin_num1 = int((xx[i][0] - mins[0])/bin_size[0]);
		bin_num2 = int((xx[i][1] - mins[1])/bin_size[1]);
		bin_num3 = int((xx[i][2] - mins[2])/bin_size[2]);
		bin_num4 = int((xx[i][3] - mins[3])/bin_size[3]);
		wx[bin_num1][bin_num2][bin_num3][bin_num4] += 1;
	}
	for (int i = 0; i < ny; i++)
	{
		bin_num1 = int((yy[i][0] - mins[0])/bin_size[0]);
		bin_num2 = int((yy[i][1] - mins[1])/bin_size[1]);
		bin_num3 = int((yy[i][2] - mins[2])/bin_size[2]);
		bin_num4 = int((yy[i][3] - mins[3])/bin_size[3]);
		wy[bin_num1][bin_num2][bin_num3][bin_num4] += 1;
	}
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					wx[i][j][k][l] /= (nx + initial_value * std::pow(ng, 4));
					wy[i][j][k][l] /= (ny + initial_value * std::pow(ng, 4));
				}
			}
		}
	}

	/* calculate the pair-wise distance for two samples */
	double res = 0;
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
			for (int k = 0; k < ng; k++) {
				for (int l = 0; l < ng; l++) {
					if (wx[i][j][k][l] + wy[i][j][k][l] != 0)
						res += abs(wx[i][j][k][l] - wy[i][j][k][l]) / (wx[i][j][k][l] + wy[i][j][k][l]);
				}
			}
		}
	}
	return(res);
}

double adtest(double * xx, double * yy, int nx, int ny)
{
	if(nx < 1000 || ny < 1000) return 0; 
	qsort((double *) xx, nx, sizeof(double), compare_doubles); 
	qsort((double *) yy, ny, sizeof(double), compare_doubles); 
	int N = nx + ny; 

//	double * xy = new double[N]; 
    double * tx = new double[N]; 
    double * ty = new double[N]; 
	for (int i =0; i < N; i++)
		tx[i] = ty[i] = 0; 

    int ix = 0; 
	int iy = 0; 
	do {
		if((xx[ix] <= yy[iy] && ix < nx) || iy == ny) { 
			tx[ix+iy] += 1; 
			if(xx[ix] == yy[iy] && iy < ny) {
				ty[ix+iy] += 1; 
				iy++;
			}
			ix++;
		}  //= cuase bias?!
		else if((xx[ix] > yy[iy] && iy < ny) || ix == nx) {
			ty[ix+iy] += 1; 
			iy++; 
		}
		else { printf("odd!!! %d  %d \n", ix,iy);} 
	} while(ix < nx || iy < ny);  


	for (int i = 1; i < N; i++)
		tx[i] += tx[i-1];  

	for (int i = 1; i < N; i++)
		ty[i] += ty[i-1];  

	double ad = 0; 
	double rmn=sqrt((double) ny / nx); 
	for (int i = 0; i < N-1; i++)
	{
		double t1 = tx[i] * rmn - ty[i] / rmn;
		ad += t1 * t1 / (i+1) / (N-i-1); 
	}                                                
	delete[] tx; 
	delete[] ty; 
	return(ad); 
}


double adtest1(double * xx, double * yy, int nx, int ny)
{
	qsort((double *) xx, nx, sizeof(double), compare_doubles); 
	qsort((double *) yy, ny, sizeof(double), compare_doubles); 
	int N = nx + ny; 

//	double * xy = new double[N]; 
    int * tx = new int[N]; 
    int * ty = new int[N]; 
	for (int i =0; i < N; i++)
		tx[i] = ty[i] = 0; 

    int ix = 0; 
	int iy = 0; 
	do {
		if((xx[ix] <= yy[iy] && ix < nx) || iy == ny) { 
//		    xy[ix+iy] = xx[ix]; 
			tx[ix+iy] ++; 
			ix++;
		}
		else if((xx[ix] > yy[iy] && iy < ny) || ix == nx) {
//		    xy[ix+iy] = yy[iy]; 
			ty[ix+iy] ++; 
			iy++; 
		}
		else { printf("odd!!! %d  %d \n", ix,iy);} 
	} while(ix < nx || iy < ny);  

	for (int i = 1; i < N; i++)
		tx[i] += tx[i-1];  

	for (int i = 1; i < N; i++)
		ty[i] += ty[i-1];  

	double ad1 = 0; 
	for (int i = 0; i < nx; i++)
	{
		double t1 = (double) (tx[i] -i); 
		ad1 += t1 * t1; 
	}                                                

	double ad2 = 0; 
	for (int i = 0; i < ny; i++)
	{
		double t2 = (double) (ty[i] -i); 
		ad2 += t2 * t2; 
	}                                                
	double ad = ad1/ny/N + ad2/nx/N; 

	delete[] tx; 
	delete[] ty; 
	return(ad); 
}

double kstest(double * xx, double * yy, int nx, int ny)
{
	qsort((double *) xx, nx, sizeof(double), compare_doubles); 
	qsort((double *) yy, ny, sizeof(double), compare_doubles); 
	int N = nx + ny; 

//	double * xy = new double[N]; 
    double * tx = new double[N]; 
    double * ty = new double[N]; 
	for (int i =0; i < N; i++)
		tx[i] = ty[i] = 0; 

    int ix = 0; 
	int iy = 0; 
	do {
		if((xx[ix] <= yy[iy] && ix < nx) || iy == ny) { 
//		    xy[ix+iy] = xx[ix]; 
			tx[ix+iy] += 1.0; 
			ix++;
		}
		else if((xx[ix] > yy[iy] && iy < ny) || ix == nx) {
//		    xy[ix+iy] = yy[iy]; 
			ty[ix+iy] += 1.0; 
			iy++; 
		}
		else { printf("odd!!! %d  %d \n", ix,iy);} 
	} while(ix < nx || iy < ny);  


	for (int i = 1; i < N; i++)
		tx[i] += tx[i-1];  

	for (int i = 1; i < N; i++)
		ty[i] += ty[i-1];  

	double ks = 0; 
	for (int i = 0; i < N; i++)
	{
		double tt = abs(tx[i]/nx - ty[i] / ny); 
		if(tt > ks) ks = tt; 
	}                                                
	delete[] tx; 
	delete[] ty; 
	return(ks); 
}

long long int filerow(string f1) 
{
		 gzFile fp = gzopen(f1.c_str(), "rb"); 
		 if(fp == NULL) {
			 printf(" can't open file, %s \n", f1.c_str()); 
			 exit(0); 
		 }
		long optimalSize = 4 * 1024 * 1024;
		char * buf = new char[optimalSize+1]; 
		long long int dim1 = 0; 
		while(1) {
			long byte_read = gzread(fp, buf, optimalSize);  
			for (int i = 0; i < byte_read; i++)
			{
				if(buf[i] =='\n') 
					dim1++; 
			}
			if(gzeof(fp)) 
				break; 
		} 
        delete[] buf; 

		 gzclose(fp); 
		 return(dim1); 
}

void regress_gc(string fin, string fout) 
{   //this reads two files twice, but significantly reduced memory footprint. 
	const int np = 3; // regressor; 
	const int ncol = 4;  //nubeam columns; 
	long long int dim1 = filerow(fin); 
	printf("File %s contains %lld lines \n", fin.c_str(), dim1); 
	fprintf(flog, "File %s contains %lld lines \n", fin.c_str(), dim1); 

	fout.append(".nogc.gz"); 

	 gzFile pout = gzopen(fout.c_str(), "wb"); 
	 if(pout == NULL) {
		 printf(" can't open file, %s \n", fout.c_str()); 
		 exit(0); 
	 }//output file ready; 
	 printf("Output file %s prepared.\n", fout.c_str());  
	 fprintf(flog, "Output file %s prepared.\n", fout.c_str());  

	 gsl_matrix * xty=gsl_matrix_alloc(np, ncol); 
	 gsl_matrix_set_zero(xty); 
	 gsl_matrix * U = gsl_matrix_alloc(np,np); 
	 gsl_matrix_set_zero(U); 

	 printf("Begin to scan the files the first time. \n");
	 fprintf(flog, "Begin to scan the files the first time. \n");

	 {
		 gzFile fp = gzopen(fin.c_str(), "rb"); 
		 if(fp == NULL) {
			 printf(" can't open file, %s \n", fin.c_str()); 
			 exit(0); 
		 }
		 // open score  and gc file; 
		 
		 char * line1 = new char[1024]; 
		 double * xx = new double[np]; 
		 double * yy = new double[ncol]; 

		 for (long long int r = 0; r < dim1; r++) {
			 gzgets(fp, line1, 1024); 
			 char * tok1 = strtok(line1, " \t"); 
			 yy[0] = atof(tok1); 
			 for (int j = 1; j < ncol; j++) 
			 {
				 tok1=strtok(NULL, " \t"); 
				 if(tok1 != NULL) 
					 yy[j] = atof(tok1); 
			 } //first four entries are scores;  
			 for (int j = 0; j < 2; j++)
			 {
				 tok1=strtok(NULL, " \t"); 
				 if(tok1 != NULL) 
					 xx[j] = log(1.0+atof(tok1)); 
			 } //remaining two entries are counts of AT and CG. 
             xx[np-1] = 1; 

			 for (int i = 0; i < np; i++)
				 for (int j = 0; j < np; j++)
				 {
					 double temp = xx[i] * xx[j] + gsl_matrix_get(U, i, j); 
					 gsl_matrix_set(U, i, j, temp); 
				 }

			 for (int i = 0; i < np; i++) 
				 for (int j = 0; j < ncol; j++)
				 {
					 double temp = xx[i] * yy[j] + gsl_matrix_get(xty, i, j); 
					 gsl_matrix_set(xty, i, j, temp); 
				 }
		 }
		 delete[] line1; 
		 delete[] xx; 
		 delete[] yy; 
		 gzclose(fp); 
	 }
	 printf("Finished scanning the files the first time.\n"); 
	 fprintf(flog, "Finished scanning the files the first time.\n"); 

	 gsl_matrix * V = gsl_matrix_alloc(np,np); 
	 gsl_vector * S = gsl_vector_alloc(np); 
	 gsl_vector * work = gsl_vector_alloc(np); 
	 for(int i = 0; i < np; i++)
		 for (int j = i+1; j < np; j++)
			 gsl_matrix_set(U, i, j, gsl_matrix_get(U, j, i)); 
	 //fill the upper; 
	 gsl_linalg_SV_decomp(U, V, S, work); 


	 gsl_matrix * mb= gsl_matrix_alloc(np,ncol); 
	 for (int c = 0; c < ncol; c++)
	 {
		 gsl_vector_view vxy = gsl_matrix_column(xty, c); 
		 gsl_vector_view gb = gsl_matrix_column(mb, c); 
		 gsl_linalg_SV_solve(U, V, S, &vxy.vector, &gb.vector); 
	 }

	 gsl_matrix_free(xty); 
	 gsl_matrix_free(U); 
	 gsl_matrix_free(V); 
	 gsl_vector_free(S); 
	 gsl_vector_free(work); 

	 printf("Obtained the regression coefficients. \n"); 
	 for(int i = 0; i < np; i++)
	 {
		 for(int c = 0; c < ncol; c++)
			printf("%lf ", gsl_matrix_get(mb, i, c)); 
	     printf("\n"); 
	 }

	 fprintf(flog, "Obtained the regression coefficients. \n"); 
	 for(int i = 0; i < np; i++)
	 {
		 for(int c = 0; c < ncol; c++)
			fprintf(flog, "%lf ", gsl_matrix_get(mb, i, c)); 
	     fprintf(flog, "\n"); 
	 }


	 printf("Begin to scan the files the second time.\n"); 
	 fprintf(flog, "Begin to scan the files the second time.\n"); 
	 {
		 gzFile fp = gzopen(fin.c_str(), "rb"); 
		 if(fp == NULL) {
			 printf(" can't open file, %s \n", fin.c_str()); 
			 exit(0); 
		 }
		 // open score  and gc contents file; 
		 
		 char * line1 = new char[1024]; 
		 double * xx = new double[np]; 
		 double * yy = new double[ncol]; 

		 for (long long int r = 0; r < dim1; r++) {
			 gzgets(fp, line1, 1024); 
			 char * tok1 = strtok(line1, " "); 
			 yy[0] = atof(tok1); 
			 for (int j = 1; j < ncol; j++) 
			 {
				 tok1=strtok(NULL, " \t"); 
				 if(tok1 != NULL) 
					 yy[j] = atof(tok1); 
			 }//read a line from fp1;  
			 for (int j = 0; j < 2; j++) 
			 {
				 tok1=strtok(NULL, " \t"); 
				 if(tok1 != NULL) 
					 xx[j] = log(1.0+atof(tok1)); 
			 }//read a line from fp1;  
             xx[np-1] = 1; 

			 for (int c = 0; c < ncol; c++)
			 {
				 double temp = 0; 
				 for (int j = 0; j < np; j++)
					 temp += xx[j] * gsl_matrix_get(mb, j, c); 
				 yy[c] -= temp; 
				 gzprintf(pout, "%lf ", yy[c]); 
			 }
			 gzprintf(pout, "\n"); 
		 }

		 delete[] line1; 
		 delete[] xx; 
		 delete[] yy; 
		 gzclose(fp); 
	 }
	 printf("Finished scanning the files the second time.\n"); 
	 fprintf(flog, "Finished scanning the files the second time.\n"); 

	 gzclose(pout); 
	 gsl_matrix_free(mb); 
}

typedef struct {
	int key; 
	double val[8]; 
} rmdup_record; 

int compare_rmdup_record(const void *a, const void *b) {
	rmdup_record * x = (rmdup_record *) a; 
	rmdup_record * y = (rmdup_record *) b; 
	double tx = x->val[0] * 1e9 + x->val[1]*1e6 + x->val[2]*1e3 + x->val[3]; 
	double ty = y->val[0] * 1e9 + y->val[1]*1e6 + y->val[2]*1e3 + y->val[3]; 
	return((int) (tx>ty)); 
}

int compare_two_record(rmdup_record *a, rmdup_record *b) {
	int ret = 0; 
	for (int i = 0; i < 4; i++)
	{
		if(fabs(a->val[i] - b->val[i]) > 1e-6) 
		{
			ret = 1; 
			break;
		}
	}
	return(ret); 
}

void rmdup_quad(string fin, string fout)
{
	int vdim = filerow(fin); 
	printf("%d \n", vdim); 
	rmdup_record * rr = new rmdup_record[vdim]; 
	{
		 gzFile fp = gzopen(fin.c_str(), "rb"); 
		 if(fp == NULL) {
			 printf(" can't open file, %s \n", fin.c_str()); 
			 exit(0); 
		 }

		 char * line = new char[1024]; 
		 for(int i = 0; i < vdim; i++) {
			 gzgets(fp, line, 1024); 
			 char * tok = strtok(line, " "); 
			 rr[i].key = i; 
			 rr[i].val[0] = atof(tok); 

			 for (int j = 1; j < 8; j++) 
			 {
				 tok=strtok(NULL, " \t"); 
				 if(tok != NULL) 
					 rr[i].val[j] = atof(tok); 
			 }  
		 }
		 gzclose(fp); 
	}
    qsort((rmdup_record *) rr, vdim, sizeof(rmdup_record), compare_rmdup_record); 

	fout.append(".rmdup.gz"); 
	gzFile fo = gzopen(fout.c_str(), "wb"); 
	if(fo == NULL) {
	 printf(" can't open file, %s \n", fout.c_str()); 
	 exit(0); 
	}
	
	int cur = 0; 
	for (int j = 0; j < 4; j++)
		gzprintf(fo, "%lf ", rr[cur].val[j]);
	for (int j = 4; j < 8; j++)
		gzprintf(fo, "%d ", (int) rr[cur].val[j]);
	gzprintf(fo, "\n"); 
    for (int next = 0; next < vdim; next++)
	{
		if(compare_two_record(&rr[next], &rr[cur]) == 1)
		{
			cur = next;
		    //then write cur into file; 
			for (int j = 0; j < 4; j++)
				gzprintf(fo, "%lf ", rr[cur].val[j]);
			for (int j = 4; j < 8; j++)
				gzprintf(fo, "%d ", (int) rr[cur].val[j]);
			gzprintf(fo, "\n"); 
		}
	}

	gzclose(fo); 
	delete[] rr; 
}

int compare_int(const void * a, const void *b)
{
	return ((int) (*(int *) a > *(int *) b));  
}

void compute_stats(vector<string> &vfn, string fout, long s, string method, unsigned long seed, int withgc, int ncol)
{
	const int nc = ncol;  
	if(withgc == 1 && strstr(vfn.at(0).c_str(), "nogc") != NULL) 
	{
		printf("### nogc file detected, toggle off withgc. \n");  
		fprintf(flog, "### nogc file detected, toggle off withgc.\n");  
		withgc = 0; 
	}
	fprintf(flog, "### cad -m %s -s %ld -R %ld -o %s ", method.c_str(), s, seed, fout.c_str());  
	for (unsigned i = 0; i < vfn.size(); i++) 
		fprintf(flog, " -i %s ", vfn.at(i).c_str());  
	fprintf(flog, " \n");  
	int nf = (int) vfn.size(); 
	double ** vdat = new double*[nf]; 
	for (int f = 0; f < nf; f++)
		vdat[f] = NULL; 
	int * vdim = new int[nf]; 
	for (int f = 0; f < nf; f++)
		vdim[f] = filerow(vfn.at(f)); 

	int size = 0; 
	double md = 1e100; 
	for (int f = 0; f < nf; f++)
		if(vdim[f] < md) md = vdim[f]; 
	if(s > 0 && s < md) 
		 size = s; 
	else  
		 size = md; 
	 for(int f = 0; f < nf; f++)
	 {
		 printf("### File %s will be sampled %d many reads out of %d.\n", vfn.at(f).c_str(), size, vdim[f]); 
		 fprintf(flog, "### File %s will be sampled %d many reads out of %d.\n", vfn.at(f).c_str(), size, vdim[f]); 
	 }

	for (int f = 0; f < nf; f++)
	{
	    int * index = new int[vdim[f]]; 
		for (int i = 0; i < vdim[f]; i++)
			index[i] = i; 
		int * chosen = new int[size]; 
		gsl_ran_choose(gsl_r, chosen, size, index, vdim[f], sizeof(int)); 
		qsort((int *) chosen, size, sizeof(int), compare_int); 
//		for(int i = 0; i < size; i++)
//			printf("%d ",chosen[i]); 
//		printf("\n"); 

		vdat[f] = new double[size * nc]; 

		 gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
		 if(fp == NULL) {
			 printf(" can't open file, %s \n", vfn.at(f).c_str()); 
			 exit(0); 
		 }

		 int i = 0; 
		 char * line = new char[1024]; 
		 for (int r = 0; r < vdim[f]; r++){
			 gzgets(fp, line, 1024); 
			 if(r == chosen[i]) {
				 if(nc == 4) 
				 	sscanf(line, "%lf %lf %lf %lf", &vdat[f][i], &vdat[f][size+i], &vdat[f][size*2+i], &vdat[f][size*3+i]);  
				 else if(nc == 2) 
				 	sscanf(line, "%lf %lf ", &vdat[f][i], &vdat[f][size+i]);  
				 else if(nc == 1) 
				 	sscanf(line, "%lf ", &vdat[f][i]);  
				 else {
					 char * tok = strtok(line, " "); 
					 vdat[f][i] = atof(tok); 
					 for (int j = 1; j < nc; j++) 
					 {
						 tok=strtok(NULL, " \t"); 
						 if(tok != NULL) 
							 vdat[f][size*j+i] = atof(tok); 
					 } //first four entries are scores;  
				 }
				 i++; 
			 }
			 if(i >= size) break; 
		 }
		 printf("Read %d random entries in file %s.\n" , i, vfn.at(f).c_str()); 
		 fprintf(flog, "### Read %d random entries in file %s.\n" , i, vfn.at(f).c_str()); 
		 delete[] line; 
		 delete[] chosen; 
		 delete[] index; 
		 gzclose(fp); 
	}
	
	double * res = new double[nf * nf]; 

	memset(res, 0, sizeof(double)*nf*nf); 
	for (int i = 0; i < nf-1; i++)
		for (int j = i+1; j < nf; j++)
		{
			double stat = 0; 
			double test = 0; 
			for (int c = 0; c < nc; c++)
			{
				if(method.compare("ad") == 0) 
					test = adtest(&vdat[i][c*size],&vdat[j][c*size], size,size);
				else if(method.compare("ks") == 0) 
					test = kstest(&vdat[i][c*size],&vdat[j][c*size], size,size);
				else if(method.compare("h2") == 0)
					test = hellinger(&vdat[i][c*size],&vdat[j][c*size], size,size);
				else if(method.compare("cos") == 0)
					test = cosine(&vdat[i][c*size],&vdat[j][c*size], size,size);
//				printf("%s\t %d\t %d\t %d\t %lf \n",method.c_str(), c, size, size, test[c]); 
//				fprintf(flog, "### %s\t %d\t %d\t %d\t %lf \n",method.c_str(), c, size, size, test[c]); 
				stat += test; 
			}
			stat /= nc; 
			//harmonic mean; 

//			{
//				if(method.compare("ad") == 0) 
//					stat = adtest(&vdat[i][0],&vdat[j][0], nc*size,nc*size);
//				else if(method.compare("ks") == 0) 
//					stat = kstest(&vdat[i][0],&vdat[j][0], nc*size,nc*size);
//				else 
//					stat = hellinger(&vdat[i][0],&vdat[j][0], nc*size,nc*size);
//			}

			printf("%s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
			fprintf(flog, "### %s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
			res[i*nf+j] = res[j*nf+i] = stat; 
		}

	for (int i = 0; i < nf; i++)
	{
		for (int j = 0; j < nf; j++)
			fprintf(flog, "%lf\t",  res[i*nf+j]); 
		fprintf(flog, "\n"); 
	}

	delete[] vdim; 
	for(int f = 0; f < nf; f++)
		delete[] vdat[f]; 
	delete[] vdat; 
	delete[] res;; 
}

// void compute_stats_unequal_datasize(vector<string> &vfn, string fout, string method, int withgc, int ncol, int n_bin, string bin_file)
// {
// 	const int nc = ncol;
// 	const int nb = n_bin;
// 	const double pseudo_count = 0;
// 	if(withgc == 1 && strstr(vfn.at(0).c_str(), "nogc") != NULL) 
// 	{
// 		printf("### nogc file detected, toggle off withgc. \n");  
// 		fprintf(flog, "### nogc file detected, toggle off withgc.\n");  
// 		withgc = 0; 
// 	}
// 	fprintf(flog, "### cad -m %s -o %s", method.c_str(), fout.c_str());  
// 	for (unsigned i = 0; i < vfn.size(); i++) 
// 		fprintf(flog, " -i %s ", vfn.at(i).c_str());  
// 	fprintf(flog, " \n");  
// 	int nf = (int) vfn.size(); // how many files are there
	
// 	int * vdim = new int[nf]; // the number of lines for each file
// 	for (int f = 0; f < nf; f++)
// 		vdim[f] = filerow(vfn.at(f)); 

// 	if(nc == 5 && method.find("_5d") != string::npos) { // if 5-dimensional joint distribution would be used
// 		/* read the data for the first time */
// 		long sum_all_size = 0;
// 		for (int f = 0; f < nf; f++) {
// 			sum_all_size += vdim[f];
// 		}
// 		double ** vdat1 = new double * [5];
// 		for (int f = 0; f < 5; f++)
// 			vdat1[f] = new double [sum_all_size];
// 		int i = 0;
// 		for (int f = 0; f < nf; f++) {
// 			int size = vdim[f];	     
// 			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 			if(fp == NULL) {
// 				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 				exit(0); 
// 			}
// 			char * line = new char[1024];
// 			for (int r = 0; r < size; r++){
// 				gzgets(fp, line, 1024); 
// 				sscanf(line, "%lf %lf %lf %lf %lf", &vdat1[0][i], &vdat1[1][i], &vdat1[2][i], &vdat1[3][i], &vdat1[4][i]);
// 				i++;
// 			}
// 			delete[] line; 
// 			gzclose(fp); 
// 		}

// 		/* calculate bins according to quantiles */
// 		double ** brk = new double*[5];
// 		for (int i = 0; i < 5; i++) {
// 			brk[i] = new double[5 + 1];
// 		}
// 		brk = select_quantiles_all(vdat1, sum_all_size, 5, 5);
// 		for ( int j = 0; j < 5; j++ )
// 			delete[] vdat1[j];
// 		delete[] vdat1;

// 		/* read the data for the second time */
// 		double *** vdat = new double ** [nf]; // store data
// 		for (int f = 0; f < nf; f++)
// 			vdat[f] = NULL; 
// 		for (int f = 0; f < nf; f++) {
// 			int size = vdim[f];	     
// 			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 			if(fp == NULL) {
// 				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 				exit(0); 
// 			}
// 			vdat[f] = new double * [size];
// 			char * line = new char[1024];
// 			int i = 0; 
// 			for (int r = 0; r < size; r++){
// 				gzgets(fp, line, 1024); 
// 				vdat[f][r] = new double [5];
// 				sscanf(line, "%lf %lf %lf %lf %lf", &vdat[f][r][0], &vdat[f][r][1], &vdat[f][r][2], &vdat[f][r][3], &vdat[f][r][4]);
// 				i++;
// 			}
// 			printf("Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			delete[] line; 
// 			gzclose(fp); 
// 		}

// 		double * res = new double[nf * nf]; 

// 		memset(res, 0, sizeof(double)*nf*nf); 
// 		for (int i = 0; i < nf-1; i++) {
// 			for (int j = i+1; j < nf; j++) {
// 				double stat = 0;
// 				if (method.compare("h2_5d") == 0)
// 					stat = hellinger_5D_quantiles(vdat[i], vdat[j], vdim[i], vdim[j], brk);
// 				else if (method.compare("cos_5d") == 0)
// 					stat = cosine_5D_quantiles(vdat[i], vdat[j], vdim[i], vdim[j], brk);
// 				printf("%s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				fprintf(flog, "### %s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				res[i*nf+j] = res[j*nf+i] = stat;
// 			}
// 		}

// 		for (int i = 0; i < nf; i++)
// 		{
// 			for (int j = 0; j < nf; j++) {
// 				fprintf(flog, "%.10f\t",  res[i*nf+j]);
// 			}
// 			fprintf(flog, "\n");
// 		}

// 		for(int f = 0; f < nf; f++) {
// 			int size = vdim[f];
// 			for (int r = 0; r < size; r++){
// 				delete[] vdat[f][r];
// 			}
// 			delete[] vdat[f]; 
// 		}
// 		delete[] vdat; 
// 		delete[] vdim;
// 		delete[] res;
// 	}
// 	else if(nc == 4 && method.find("_4d") != string::npos) { // if 4-dimensional joint distribution would be used
// 		/* calculate bins according to quantiles */
// 		double * bins_1st = new double[nb + 1];

// 		double ** bins_2nd_all = new double * [nb];
// 		for (int f = 0; f < nb; f++) {
// 			bins_2nd_all[f] = new double[nb + 1];
// 		}

// 		double *** bins_3rd_all = new double ** [nb];
// 		for (int f = 0; f < nb; f++) {
// 			bins_3rd_all[f] = new double * [nb];
// 			for (int g = 0; g < nb; g++)
// 				bins_3rd_all[f][g] = new double[nb + 1];
// 		}

// 		double **** bins_4th_all = new double *** [nb];
// 		for (int f = 0; f < nb; f++) {
// 			bins_4th_all[f] = new double ** [nb];
// 			for (int g = 0; g < nb; g++) {
// 				bins_4th_all[f][g] = new double * [nb];
// 				for (int h = 0; h < nb; h++)
// 					bins_4th_all[f][g][h] = new double[nb + 1];
// 			}
// 		}

// 		gzFile myfile = gzopen(bin_file.c_str(), "rb");
// 		if (myfile == NULL) {
// 			printf(" No bin file detected, calculate bin partitions\n");
// 			/* read the data for the first time */
// 			long sum_all_size = 0;
// 			for (int f = 0; f < nf; f++) {
// 				sum_all_size += vdim[f];
// 			}
// 			double ** vdat1 = new double * [4];
// 			for (int f = 0; f < 4; f++)
// 				vdat1[f] = new double [sum_all_size];
// 			int i = 0;
// 			for (int f = 0; f < nf; f++) {
// 				int size = vdim[f];	     
// 				gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 				if(fp == NULL) {
// 					printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 					exit(0); 
// 				}
// 				char * line = new char[1024];
// 				for (int r = 0; r < size; r++){
// 					gzgets(fp, line, 1024); 
// 					sscanf(line, "%lf %lf %lf %lf", &vdat1[0][i], &vdat1[1][i], &vdat1[2][i], &vdat1[3][i]);
// 					i++;
// 				}
// 				delete[] line; 
// 				gzclose(fp); 
// 			}

// 			bins_1st = select_quantiles(vdat1[0], sum_all_size, nb);
// 			for (int i = 0; i < nb + 1; i++)
// 				printf("%f\n", bins_1st[i]);

// 			for (int i = 0; i < nb; i++) {
// 				double ** vdat1_2 = new double * [3];
// 				for (int f = 0; f < 3; f++)
// 					vdat1_2[f] = new double [sum_all_size]; // excessively large size
// 				int count1 = 0;
// 				for (int l = 0; l < sum_all_size; l++) {
// 					if (vdat1[0][l] < bins_1st[i + 1] && vdat1[0][l] >= bins_1st[i]) {
// 						vdat1_2[0][count1] = vdat1[1][l]; 
// 						vdat1_2[1][count1] = vdat1[2][l]; 
// 						vdat1_2[2][count1] = vdat1[3][l];
// 						count1 ++; 
// 					}
// 				}
// 				printf("%d\n", count1);
// 				double ** vdat1_2_copy = new double * [3];
// 				for (int f = 0; f < 3; f++)
// 					vdat1_2_copy[f] = new double [count1]; // correct size
// 				for (int l = 0; l < count1; l++) {
// 					vdat1_2_copy[0][l] = vdat1_2[0][l]; 
// 					vdat1_2_copy[1][l] = vdat1_2[1][l];
// 					vdat1_2_copy[2][l] = vdat1_2[2][l];
// 				}
// 				for (int f = 0; f < 3; f++)
// 					delete[] vdat1_2[f];
// 				delete[] vdat1_2;
// 				printf("i%d, finish copy and delete\n", i);
// 				bins_2nd_all[i] = select_quantiles(vdat1_2_copy[0], count1, nb);
// 				for (int hehe = 0; hehe < nb + 1; hehe++)
// 					printf("%f\n", bins_2nd_all[i][hehe]);
// 				printf("i%d, finish select quantiles\n", i);

// 				for (int j = 0; j < nb; j++) {
// 					double ** vdat1_3 = new double * [2];
// 					for (int f = 0; f < 2; f++)
// 						vdat1_3[f] = new double [count1]; // excessively large size
// 					int count2 = 0;
// 					for (int l = 0; l < count1; l++) {
// 						if (vdat1_2_copy[0][l] < bins_2nd_all[i][j + 1] && vdat1_2_copy[0][l] >= bins_2nd_all[i][j]) {
// 							vdat1_3[0][count2] = vdat1_2_copy[1][l]; 
// 							vdat1_3[1][count2] = vdat1_2_copy[2][l]; 
// 							count2 ++;
// 						}
// 					}
// 					printf("%d\n", count2);
// 					double ** vdat1_3_copy = new double * [2];
// 					for (int f = 0; f < 2; f++)
// 						vdat1_3_copy[f] = new double [count2]; // correct size
// 					for (int l = 0; l < count2; l++) {
// 						vdat1_3_copy[0][l] = vdat1_3[0][l];
// 						vdat1_3_copy[1][l] = vdat1_3[1][l];
// 					}
// 					for (int f = 0; f < 2; f++)
// 						delete[] vdat1_3[f];
// 					delete[] vdat1_3;
// 					printf("j%d, finish copy and delete\n", j);
// 					bins_3rd_all[i][j] = select_quantiles(vdat1_3_copy[0], count2, nb);
// 					for (int hehe = 0; hehe < nb + 1; hehe++)
// 						printf("%f\n", bins_3rd_all[i][j][hehe]);
// 					printf("j%d, finish select quantiles\n", j);

// 					for (int k = 0; k < nb; k++) {
// 						double * vdat1_4 = new double [count2]; // excessively large size
// 						int count3 = 0;
// 						for (int l = 0; l < count2; l++) {
// 							if (vdat1_3_copy[0][l] < bins_3rd_all[i][j][k + 1] && vdat1_3_copy[0][l] >= bins_3rd_all[i][j][k]) {
// 								vdat1_4[count3] = vdat1_3_copy[1][l];
// 								count3++;
// 							}
// 						}
// 						printf("%d\n", count3);
// 						double * vdat1_4_copy = new double [count3]; // correct size
// 						for (int l = 0; l < count3; l++) {
// 							vdat1_4_copy[l] = vdat1_4[l];
// 						}
// 						delete[] vdat1_4;
// 						printf("k%d, finish copy and delete\n", k);
// 						bins_4th_all[i][j][k] = select_quantiles(vdat1_4_copy, count3, nb);
// 						for (int hehe = 0; hehe < nb + 1; hehe++)
// 							printf("%f\n", bins_4th_all[i][j][k][hehe]);
// 						printf("k%d, finish select quantiles\n", k);
// 						delete[] vdat1_4_copy;
// 					}
// 					for (int f = 0; f < 2; f++)
// 						delete[] vdat1_3_copy[f];
// 					delete[] vdat1_3_copy;
// 				}
// 				for (int f = 0; f < 3; f++)
// 					delete[] vdat1_2_copy[f];
// 				delete[] vdat1_2_copy;
// 			}
// 			for ( int j = 0; j < 4; j++ )
// 				delete[] vdat1[j];
// 			delete[] vdat1;

// 			/* write bin files for future use */
// 			FILE * fbin = fopen("bin.txt", "w");
// 			printf("Bin file name bin.txt\n");
// 			for (int i = 0; i < nb + 1; i++) {
// 				fprintf(fbin, "%.10f ", bins_1st[i]);
// 			}
// 			fprintf(fbin, "\n");

// 			for (int i = 0; i < nb; i++) {
// 				for (int j = 0; j < nb + 1; j++) {
// 					fprintf(fbin, "%.10f ", bins_2nd_all[i][j]);
// 				}
// 				fprintf(fbin, "\n");
// 			}

// 			for (int i = 0; i < nb; i++) {
// 				for (int j = 0; j < nb; j++) {
// 					for (int k = 0; k < nb + 1; k++) {
// 						fprintf(fbin, "%.10f ", bins_3rd_all[i][j][k]);
// 					}
// 					fprintf(fbin, "\n");
// 				}
// 			}

// 			for (int i = 0; i < nb; i++) {
// 				for (int j = 0; j < nb; j++) {
// 					for (int k = 0; k < nb; k++) {
// 						for (int l = 0; l < nb + 1; l++) {
// 							fprintf(fbin, "%.10f ", bins_4th_all[i][j][k][l]);
// 						}
// 						fprintf(fbin, "\n");
// 					}
// 				}
// 			}
// 			fclose(fbin);
// 		}
// 		else {
// 			char * line = new char[1024000];
// 			for (int i = 0; i < (std::pow(nb, 3) + std::pow(nb, 2) + std::pow(nb, 1) + std::pow(nb, 0)); i++) {
// 				// printf("%d\n", i);
// 				gzgets(myfile, line, 1024000);
// 				char * tok = strtok(line, " ");
// 				if (i < std::pow(nb, 0)) {
// 					int j = 0;
// 			  		while (tok != NULL && j < nb + 1)
// 			  		{
// 			    		bins_1st[j] = atof(tok);
// 			    		tok = strtok(NULL, " ");
// 			    		j++;
// 			  		}
// 				}
// 				else if (i < (std::pow(nb, 1) + std::pow(nb, 0)))
// 				{
// 					int j = 0;
// 			  		while (tok != NULL && j < nb + 1)
// 			  		{
// 			    		bins_2nd_all[i - 1][j] = atof(tok);
// 			    		tok = strtok(NULL, " ");
// 			    		j++;
// 			  		}
// 			  		// printf("finish line %d\n", i);
// 				}
// 				else if (i < (std::pow(nb, 2) + std::pow(nb, 1) + std::pow(nb, 0)))
// 				{
// 					int j = 0;
// 			  		while (tok != NULL && j < nb + 1)
// 			  		{
// 			    		bins_3rd_all[static_cast<int>(ceil(i / std::pow(nb, 1)) - 2)][static_cast<int>(i - (ceil(i / std::pow(nb, 1)) - 1) * std::pow(nb, 1) - 1)][j] = atof(tok);
// 			    		tok = strtok(NULL, " ");
// 			    		j++;
// 			  		}
// 			  		// printf("finish line %d\n", i);
// 				}
// 				else
// 				{
// 					int l = nb - 1;
// 					if (i % nb != 0) {
// 						l = i % nb - 1;
// 					}
// 					// printf("%d\n", static_cast<int>(ceil((i - 10) / std::pow(nb, 2)) - 2));
// 					// printf("%d\n", static_cast<int>(ceil((i - (ceil((i - 10) / std::pow(nb, 2)) - 1) * std::pow(nb, 2)) / std::pow(nb, 1)) - 2));
// 					// printf("%d\n", l);
// 					int j = 0;
// 			  		while (tok != NULL && j < nb + 1)
// 			  		{
// 			    		bins_4th_all[static_cast<int>(ceil((i - 10) / std::pow(nb, 2)) - 2)][static_cast<int>(ceil((i - (ceil((i - 10) / std::pow(nb, 2)) - 1) * std::pow(nb, 2)) / std::pow(nb, 1)) - 2)][l][j] = atof(tok);
// 			    		tok = strtok(NULL, " ");
// 			    		j++;
// 			  		}
// 			  		// printf("finish line %d\n", i);
// 				}
// 			}
// 			delete[] line;
// 			gzclose(myfile);
// 			// for (int i = 0; i < 4; i++) {
// 			// 	brk[i][0] = -1 * INFINITY;
// 			// 	brk[i][nb] = INFINITY;
// 			// }
// 		}

// 		/* read the data for the second time */
// 		unordered_map<string, double> * vdat = new unordered_map<string, double>[nf];
// 		long int pseudo_size = 1000000 + pseudo_count * std::pow(nb, 4);
// 		double pseudo_freq = pseudo_count / pseudo_size;
// 		// printf("%.15lf\n", pseudo_freq);
// 		for (int f = 0; f < nf; f++) {
// 			int size = vdim[f];	     
// 			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 			if(fp == NULL) {
// 				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 				exit(0); 
// 			}
// 			char * line = new char[1024];
// 			double * xx = new double [4];
// 			int i = 0; 
// 			for (int r = 0; r < size; r++){
// 				gzgets(fp, line, 1024); 
// 				sscanf(line, "%lf %lf %lf %lf", &xx[0], &xx[1], &xx[2], &xx[3]);
// 				string bin_num; //bin_num track which bin should the data point go in for each of the 4 columns, like "3,76,24,18"
// 				int cursor1 = quickselect(bins_1st, nb, xx[0]);
// 				bin_num += to_string((long long) cursor1);
// 				bin_num += ",";
// 				int cursor2 = quickselect(bins_2nd_all[cursor1], nb, xx[1]);
// 				bin_num += to_string((long long) cursor2);
// 				bin_num += ",";
// 				int cursor3 = quickselect(bins_3rd_all[cursor1][cursor2], nb, xx[2]);
// 				bin_num += to_string((long long) cursor3);
// 				bin_num += ",";
// 				int cursor4 = quickselect(bins_4th_all[cursor1][cursor2][cursor3], nb, xx[3]);
// 				bin_num += to_string((long long) cursor4);
// 				vdat[f][bin_num] += 1;
// 				i++;
// 			}
// 			// printf("%zu\n", vdat[f].size());
// 			/* get the density for each 'volumn' */
// 			unordered_map<string, double>:: iterator itr;
// 			for (itr = vdat[f].begin(); itr != vdat[f].end(); itr++) { 
// 				// itr works as a pointer to pair<string, double> 
// 				// type itr->first stores the key part  and 
// 				// itr->second stroes the value part
// 				string key = itr->first;
// 				// cout << key << "\t" << itr->second << "\n";
// 				vdat[f][key] = ((1000000 * vdat[f][key] / size) + pseudo_count) / pseudo_size;
// 			}

// 			printf("Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			delete[] line; 
// 			delete[] xx;
// 			gzclose(fp); 
// 		}

// 		double * res = new double[nf * nf]; 

// 		memset(res, 0, sizeof(double)*nf*nf); 
// 		for (int i = 0; i < nf-1; i++) {
// 			for (int j = i+1; j < nf; j++) {
// 				double stat = 0;
// 				if (method.compare("h2_4d") == 0)
// 					stat = hellinger_4D_quantiles(vdat[i], vdat[j], pseudo_freq);
// 				else if (method.compare("cos_4d") == 0)
// 					stat = cosine_4D_quantiles(vdat[i], vdat[j]);
// 				else if (method.compare("ham_4d") == 0)
// 					stat = hamming_4D_quantiles(vdat[i], vdat[j]);
// 				else if (method.compare("l1_4d") == 0)
// 					stat = l1_4D_quantiles(vdat[i], vdat[j]);
// 				printf("%s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				fprintf(flog, "### %s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				res[i*nf+j] = res[j*nf+i] = stat;
// 			}
// 		}

// 		for (int i = 0; i < nf; i++)
// 		{
// 			for (int j = 0; j < nf; j++) {
// 				fprintf(flog, "%.10f\t",  res[i*nf+j]);
// 			}
// 			fprintf(flog, "\n");
// 		}

// 		// for(int f = 0; f < nf; f++) {
// 		// 	delete[] vdat[f]; 
// 		// }
// 		delete[] vdat; 
// 		delete[] vdim;
// 		delete[] res;
// 		delete[] bins_1st;

// 		for (int i = 0; i < nb; i++) {
// 			delete[] bins_2nd_all[i];
// 		}
// 		delete[] bins_2nd_all;

// 		for (int i = 0; i < nb; i++) {
// 			for (int j = 0; j < nb; j++) {
// 				delete[] bins_3rd_all[i][j];
// 			}
// 			delete[] bins_3rd_all[i];
// 		}
// 		delete[] bins_3rd_all;

// 		for (int i = 0; i < nb; i++) {
// 			for (int j = 0; j < nb; j++) {
// 				for (int k = 0; k < nb; k++) {
// 					delete[] bins_4th_all[i][j][k];
// 				}
// 				delete[] bins_4th_all[i][j];
// 			}
// 			delete[] bins_4th_all[i];
// 		}
// 		delete[] bins_4th_all;
// 	}
// 	else if(nc == 3 && method.find("_3d") != string::npos) { // if 3-dimensional joint distribution would be used
// 		/* read the data for the first time */
// 		long sum_all_size = 0;
// 		for (int f = 0; f < nf; f++) {
// 			sum_all_size += vdim[f];
// 		}
// 		double ** vdat1 = new double * [3];
// 		for (int f = 0; f < 3; f++)
// 			vdat1[f] = new double [sum_all_size];
// 		int i = 0;
// 		for (int f = 0; f < nf; f++) {
// 			int size = vdim[f];	     
// 			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 			if(fp == NULL) {
// 				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 				exit(0); 
// 			}
// 			char * line = new char[1024];
// 			for (int r = 0; r < size; r++){
// 				gzgets(fp, line, 1024); 
// 				sscanf(line, "%lf %lf %lf", &vdat1[0][i], &vdat1[1][i], &vdat1[2][i]);
// 				i++;
// 			}
// 			delete[] line; 
// 			gzclose(fp); 
// 		}

// 		/* calculate bins according to quantiles */
// 		double ** brk = new double*[3];
// 		for (int i = 0; i < 3; i++) {
// 			brk[i] = new double[nb + 1];
// 		}
// 		brk = select_quantiles_all(vdat1, sum_all_size, nb, 3);
// 		for ( int j = 0; j < 3; j++ )
// 			delete[] vdat1[j];
// 		delete[] vdat1;

// 		/* read the data for the second time */
// 		double *** vdat = new double ** [nf]; // store data
// 		for (int f = 0; f < nf; f++)
// 			vdat[f] = NULL; 
// 		for (int f = 0; f < nf; f++) {
// 			int size = vdim[f];	     
// 			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 			if(fp == NULL) {
// 				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 				exit(0); 
// 			}
// 			vdat[f] = new double * [size];
// 			char * line = new char[1024];
// 			int i = 0; 
// 			for (int r = 0; r < size; r++){
// 				gzgets(fp, line, 1024); 
// 				vdat[f][r] = new double [3];
// 				sscanf(line, "%lf %lf %lf", &vdat[f][r][0], &vdat[f][r][1], &vdat[f][r][2]);
// 				i++;
// 			}
// 			printf("Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			delete[] line; 
// 			gzclose(fp); 
// 		}

// 		double * res = new double[nf * nf]; 

// 		memset(res, 0, sizeof(double)*nf*nf); 
// 		for (int i = 0; i < nf-1; i++) {
// 			for (int j = i+1; j < nf; j++) {
// 				double stat = 0;
// 				if (method.compare("h2_3d") == 0)
// 					stat = hellinger_3D_quantiles(vdat[i], vdat[j], vdim[i], vdim[j], brk, nb);
// 				else if (method.compare("cos_3d") == 0)
// 					stat = cosine_3D_quantiles(vdat[i], vdat[j], vdim[i], vdim[j], brk, nb);
// 				printf("%s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				fprintf(flog, "### %s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				res[i*nf+j] = res[j*nf+i] = stat;
// 			}
// 		}

// 		for (int i = 0; i < nf; i++)
// 		{
// 			for (int j = 0; j < nf; j++) {
// 				fprintf(flog, "%.10f\t",  res[i*nf+j]);
// 			}
// 			fprintf(flog, "\n");
// 		}

// 		for(int f = 0; f < nf; f++) {
// 			int size = vdim[f];
// 			for (int r = 0; r < size; r++){
// 				delete[] vdat[f][r];
// 			}
// 			delete[] vdat[f]; 
// 		}
// 		delete[] vdat; 
// 		delete[] vdim;
// 		delete[] res;
// 	}
// 	else if(nc == 2 && method.find("_2d") != string::npos) { // if 2-dimensional joint distribution would be used
// 		/* read the data for the first time */
// 		long sum_all_size = 0;
// 		for (int f = 0; f < nf; f++) {
// 			sum_all_size += vdim[f];
// 		}
// 		double ** vdat1 = new double * [2];
// 		for (int f = 0; f < 2; f++)
// 			vdat1[f] = new double [sum_all_size];
// 		int i = 0;
// 		for (int f = 0; f < nf; f++) {
// 			int size = vdim[f];	     
// 			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 			if(fp == NULL) {
// 				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 				exit(0); 
// 			}
// 			char * line = new char[1024];
// 			for (int r = 0; r < size; r++){
// 				gzgets(fp, line, 1024); 
// 				sscanf(line, "%lf %lf", &vdat1[0][i], &vdat1[1][i]);
// 				i++;
// 			}
// 			delete[] line; 
// 			gzclose(fp); 
// 		}

// 		/* calculate bins according to quantiles */
// 		double ** brk = new double*[2];
// 		for (int i = 0; i < 2; i++) {
// 			brk[i] = new double[nb + 1];
// 		}
// 		brk = select_quantiles_all(vdat1, sum_all_size, nb, 2);
// 		for ( int j = 0; j < 2; j++ )
// 			delete[] vdat1[j];
// 		delete[] vdat1;

// 		/* read the data for the second time */
// 		double *** vdat = new double ** [nf]; // store data
// 		for (int f = 0; f < nf; f++)
// 			vdat[f] = NULL; 
// 		for (int f = 0; f < nf; f++) {
// 			int size = vdim[f];	     
// 			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 			if(fp == NULL) {
// 				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 				exit(0); 
// 			}
// 			vdat[f] = new double * [size];
// 			char * line = new char[1024];
// 			int i = 0; 
// 			for (int r = 0; r < size; r++){
// 				gzgets(fp, line, 1024); 
// 				vdat[f][r] = new double [2];
// 				sscanf(line, "%lf %lf", &vdat[f][r][0], &vdat[f][r][1]);
// 				i++;
// 			}
// 			printf("Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			delete[] line; 
// 			gzclose(fp); 
// 		}

// 		double * res = new double[nf * nf]; 

// 		memset(res, 0, sizeof(double)*nf*nf); 
// 		for (int i = 0; i < nf-1; i++) {
// 			for (int j = i+1; j < nf; j++) {
// 				double stat = 0;
// 				if (method.compare("h2_2d") == 0)
// 					stat = hellinger_2D_quantiles(vdat[i], vdat[j], vdim[i], vdim[j], brk, nb);
// 				else if (method.compare("cos_2d") == 0)
// 					stat = cosine_2D_quantiles(vdat[i], vdat[j], vdim[i], vdim[j], brk, nb);
// 				printf("%s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				fprintf(flog, "### %s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				res[i*nf+j] = res[j*nf+i] = stat;
// 			}
// 		}

// 		for (int i = 0; i < nf; i++)
// 		{
// 			for (int j = 0; j < nf; j++) {
// 				fprintf(flog, "%.10f\t",  res[i*nf+j]);
// 			}
// 			fprintf(flog, "\n");
// 		}

// 		for(int f = 0; f < nf; f++) {
// 			int size = vdim[f];
// 			for (int r = 0; r < size; r++){
// 				delete[] vdat[f][r];
// 			}
// 			delete[] vdat[f]; 
// 		}
// 		delete[] vdat; 
// 		delete[] vdim;
// 		delete[] res;
// 		for ( int j = 0; j < 2; j++ )
// 			delete[] brk[j];
// 		delete[] brk;
// 	}
// 	else if(nc == 1 && method.find("_1d") != string::npos) { // if 1-dimensional distribution would be used
// 		double ** brk = new double*[1];
// 		for (int i = 0; i < 1; i++) {
// 			brk[i] = new double[nb + 1];
// 		}
// 		gzFile myfile = gzopen(bin_file.c_str(), "rb");
// 		if (myfile == NULL) {
// 			printf(" No bin file detected, calculate bin partitions\n");
// 			/* read the data for the first time */
// 			long sum_all_size = 0;
// 			for (int f = 0; f < nf; f++) {
// 				sum_all_size += vdim[f];
// 			}
// 			double ** vdat1 = new double * [1];
// 			for (int f = 0; f < 1; f++)
// 				vdat1[f] = new double [sum_all_size];
// 			int i = 0;
// 			for (int f = 0; f < nf; f++) {
// 				int size = vdim[f];	     
// 				gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 				if(fp == NULL) {
// 					printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 					exit(0); 
// 				}
// 				char * line = new char[1024];
// 				for (int r = 0; r < size; r++){
// 					gzgets(fp, line, 1024); 
// 					sscanf(line, "%lf", &vdat1[0][i]);
// 					i++;
// 				}
// 				delete[] line; 
// 				gzclose(fp); 
// 			}

// 			/* calculate bins according to quantiles */
// 			brk = select_quantiles_all(vdat1, sum_all_size, nb, 1);
// 			for ( int j = 0; j < 1; j++ )
// 				delete[] vdat1[j];
// 			delete[] vdat1;

// 			/* write a bin file for future use */
// 			FILE * fbin = fopen("bin.txt", "w");
// 			printf("Bin file name bin.txt\n");
// 			for (int i = 0; i < 1; i++) {
// 				for (int j = 0; j < nb + 1; j++) {
// 					fprintf(fbin, "%.10f ", brk[i][j]);
// 				}
// 				fprintf(fbin, "\n");
// 			}
// 			fclose(fbin);
// 		}
// 		else {
// 			char * line = new char[1024000];
// 			for (int i = 0; i < 1; i++) {
// 				gzgets(myfile, line, 1024000);
// 				char * tok = strtok(line, " ");
// 				int j = 0;
// 		  		while (tok != NULL && j < nb + 1)
// 		  		{
// 		    		brk[i][j] = atof(tok);
// 		    		tok = strtok(NULL, " ");
// 		    		j++;
// 		  		}
// 			}
// 			delete[] line;
// 			gzclose(myfile);
// 			for (int i = 0; i < 1; i++) {
// 				brk[i][0] = -1 * INFINITY;
// 				brk[i][nb] = INFINITY;
// 			}
// 		}

// 		/* read the data for the second time */
// 		unordered_map<string, double> * vdat = new unordered_map<string, double>[nf];
// 		for (int f = 0; f < nf; f++) {
// 			int size = vdim[f];	     
// 			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 			if(fp == NULL) {
// 				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 				exit(0); 
// 			}
// 			char * line = new char[1024];
// 			double * xx = new double [1];
// 			int i = 0; 
// 			for (int r = 0; r < size; r++){
// 				gzgets(fp, line, 1024); 
// 				sscanf(line, "%lf", &xx[0]);
// 				string bin_num; //bin_num track which bin should the data point go in for the column, like "3"
// 				for ( int j = 0; j < 1; j++ ) { // for each column of only 1 column
// 					if (j != 0) {
// 	        			bin_num += ",";
// 	    			}
// 	    			bin_num += to_string((long long) quickselect(brk[j], nb, xx[j]));
// 				}
// 				vdat[f][bin_num] += 1;
// 				i++;
// 			}
// 			printf("%zu\n", vdat[f].size());
// 			/* get the density for each 'volumn' */
// 			unordered_map<string, double>:: iterator itr;
// 			for (itr = vdat[f].begin(); itr != vdat[f].end(); itr++) { 
// 				// itr works as a pointer to pair<string, double> 
// 				// type itr->first stores the key part  and 
// 				// itr->second stroes the value part
// 				string key = itr->first;
// 				vdat[f][key] /= size;
// 			}

// 			printf("Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			delete[] line; 
// 			delete[] xx;
// 			gzclose(fp); 
// 		}

// 		double * res = new double[nf * nf]; 

// 		memset(res, 0, sizeof(double)*nf*nf); 
// 		for (int i = 0; i < nf-1; i++) {
// 			for (int j = i+1; j < nf; j++) {
// 				double stat = 0;
// 				if (method.compare("h2_1d") == 0)
// 					stat = hellinger_4D_quantiles(vdat[i], vdat[j], 0);
// 				else if (method.compare("cos_1d") == 0)
// 					stat = cosine_4D_quantiles(vdat[i], vdat[j]);
// 				else if (method.compare("l1_1d") == 0)
// 					stat = l1_4D_quantiles(vdat[i], vdat[j]);
// 				printf("%s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				fprintf(flog, "### %s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				res[i*nf+j] = res[j*nf+i] = stat;
// 			}
// 		}

// 		for (int i = 0; i < nf; i++)
// 		{
// 			for (int j = 0; j < nf; j++) {
// 				fprintf(flog, "%.10f\t",  res[i*nf+j]);
// 			}
// 			fprintf(flog, "\n");
// 		}

// 		// for(int f = 0; f < nf; f++) {
// 		// 	delete[] vdat[f]; 
// 		// }
// 		delete[] vdat; 
// 		delete[] vdim;
// 		delete[] res;
// 		for (int i = 0; i < 1; i++) {
// 			delete[] brk[i];
// 		}
// 		delete[] brk;
// 	}
// 	else {
// 		double ** vdat = new double*[nf]; // store data
// 		for (int f = 0; f < nf; f++)
// 			vdat[f] = NULL; 
// 		for (int f = 0; f < nf; f++) {
// 			int size = vdim[f];	     
// 			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
// 			if(fp == NULL) {
// 				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
// 				exit(0); 
// 			}
// 			vdat[f] = new double[size * nc];
// 			int i = 0; 
// 			char * line = new char[1024]; 
// 			for (int r = 0; r < size; r++){
// 				gzgets(fp, line, 1024); 
// 				if(nc == 4) 
// 					sscanf(line, "%lf %lf %lf %lf", &vdat[f][i], &vdat[f][size+i], &vdat[f][size*2+i], &vdat[f][size*3+i]);  
// 				else if(nc == 2) 
// 					sscanf(line, "%lf %lf ", &vdat[f][i], &vdat[f][size+i]);  
// 				else if(nc == 1) 
// 					sscanf(line, "%lf ", &vdat[f][i]);  
// 				else {
// 					char * tok = strtok(line, " "); 
// 					vdat[f][i] = atof(tok); 
// 					for (int j = 1; j < nc; j++) 
// 					{
// 						tok=strtok(NULL, " \t"); 
// 						if(tok != NULL) 
// 							vdat[f][size*j+i] = atof(tok); 
// 					} //first four entries are scores;  
// 					}
// 				i++; 
// 			}
// 			printf("Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
// 			delete[] line; 
// 			gzclose(fp); 
// 		}

// 		double * res = new double[nf * nf]; 

// 		memset(res, 0, sizeof(double)*nf*nf); 
// 		for (int i = 0; i < nf-1; i++) {
// 			for (int j = i+1; j < nf; j++) {
// 				double test = 0; 
// 				double stat = 0;
// 				for (int c = 0; c < nc; c++)
// 				{
// 					if(method.compare("ad") == 0) 
// 						test = adtest(&vdat[i][c*vdim[i]],&vdat[j][c*vdim[j]], vdim[i],vdim[j]);
// 					if(method.compare("ad1") == 0) 
// 						test = adtest1(&vdat[i][c*vdim[i]],&vdat[j][c*vdim[j]], vdim[i],vdim[j]);
// 					else if(method.compare("ks") == 0) 
// 						test = kstest(&vdat[i][c*vdim[i]],&vdat[j][c*vdim[j]], vdim[i],vdim[j]);
// 					else if(method.compare("h2") == 0)
// 						test = hellinger_1D_quantiles_old(&vdat[i][c*vdim[i]],&vdat[j][c*vdim[j]], vdim[i],vdim[j]);
// 					else if(method.compare("cos") == 0)
// 						test = cosine_1D_quantiles_old(&vdat[i][c*vdim[i]],&vdat[j][c*vdim[j]], vdim[i],vdim[j]);
// 					stat += test;
// 				}
// 				stat /= nc;
// 				printf("%s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				fprintf(flog, "### %s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
// 				res[i*nf+j] = res[j*nf+i] = stat; 
// 			}
// 		}

// 		for (int i = 0; i < nf; i++)
// 		{
// 			for (int j = 0; j < nf; j++)
// 				fprintf(flog, "%.10f\t",  res[i*nf+j]); 
// 			fprintf(flog, "\n"); 
// 		}

// 		delete[] vdim; 
// 		for(int f = 0; f < nf; f++)
// 			delete[] vdat[f]; 
// 		delete[] vdat; 
// 		delete[] res;
// 	}
// }

void compute_stats2(vector<string> &vfn1, vector<string> &vfn2, string fout, long s, string method, unsigned long seed, int ncol)
{
	const int nc = ncol;  
	fprintf(flog, "### cad2 -m %s -s %ld -R %ld -o %s ", method.c_str(), s, seed, fout.c_str());  
	for (unsigned i = 0; i < vfn1.size(); i++) 
		fprintf(flog, " -i %s ", vfn1.at(i).c_str());  
	for (unsigned i = 0; i < vfn2.size(); i++) 
		fprintf(flog, " -j %s ", vfn2.at(i).c_str());  
	fprintf(flog, " \n");  
	int nf1 = (int) vfn1.size(); 
	double ** vdat1 = new double*[nf1]; 
	for (int f = 0; f < nf1; f++)
		vdat1[f] = NULL; 
	int * vdim1 = new int[nf1]; 
	for (int f = 0; f < nf1; f++)
		vdim1[f] = filerow(vfn1.at(f)); 

	int nf2 = (int) vfn2.size(); 
	double * vdat2 = NULL; 
	int * vdim2 = new int[nf2]; 
	for (int f = 0; f < nf2; f++)
		vdim2[f] = filerow(vfn2.at(f)); 

	int size = 0; 
	double md = 1e100; 
	for (int f = 0; f < nf1; f++)
		if(vdim1[f] < md) md = vdim1[f]; 
	for (int f = 0; f < nf2; f++)
		if(vdim2[f] < md) md = vdim2[f]; 
	if(s > 0 && s < md) 
		 size = s; 
	else  
		 size = md; 
	 for(int f = 0; f < nf1; f++)
	 {
		 printf("### File %s will be sampled %d many reads out of %d.\n", vfn1.at(f).c_str(), size, vdim1[f]); 
		 fprintf(flog, "### File %s will be sampled %d many reads out of %d.\n", vfn1.at(f).c_str(), size, vdim1[f]); 
	 }
	 for(int f = 0; f < nf2; f++)
	 {
		 printf("### File %s will be sampled %d many reads out of %d.\n", vfn2.at(f).c_str(), size, vdim2[f]); 
		 fprintf(flog, "### File %s will be sampled %d many reads out of %d.\n", vfn2.at(f).c_str(), size, vdim2[f]); 
	 }

	for (int f = 0; f < nf1; f++)
	{
	    int * index = new int[vdim1[f]]; 
		for (int i = 0; i < vdim1[f]; i++)
			index[i] = i; 
		int * chosen = new int[size]; 
		gsl_ran_choose(gsl_r, chosen, size, index, vdim1[f], sizeof(int)); 
		qsort((int *) chosen, size, sizeof(int), compare_int); 

		vdat1[f] = new double[size * nc]; 

		 gzFile fp = gzopen(vfn1.at(f).c_str(), "rb"); 
		 if(fp == NULL) {
			 printf(" can't open file, %s \n", vfn1.at(f).c_str()); 
			 exit(0); 
		 }

		 int i = 0; 
		 char * line = new char[1024]; 
		 for (int r = 0; r < vdim1[f]; r++){
			 gzgets(fp, line, 1024); 
			 if(r == chosen[i]) {
//				 if(nc == 4) 
//				 	sscanf(line, "%lf %lf %lf %lf", &vdat1[f][i], &vdat1[f][size+i], &vdat1[f][size*2+i], &vdat1[f][size*3+i]);  
//				 else if(nc == 2) 
//				 	sscanf(line, "%lf %lf", &vdat1[f][i], &vdat1[f][size+i]);  
//				 else if(nc == 1) 
//				 	sscanf(line, "%lf ", &vdat1[f][i]);  
//				 else 
				 {
					 char * tok = strtok(line, " "); 
					 vdat1[f][i] = atof(tok); 
					 for (int j = 1; j < nc; j++) 
					 {
						 tok=strtok(NULL, " \t"); 
						 if(tok != NULL) 
							 vdat1[f][size*j+i] = atof(tok); 
					 } 
				 }
				 i++; 
			 }
		 }
		 printf("Read %d random entries in file %s.\n" , i, vfn1.at(f).c_str()); 
		 fprintf(flog, "### Read %d random entries in file %s.\n" , i, vfn1.at(f).c_str()); 
		 delete[] line; 
		 delete[] chosen; 
		 delete[] index; 
		 gzclose(fp); 
	}
	
	double * res = new double[nf1 * nf2]; 
	memset(res, 0, sizeof(double)*nf1*nf2); 
	vdat2 = new double[size * nc]; 
	for (int f = 0; f < nf2; f++)
	{
	    int * index = new int[vdim2[f]]; 
		for (int i = 0; i < vdim2[f]; i++)
			index[i] = i; 
		int * chosen = new int[size]; 
		gsl_ran_choose(gsl_r, chosen, size, index, vdim2[f], sizeof(int)); 
		qsort((int *) chosen, size, sizeof(int), compare_int); 

		 gzFile fp = gzopen(vfn2.at(f).c_str(), "rb"); 
		 if(fp == NULL) {
			 printf(" can't open file, %s \n", vfn2.at(f).c_str()); 
			 exit(0); 
		 }

		 int i = 0; 
		 char * line = new char[1024]; 
		 for (int r = 0; r < vdim2[f]; r++){
			 gzgets(fp, line, 1024); 
			 if(r == chosen[i]) {
//				 if(nc == 4) 
//				 	sscanf(line, "%lf %lf %lf %lf", &vdat2[i], &vdat2[size+i], &vdat2[size*2+i], &vdat2[size*3+i]);  
//				 else if (nc == 2) 
//				 	sscanf(line, "%lf %lf", &vdat2[i], &vdat2[size+i]);  
//				 else if (nc == 1) 
//				 	sscanf(line, "%lf ", &vdat2[i]);  
//				 else
				 {
					 char * tok = strtok(line, " "); 
					 vdat2[i] = atof(tok); 
					 for (int j = 1; j < nc; j++) 
					 {
						 tok=strtok(NULL, " \t"); 
						 if(tok != NULL) 
							 vdat2[size*j+i] = atof(tok); 
					 } 
				 }
				 i++; 
			 }
		 }
		 printf("Read %d random entries in file %s.\n" , i, vfn2.at(f).c_str()); 
		 fprintf(flog, "### Read %d random entries in file %s.\n" , i, vfn2.at(f).c_str()); 

		for (int i = 0; i < nf1; i++)
		{
			double stat = 0; 
			double test = 0; 
			for (int c = 0; c < nc; c++)
			{
				if(method.compare("ad") == 0) 
					test = adtest(&vdat1[i][c*size], &vdat2[c*size], size, size);
				else if(method.compare("ks") == 0) 
					test = kstest(&vdat1[i][c*size], &vdat2[c*size], size, size);
				else if(method.compare("h2") == 0)
					test = hellinger(&vdat1[i][c*size], &vdat2[c*size], size, size);
				else if(method.compare("cos") == 0)
					test = cosine(&vdat1[i][c*size], &vdat2[c*size], size, size);
//			printf("%s\t %d\t %d\t %lf \n",method.c_str(), size*4, size*4, stat); 
//			fprintf(flog, "### %s\t %d\t %d\t %lf \n",method.c_str(), size*4, size*4, stat); 
                stat += test; 
			}
			stat /= nc; 

//			{
//				if(method.compare("ad") == 0) 
//					stat = adtest(&vdat1[i][0], &vdat2[0], nc*size, nc*size);
//				else if(method.compare("ks") == 0) 
//					stat = kstest(&vdat1[i][0], &vdat2[0], nc*size, nc*size);
//				else 
//					stat = hellinger(&vdat1[i][0], &vdat2[0], nc*size, nc*size);
//			}

			printf("%s\t -j %d\t -i %d\t %lf \n",method.c_str(), f+1, i+1, stat); 
			fprintf(flog, "### %s\t -j %d\t -i %d\t %lf \n",method.c_str(), f+1, i+1, stat); 
			res[i*nf2+f] = stat; 
		}

		 delete[] line; 
		 delete[] chosen; 
		 delete[] index; 
		 gzclose(fp); 
	}
	delete[] vdat2; 

	for (int i = 0; i < nf1; i++)
	{
		for (int j = 0; j < nf2; j++)
			fprintf(flog, "%lf\t",  res[i*nf2+j]); 
		fprintf(flog, "\n"); 
	}

	delete[] vdim1; 
	delete[] vdim2; 
	for(int f = 0; f < nf1; f++)
		delete[] vdat1[f]; 
	delete[] vdat1; 
	delete[] res;; 
}

void compute_stats2_unequal_datasize(vector<string> &vfn1, vector<string> &vfn2, string fout, string method, int ncol, int n_bin, string bin_file)
{
	const int nc = ncol;
	const int nb = n_bin; 
	const double pseudo_count = 0;
	fprintf(flog, "### cad2 -m %s -o %s ", method.c_str(), fout.c_str());  
	for (unsigned i = 0; i < vfn1.size(); i++) 
		fprintf(flog, " -i %s ", vfn1.at(i).c_str());  
	for (unsigned i = 0; i < vfn2.size(); i++) 
		fprintf(flog, " -j %s ", vfn2.at(i).c_str());  
	fprintf(flog, " \n");  
	int nf1 = (int) vfn1.size(); // how many files are there
	int * vdim1 = new int[nf1]; // the number of lines for each file
	for (int f = 0; f < nf1; f++)
		vdim1[f] = filerow(vfn1.at(f)); 

	int nf2 = (int) vfn2.size(); // how many files are there
	int * vdim2 = new int[nf2]; // the number of lines for each file
	for (int f = 0; f < nf2; f++)
		vdim2[f] = filerow(vfn2.at(f)); 

	if(nc == 4 && method.find("_4d") != string::npos) { // if 4-dimensional joint distribution would be used
		/* get bins according to input bin_file */
		double * bins_1st = new double[nb + 1];

		double ** bins_2nd_all = new double * [nb];
		for (int f = 0; f < nb; f++) {
			bins_2nd_all[f] = new double[nb + 1];
		}

		double *** bins_3rd_all = new double ** [nb];
		for (int f = 0; f < nb; f++) {
			bins_3rd_all[f] = new double * [nb];
			for (int g = 0; g < nb; g++)
				bins_3rd_all[f][g] = new double[nb + 1];
		}

		double **** bins_4th_all = new double *** [nb];
		for (int f = 0; f < nb; f++) {
			bins_4th_all[f] = new double ** [nb];
			for (int g = 0; g < nb; g++) {
				bins_4th_all[f][g] = new double * [nb];
				for (int h = 0; h < nb; h++)
					bins_4th_all[f][g][h] = new double[nb + 1];
			}
		}

		gzFile myfile = gzopen(bin_file.c_str(), "rb");
		if (myfile == NULL)  {
			printf(" No bin file detected, calculate bin partitions\n");
			exit(0);
		}
		else {
			char * line = new char[1024000];
			for (int i = 0; i < (std::pow(nb, 3) + std::pow(nb, 2) + std::pow(nb, 1) + std::pow(nb, 0)); i++) {
				// printf("%d\n", i);
				gzgets(myfile, line, 1024000);
				char * tok = strtok(line, " ");
				if (i < std::pow(nb, 0)) {
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_1st[j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
				}
				else if (i < (std::pow(nb, 1) + std::pow(nb, 0)))
				{
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_2nd_all[i - 1][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
			  		// printf("finish line %d\n", i);
				}
				else if (i < (std::pow(nb, 2) + std::pow(nb, 1) + std::pow(nb, 0)))
				{
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_3rd_all[static_cast<int>(ceil(i / std::pow(nb, 1)) - 2)][static_cast<int>(i - (ceil(i / std::pow(nb, 1)) - 1) * std::pow(nb, 1) - 1)][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
			  		// printf("finish line %d\n", i);
				}
				else
				{
					int l = nb - 1;
					if (i % nb != 0) {
						l = i % nb - 1;
					}
					// printf("%d\n", static_cast<int>(ceil((i - 10) / std::pow(nb, 2)) - 2));
					// printf("%d\n", static_cast<int>(ceil((i - (ceil((i - 10) / std::pow(nb, 2)) - 1) * std::pow(nb, 2)) / std::pow(nb, 1)) - 2));
					// printf("%d\n", l);
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_4th_all[static_cast<int>(ceil((i - 10) / std::pow(nb, 2)) - 2)][static_cast<int>(ceil((i - (ceil((i - 10) / std::pow(nb, 2)) - 1) * std::pow(nb, 2)) / std::pow(nb, 1)) - 2)][l][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
			  		// printf("finish line %d\n", i);
				}
			}
			delete[] line;
			gzclose(myfile);
			// for (int i = 0; i < 4; i++) {
			// 	brk[i][0] = -1 * INFINITY;
			// 	brk[i][nb] = INFINITY;
			// }
		}

		unordered_map<string, double> * vdat1 = new unordered_map<string, double>[nf1];
		long int pseudo_size = 1000000 + pseudo_count * std::pow(nb, 4);
		double pseudo_freq = pseudo_count / pseudo_size;
		// printf("%.15lf\n", pseudo_freq);
		for (int f = 0; f < nf1; f++) {
			int size = vdim1[f];    
			gzFile fp = gzopen(vfn1.at(f).c_str(), "rb"); 
			if(fp == NULL) {
				printf(" can't open file, %s \n", vfn1.at(f).c_str()); 
				exit(0); 
			}
			char * line = new char[1024];
			double * xx = new double [4];
			int i = 0; 
			for (int r = 0; r < size; r++){
				gzgets(fp, line, 1024); 
				sscanf(line, "%lf %lf %lf %lf", &xx[0], &xx[1], &xx[2], &xx[3]);
				string bin_num; //bin_num track which bin should the data point go in for each of the 4 columns, like "3,76,24,18"
				int cursor1 = quickselect(bins_1st, nb, xx[0]);
				bin_num += to_string((long long) cursor1);
				bin_num += ",";
				int cursor2 = quickselect(bins_2nd_all[cursor1], nb, xx[1]);
				bin_num += to_string((long long) cursor2);
				bin_num += ",";
				int cursor3 = quickselect(bins_3rd_all[cursor1][cursor2], nb, xx[2]);
				bin_num += to_string((long long) cursor3);
				bin_num += ",";
				int cursor4 = quickselect(bins_4th_all[cursor1][cursor2][cursor3], nb, xx[3]);
				bin_num += to_string((long long) cursor4);
				vdat1[f][bin_num] += 1;
				i++;
			}
			// printf("%zu\n", vdat[f].size());
			/* get the density for each 'volumn' */
			unordered_map<string, double>:: iterator itr;
			for (itr = vdat1[f].begin(); itr != vdat1[f].end(); itr++) { 
				// itr works as a pointer to pair<string, double> 
				// type itr->first stores the key part  and 
				// itr->second stroes the value part
				string key = itr->first;
				// cout << key << "\t" << itr->second << "\n";
				vdat1[f][key] = ((1000000 * vdat1[f][key] / size) + pseudo_count) / pseudo_size;
			}
			printf("Read %d entries in file %s.\n" , i, vfn1.at(f).c_str()); 
			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn1.at(f).c_str()); 
			delete[] line; 
			delete[] xx;
			gzclose(fp); 
		}

		double * res = new double[nf1 * nf2]; 
		memset(res, 0, sizeof(double)*nf1*nf2);

		for (int f = 0; f < nf2; f++)
		{
			unordered_map<string, double> vdat2;

			gzFile fp = gzopen(vfn2.at(f).c_str(), "rb"); 
			if(fp == NULL) {
				printf(" can't open file, %s \n", vfn2.at(f).c_str()); 
				exit(0); 
			}
			char * line = new char[1024]; 
			double * xx = new double [4];
			int i = 0;
			for (int r = 0; r < vdim2[f]; r++){
				gzgets(fp, line, 1024); 
				sscanf(line, "%lf %lf %lf %lf", &xx[0], &xx[1], &xx[2], &xx[3]); 
				string bin_num; //bin_num track which bin should the data point go in for each of the 4 columns, like "3,76,24,18"
				int cursor1 = quickselect(bins_1st, nb, xx[0]);
				bin_num += to_string((long long) cursor1);
				bin_num += ",";
				int cursor2 = quickselect(bins_2nd_all[cursor1], nb, xx[1]);
				bin_num += to_string((long long) cursor2);
				bin_num += ",";
				int cursor3 = quickselect(bins_3rd_all[cursor1][cursor2], nb, xx[2]);
				bin_num += to_string((long long) cursor3);
				bin_num += ",";
				int cursor4 = quickselect(bins_4th_all[cursor1][cursor2][cursor3], nb, xx[3]);
				bin_num += to_string((long long) cursor4);
				vdat2[bin_num] += 1;
				i++;
			}
			/* get the density for each 'volumn' */ 
			unordered_map<string, double>:: iterator itr;
			for (itr = vdat2.begin(); itr != vdat2.end(); itr++) { 
				// itr works as a pointer to pair<string, double> 
				// type itr->first stores the key part  and 
				// itr->second stroes the value part
				string key = itr->first;
				// cout << key << "\t" << itr->second << "\n";
				vdat2[key] = ((1000000 * vdat2[key] / vdim2[f]) + pseudo_count) / pseudo_size;
			}
			printf("Read %d entries in file %s.\n" , i, vfn2.at(f).c_str()); 
			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn2.at(f).c_str()); 

			for (int i = 0; i < nf1; i++)
			{
				double stat = 0;  
				if(method.compare("h2_4d") == 0)
					stat = hellinger_4D_quantiles(vdat1[i], vdat2, pseudo_freq);
				else if(method.compare("cos_4d") == 0)
					stat = cosine_4D_quantiles(vdat1[i], vdat2);
				else if (method.compare("ham_4d") == 0)
					stat = hamming_4D_quantiles(vdat1[i], vdat2);
				else if (method.compare("l1_4d") == 0)
					stat = l1_4D_quantiles(vdat1[i], vdat2);
				printf("%s\t -j %d\t -i %d\t %lf \n",method.c_str(), f+1, i+1, stat); 
				fprintf(flog, "### %s\t -j %d\t -i %d\t %lf \n",method.c_str(), f+1, i+1, stat); 
				res[i*nf2+f] = stat; 
			}

			delete[] line; 
			delete[] xx;
			gzclose(fp); 
		}
		 

		for (int i = 0; i < nf1; i++)
		{
			for (int j = 0; j < nf2; j++)
				fprintf(flog, "%.10f\t",  res[i*nf2+j]); 
			fprintf(flog, "\n"); 
		}

		// for(int f = 0; f < nf1; f++) {
		// 	delete[] vdat1[f]; 
		// }
		delete[] vdat1; 
		delete[] vdim1; 
		delete[] vdim2; 
		delete[] res;
		delete[] bins_1st;

		for (int i = 0; i < nb; i++) {
			delete[] bins_2nd_all[i];
		}
		delete[] bins_2nd_all;

		for (int i = 0; i < nb; i++) {
			for (int j = 0; j < nb; j++) {
				delete[] bins_3rd_all[i][j];
			}
			delete[] bins_3rd_all[i];
		}
		delete[] bins_3rd_all;

		for (int i = 0; i < nb; i++) {
			for (int j = 0; j < nb; j++) {
				for (int k = 0; k < nb; k++) {
					delete[] bins_4th_all[i][j][k];
				}
				delete[] bins_4th_all[i][j];
			}
			delete[] bins_4th_all[i];
		}
		delete[] bins_4th_all;
	}
}

void twoalleles(void) 
{
	printf("%d \n", INT_MAX); 
	exit(0); 
	 gzFile fi = gzopen("test.input.txt", "rb"); 
	 if(fi == NULL) {
		 printf(" can't open file\n"); 
		 exit(0); 
	 }
	 long bufsize = 100; 
	 char * buf = new char[bufsize+1]; 
	 int d = 0; 
	 while(1) 
	 {
		 long len=gzread(fi, buf+d, (bufsize-d)*sizeof(char)); 
		 printf("%s\n", buf); 
		 if(len == 0) break; 
		 char * next = buf; 
		 long curlen = 0; 
		 while(1) {
		     d=strcspn(next, "\n"); 
			 curlen += d; 
			 buf[d] = '\0'; 
			 printf("%ld %d \t %s \n", curlen, d, next); 
			 if(curlen >= bufsize) 
				 break; 
			 next += d+1;  
		 } 
		 memcpy(buf, next, d* sizeof(char)); 
	 }
     
	 gzclose(fi); 


//	int len = 1000000; 
//	int rl = 100; 
//	int nr = len/2; 
//	int mut = len/100; 
//	int indel = mut/10; 
//	double af, at, ac, ag; 
//	af=at= 0.2;
//	ac=ag=0.3; 
//	
//    char * a1 = new char[len]; 
//	char * a2 = new char[len];
//	double * xx = new double[nr*4]; 
//	double * yy = new double[nr*4]; 
//	char * read = new char[rl]; 
//
//	for(int rep = 0; rep<1000; rep++)
//	{
//		for (int i = 0; i < len; i++)
//		{
//			double rr = gsl_rng_uniform(gsl_r); 
//			char nu;  
//			if(rr < af) nu = 'A';  
//			else if(rr < (af + at) ) nu = 'T'; 
//			else if (rr < (af + at + ac)) nu = 'C'; 
//			else nu = 'G'; 
//			a1[i] = a2[i] = nu; 
//		}
//	//    for (int i = 0; i < len; i++)
//	//	{
//	//    	double rr = gsl_rng_uniform(gsl_r); 
//	//		char nu;  
//	//        if(rr < af) nu = 'A';  
//	//		else if(rr < (af + at) ) nu = 'T'; 
//	//		else if (rr < (af + at + ac)) nu = 'C'; 
//	//		else nu = 'G'; 
//	//		a2[i] = nu; 
//	//	}
//		//simulate an allele; need to read from the reference genome; 
//		
//		char nuc[4] = {'A','T','C','G'}; 
//		for (int j = 0; j < mut; j++)
//		{
//			int loc = gsl_rng_uniform_int(gsl_r, len); 
//			int rr = gsl_rng_uniform_int(gsl_r, 4); 
//			rr = (rr+1) % 4; 
//			a2[loc] = nuc[rr]; 
//		} //random mutations; 
//
//		for (int i = 0; i < indel; i++)
//		{
//			int loc = gsl_rng_uniform_int(gsl_r, len-20); 
//			int size = gsl_rng_uniform_int(gsl_r, 10)+2;
//			for (int k = 0; k < size; k++) 
//			{
//				int rr = gsl_rng_uniform_int(gsl_r, 4); 
//				rr = (rr+1) % 4; 
//				a2[loc+k] = 'N'; //nuc[rr]; 
//			} 
//		} //mask a consecutive segment to approximate indel; 
//
//		map<int, char> i2c; 
//		i2c[0] = 'A'; 
//		i2c[1] = 'T'; 
//		i2c[2] = 'C'; 
//		i2c[3] = 'G'; 
//		for(int i = 0; i < nr; i++)
//		{
//			int beg = gsl_rng_uniform_int(gsl_r, len-rl); 
//			for (int j = 0; j < rl; j++)
//				read[j] = a1[beg+j]; 
//
//			for(int nu = 0; nu < 4; nu++) {
//				double prod[dd][dd] = {{1,0},{0,1}}; 
//				for(int s = 0; s < rl; s++)
//					operation(prod, (int) (toupper(read[s]) == i2c[nu])); 
//				xx[nu*nr+i] = log(prod[0][0]+prod[1][1]);
//			}
//		}
//
//		for(int i = 0; i < nr; i++)
//		{
//			int beg = gsl_rng_uniform_int(gsl_r, len-rl); 
//			for (int j = 0; j < rl; j++)
//				read[j] = a2[beg+j]; 
//
//			for(int nu = 0; nu < 4; nu++) {
//				double prod[dd][dd] = {{1,0},{0,1}}; 
//				for(int s = 0; s < rl; s++)
//					operation(prod, (int) (toupper(read[s]) == i2c[nu])); 
//				yy[nu*nr+i] = log(prod[0][0]+prod[1][1]);
//			}
//		}
//
//		double res[4]; 
//		double mean = 0; 
//		for (int i = 0; i < 4; i++)
//		{
//			res[i] = adtest(&xx[i*nr], &yy[i*nr], nr, nr); 
//			mean += res[i]; 
//		}
//		mean /= 4.0; 
//		printf("%lf ", mean); 
//		for (int i = 0; i < 4; i++)
//			printf("%lf ", res[i]); 
//		printf("\n "); 
//
//	}
//	delete[] xx; 
//	delete[] yy; 
//	delete[] read; 
//	delete[] a1; 
//	delete[] a2; 
}

void compute_alpha_unequal_datasize(vector<string> &vfn, string fout, string method, int withgc, int ncol, int n_bin, string bin_file)
{
	const int nc = ncol;
	const int nb = n_bin;
	if(withgc == 1 && strstr(vfn.at(0).c_str(), "nogc") != NULL) 
	{
		printf("### nogc file detected, toggle off withgc. \n");  
		fprintf(flog, "### nogc file detected, toggle off withgc.\n");  
		withgc = 0; 
	}
	fprintf(flog, "### cad -m %s -o %s", method.c_str(), fout.c_str());  
	for (unsigned i = 0; i < vfn.size(); i++) 
		fprintf(flog, " -i %s ", vfn.at(i).c_str());  
	fprintf(flog, " \n");  
	int nf = (int) vfn.size(); // how many files are there
	
	int * vdim = new int[nf]; // the number of lines for each file
	for (int f = 0; f < nf; f++)
		vdim[f] = filerow(vfn.at(f)); 

	if(nc == 4 && method.find("_4d") != string::npos) { // if 4-dimensional joint distribution would be used
		/* calculate bins according to quantiles */
		double * bins_1st = new double[nb + 1];

		double ** bins_2nd_all = new double * [nb];
		for (int f = 0; f < nb; f++) {
			bins_2nd_all[f] = new double[nb + 1];
		}

		double *** bins_3rd_all = new double ** [nb];
		for (int f = 0; f < nb; f++) {
			bins_3rd_all[f] = new double * [nb];
			for (int g = 0; g < nb; g++)
				bins_3rd_all[f][g] = new double[nb + 1];
		}

		double **** bins_4th_all = new double *** [nb];
		for (int f = 0; f < nb; f++) {
			bins_4th_all[f] = new double ** [nb];
			for (int g = 0; g < nb; g++) {
				bins_4th_all[f][g] = new double * [nb];
				for (int h = 0; h < nb; h++)
					bins_4th_all[f][g][h] = new double[nb + 1];
			}
		}

		gzFile myfile = gzopen(bin_file.c_str(), "rb");
		if (myfile == NULL) {
			printf(" No bin file detected, calculate bin partitions\n");
			/* read the data for the first time */
			long sum_all_size = 0;
			for (int f = 0; f < nf; f++) {
				sum_all_size += vdim[f];
			}
			double ** vdat1 = new double * [4];
			for (int f = 0; f < 4; f++)
				vdat1[f] = new double [sum_all_size];
			int i = 0;
			for (int f = 0; f < nf; f++) {
				int size = vdim[f];	     
				gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
				if(fp == NULL) {
					printf(" can't open file, %s \n", vfn.at(f).c_str()); 
					exit(0); 
				}
				char * line = new char[1024];
				for (int r = 0; r < size; r++){
					gzgets(fp, line, 1024); 
					sscanf(line, "%lf %lf %lf %lf", &vdat1[0][i], &vdat1[1][i], &vdat1[2][i], &vdat1[3][i]);
					i++;
				}
				delete[] line; 
				gzclose(fp); 
			}

			bins_1st = select_quantiles(vdat1[0], sum_all_size, nb);
			for (int i = 0; i < nb + 1; i++)
				printf("%f\n", bins_1st[i]);

			for (int i = 0; i < nb; i++) {
				double ** vdat1_2 = new double * [3];
				for (int f = 0; f < 3; f++)
					vdat1_2[f] = new double [sum_all_size]; // excessively large size
				int count1 = 0;
				for (int l = 0; l < sum_all_size; l++) {
					if (vdat1[0][l] < bins_1st[i + 1] && vdat1[0][l] >= bins_1st[i]) {
						vdat1_2[0][count1] = vdat1[1][l]; 
						vdat1_2[1][count1] = vdat1[2][l]; 
						vdat1_2[2][count1] = vdat1[3][l];
						count1 ++; 
					}
				}
				printf("%d\n", count1);
				double ** vdat1_2_copy = new double * [3];
				for (int f = 0; f < 3; f++)
					vdat1_2_copy[f] = new double [count1]; // correct size
				for (int l = 0; l < count1; l++) {
					vdat1_2_copy[0][l] = vdat1_2[0][l]; 
					vdat1_2_copy[1][l] = vdat1_2[1][l];
					vdat1_2_copy[2][l] = vdat1_2[2][l];
				}
				for (int f = 0; f < 3; f++)
					delete[] vdat1_2[f];
				delete[] vdat1_2;
				printf("i%d, finish copy and delete\n", i);
				bins_2nd_all[i] = select_quantiles(vdat1_2_copy[0], count1, nb);
				for (int hehe = 0; hehe < nb + 1; hehe++)
					printf("%f\n", bins_2nd_all[i][hehe]);
				printf("i%d, finish select quantiles\n", i);

				for (int j = 0; j < nb; j++) {
					double ** vdat1_3 = new double * [2];
					for (int f = 0; f < 2; f++)
						vdat1_3[f] = new double [count1]; // excessively large size
					int count2 = 0;
					for (int l = 0; l < count1; l++) {
						if (vdat1_2_copy[0][l] < bins_2nd_all[i][j + 1] && vdat1_2_copy[0][l] >= bins_2nd_all[i][j]) {
							vdat1_3[0][count2] = vdat1_2_copy[1][l]; 
							vdat1_3[1][count2] = vdat1_2_copy[2][l]; 
							count2 ++;
						}
					}
					printf("%d\n", count2);
					double ** vdat1_3_copy = new double * [2];
					for (int f = 0; f < 2; f++)
						vdat1_3_copy[f] = new double [count2]; // correct size
					for (int l = 0; l < count2; l++) {
						vdat1_3_copy[0][l] = vdat1_3[0][l];
						vdat1_3_copy[1][l] = vdat1_3[1][l];
					}
					for (int f = 0; f < 2; f++)
						delete[] vdat1_3[f];
					delete[] vdat1_3;
					printf("j%d, finish copy and delete\n", j);
					bins_3rd_all[i][j] = select_quantiles(vdat1_3_copy[0], count2, nb);
					for (int hehe = 0; hehe < nb + 1; hehe++)
						printf("%f\n", bins_3rd_all[i][j][hehe]);
					printf("j%d, finish select quantiles\n", j);

					for (int k = 0; k < nb; k++) {
						double * vdat1_4 = new double [count2]; // excessively large size
						int count3 = 0;
						for (int l = 0; l < count2; l++) {
							if (vdat1_3_copy[0][l] < bins_3rd_all[i][j][k + 1] && vdat1_3_copy[0][l] >= bins_3rd_all[i][j][k]) {
								vdat1_4[count3] = vdat1_3_copy[1][l];
								count3++;
							}
						}
						printf("%d\n", count3);
						double * vdat1_4_copy = new double [count3]; // correct size
						for (int l = 0; l < count3; l++) {
							vdat1_4_copy[l] = vdat1_4[l];
						}
						delete[] vdat1_4;
						printf("k%d, finish copy and delete\n", k);
						bins_4th_all[i][j][k] = select_quantiles(vdat1_4_copy, count3, nb);
						for (int hehe = 0; hehe < nb + 1; hehe++)
							printf("%f\n", bins_4th_all[i][j][k][hehe]);
						printf("k%d, finish select quantiles\n", k);
						delete[] vdat1_4_copy;
					}
					for (int f = 0; f < 2; f++)
						delete[] vdat1_3_copy[f];
					delete[] vdat1_3_copy;
				}
				for (int f = 0; f < 3; f++)
					delete[] vdat1_2_copy[f];
				delete[] vdat1_2_copy;
			}
			for ( int j = 0; j < 4; j++ )
				delete[] vdat1[j];
			delete[] vdat1;

			/* write bin files for future use */
			FILE * fbin = fopen("bin.txt", "w");
			printf("Bin file name bin.txt\n");
			for (int i = 0; i < nb + 1; i++) {
				fprintf(fbin, "%.10f ", bins_1st[i]);
			}
			fprintf(fbin, "\n");

			for (int i = 0; i < nb; i++) {
				for (int j = 0; j < nb + 1; j++) {
					fprintf(fbin, "%.10f ", bins_2nd_all[i][j]);
				}
				fprintf(fbin, "\n");
			}

			for (int i = 0; i < nb; i++) {
				for (int j = 0; j < nb; j++) {
					for (int k = 0; k < nb + 1; k++) {
						fprintf(fbin, "%.10f ", bins_3rd_all[i][j][k]);
					}
					fprintf(fbin, "\n");
				}
			}

			for (int i = 0; i < nb; i++) {
				for (int j = 0; j < nb; j++) {
					for (int k = 0; k < nb; k++) {
						for (int l = 0; l < nb + 1; l++) {
							fprintf(fbin, "%.10f ", bins_4th_all[i][j][k][l]);
						}
						fprintf(fbin, "\n");
					}
				}
			}
			fclose(fbin);
		}
		else {
			char * line = new char[1024000];
			for (int i = 0; i < (std::pow(nb, 3) + std::pow(nb, 2) + std::pow(nb, 1) + std::pow(nb, 0)); i++) {
				// printf("%d\n", i);
				gzgets(myfile, line, 1024000);
				char * tok = strtok(line, " ");
				if (i < std::pow(nb, 0)) {
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_1st[j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
				}
				else if (i < (std::pow(nb, 1) + std::pow(nb, 0)))
				{
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_2nd_all[i - 1][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
			  		// printf("finish line %d\n", i);
				}
				else if (i < (std::pow(nb, 2) + std::pow(nb, 1) + std::pow(nb, 0)))
				{
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_3rd_all[static_cast<int>(ceil(i / std::pow(nb, 1)) - 2)][static_cast<int>(i - (ceil(i / std::pow(nb, 1)) - 1) * std::pow(nb, 1) - 1)][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
			  		// printf("finish line %d\n", i);
				}
				else
				{
					int l = nb - 1;
					if (i % nb != 0) {
						l = i % nb - 1;
					}
					// printf("%d\n", static_cast<int>(ceil((i - 10) / std::pow(nb, 2)) - 2));
					// printf("%d\n", static_cast<int>(ceil((i - (ceil((i - 10) / std::pow(nb, 2)) - 1) * std::pow(nb, 2)) / std::pow(nb, 1)) - 2));
					// printf("%d\n", l);
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_4th_all[static_cast<int>(ceil((i - 10) / std::pow(nb, 2)) - 2)][static_cast<int>(ceil((i - (ceil((i - 10) / std::pow(nb, 2)) - 1) * std::pow(nb, 2)) / std::pow(nb, 1)) - 2)][l][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
			  		// printf("finish line %d\n", i);
				}
			}
			delete[] line;
			gzclose(myfile);
			// for (int i = 0; i < 4; i++) {
			// 	brk[i][0] = -1 * INFINITY;
			// 	brk[i][nb] = INFINITY;
			// }
		}

		/* read the data for the second time */
		unordered_map<string, double> * vdat = new unordered_map<string, double>[nf];
		for (int f = 0; f < nf; f++) {
			int size = vdim[f];	     
			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
			if(fp == NULL) {
				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
				exit(0); 
			}
			char * line = new char[1024];
			double * xx = new double [4];
			int i = 0; 
			for (int r = 0; r < size; r++){
				gzgets(fp, line, 1024); 
				sscanf(line, "%lf %lf %lf %lf", &xx[0], &xx[1], &xx[2], &xx[3]);
				string bin_num; //bin_num track which bin should the data point go in for each of the 4 columns, like "3,76,24,18"
				int cursor1 = quickselect(bins_1st, nb, xx[0]);
				bin_num += to_string((long long) cursor1);
				bin_num += ",";
				int cursor2 = quickselect(bins_2nd_all[cursor1], nb, xx[1]);
				bin_num += to_string((long long) cursor2);
				bin_num += ",";
				int cursor3 = quickselect(bins_3rd_all[cursor1][cursor2], nb, xx[2]);
				bin_num += to_string((long long) cursor3);
				bin_num += ",";
				int cursor4 = quickselect(bins_4th_all[cursor1][cursor2][cursor3], nb, xx[3]);
				bin_num += to_string((long long) cursor4);
				vdat[f][bin_num] += 1;
				i++;
			}
			printf("%zu\n", vdat[f].size());
			/* get the density for each 'volumn' */
			unordered_map<string, double>:: iterator itr;
			for (itr = vdat[f].begin(); itr != vdat[f].end(); itr++) { 
				// itr works as a pointer to pair<string, double> 
				// type itr->first stores the key part  and 
				// itr->second stroes the value part
				string key = itr->first;
				vdat[f][key] = vdat[f][key] / size;
			}

			printf("Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
			delete[] line; 
			delete[] xx;
			gzclose(fp); 
		}

		double * res = new double[nf]; 

		memset(res, 0, sizeof(double)*nf); 
		for (int i = 0; i < nf; i++) {
			double stat = 0;
			if (method.compare("simpson_4d") == 0)
				stat = simpson(vdat[i]);
			if (method.compare("shannon_4d") == 0)
				stat = shannon(vdat[i]);
			printf("%s\t %d\t %lf \n",method.c_str(), i+1, stat); 
			res[i] = stat;
		}

		for (int i = 0; i < nf; i++)
		{
			fprintf(flog, "%.10f\n",  res[i]);
		}

		// for(int f = 0; f < nf; f++) {
		// 	delete[] vdat[f]; 
		// }
		delete[] vdat; 
		delete[] vdim;
		delete[] res;
		delete[] bins_1st;

		for (int i = 0; i < nb; i++) {
			delete[] bins_2nd_all[i];
		}
		delete[] bins_2nd_all;

		for (int i = 0; i < nb; i++) {
			for (int j = 0; j < nb; j++) {
				delete[] bins_3rd_all[i][j];
			}
			delete[] bins_3rd_all[i];
		}
		delete[] bins_3rd_all;

		for (int i = 0; i < nb; i++) {
			for (int j = 0; j < nb; j++) {
				for (int k = 0; k < nb; k++) {
					delete[] bins_4th_all[i][j][k];
				}
				delete[] bins_4th_all[i][j];
			}
			delete[] bins_4th_all[i];
		}
		delete[] bins_4th_all;
	}
}

double simpson(unordered_map<string, double> wx) // wx is unordered hash maps
{
	/* calculate the pair-wise distance for two samples */
	double res = 0;
	double temp = 0;
	unordered_map<string, double>:: iterator itr; 
	for (itr = wx.begin(); itr != wx.end(); itr++)
	{
		// itr works as a pointer to pair<string, double> 
		// type itr->first stores the key part  and 
		// itr->second stroes the value part
		// string key = itr->first;
		temp = itr->second;
		res += temp * temp;
	}

	res = 1 / res;
	return(res);
}

double shannon(unordered_map<string, double> wx) // wx is unordered hash maps
{
	/* calculate the pair-wise distance for two samples */
	double res = 0;
	double temp = 0;
	unordered_map<string, double>:: iterator itr; 
	for (itr = wx.begin(); itr != wx.end(); itr++)
	{
		// itr works as a pointer to pair<string, double> 
		// type itr->first stores the key part  and 
		// itr->second stroes the value part
		// string key = itr->first;
		temp = itr->second;
		res += temp * log2(temp);
	}

	res = (-1) * res;
	return(res);
}

void compute_stats_unequal_datasize(vector<string> &vfn, string fout, string method, int withgc, int ncol, int n_bin, string bin_file)
{
	const int nc = ncol;
	const int nb = n_bin;
	const double pseudo_count = 0;
	if(withgc == 1 && strstr(vfn.at(0).c_str(), "nogc") != NULL) 
	{
		printf("### nogc file detected, toggle off withgc. \n");  
		fprintf(flog, "### nogc file detected, toggle off withgc.\n");  
		withgc = 0; 
	}
	fprintf(flog, "### cad -m %s -o %s", method.c_str(), fout.c_str());  
	for (unsigned i = 0; i < vfn.size(); i++) 
		fprintf(flog, " -i %s ", vfn.at(i).c_str());  
	fprintf(flog, " \n");  
	int nf = (int) vfn.size(); // how many files are there
	
	int * vdim = new int[nf]; // the number of lines for each file
	for (int f = 0; f < nf; f++)
		vdim[f] = filerow(vfn.at(f)); 

	if(nc == 4 && method.find("_4d") != string::npos) { // if 4-dimensional joint distribution would be used
		/* calculate bins according to quantiles */
		double * bins_1st = new double[nb + 1];

		double ** bins_2nd_all = new double * [nb];
		for (int f = 0; f < nb; f++) {
			bins_2nd_all[f] = new double[nb + 1];
		}

		double *** bins_3rd_all = new double ** [nb];
		for (int f = 0; f < nb; f++) {
			bins_3rd_all[f] = new double * [nb];
			for (int g = 0; g < nb; g++)
				bins_3rd_all[f][g] = new double[nb + 1];
		}

		double **** bins_4th_all = new double *** [nb];
		for (int f = 0; f < nb; f++) {
			bins_4th_all[f] = new double ** [nb];
			for (int g = 0; g < nb; g++) {
				bins_4th_all[f][g] = new double * [nb];
				for (int h = 0; h < nb; h++)
					bins_4th_all[f][g][h] = new double[nb + 1];
			}
		}

		gzFile myfile = gzopen(bin_file.c_str(), "rb");
		if (myfile == NULL) {
			printf(" No bin file detected, calculate bin partitions\n");
			/* read the data for the first time */
			long sum_all_size = 0;
			for (int f = 0; f < nf; f++) {
				sum_all_size += vdim[f];
			}
			double ** vdat1 = new double * [4];
			for (int f = 0; f < 4; f++)
				vdat1[f] = new double [sum_all_size];
			int i = 0;
			for (int f = 0; f < nf; f++) {
				int size = vdim[f];	     
				gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
				if(fp == NULL) {
					printf(" can't open file, %s \n", vfn.at(f).c_str()); 
					exit(0); 
				}
				char * line = new char[1024];
				for (int r = 0; r < size; r++){
					gzgets(fp, line, 1024); 
					sscanf(line, "%lf %lf %lf %lf", &vdat1[0][i], &vdat1[1][i], &vdat1[2][i], &vdat1[3][i]);
					i++;
				}
				delete[] line; 
				gzclose(fp); 
			}

			bins_1st = select_quantiles(vdat1[0], sum_all_size, nb);
			for (int i = 0; i < nb + 1; i++)
				printf("%f\n", bins_1st[i]);

			for (int i = 0; i < nb; i++) {
				double ** vdat1_2 = new double * [3];
				for (int f = 0; f < 3; f++)
					vdat1_2[f] = new double [sum_all_size]; // excessively large size
				int count1 = 0;
				for (int l = 0; l < sum_all_size; l++) {
					if (vdat1[0][l] < bins_1st[i + 1] && vdat1[0][l] >= bins_1st[i]) {
						vdat1_2[0][count1] = vdat1[1][l]; 
						vdat1_2[1][count1] = vdat1[2][l]; 
						vdat1_2[2][count1] = vdat1[3][l];
						count1 ++; 
					}
				}
				printf("%d\n", count1);
				double ** vdat1_2_copy = new double * [3];
				for (int f = 0; f < 3; f++)
					vdat1_2_copy[f] = new double [count1]; // correct size
				for (int l = 0; l < count1; l++) {
					vdat1_2_copy[0][l] = vdat1_2[0][l]; 
					vdat1_2_copy[1][l] = vdat1_2[1][l];
					vdat1_2_copy[2][l] = vdat1_2[2][l];
				}
				for (int f = 0; f < 3; f++)
					delete[] vdat1_2[f];
				delete[] vdat1_2;
				printf("i%d, finish copy and delete\n", i);
				bins_2nd_all[i] = select_quantiles(vdat1_2_copy[0], count1, nb);
				for (int hehe = 0; hehe < nb + 1; hehe++)
					printf("%f\n", bins_2nd_all[i][hehe]);
				printf("i%d, finish select quantiles\n", i);

				for (int j = 0; j < nb; j++) {
					double ** vdat1_3 = new double * [2];
					for (int f = 0; f < 2; f++)
						vdat1_3[f] = new double [count1]; // excessively large size
					int count2 = 0;
					for (int l = 0; l < count1; l++) {
						if (vdat1_2_copy[0][l] < bins_2nd_all[i][j + 1] && vdat1_2_copy[0][l] >= bins_2nd_all[i][j]) {
							vdat1_3[0][count2] = vdat1_2_copy[1][l]; 
							vdat1_3[1][count2] = vdat1_2_copy[2][l]; 
							count2 ++;
						}
					}
					printf("%d\n", count2);
					double ** vdat1_3_copy = new double * [2];
					for (int f = 0; f < 2; f++)
						vdat1_3_copy[f] = new double [count2]; // correct size
					for (int l = 0; l < count2; l++) {
						vdat1_3_copy[0][l] = vdat1_3[0][l];
						vdat1_3_copy[1][l] = vdat1_3[1][l];
					}
					for (int f = 0; f < 2; f++)
						delete[] vdat1_3[f];
					delete[] vdat1_3;
					printf("j%d, finish copy and delete\n", j);
					bins_3rd_all[i][j] = select_quantiles(vdat1_3_copy[0], count2, nb);
					for (int hehe = 0; hehe < nb + 1; hehe++)
						printf("%f\n", bins_3rd_all[i][j][hehe]);
					printf("j%d, finish select quantiles\n", j);

					for (int k = 0; k < nb; k++) {
						double * vdat1_4 = new double [count2]; // excessively large size
						int count3 = 0;
						for (int l = 0; l < count2; l++) {
							if (vdat1_3_copy[0][l] < bins_3rd_all[i][j][k + 1] && vdat1_3_copy[0][l] >= bins_3rd_all[i][j][k]) {
								vdat1_4[count3] = vdat1_3_copy[1][l];
								count3++;
							}
						}
						printf("%d\n", count3);
						double * vdat1_4_copy = new double [count3]; // correct size
						for (int l = 0; l < count3; l++) {
							vdat1_4_copy[l] = vdat1_4[l];
						}
						delete[] vdat1_4;
						printf("k%d, finish copy and delete\n", k);
						bins_4th_all[i][j][k] = select_quantiles(vdat1_4_copy, count3, nb);
						for (int hehe = 0; hehe < nb + 1; hehe++)
							printf("%f\n", bins_4th_all[i][j][k][hehe]);
						printf("k%d, finish select quantiles\n", k);
						delete[] vdat1_4_copy;
					}
					for (int f = 0; f < 2; f++)
						delete[] vdat1_3_copy[f];
					delete[] vdat1_3_copy;
				}
				for (int f = 0; f < 3; f++)
					delete[] vdat1_2_copy[f];
				delete[] vdat1_2_copy;
			}
			for ( int j = 0; j < 4; j++ )
				delete[] vdat1[j];
			delete[] vdat1;

			/* write bin files for future use */
			FILE * fbin = fopen("bin.txt", "w");
			printf("Bin file name bin.txt\n");
			for (int i = 0; i < nb + 1; i++) {
				fprintf(fbin, "%.10f ", bins_1st[i]);
			}
			fprintf(fbin, "\n");

			for (int i = 0; i < nb; i++) {
				for (int j = 0; j < nb + 1; j++) {
					fprintf(fbin, "%.10f ", bins_2nd_all[i][j]);
				}
				fprintf(fbin, "\n");
			}

			for (int i = 0; i < nb; i++) {
				for (int j = 0; j < nb; j++) {
					for (int k = 0; k < nb + 1; k++) {
						fprintf(fbin, "%.10f ", bins_3rd_all[i][j][k]);
					}
					fprintf(fbin, "\n");
				}
			}

			for (int i = 0; i < nb; i++) {
				for (int j = 0; j < nb; j++) {
					for (int k = 0; k < nb; k++) {
						for (int l = 0; l < nb + 1; l++) {
							fprintf(fbin, "%.10f ", bins_4th_all[i][j][k][l]);
						}
						fprintf(fbin, "\n");
					}
				}
			}
			fclose(fbin);
		}
		else {
			char * line = new char[1024000];
			for (int i = 0; i < (std::pow(nb, 3) + std::pow(nb, 2) + std::pow(nb, 1) + std::pow(nb, 0)); i++) {
				// printf("%d\n", i);
				gzgets(myfile, line, 1024000);
				char * tok = strtok(line, " ");
				if (i < std::pow(nb, 0)) {
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_1st[j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
				}
				else if (i < (std::pow(nb, 1) + std::pow(nb, 0)))
				{
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_2nd_all[i - 1][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
			  		// printf("finish line %d\n", i);
				}
				else if (i < (std::pow(nb, 2) + std::pow(nb, 1) + std::pow(nb, 0)))
				{
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_3rd_all[static_cast<int>(ceil(i / std::pow(nb, 1)) - 2)][static_cast<int>(i - (ceil(i / std::pow(nb, 1)) - 1) * std::pow(nb, 1) - 1)][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
			  		// printf("finish line %d\n", i);
				}
				else
				{
					int l = nb - 1;
					if (i % nb != 0) {
						l = i % nb - 1;
					}
					// printf("%d\n", static_cast<int>(ceil((i - 10) / std::pow(nb, 2)) - 2));
					// printf("%d\n", static_cast<int>(ceil((i - (ceil((i - 10) / std::pow(nb, 2)) - 1) * std::pow(nb, 2)) / std::pow(nb, 1)) - 2));
					// printf("%d\n", l);
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_4th_all[static_cast<int>(ceil((i - 10) / std::pow(nb, 2)) - 2)][static_cast<int>(ceil((i - (ceil((i - 10) / std::pow(nb, 2)) - 1) * std::pow(nb, 2)) / std::pow(nb, 1)) - 2)][l][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
			  		// printf("finish line %d\n", i);
				}
			}
			delete[] line;
			gzclose(myfile);
			// for (int i = 0; i < 4; i++) {
			// 	brk[i][0] = -1 * INFINITY;
			// 	brk[i][nb] = INFINITY;
			// }
		}

		/* read the data for the second time */
		unordered_map<string, double> ** vdat = new unordered_map<string, double> * [nf];
		for (int i = 0; i < nf; i++)
			vdat[i] = new unordered_map<string, double>[76];
		double ** gc_count = new double * [nf];
		for (int i = 0; i < nf; i++) {
			gc_count[i] = new double[76];
			for (int j = 0; j < 76; j++)
				gc_count[i][j] = 0;
		}
		long int pseudo_size = 1000000 + pseudo_count * std::pow(nb, 4);
		double pseudo_freq = pseudo_count / pseudo_size;
		int size_total = 0; 
		for (int f = 0; f < nf; f++) {
			int size = vdim[f];	
			size_total += size;     
			gzFile fp = gzopen(vfn.at(f).c_str(), "rb"); 
			if(fp == NULL) {
				printf(" can't open file, %s \n", vfn.at(f).c_str()); 
				exit(0); 
			}
			char * line = new char[1024];
			double * xx = new double [6];
			for (int r = 0; r < size; r++){
				gzgets(fp, line, 1024); 
				sscanf(line, "%lf %lf %lf %lf %lf %lf", &xx[0], &xx[1], &xx[2], &xx[3], &xx[4], &xx[5]);
				int gc = static_cast<int>(xx[5]);
				gc_count[f][gc]++;
				string bin_num; //bin_num track which bin should the data point go in for each of the 4 columns, like "3,76,24,18"
				int cursor1 = quickselect(bins_1st, nb, xx[0]);
				bin_num += to_string((long long) cursor1);
				bin_num += ",";
				int cursor2 = quickselect(bins_2nd_all[cursor1], nb, xx[1]);
				bin_num += to_string((long long) cursor2);
				bin_num += ",";
				int cursor3 = quickselect(bins_3rd_all[cursor1][cursor2], nb, xx[2]);
				bin_num += to_string((long long) cursor3);
				bin_num += ",";
				int cursor4 = quickselect(bins_4th_all[cursor1][cursor2][cursor3], nb, xx[3]);
				bin_num += to_string((long long) cursor4);
				vdat[f][gc][bin_num] += 1;
			}
			// printf("%zu\n", vdat[f].size());
			/* get the density for each 'volumn' */
			for (int i = 0; i < 76; i++) {
				if (gc_count[f][i] != 0) {
					unordered_map<string, double>:: iterator itr;
					for (itr = vdat[f][i].begin(); itr != vdat[f][i].end(); itr++) { 
						// itr works as a pointer to pair<string, double> 
						// type itr->first stores the key part  and 
						// itr->second stroes the value part
						string key = itr->first;
						// cout << key << "\t" << itr->second << "\n";
						// vdat[f][key] = ((1000000 * vdat[f][key] / size) + pseudo_count) / pseudo_size;
						vdat[f][i][key] /= gc_count[f][i];
					}
				}
			}

			// printf("Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
			// fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
			delete[] line; 
			delete[] xx;
			gzclose(fp); 
		}

		double * gc = new double[76];
		for (int i = 0; i < 76; i++) {
			for (int j = 0; j < nf; j++) {
				gc[i] += gc_count[j][i];
			}
			gc[i] /= size_total;
		}

		double * res = new double[nf * nf]; 

		memset(res, 0, sizeof(double)*nf*nf); 
		for (int i = 0; i < nf-1; i++) {
			for (int j = i+1; j < nf; j++) {
				double stat = 0;
				for (int g = 0; g < 76; g++) {
					if (method.compare("h2_4d") == 0) {
						printf("%lf\t%10lf\t%lf\t%lf\n", hellinger_4D_quantiles(vdat[i][g], vdat[j][g], pseudo_freq), gc[g], gc_count[i][g], gc_count[j][g]);
						stat += hellinger_4D_quantiles(vdat[i][g], vdat[j][g], pseudo_freq) * gc[g];
					}
					else if (method.compare("cos_4d") == 0)
						stat += cosine_4D_quantiles(vdat[i][g], vdat[j][g]) * gc[g];
					else if (method.compare("ham_4d") == 0)
						stat += hamming_4D_quantiles(vdat[i][g], vdat[j][g]) * gc[g];
					else if (method.compare("l1_4d") == 0)
						stat += l1_4D_quantiles(vdat[i][g], vdat[j][g]) * gc[g];
				}
				printf("%s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
				fprintf(flog, "### %s\t %d\t %d\t %lf \n",method.c_str(), i+1, j+1, stat); 
				res[i*nf+j] = res[j*nf+i] = stat;
			}
		}

		for (int i = 0; i < nf; i++)
		{
			for (int j = 0; j < nf; j++) {
				fprintf(flog, "%.10f\t",  res[i*nf+j]);
			}
			fprintf(flog, "\n");
		}

		for(int f = 0; f < nf; f++) {
			delete[] vdat[f]; 
		}
		delete[] vdat; 
		delete[] vdim;
		delete[] res;
		delete[] bins_1st;
		delete[] gc;
		for (int i = 0; i < nf; i++) {
			delete[] gc_count[i];
		}
		delete[] gc_count;

		for (int i = 0; i < nb; i++) {
			delete[] bins_2nd_all[i];
		}
		delete[] bins_2nd_all;

		for (int i = 0; i < nb; i++) {
			for (int j = 0; j < nb; j++) {
				delete[] bins_3rd_all[i][j];
			}
			delete[] bins_3rd_all[i];
		}
		delete[] bins_3rd_all;

		for (int i = 0; i < nb; i++) {
			for (int j = 0; j < nb; j++) {
				for (int k = 0; k < nb; k++) {
					delete[] bins_4th_all[i][j][k];
				}
				delete[] bins_4th_all[i][j];
			}
			delete[] bins_4th_all[i];
		}
		delete[] bins_4th_all;
	}
}
