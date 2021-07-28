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

#include <fstream> 
#include <iostream> 
#include <array>
#include <unordered_map>
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
     
using namespace std; 
                          
#define dd 2  
void operation(double a[][dd], int, double); 
void compute_stats_unequal_datasize(vector<string> &vfn, string fout, string method, int withgc, int ncol, int n_bin, string bin_file); 
void compute_stats2_unequal_datasize(vector<string> &vfn1, vector<string> &vfn2, string fout, string method, int ncol, int n_bin, string bin_file);
void adio(string, double*, int *); 
int compare_doubles(const void * a, const void *b); 
int compare_doubles_2D(const void *pa, const void *pb);
double * select_quantiles(double * xx, long nx, int ng);
double ** select_quantiles_all(double ** xx, long nx, int ng, int ncol);
int quickselect(double * brk, int ng, double num);
double hellinger(unordered_map<string, double> wx, unordered_map<string, double> wy);
double cosine(unordered_map<string, double> wx, unordered_map<string, double> wy);
void quantify_reads(string fin, string fout, int dim, int w, int step, int nm, int, int, double alpha); 
void regress_gc(string fin, string fout); 
void regress_gc_residual(string fin, string beta_file, string fout);
void regress_gc_beta(string fin);

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
	string method("h2"); 
	vector<string> vfn; 
	vector<string> vfn2; 
	int dim = 75; 
	int w = dim; //sliding window size
	int step = w; //sliding window step
	int nm = 1; 
	double alpha = 1.0;
	int adaptor_size = 0; 
	int cad = 0; 
	int cad2 = 0; 
	int qtf = 0; 
	int rgc = 0;  
	int rgc_beta = 0;
	int rgc_res = 0;
	int tphred = 20; 
	int withgc = 0; //toggle for whether to compute ad distance for gc contents.
	int ncol=4; 
	int n_bin=15; 
	string bin_file = "";
	string beta_file = "";
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
		else if(str.compare("qtf") == 0) {			
			qtf = 1; 
		}
		else if(str.compare("rgc") == 0) {			
			rgc = 1; 
		}
		else if(str.compare("rgc_beta") == 0) {			
			rgc_beta = 1; 
		}
		else if(str.compare("rgc_res") == 0) {			
			rgc_res = 1; 
		}
		else if(str.compare("-h") == 0) {			
			if(qtf == 1)
			{
				printf("./nubeam quad [-idonwh]\n"); 
				printf("compute quadriples for reads in (gzipped) fastq format.\n");
				printf("produces prefix.quad.gz (gc content is within) and prefix.quad.log.\n");  
				printf("-a : adaptor size (default 0)\n"); 
				printf("-i : input filename\n"); 
				printf("-o : output prefix\n"); 
				printf("-d : dimension, or length of the reads (default d=75).\n"); 
				printf("-n : number of missing nucleotide allowed.\n"); 
				printf("-w : sliding window size (default w=d).\n"); 
				printf("-S : sliding window step (default S=w).\n"); 
				printf("-f : value, plus 33 is the PHRED quality value of fastq reads.\n");
				printf("-h : print this help\n"); 
			}
			else if(cad == 1)
			{
				printf("./nubeam cad [-iomsh]\n"); 
				printf("compute pariwise distances of a set; the inputs are nubeam qtf outputs.\n"); 
				printf("produces prefix.cad.log.\n");  
				printf("-i : specifies input file which is output of nubeam.\n"); 
				printf("-o : output prefix (prefix.log contains pairwise distance matrix)\n"); 
				printf("-m : choice of methods: h2, cos.\n"); 
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
			else if(rgc_beta == 1) 
			{
                printf("./nubeam rgc_beta [-ioh]\n"); 
				printf("perform regression on gc contents from read quantification and output regression coefficients.\n");  
				printf("produces prefix.beta.log.\n");  
				printf("-i : input file name.\n"); 
				printf("-o : output prefix.\n"); 
				printf("-h : print this help\n"); 
			}
			else if(rgc_res == 1) 
			{
                printf("./nubeam rgc_res [-ioh beta]\n"); 
				printf("regress out gc contents from read quantification and output residuals.\n");  
				printf("produces prefix.nogc.gz and prefix.nogc.log.\n");  
				printf("-i : input file name.\n");
				printf("-beta : beta file name.\n");
				printf("-o : output prefix.\n"); 
				printf("-h : print this help\n"); 
			}
			else {
                printf("./nubeam [qtf, rgc_beta, rgc_res, cad, cad2]\n"); 
				printf("For example, use ./nubeam qtf -h for more options.\n"); 
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
		}	
		else if (str.compare("-method") == 0 || str.compare("-m") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				method.clear(); 
                method.assign(argv[i+1]);
				if(method.compare("h2") != 0 && method.compare("cos") !=0)
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
		else if (str.compare("-R") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    printf("wrong augument after option.\n");
                    exit(0);
                }
                seed = atol(argv[i+1]);
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
		else if (str.compare("-beta") == 0) {
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue; 
                beta_file.assign(argv[i+1]);
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
	if(rgc_beta) {
		fnlog.append(".beta.log");
		printf("Log file name: %s \n", fnlog.c_str()); 
		flog = fopen(fnlog.c_str(), "w");
		if(flog == NULL) {
			printf("can't open log file %s to write.", fnlog.c_str());
			exit(0);
		}
		regress_gc_beta(fin); 
		fclose(flog); 
		return 1; 
	}
	if(rgc_res) {
		fnlog.append(".nogc.log");
		printf("Log file name: %s \n", fnlog.c_str()); 
	   flog = fopen(fnlog.c_str(), "w");
	   if(flog == NULL) 
		   printf("can't open log file %s to write, proceed anyway.", fnlog.c_str());
		regress_gc_residual(fin, beta_file, fout); 
		fclose(flog); 
		return 1; 
	}
	if(cad) {        
		fnlog.append(".cad.log");
		printf("Log file name: %s \n", fnlog.c_str()); 
	   flog = fopen(fnlog.c_str(), "w");
	   if(flog == NULL) 
		   printf("can't open log file %s to write, proceed anyway.", fnlog.c_str());
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
		compute_stats2_unequal_datasize(vfn, vfn2, fout, method, ncol, n_bin, bin_file); 
		fclose(flog); 
		return 1;  
	}
	return 0; 
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

	string f1(fout); 
	f1.append(".quad.gz"); 
	gzFile fp1 = gzopen(f1.c_str(), "wb"); 
	if(fp1 == NULL) 
	{
		printf("can't open %s file to write\n", fout.c_str()); 
		exit(0); 
    }

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
        dim = strlen(line) - 1;
		if((dim + 1) < (w + adaptor)) 
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
			for(int s = 0; s <= (dim-w); s+=step)
			{
				string zz(""); 
				zz.append(line, s, w); 
				// write gc contents.
				int nc[4]; 
				int na = 0; 
				memset(nc, 0, 4*sizeof(int)); 
				// start the for loop
				for (int k = 0; k < w; k++) {
					// all capital letters
					zz[k] = toupper(zz[k]);
					// get GC contents
					int t1=c2i[zz[k]];
					if(t1>=0) nc[t1]++;
					else na++;
				}
				// reverse complement of zz
				string zz_rc(w, 'N');
				for (int k = 0; k < w; k++) {
					// get reverse complement
					zz_rc[k] = c2c[zz[w - 1 - k]];
				}

 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 						operation(prod, (int) (zz.at(k) == 'A'), alpha); 
 					}
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1]));
 				}
 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 						operation(prod, (int) (zz.at(k) == 'T'), alpha); 
 					}
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1]));
 				}
 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 						operation(prod, (int) (zz.at(k) == 'C'), alpha); 
 					}
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1]));
 				}
 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 						operation(prod, (int) (zz.at(k) == 'G'), alpha); 
 					}
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1]));
 				}
 				gzprintf(fp1, "%d %d\n", nc[0]+nc[1], nc[2]+nc[3]);
 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 						operation(prod, (int) (zz_rc.at(k) == 'A'), alpha); 
 					}
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1]));
 				}
 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 						operation(prod, (int) (zz_rc.at(k) == 'T'), alpha); 
 					}
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1]));
 				}
 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) {
 						operation(prod, (int) (zz_rc.at(k) == 'C'), alpha); 
 					}
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1]));
 				}
 				{
 					double prod[dd][dd] = {{1,0},{0,1}}; 
 					for(unsigned k = 0; k < zz.size(); k++) { 
 						operation(prod, (int) (zz_rc.at(k) == 'G'), alpha); 
 					}
 					gzprintf(fp1, "%lf ",  log(prod[0][0] + sqrt(3) * prod[0][1] + M_SQRT2 * prod[1][0] + sqrt(5) * prod[1][1]));
 				}
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
	gzclose(fgz); 
	delete[] tline; 
	delete[] first; 
	delete[] tscore; 
}

void operation(double prod[][dd], int yes, double aa) 
{
	if(yes == 1) {
		prod[0][1] += aa*prod[0][0]; 
		prod[1][1] += aa*prod[1][0]; 
	}
	else {
		prod[0][0] += aa*prod[0][1]; 
		prod[1][0] += aa*prod[1][1]; 
	}
}

int compare_doubles(const void * a, const void * b)
{
	return ((int) (*(double *) a > *(double *) b));  
}

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

double hellinger(unordered_map<string, double> wx, unordered_map<string, double> wy) // wx and wy are unordered hash maps
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
		if (wy.find(itr->first) != wy.end()) {
			temp = sqrt((itr->second) * wy[itr->first]);
			res += temp;
		}
	}

	res = sqrt(1 - res);
	return(res);
}

double cosine(unordered_map<string, double> wx, unordered_map<string, double> wy) // wx and wy are unordered hash maps
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

void regress_gc_beta(string fin) 
{   const int np = 3; // regressor; 
	const int ncol = 4;  //nubeam columns; 
	long long int dim1 = filerow(fin); 
	printf("File %s contains %lld lines \n", fin.c_str(), dim1);     

	 gsl_matrix * xty=gsl_matrix_alloc(np, ncol); 
	 gsl_matrix_set_zero(xty); 
	 gsl_matrix * U = gsl_matrix_alloc(np,np); 
	 gsl_matrix_set_zero(U); 

	 printf("Begin to scan the file. \n");

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
					// xx[j] = log(100 * atof(tok1)); 
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
	 printf("Finished scanning the file. Begin to calculate regression coefficients.\n"); 

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

	 printf("Obtained regression coefficients.\n"); 
	 for(int i = 0; i < np; i++)
	 {
		 for(int c = 0; c < ncol; c++)
			fprintf(flog, "%.10lf ", gsl_matrix_get(mb, i, c)); 
	     fprintf(flog, "\n"); 
	 }

	 gsl_matrix_free(mb); 
}

void regress_gc_residual(string fin, string beta_file, string fout) 
{   const int np = 3; // regressor; 
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

	gsl_matrix * mb= gsl_matrix_alloc(np,ncol); 
	gzFile myfile = gzopen(beta_file.c_str(), "rb");
	char * line = new char[1024];
	for (int i = 0; i < np; i++) {
		gzgets(myfile, line, 1024);
		char * tok = strtok(line, " ");
		int j = 0;
  		while (tok != NULL && j < ncol)
  		{
    		gsl_matrix_set(mb, i, j, atof(tok));
    		tok = strtok(NULL, " ");
    		j++;
  		}
	}
	delete[] line;
	gzclose(myfile);

	printf("Obtained the regression coefficients.\n"); 
	for(int i = 0; i < np; i++)
	{
		for(int c = 0; c < ncol; c++)
			printf("%lf ", gsl_matrix_get(mb, i, c)); 
	    printf("\n"); 
	}


	 printf("Begin to scan the file and obtain residuals.\n"); 
	 fprintf(flog, "Begin to scan the file and obtain residuals.\n"); 
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
					// xx[j] = log(100 * atof(tok1));
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
	 printf("Obtained residuals.\n"); 
	 fprintf(flog, "Obtained residuals.\n"); 

	 gzclose(pout); 
	 gsl_matrix_free(mb); 
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

int compare_int(const void * a, const void *b)
{
	return ((int) (*(int *) a > *(int *) b));  
}

void compute_stats_unequal_datasize(vector<string> &vfn, string fout, string method, int withgc, int ncol, int n_bin, string bin_file)
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

	if(nc == 4) { // if 4-dimensional joint distribution would be used
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
				}
				else
				{
					int l = nb - 1;
					if (i % nb != 0) {
						l = i % nb - 1;
					}
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_4th_all[static_cast<int>(ceil((i - 10) / std::pow(nb, 2)) - 2)][static_cast<int>(ceil((i - (ceil((i - 10) / std::pow(nb, 2)) - 1) * std::pow(nb, 2)) / std::pow(nb, 1)) - 2)][l][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
				}
			}
			delete[] line;
			gzclose(myfile);
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
			/* get the density for each 'volumn' */
			unordered_map<string, double>:: iterator itr;
			for (itr = vdat[f].begin(); itr != vdat[f].end(); itr++) { 
				// itr works as a pointer to pair<string, double> 
				// type itr->first stores the key part  and 
				// itr->second stroes the value part
				string key = itr->first;
				vdat[f][key] /= size;
			}

			printf("Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn.at(f).c_str()); 
			delete[] line; 
			delete[] xx;
			gzclose(fp); 
		}

		double * res = new double[nf * nf]; 

		memset(res, 0, sizeof(double)*nf*nf); 
		for (int i = 0; i < nf-1; i++) {
			for (int j = i+1; j < nf; j++) {
				double stat = 0;
				if (method.compare("h2") == 0)
					stat = hellinger(vdat[i], vdat[j]);
				else if (method.compare("cos") == 0)
					stat = cosine(vdat[i], vdat[j]);
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

void compute_stats2_unequal_datasize(vector<string> &vfn1, vector<string> &vfn2, string fout, string method, int ncol, int n_bin, string bin_file)
{
	const int nc = ncol;
	const int nb = n_bin; 
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

	if(nc == 4) { // if 4-dimensional joint distribution would be used
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
				}
				else
				{
					int l = nb - 1;
					if (i % nb != 0) {
						l = i % nb - 1;
					}
					int j = 0;
			  		while (tok != NULL && j < nb + 1)
			  		{
			    		bins_4th_all[static_cast<int>(ceil((i - 10) / std::pow(nb, 2)) - 2)][static_cast<int>(ceil((i - (ceil((i - 10) / std::pow(nb, 2)) - 1) * std::pow(nb, 2)) / std::pow(nb, 1)) - 2)][l][j] = atof(tok);
			    		tok = strtok(NULL, " ");
			    		j++;
			  		}
				}
			}
			delete[] line;
			gzclose(myfile);
		}

		unordered_map<string, double> * vdat1 = new unordered_map<string, double>[nf1];
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
			/* get the density for each 'volumn' */
			unordered_map<string, double>:: iterator itr;
			for (itr = vdat1[f].begin(); itr != vdat1[f].end(); itr++) { 
				// itr works as a pointer to pair<string, double> 
				// type itr->first stores the key part  and 
				// itr->second stroes the value part
				string key = itr->first;
				vdat1[f][key] /= size;
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
				vdat2[key] /= vdim2[f];
			}
			printf("Read %d entries in file %s.\n" , i, vfn2.at(f).c_str()); 
			fprintf(flog, "### Read %d entries in file %s.\n" , i, vfn2.at(f).c_str()); 

			for (int i = 0; i < nf1; i++)
			{
				double stat = 0;  
				if(method.compare("h2") == 0)
					stat = hellinger(vdat1[i], vdat2);
				else if(method.compare("cos") == 0)
					stat = cosine(vdat1[i], vdat2);
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