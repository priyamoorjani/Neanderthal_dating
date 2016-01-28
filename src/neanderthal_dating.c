/* Author: Priya Moorjani, Columbia University 
   Method described in Fu et al. 2014 and Moorjani et al. 2015 

   Copyright (C) Priya Moorjani (moorjani@gmail.com)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This software is distributed in the hope that it will be useful without any warranty or 
   guaranteed support whatsoever. The author is not responsible for its use, misuse, or 
   functionality. The software may be freely copied for non-commercial purposes, 
   provided this copyright notice is retained.

*/

#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <globals.h>
#include <mcmcpars.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  

#define WVERSION   "5" 
#define MAXDIS 0.01
#define BINSIZE 0.00001
#define MAXFL  50   
#define MAXSTR  512
#define NUMCH  22
extern int packmode ;

char *trashdir = "/var/tmp" ;
extern int verbose  ;
int qtmode = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 
int jackknife = YES ;

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *genooutfilename = NULL ;
char  *indoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *goodsnpname = NULL ;
char *badpairsname = NULL ;
char *markername = NULL ;
char *idname = NULL ;
double max_dist=MAXDIS, bin_size=BINSIZE;
int xchrom = -1 ;
int xnochrom = -1;
int nchrom=NUMCH;

char *outputname = NULL ;
FILE *ofile ;

double fakespacing = 0.0 ;

char  unknowngender = 'U' ;

void readcommands(int argc, char **argv) ;
void dophyscheck(SNP **snpm, int numsnps) ;
static void free_3d(int ***a3d, int levels, int rows);

int main(int argc, char **argv)
{
	int i, j, k, g, chrom ; 
	char string1[MAXSTR], string2[MAXSTR] ;
		
	// Data structure SNP and Indiv store SNP and ind info
  	SNP *cupt ;
  	Indiv *indx ;
  	int ch1, ch2 ;

	int numvind, nignore, numrisks = 1 ;
  	int markernum, idnum ;

  	ofile = stdout; 
  	packmode = YES ;

	// read par file
  	readcommands(argc, argv) ;
	printf("\n## indroll version: %s\n", WVERSION) ;

	//Check if test_indivname set
        if (idname == NULL) {
                fatalx("Fatalx: Test individual not entered\n");
	}

	if (outputname == NULL) {
		fatalx("Fatalx: Enter output file name\n");
	}
	
  	numsnps = 
    	getsnps(snpname, &snpmarkers, fakespacing,  badsnpname, &nignore, numrisks) ;

	// fakespacing 0.0 (default)

  	numindivs = getindivs(indivname, &indivmarkers) ;
  	setstatus(indivmarkers, numindivs, "Case") ;

  	setgenotypename(&genotypename, indivname) ;

  	printf("genotypename:  %s\n", genotypename) ;

  	if (genotypename != NULL)  {
   		getgenos(genotypename, snpmarkers, indivmarkers, numsnps, numindivs, nignore) ;

		// badpairsname and dophyscheck not used
   		if (badpairsname != NULL) {
   	 		loadbadpsc(snpmarkers, numsnps, NO, goodsnpname) ;
    			dobadpairs(badpairsname, snpmarkers, numsnps) ;
   		}
  	}
  	dophyscheck(snpmarkers,  numsnps) ;

  	numvind = numvalidind(indivmarkers, numindivs) ;
  	printf("\n\n") ;
  	printf("numindivs: %d numsnps: %d nignore: %d\n" , numindivs, numsnps, nignore) ; 

	// Check if any SNPs left for analysis
	if (numsnps <= nignore) { 
		fatalx("Error: No SNPs left for analysis, numsnps <= nignore\n");
	}

	// If verbose mode chosen
  	if (verbose)  {  
    		for (i=0; i<numindivs; ++i)  { 
     			indx = indivmarkers[i] ;
     			printf("%20s ", indx -> ID) ;
     			for (j=0; j<numsnps; ++j) { 
      				cupt = snpmarkers[j] ;
      				if (cupt -> ignore) continue ;
      				g = getgtypes(cupt, i) ; 
      				if (g<0) g = 9 ;  
      				printf("%1d", g) ;
     			}
     			printf("  %20s", indx -> egroup) ;
     			printnl() ;
    		}
  	}
	// Print marker information for selected markers
  	if (markername != NULL) {  
   		markernum = snpindex(snpmarkers, numsnps, markername) ;
   		if (markernum < 0) fatalx("markername %s not found\n", markername) ;
   		cupt = snpmarkers[markernum] ;
   		printf("markername: %s  %d   %9.3f %12.0f\n",cupt -> ID, cupt -> chrom, cupt -> genpos, cupt -> physpos) ;
   		for (i=0; i<numindivs; ++i) { 
    			indx = indivmarkers[i] ;
    			g = getgtypes(cupt, i) ;
    			printf("%20s %20s %2d\n", cupt -> ID, indx -> ID, g) ;
  		 }
  	}

	// Print  ind information for selected ind
   	idnum = indindex(indivmarkers, numindivs, idname) ;
   	if (idnum < 0) fatalx("Test Individual %s not found\n", idname) ;
   	indx = indivmarkers[idnum] ;

/*Start printing ind
           // ID = indname, gender = gender, egroup = popgroup
                printf("idname: %20s  %c   %20s\n", indx -> ID, indx -> gender, indx -> egroup) ;

                // Print genotype and marker info
                for (j=0; j<numsnps; ++j) {
                        cupt = snpmarkers[j] ;
                        if (cupt -> ignore) continue ;
			chrom = cupt -> chrom;
			if (chrom == xnochrom) continue;
                        g = getgtypes(cupt, idnum) ;
                        // cupt-> ID = snpname, indx-> ID = indname, g = genotype = 0/1/2 or -1 (instead of 9), cupt-> chrom = chrom, cupt->physpos = physpos
                        printf("%20s %20s %2d", cupt -> ID, indx -> ID, g) ;
                        printf("  %3d %12.0f", cupt -> chrom, cupt -> physpos) ;
                        printf(" %c %c", cupt -> alleles[0], cupt -> alleles[1]) ;
                        printnl() ;
                } 
 End printing ind */
	
	// jackknife
	int ychrom = 0; // index for chromosome to remove
        int max_bin = ceil(max_dist/bin_size);
	//printf("maxbin: %d\n", max_bin);
        int snp1, snp2,  bin; 
	int *count, *sum_ab, *sum_a, *sum_b;
        double dis;
	int t;
	SNP *cupt1, *cupt2 ;
	int geno1, geno2;
        int bin_index, check=0;
	double *cov_bin, obin;
	int count_bin;
	double mean_a=0.0, mean_b=0.0;

	for (ychrom = 0; ychrom <= nchrom; ychrom++)
	{
		if (ychrom > 0) { if (jackknife == NO) break ; }
		ZALLOC(count, max_bin, int);
		ZALLOC(sum_ab, max_bin, int);
		ZALLOC(sum_a, max_bin, int);
		ZALLOC(sum_b, max_bin, int);

        	//Initialize count, sum_ab, sum_a, sum_b
       	 	for (i = 0; i < max_bin; i++)
        	{
                	count[i]= 0;
			sum_ab[i] = 0; sum_a[i] = 0; sum_b[i] = 0;
        	}

        	for (snp1 = 0; snp1 < numsnps-1; snp1++)
        	{
			cupt1 = snpmarkers[snp1] ;
			if (cupt1 -> ignore) continue ;
			chrom = cupt1 -> chrom ;
			if (chrom == xnochrom) continue;
			if (chrom == ychrom) continue; 
   			if (chrom<1) { cupt1 -> ignore = YES; continue; }
    			if (chrom>nchrom) { cupt1 -> ignore = YES; continue; }
    			//if ((xchrom>0) && (chrom != xchrom)) { cupt1 -> ignore = YES; continue ; }
			if (cupt1 -> ignore) continue ;
			geno1 = getgtypes(cupt1, idnum);
		//	printf("geno1: %d\n", geno1);
                	if (geno1 == -1) continue;

                	for(snp2 = snp1+1; snp2 < numsnps; snp2++)
                	{
        			cupt2 = snpmarkers[snp2] ;
				if (cupt2 -> ignore) continue ;
	               		chrom = cupt2 -> chrom ;
				if (chrom == xnochrom) continue;
				if (chrom == ychrom) continue;
        	        	if (chrom<1) { cupt2 -> ignore = YES; continue; }
                		if (chrom>nchrom) { cupt2 -> ignore = YES; continue; }
                	//	if ((xchrom>0) && (chrom != xchrom)) { cupt2 -> ignore = YES; continue ; }
				if (cupt2 -> ignore) continue ;
  				geno2 = getgtypes(cupt2, idnum);
	                	if (geno2 == -1) continue;
			
				t = cupt2 -> chrom - cupt1 -> chrom ;
				if (t != 0) break ;		
		        	dis = cupt2 -> genpos - cupt1 -> genpos ; 
                        	if (dis >= max_dist) break ;
                        	bin = dis/bin_size;
				//printf("bin: %d\n", bin);
			
				// Print the pair of snps for which we are computing covariance
				//printf("bin, snp1, snp2: %d %20s %20s\n", bin, cupt1 -> ID, cupt2 -> ID);
			
				// compute covariance
                       		count[bin]++;
				sum_ab[bin] += (geno1*geno2);
				sum_a[bin] += geno1;
				sum_b[bin] += geno2;
                	}
        	}

		check =0;
        	ZALLOC(cov_bin, max_bin, double);

		// Print output
		if (ychrom == 0) {
			openit(outputname, &ofile, "w") ;
        		fprintf(ofile,"## Bin(cM)\tCovariance\n");
		}
		else {
			sprintf(string1, "%s:%d", outputname, ychrom) ;
			openit(string1, &ofile, "w");
			fprintf(ofile, "## Jackknife output: chrom %d\n", ychrom);
			fprintf(ofile,"## Bin(cM)\tCovariance\n");			
		}		
		mean_a=0.0, mean_b=0.0;
 
        	for (bin_index =0; bin_index < max_bin; bin_index++)
        	{
			count_bin = count[bin_index];
                	if (count_bin < 2) continue;
                	//printf("bin_index, pairs: %d, %d\n", bin_index, count_bin);
	
			mean_a = (double) sum_a[bin_index]/(double) count_bin;
			mean_b = (double) sum_b[bin_index]/ (double) count_bin;
	
        		//printf("bin, sum_ab, mean_a, mean_b, count: %d, %d, %3.2f, %3.2f, %d\n", bin_index, sum_ab[bin_index], mean_a, mean_b, count[bin_index]); 
			// sample covariance = 1/(n-1)* [ sum_xy - n*mean_x*mean_y]
	        	cov_bin[bin_index] = (double)(sum_ab[bin_index] - (count_bin*mean_a*mean_b))/(double) (count_bin-1);

	        	check = isnan(cov_bin[bin_index]);
                	if (check ==0)
               		{
                        	obin = (bin_index+1)*bin_size*100; //convert to cM
                        	fprintf(ofile, "%9.4f\t%9.5f\n", obin, cov_bin[bin_index]);
                	}

        	}

        	fclose(ofile);
		free(count);
		free(sum_ab);
		free(sum_a);;
		free(sum_b);
		free(cov_bin);
	}

  	printf("##end of run\n") ;
  	return 0 ;
}

// function reads in parfile
void readcommands(int argc, char **argv) 
{
  	int i,haploid=0;
  	char *parname = NULL ;
  	phandle *ph ;
  	char str[5000]  ;
  	char *tempname ;
 	int n ;

  	while ((i = getopt (argc, argv, "p:vV")) != -1) {
    		switch (i)
      		{
      			case 'p':
			parname = strdup(optarg) ;
			break;

      			case 'v':
			printf("version: %s\n", WVERSION) ; 
			break; 

      			case 'V':
			verbose = YES ;
			break; 

      			case '?':
			printf ("Usage: bad params.... \n") ;
			fatalx("bad params\n") ;
      		}
  	}

         
   	pcheck(parname,'p') ;
   	printf("parameter file: %s\n", parname) ;
   	ph = openpars(parname) ;
   	dostrsub(ph) ;

	// Read in parameters from parfile
	getstring(ph, "genotypename:", &genotypename) ;
   	getstring(ph, "snpname:", &snpname) ;
   	getstring(ph, "indivname:", &indivname) ;
   	getstring(ph, "output:", &outputname) ;
   	getstring(ph, "badsnpname:", &badsnpname) ;
   	getstring(ph, "goodsnpname:", &goodsnpname) ;
   	getstring(ph, "idname:", &idname) ;
	getdbl(ph, "maxdis:", &max_dist);
	getdbl(ph, "binsize:", &bin_size) ;
	getint(ph, "nochrom:", &xnochrom) ;   
	getint(ph, "numchrom:", &nchrom);
	getint(ph, "jackknife:", &jackknife) ;
	writepars(ph) ;
   	closepars(ph) ;
}

void dophyscheck(SNP **snpm, int numsnps) 
{
	//NP: catch places where physpos genpos are in opposite order
  	SNP *cupt, *cuptold ;
  	int i ;

  	for (i=0; i<numsnps; i++) {   
   		cupt = snpm[i] ;
   		if (i==0) cuptold = cupt ;
   		if (cupt -> isfake) continue ;
   		if (cupt -> ignore) continue ;
   		if (cupt -> chrom == cuptold -> chrom)  {
    			if (cupt -> physpos < cuptold -> physpos) {  
     				printf("physcheck %20s %15s %12.3f %12.3f %13.0f %13.0f\n", cuptold->ID, cupt -> ID, cuptold -> genpos, cupt -> genpos, cuptold -> physpos, cupt -> physpos);
    			}
   		}
   		cuptold = cupt ;
  	}
}

