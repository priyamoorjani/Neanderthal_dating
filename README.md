# Neanderthal_dating
This is a method for estimating the date of Neanderthal gene flow using a single diploid genome. The main idea is to measure the extent of covariance in Neanderthal ancestry present in modern human genome to estimate the time of the Neanderthal gene flow. Details of the method and algorithm can be found in the following papers:

*Genome sequence of a 45,000-year-old modern human from western Siberia* <br />
http://www.nature.com/nature/journal/v514/n7523/full/nature13810.html

*A genetic method for dating ancient genomes provides a direct estimate of human generation interval in the last 45,000 years* <br />
Priya Moorjani, Sriram Sankararaman, Qiaomei Fu, Molly Przeworski, Nick J Patterson, David E. Reich
http://www.pnas.org/content/113/20/5652.long

### Installation
```
tar xvf neandertal_dating.tar
cd neandertal_dating/src
make
make install
```
This will create all executables in the neandertal_dating/bin/ folder. 


### Input
This method requires that the input data is available in EIGENSTRAT format (See https://reich.hms.harvard.edu/software/InputFileFormats).  

### Parameter file
```
genotypename:   input genotype file (in eigenstrat format).
snpname:   input snp file      (in eigenstrat format).
indivname:   input indiv file    (in eigenstrat format).
idname: individual_ID for the sample of interest.
binsize:   binsize (in Morgans). Range is from 0-1. Optimal binsize of 0.00001 is recommended.
maxdis:   maximum_distance (in Morgans). Range is 0-1. 
jackknife: YES/ NO (If YES, selected then program performs chromosome jackknife where one chromosome is removed in each run. Hence there will be 23 outputs (one for the whole genome and 22 for each chromosome removed).
output: output file name.
```
### Output

Output file contains the covariance values across various bins of genetic distances (2 columns: Bin_distance and Covariance (in cM)). You can use least squares to fit an exponential to this output to estimate the date of Neanderthal admixture in generations. R code for fitting a single exponential and visualizing the the output is available at:
https://github.com/DReichLab/AdmixTools/blob/master/src/rexpfit.r
and documentation for running this script can be found at: https://github.com/DReichLab/AdmixTools/blob/master/README.REXPFIT

### Example

Test data is available in the neandertal_dating/examples/ directory. To ensure that program is working, run:
```
./bin/neanderthal_dating -p parExample1
```

# Support
Send queries to Priya Moorjani <moorjani@berkeley.edu>
