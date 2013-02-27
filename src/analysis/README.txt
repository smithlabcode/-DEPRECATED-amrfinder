amrfinder
============================================================================
The following command can find allelically methylated regions (AMRs) from the
input mapped read file "Human_test.mr":

$$ amrfinder -o Human_test_AMR.bed -i 10 -w 10 -m 1 -b -c hg18 Human_test.mr

Option '-i' indicates the maximum iterations in EM algorithm to partition the
reads. '-w' is the size of the sliding window to scan AMRs. '-m' is the
minimum reads for each CpG. '-b' means using the BIC criterion to compare the
two likelihood models. 'c' indicates the chromosome directory. The output in
"Human_test_AMR.bed" is in the BED format and the score column is either the
p-value from likelihood ratio test of the two models (without '-b') or the
difference of the BIC values of the two models (with '-b').

amrfinder_BAM is specially designed to take BAM input files. It is required to
intall the "BAMTOOLS"(https://github.com/pezmaster31/bamtools) to use the
program.

amrrefiner
============================================================================
The following command detects AMRs with accurate boundaires in several genomic
intervals around the gene GNAS defined by "GNAS_cpg.bed", which denotes the
numbering orders of CpGs in these intervals:

$$ ./amrrefiner -o Human_test_AMR_refined.bed -s 10 -m 100 -M 400 -e 1 \
   -c hg18 -E Human_test.txt GNAS_cpg.bed 

Option '-s' is the minimum size of AMR, '-m' is the mean size of AMR, and '-M'
is the maximum size of AMR. '-e' is the expected number of AMRs in each
interval. '-E' indicates the input read file is in the 'Epiread' format. 

amrtester
============================================================================
The following command test the regions in the input BED file
"GNAS_cpg.bed" are AMRs or not:

$$ ./amrtester -o Human_test_AMRs.bed -c hg18 -E GNAS_cpg.bed Human_test.txt

The output file is in BED format and contains the score column the same
meaning as in 'amrfinder'.
