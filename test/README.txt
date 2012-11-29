Human_test.mr
===================================================================
This file contains reads in "MappedRead' format.

Human_test.txt
===================================================================
This file contains reads in "Epiread" format.

GNAS_cpg.bed
===================================================================
This file contains 4 regions around the gene GNAS in human in the 
format of numbering order of CpGs.

Examples
===================================================================
methstates -c hg18 Human_test.mr -o Human_test.txt

amrfinder -c hg18 Human_test.mr

cpg2base -c hg18 GNAS_cpg.bed
