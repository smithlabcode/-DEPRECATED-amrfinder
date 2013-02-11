methstates
============================================================================
The following command converts the reads from "MappedRead" format into
"Epiread" format:

$$ ./methstates -c hg18 Human_test.mr -o Human_test.txt

Option '-c' indicates the chromosome directory. 
The input file "Human_test.mr" is in "MappedRead" format which consists of 8
columns with the first 6 the same as in BED format and the last two columns  
being the sequence and quality score of the read.
The output file "Human_test.txt" is in 'Epiread format which consists of 3 
columns: chromosome, numbering order of the first CpG in the read, and the 
sequence of CpG sites in the read.

cpg2base
============================================================================
The following command converts the locations from numbering order of CpGs to
genomic location:

$$ ./cpg2base -o GNAS_region.bed -c hg18 GNAS_cpg.bed
