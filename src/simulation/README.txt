simasm
=============================================================================
The following command simulates an AMR in with 40 CpGs, reads of width
100, 250 reads and inter-CpG distance distribution take from the
histogram file "inter_cpg_dist_prom.txt":

$$ ./simasm -c 40 -D inter_cpg_dist_prom.txt -w 100 -n 250 -a alleles.txt -r reads.txt

The output in the "reads.txt" file will be in the form of 3 columns:
(1) The chromosome, which will always be chr1; (2) The CpG where the read
starts; and (3) the C/T sequence indicating methylation states at covered reads.


simasmseg
=============================================================================
The following simulates a genomic region which may have multiple
AMRs. The mean size (-m) is 20 CpGs, the number of CpGs to simulate is
(-c) 200, there are (-n) 400 reads having a width (-w) of 100:

$$ ./simasmseg -v -m 20 -c 200 -w 100 -D inter_cpg_dist_prom.txt -n 400 -a alleles.txt -r reads.txt


simasmreal
=============================================================================
The following simulates AMRs in one chromosome using real sequenced reads locations
in that chromosome.

$$ ./simasmreal -m meth.bed -r reads.mr chr10.fa chr10_reads.mr AMRs.bed

The input "AMRs.bed" indicates the locations of simulated AMRs. The output file
"meth.bed" contains the simulated methylation states of all CpGs. The output
file "reads.mr" contain all simulated reads.
