/*    simasmreal: a program for simulating allele-specific methylation
 *    using real chromosomes and real mapped read locations
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith and Fang Fang
 *
 *    Authors: Fang Fang and Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <iterator>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

using std::tr1::unordered_map;


inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}


static bool
check_sorted(const vector<MappedRead> &mr) {
  for (size_t i = 1; i < mr.size(); ++i)
    if (mr[i].r < mr[i-1].r) return false;
  return true;
}


static void
collect_cpgs(const string &chrom, const string &seq,
	     unordered_map<size_t, size_t> &cpg_lookup,
	     vector<GenomicRegion> &cpgs) {
  size_t count = 0;
  for (size_t i = 0; i < seq.length() - 1; ++i)
    if (is_cpg(seq, i)) {
      cpgs.push_back(GenomicRegion(chrom, i, i+1, "CpG", 0, '+'));
      cpg_lookup[i] = count;
      ++count;
    }
}


static void
simulate_methylation(const gsl_rng *rng, const size_t n_cpgs, 
		     const double alpha, const double beta, 
		     const double meth_var, vector<double> &meth) {
  const double adjusted_alpha = alpha*meth_var;
  const double adjusted_beta = beta*meth_var;
  for (size_t i = 0; i < n_cpgs; ++i)
    meth.push_back(gsl_ran_beta(rng, adjusted_alpha, adjusted_beta));
}


static void
fix_non_allelic_meth(const gsl_rng *rng, const vector<GenomicRegion> &amrs,
		     const vector<GenomicRegion> &cpgs,
		     vector<double> &meth1, vector<double> &meth2,
		     vector<bool> &allelic_id) {
  allelic_id.resize(cpgs.size(), true);
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    while (j < amrs.size() && amrs[j].get_end() < cpgs[i].get_end()) ++j;
    if (j < amrs.size() && !amrs[j].overlaps(cpgs[i])) {
      allelic_id[i] = false;
      if (gsl_ran_flat(rng, 0, 1) > 0.5) meth1[i] = meth2[i];
      else meth2[i] = meth1[i];
    }
  }
}


static void
fake_the_allelic_meth(const gsl_rng *rng, const vector<GenomicRegion> &amrs,
		      const vector<GenomicRegion> &cpgs,
		      vector<double> &meth1, vector<double> &meth2,
		      vector<bool> &allelic_id) {
  allelic_id.resize(cpgs.size(), true);
  for (size_t i = 0; i < cpgs.size(); ++i) {
    allelic_id[i] = false;
    if (gsl_ran_flat(rng, 0, 1) > 0.5) meth1[i] = meth2[i];
    else meth2[i] = meth1[i];
  }
}


static void
mix_the_patterns(const gsl_rng *rng, 
		 vector<double> &meth1, vector<double> &meth2) {
  for (size_t i = 0; i < meth1.size(); ++i)
    if (gsl_ran_flat(rng, 0, 1) > 0.5) std::swap(meth1[i],  meth2[i]);
}


static char
sample_from_methylation(const gsl_rng *rng,
			const vector<double> &meth, const size_t cpg_idx) {
  return (gsl_ran_flat(rng, 0.0, 1.0) < meth[cpg_idx]) ? 'C' : 'T';
}


// adjust the methylation indicated in the reads
static void
simulate_reads(const gsl_rng *rng, const double allelic_bias,
	       const string &seq, 
	       const unordered_map<size_t, size_t> &cpg_lookup,
	       const vector<double> &meth1, const vector<double> &meth2, 
	       vector<MappedRead> &reads, vector<size_t> &coverage) {
  
  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t readlen = reads[i].seq.length();
    reads[i].seq = string(readlen, 'N');
    reads[i].r.set_strand('+');
    const bool is_allele1 = gsl_ran_flat(rng, 0, 1) > allelic_bias;
    reads[i].r.set_name(reads[i].r.get_name() + ":" + toa(is_allele1));
    
    const size_t seq_idx = reads[i].r.get_start();
    size_t cpg_count = 0;
     for (size_t j = 0; j < reads[i].r.get_width(); ++j)
       if (is_cpg(seq, seq_idx + j)) {
	 cpg_count += 1;
	 const unordered_map<size_t, size_t>::const_iterator l = 
	   cpg_lookup.find(seq_idx + j);
	 assert(l != cpg_lookup.end());
	 const size_t cpg_idx = l->second;
	 reads[i].seq[j] = sample_from_methylation(rng, is_allele1 ? meth1 : meth2, cpg_idx);
	 coverage[cpg_idx]++;
       }
    reads[i].r.set_score(cpg_count);
  }
}


int 
main(int argc, const char **argv) {
  
  try {
    
    //   Main Inputs:
    // 	(1) Chromosome file (only a single one is needed for this program)
    // 	(2) Reads file (mapped reads from the chromosome listed above)
    // 	(3) Set of simulated AMR locations
    
    // Outputs:
    // 	(1) Full set of CpGs from the chromosome (their locations are
    // 	    determined by the chromosome sequence from the input) with names:
    // 	    (a) Indicating if the CpG is in a simulated AMR or not
    // 	    (b) 2 methylation levels, one for each allele, and the two levels identical if the CpG is not in an AMR
    // 	    (c) The coverage at that CpG according to the reads passed in as input
    // 	(2) Set of reads modified from the input reads with:
    // 	    (a) [**** NOT IMPLEMENTED] Read names indicating whether or not the read falls within a simulated AMR
    // 	    (b) Name also indicating which allele the read belongs to (and always allele 1 if not in AMR)
    // 	    (c) Every position in the read sequences has its CpG
    // 	        methylation state re-sampled from the simulated
    // 		methylation level of its allele
    
    
    bool VERBOSE = false;
    bool GENERATE_NON_ASM_DATA = false;
    bool GENERATE_PATTERNED_ASM = false;
    
    string methoutfile;
    string readsoutfile;
    
    double alpha1 = 0.25, beta1 = 0.75;
    double alpha2 = 0.75, beta2 = 0.25;
    double meth_var = 1.0;
    
    double allelic_bias = 0.5;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "simulate ASMs with real reads and sequences", 
			   "<chrom> <reads> <amrs>");
    opt_parse.add_opt("methout", 'm', "meth output file", true, methoutfile);
    opt_parse.add_opt("readsout", 'r', "reads output file", true, readsoutfile);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("fake", 'F', "make fake amrs", false, GENERATE_NON_ASM_DATA);
    opt_parse.add_opt("pattern", 'P', "make patterned amrs", false, GENERATE_PATTERNED_ASM);
    opt_parse.add_opt("bias", 'b', "allelic bias [0,1]", false, allelic_bias);
    opt_parse.add_opt("var", 'V', "methylation variance (lower means more variance)", 
		      false, meth_var);
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 3) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string chromfile(leftover_args[0]);
    const string readsfile(leftover_args[1]);
    const string amrsfile(leftover_args[2]);
    /****************** END COMMAND LINE OPTIONS *****************/
    
    // load the mapped reads (assuming only CpG locations used)
    if (VERBOSE)
      cerr << "LOADING MAPPED READS" << endl;
    vector<MappedRead> reads;
    LoadMappedReadsFile(readsfile, reads);
    check_sorted(reads);
    if (VERBOSE)
      cerr << "READS:\t" << reads.size() << endl;
    
    // load the AMRs to simulate
    vector<GenomicRegion> amrs;
    ReadBEDFile(amrsfile, amrs);
    check_sorted(amrs);
    if (VERBOSE)
      cerr << "AMRS:\t" << amrs.size() << endl;
    
    // load the chromosome
    vector<string> chrom_names, sequences;
    read_fasta_file(chromfile.c_str(), chrom_names, sequences);
    assert(chrom_names.size() == 1);
    if (VERBOSE)
      cerr << "CHROM:\t" << chrom_names.front() << endl
	   << "LENGTH:\t" << sequences.front().length() << endl;

    // determine CpG locations and make a lookup table
    unordered_map<size_t, size_t> cpg_lookup;
    vector<GenomicRegion> cpgs;
    collect_cpgs(chrom_names.front(), sequences.front(), cpg_lookup, cpgs);
    const size_t n_cpgs = cpgs.size();
    if (VERBOSE)
      cerr << "CPGS:\t" << n_cpgs << endl;
    
    // setup the random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(0) + getpid());
    
    // simulate the methylation profiles along the full chromosomes
    vector<double> meth1, meth2;
    simulate_methylation(rng, n_cpgs, alpha1, beta1, meth_var, meth1);
    simulate_methylation(rng, n_cpgs, alpha2, beta2, meth_var, meth2);

    if (VERBOSE)
      cerr << "SIMUALTED METHYLATION" << endl;
    
    // make the non-AMR parts have a single shared methylation pattern
    vector<bool> allelic_id;
    fix_non_allelic_meth(rng, amrs, cpgs, meth1, meth2, allelic_id);
    
    if (VERBOSE)
      cerr << "FIXED NON-ALLELIC" << endl;
    
    if (GENERATE_NON_ASM_DATA)
      fake_the_allelic_meth(rng, amrs, cpgs, meth1, meth2, allelic_id);
    else if (GENERATE_PATTERNED_ASM) 
      mix_the_patterns(rng, meth1, meth2);
    
    // adjust the methylation indicated in the reads
    vector<size_t> coverage(n_cpgs, 0ul);
    simulate_reads(rng, allelic_bias, sequences.front(), cpg_lookup, 
		   meth1, meth2, reads, coverage);
    
    if (VERBOSE)
      cerr << "WRITING SIMULATED METH PROFILES" << endl;
    std::ofstream methout(methoutfile.c_str());
    for (size_t i = 0; i < cpgs.size(); ++i) {
      if (coverage[i] > 0) {
	cpgs[i].set_name("CpG:" + toa(meth1[i]) + ":" + 
			 toa(meth2[i]) + ":" + toa(allelic_id[i]));
	cpgs[i].set_score(coverage[i]);
	methout << cpgs[i] << endl;
      }
    }
    
    if (VERBOSE)
      cerr << "WRITING SIMULATED READS" << endl;
    std::ofstream readsout(readsoutfile.c_str());
    copy(reads.begin(), reads.end(), 
	 std::ostream_iterator<MappedRead>(readsout, "\n"));
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
