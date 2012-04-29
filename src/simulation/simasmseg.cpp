/*    simasmseg: a program for simulating a segmented sequence of AMRs
 *    and non-AMRs
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
#include <numeric>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

static size_t
sample_from_histogram(const gsl_rng *rng, vector<double> cumul) {
  return (lower_bound(cumul.begin(), cumul.end(), gsl_ran_flat(rng, 0, 1)) - 
	  cumul.begin());
}

static void
read_distribution(const string &filename, vector<double> &probs) {
  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("cannot open input file " + filename);
  double val = 0.0, count = 0.0;
  while (in >> val >> count)
    probs.push_back(count);
  const double total = accumulate(probs.begin(), probs.end(), 0.0);
  transform(probs.begin(), probs.end(), probs.begin(),
	    std::bind2nd(std::divides<double>(), total));
  std::partial_sum(probs.begin(), probs.end(), probs.begin());
}


static size_t
sample_cpg_distance(const gsl_rng *rng, const double mean_cpg_dist) {
  return static_cast<size_t>(gsl_ran_gamma(rng, 2.0, mean_cpg_dist));
}

static size_t
simulate_cpg_locations(const gsl_rng *rng, const size_t readlen,
		       const size_t n_cpgs, const double mean_cpg_dist, 
		       vector<SimpleGenomicRegion> &cpgs) {
  const string chrom("chr1");
  size_t offset = readlen;
  cpgs.push_back(SimpleGenomicRegion(chrom, offset, offset + 1));
  for (size_t i = 1; i < n_cpgs; ++i) {
    const size_t location = offset + sample_cpg_distance(rng, mean_cpg_dist);
    cpgs.push_back(SimpleGenomicRegion(chrom, location, location + 1));
    offset = location;
  }
  return offset + readlen;
}

static size_t
simulate_cpg_locations(const gsl_rng *rng, const size_t readlen,
		       const size_t n_cpgs, const vector<double> &dist_probs,
		       vector<SimpleGenomicRegion> &cpgs) {
  const string chrom("chr1");
  size_t offset = readlen;
  cpgs.push_back(SimpleGenomicRegion(chrom, offset, offset + 1));
  for (size_t i = 0; i < n_cpgs; ++i) {
    const size_t location = offset + sample_from_histogram(rng, dist_probs);
    cpgs.push_back(SimpleGenomicRegion(chrom, location, location + 1));
    offset = location;
  }
  return offset + readlen;
}


static void
simulate_methylation(const gsl_rng *rng, const size_t n_cpgs, 
		     const double alpha, const double beta, 
		     const double meth_var, vector<double> &meth) {
  const double adjusted_alpha = alpha*100.0/meth_var;
  const double adjusted_beta = beta*100.0/meth_var;
  for (size_t i = 0; i < n_cpgs; ++i)
    meth.push_back(gsl_ran_beta(rng, adjusted_alpha, adjusted_beta));
}


static MappedRead
make_read(const string &chrom, const size_t st, 
	  const bool allele, const string &seq) {
  MappedRead mr;
  mr.r = GenomicRegion(chrom, st, st + seq.length(), "A" + toa(allele), 0.0, '+');
  mr.seq = seq;
  mr.scr = string(seq.length(), 'B');
  return mr;
}


static void
simulate_read_locations(const gsl_rng *rng, const size_t chrom_len, 
			const size_t n_reads, const size_t readlen, 
			vector<size_t> &starts) {
  const size_t max_pos = chrom_len - readlen + 1;
  starts.resize(n_reads);
  for (size_t i = 0; i < n_reads; ++i) 
    starts[i] = gsl_ran_flat(rng, 0, max_pos);
  sort(starts.begin(), starts.end());
}


static void
simulate_reads(const gsl_rng *rng, 
	       const size_t chrom_len, const size_t readlen, 
	       const size_t n_reads, const vector<SimpleGenomicRegion> &cpgs, 
	       const vector<double> &meth1, const vector<double> &meth2, 
	       vector<MappedRead> &reads, vector<size_t> &coverage) {
  
  // this variable below assumes only 1 chromosome...
  const string chrom(cpgs.front().get_chrom());
  
  // simulate start positions
  vector<size_t> starts;
  simulate_read_locations(rng, chrom_len, n_reads, readlen, starts);
  
  // initialize the vector of coverage values to 0 coverage per CpG
  coverage.resize(cpgs.size(), 0ul);
  
  size_t cpg_idx = 0;
  for (size_t i = 0; i < starts.size(); ++i) {
    
    // move to the first CpG after the start of the current read
    while (cpg_idx < cpgs.size() && cpgs[cpg_idx].get_start() < starts[i])
      ++cpg_idx;
    
    // determine which allele to simulate
    const bool allele(gsl_ran_flat(rng, 0, 1) > 0.5);

    // flag to indicate if a CpG is covered by this read
    bool contains_cpg = false;
    
    string seq; 

    // process each CpG covered by the current read
    const size_t lim = starts[i] + readlen;
    for (size_t j = cpg_idx; j < cpgs.size() && cpgs[j].get_start() < lim; ++j) {
      const double coin_flip = gsl_ran_flat(rng, 0.0, 1.0);
      const char the_base = (coin_flip < (allele ? meth1[j] : meth2[j])) ? 'C' : 'T';
      seq += the_base;
      coverage[j] += 1;
      contains_cpg = true;
    }
    const size_t current_start = cpg_idx;
    if (contains_cpg)
      reads.push_back(make_read(chrom, current_start, allele, seq));
  }
}


int 
main(int argc, const char **argv) {
  
  try {
    
    size_t n_cpgs = 100;
    size_t n_reads = 100;
    size_t readlen = 100;
    size_t mean_cpg_dist = 30;
    
    double alpha1 = 0.15, beta1 = 0.85;
    double alpha2 = 0.85, beta2 = 0.15;
    double meth_var = 1.0;
    
    bool VERBOSE = false;

    string regions_out_name;
    string alleles_out_name;
    string reads_out_name;
    string inter_cpg_dist_file;

    size_t exp_amrs = 1;
    size_t mean_amr_size = 20;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "sim AMR and corresponding reads");
    opt_parse.add_opt("regions", 'R', "regions out file", true, regions_out_name);
    opt_parse.add_opt("alleles", 'a', "alleles out file", true, alleles_out_name);
    opt_parse.add_opt("alpha", 'l', "alpha value for beta distro", false, alpha1);
    opt_parse.add_opt("reads", 'r', "reads out file", true, reads_out_name);
    opt_parse.add_opt("n_cpgs", 'c', "number of cpgs", false, n_cpgs);
    opt_parse.add_opt("dist", 'd', "mean inter-cpg distance", false, 
		      mean_cpg_dist);
    opt_parse.add_opt("distfile", 'D', "inter-cpg distances", false, 
		      inter_cpg_dist_file);
    opt_parse.add_opt("width", 'w', "width of reads", false , readlen);
    opt_parse.add_opt("n_reads", 'n', "number of reads", false, n_reads);
    opt_parse.add_opt("var", 'V', "variance of methylation", false, meth_var);
    opt_parse.add_opt("mean-size", 'm', "mean AMR size", false, mean_amr_size);
    opt_parse.add_opt("expected", 'e', "expected number of AMRs", false, exp_amrs);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    
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
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "CPGs:    \t" << n_cpgs << endl
	   << "CpG dist:\t" << mean_cpg_dist << endl
	   << "width:   \t" << readlen << endl
	   << "reads:   \t" << n_reads << endl;
    
    // setup the random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    srand(time(0) + getpid());
    gsl_rng_set(rng, rand());
    
    // get the distribution for inter-CpG distances
    vector<double> dist_probs;
    if (!inter_cpg_dist_file.empty())
      read_distribution(inter_cpg_dist_file, dist_probs);
    
    // simulate the CpG locations
    vector<SimpleGenomicRegion> cpgs;
    const size_t chrom_len = (dist_probs.empty()) ?
      simulate_cpg_locations(rng, readlen, n_cpgs, mean_cpg_dist, cpgs) :
      simulate_cpg_locations(rng, readlen, n_cpgs, dist_probs, cpgs);
    if (VERBOSE)
      cerr << "chrom len:\t" << chrom_len << endl;
    
    // setup the parameters for simulation
    alpha2 = beta1 = 1.0 - alpha1;
    beta2 = alpha1;
    
    // simulate the methylation profiles
    vector<double> meth1, meth2;
    simulate_methylation(rng, n_cpgs, alpha1, beta1, meth_var, meth1);
    simulate_methylation(rng, n_cpgs, alpha2, beta2, meth_var, meth2);

    // sample a number of AMRs
    const size_t n_amrs = std::max(1u, gsl_ran_poisson(rng, exp_amrs));
    if (VERBOSE)
      cerr << "amrs:\t" << n_amrs << endl;
    
    // sample AMR sizes
    vector<size_t> amr_sizes;
    for (size_t i = 0; i < n_amrs; ++i)
      amr_sizes.push_back(gsl_ran_poisson(rng, mean_amr_size));
    
    // sample AMR positions
    vector<size_t> amr_positions;
    for (size_t i = 0; i < n_amrs; ++i)
      amr_positions.push_back(gsl_ran_flat(rng, 0, n_cpgs - mean_amr_size));
    sort(amr_positions.begin(), amr_positions.end());
    
    // create the inter-AMR methylation profiles
    size_t curr_amr = 0;
    for (size_t i = 0; i < n_cpgs; ++i) {
      if (curr_amr == n_amrs || i < amr_positions[curr_amr]) {
	if (gsl_ran_flat(rng, 0, 1) > 0.5) meth1[i] = meth2[i];
	else meth2[i] = meth1[i];
      }
      else if (i == amr_positions[curr_amr] + amr_sizes[curr_amr] - 1)
	++curr_amr;
    }
    
    
    // simulate the reads
    vector<MappedRead> reads;
    vector<size_t> coverage;
    simulate_reads(rng, chrom_len, readlen, n_reads, cpgs, 
		   meth1, meth2, reads, coverage);
    
    if (VERBOSE)
      cerr << "WRITING READS OUTPUT" << endl;
    std::ofstream reads_out(reads_out_name.c_str());
    for (size_t i = 0; i < reads.size(); ++i) {
      reads[i].r.set_name(reads[i].r.get_name() + ":" + reads[i].seq);
      reads_out << reads[i].r.get_chrom() << '\t'
		<< reads[i].r.get_start() << '\t'
		<< reads[i].seq << endl;
    }
    
    if (VERBOSE)
      cerr << "WRITING ALLELES OUTPUT" << endl;
    std::ofstream alleles_out(alleles_out_name.c_str());
    for (size_t i = 0; i < cpgs.size(); ++i) {
      cpgs[i].set_start(i);
      cpgs[i].set_end(i+1);
      alleles_out << std::setprecision(3) << std::fixed 
		  << cpgs[i] << "\tCpG:" 
 		  << meth1[i] << ":" << meth2[i] << "\t"
		  << coverage[i] << "\t+" << endl;
    }
    
    if (VERBOSE)
      cerr << "WRITING REGIONS OUTPUT" << endl;
    std::ofstream regions_out(regions_out_name.c_str());
    regions_out << reads.front().r.get_chrom() << '\t'
		<< reads.front().r.get_start() << '\t'
		<< reads.back().r.get_start() 
		<< "\tX\t0\t+" << endl;
    
    for (size_t i = 0; i < n_amrs; ++i)
      cout << "chr1" << "\t" 
	   << amr_positions[i] << "\t"
	   << amr_positions[i] + amr_sizes[i] << endl;
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
