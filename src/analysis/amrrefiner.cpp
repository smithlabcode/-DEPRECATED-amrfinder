/*    amrrefiner: A program for refining AMRs identified with
 *    amrfinder
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

#include "Epiread.hpp"
#include "EpireadStats.hpp"
#include "EpireadIO.hpp"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::pair;
using std::numeric_limits;

////////////////////////////////////////////////////////////////////////
// EPIREAD CODE

static void
get_read_interval_endpoints(const vector<epiread> &reads, 
			    const size_t n_cpgs,
			    vector<size_t> &start_read, 
			    vector<size_t> &end_read) {
  for (size_t i = 0; i < n_cpgs; ++i) {
    size_t j = (start_read.empty()) ? 0 : start_read.back();
    while (j < reads.size() && reads[j].pos + reads[j].seq.length() <= i) 
      ++j;
    start_read.push_back(j);
  }
  for (size_t i = n_cpgs; i > 0; --i) {
    size_t j = (end_read.empty()) ? reads.size() : end_read.back();
    while (j > 0 && reads[j-1].pos >= i) 
      --j;
    end_read.push_back(j);
  }
  reverse(end_read.begin(), end_read.end());
}


static double
get_single_allele_score(const size_t start, const size_t end,
			const vector<epiread> &reads) {
  vector<double> a0(end - start, 0.5);
 return fit_single_epiallele(0, reads.size(), start, end, reads, a0);
}


static double
get_pair_score(const size_t max_itr, 
	       const size_t read_start, const size_t read_end,
	       const size_t start, const size_t end, const vector<epiread> &reads,
	       vector<double> &a1, vector<double> &a2, vector<double> &indicators) {

  static const double MIN_INDICATOR_SUM = 1e-10;
  resolve_epialleles(read_start, read_end, start, end, 
		     max_itr, reads, indicators, a1, a2);
  const double indic_sum = accumulate(indicators.begin() + read_start, 
					   indicators.begin() + read_end, 
					   MIN_INDICATOR_SUM);
  
  assert(indic_sum < (read_end - read_start) && indic_sum > 0.0);
  const double balance_correction_observed = 
    gsl_sf_lnbeta(indic_sum, float(read_end - read_start) - indic_sum);
  
  const double balance_correction_optimal = 
    gsl_sf_lnbeta((read_end - read_start)/2.0, (read_end - read_start)/2.0);
  
  const double balance_correction = read_end > read_start ? 
    -(balance_correction_observed - balance_correction_optimal) : 0;
  
  const double score = 
    log_likelihood(read_start, read_end, start, end, reads, indicators, a1, a2) 
    // BELOW IS AN APPROXIMATION -- NOT ALL READS BETWEEN "LAST" AND "FIRST" ARE USED
    + 0.05*balance_correction
    + (read_end - read_start)*log(0.5);
  
  /* WE COULD ALTERNATIVELY ACTUALLY COUNT THE READS, WHICH WOULD BE
     SLOWER... */
  // get_reads_count(start, end, reads)*log(0.5);
  return score;
}


static void
update_single_tables(const double inc_score, const double trans_score, 
		     const vector<double> &single_cpg_scores,
		     const vector<double> &pair_scores, const size_t i,
		     vector<double> &single_scores, vector<size_t> &single_lookback) {
  const double using_prev_single = single_scores[i - 1] + inc_score;
  const double using_prev_pair = pair_scores[i - 1] + trans_score;
  if (using_prev_single > using_prev_pair) {
    single_scores[i] = single_cpg_scores[i] + using_prev_single;
    single_lookback[i] = single_lookback[i - 1];
  }
  else {
    single_scores[i] = single_cpg_scores[i] + using_prev_pair;
    single_lookback[i] = i - 1;
  }
}


static void
update_pair_tables(const vector<double> &duration_probs,
		   const double trans_score,  
		   const size_t min_amr_size, const double mean_amr_size,
		   const size_t max_amr_size,
		   const size_t max_itr,
		   const vector<size_t> &start_read,
		   const vector<size_t> &end_read,
		   const vector<epiread> &reads,
		   const vector<double> &single_scores,
		   const size_t current_position,
		   vector<double> a1, vector<double> a2, 
		   vector<double> &indicators, 
		   vector<double> &pair_scores,
		   vector<size_t> &pair_lookback) {
  
  for (size_t amr_size = std::min(current_position, max_amr_size); 
       amr_size > min_amr_size; --amr_size) {
    const size_t amr_start = current_position - amr_size;
    
    copy(a1.begin() + 1, a1.end(), a1.begin());
    copy(a2.begin() + 1, a2.end(), a2.begin());
    
    const double region_score = (start_read[amr_start] < end_read[current_position]) ? 
      get_pair_score(max_itr, start_read[amr_start], 
		     end_read[current_position], amr_start, current_position, 
		     reads, a1, a2, indicators) : 
      -std::numeric_limits<double>::max();
    
    const double prev_single = (current_position > amr_size) ? 
      single_scores[amr_start - 1] : 0;
    const double pair_score = region_score + trans_score + 
      prev_single + duration_probs[amr_size];
    
    //     const double bic_pair = 
    //       2*(current_position - amr_start)*log(end_read[current_position] - start_read[amr_start]) - 
    //       2*(region_score + (current_position - amr_start)*log(0.5));
    
    //     const double region_score_single =
    //       get_single_allele_score(amr_start, current_position, reads);
    
    //     const double bic_single = 
    //       (current_position - amr_start)*log(end_read[current_position] - 
    //     start_read[amr_start]) - 
    //       2*region_score_single;
    
    if (pair_score > pair_scores[current_position]) { // && bic_pair < bic_single) {
      pair_scores[current_position] = pair_score;
      pair_lookback[current_position] = amr_start - 1;
    }
  }
}


////////////////////////////////////////////////////////////////////////
// CODE TO GET THE SEGMENTATION AFTER THE VALUES HAVE BEEN COMPUTED


static void
traceback(const vector<double> &single_scores, const vector<double> &pair_scores,
	  const vector<size_t> &single_lookback, const vector<size_t> &pair_lookback,
	  vector<pair<size_t, size_t> > &amrs) {

  int pos = single_scores.size() - 1;
  cerr << single_scores[pos] << "\t" << pair_scores[pos] << endl;
  while (pos >= 0) { 
    if (single_scores[pos] > pair_scores[pos])
      pos = single_lookback[pos];
    else {
      const size_t curr = pos;
      pos = pair_lookback[pos];
      amrs.push_back(std::make_pair(pos, curr - pos));
    }
  }
  reverse(amrs.begin(), amrs.end());
}


static void
dynamic_programming_segmentation(const bool VERBOSE, const bool PROGRESS,
				 const vector<epiread> &reads, const size_t max_itr,
				 const double low_prob, const double high_prob,
				 const vector<double> &duration_probs,
				 const size_t min_amr_size,
				 const double mean_amr_size,
				 const size_t max_amr_size, 
				 const double exp_amrs,
				 vector<double> &single_scores, 
				 vector<size_t> &single_lookback,
				 vector<double> &pair_scores, 
				 vector<size_t> &pair_lookback) {
  
  const size_t n_cpgs = get_n_cpgs(reads);
  
  single_scores = pair_scores = 
    vector<double>(n_cpgs, -numeric_limits<double>::max());
  single_lookback = pair_lookback = 
    vector<size_t>(n_cpgs, numeric_limits<size_t>::max());
  
  if (VERBOSE)
    cerr << "[STATUS] finding read interval endpoints" << endl;
  vector<size_t> start_read, end_read;
  get_read_interval_endpoints(reads, n_cpgs, start_read, end_read);
  
  if (VERBOSE)
    cerr << "[STATUS] computing single allele scores" << endl;
  vector<double> single_cpg_scores;
  for (size_t i = 0; i < n_cpgs; ++i)
    single_cpg_scores.push_back(get_single_allele_score(i, i + 1, reads));
  
  // These below are the log values used in computations related to
  // the exponential duration of the background:
  const double inc_score = log(1.0 - exp_amrs); // score for incrementing
  const double trans_score = log(exp_amrs); // transition score
  
  // Allocate the space for vectors needed in much of the
  // computation -- don't want to have to keep re-allocating
  // these...
  vector<double> indicators(reads.size(), 0.0);
  vector<double> a1(n_cpgs, low_prob), a2(n_cpgs, high_prob);
  
  // Initialize all the tables for the base-case
  single_scores[0] = single_cpg_scores.front() + trans_score;
  pair_scores[0] = duration_probs[1] +
    get_pair_score(max_itr, start_read[0], end_read[1], 0, 1, 
		   reads, a1, a2, indicators);
  
  if (VERBOSE)
    cerr << "[STATUS] filling dynamic programming tables" << endl;
  for (size_t i = 1; i < n_cpgs; ++i) {
    if (PROGRESS) 
      cerr << '\r' << percent(i, n_cpgs) << "%\r";
    
    update_single_tables(inc_score, trans_score, single_cpg_scores,
			 pair_scores, i, single_scores, single_lookback);
    update_pair_tables(duration_probs, trans_score, min_amr_size, mean_amr_size,
		       max_amr_size, max_itr, start_read, end_read, reads, 
		       single_scores, i, a1, a2, indicators, pair_scores, 
		       pair_lookback);
  }
  if (PROGRESS) cerr << "\r100%" << endl;
}


static void
get_amrs(const GenomicRegion &region, const size_t first_read_offset, 
	 const vector<pair<size_t, size_t> > &segments, 
	 vector<GenomicRegion> &amrs) {
  const string chrom(region.get_chrom());
  for (size_t i = 0; i < segments.size(); ++i) {
    const size_t start = segments[i].first + first_read_offset;
    const size_t end = start + segments[i].second;
    amrs.push_back(GenomicRegion(chrom, start, end));
  }
}


static void
expand_regions(const size_t expansion_size, vector<GenomicRegion> &regions) {
  for (size_t i = 0; i < regions.size(); ++i) {
    regions[i].set_start((regions[i].get_start() > expansion_size) ? 
			 regions[i].get_start() - expansion_size : 0);
    regions[i].set_end(regions[i].get_end() + expansion_size);
  }
}


int
main(int argc, const char **argv) {
  
  try {

    bool EPIREAD_FORMAT = false;
    bool VERBOSE = false;
    bool PROGRESS = false;
    string outfile;
    string chroms_dir;

    size_t max_itr = 5;
    double high_prob = 0.75, low_prob = 0.25;
    
    size_t max_amr_size = 400;
    size_t min_amr_size = 10;
    double mean_amr_size = 100.0;
    double exp_amrs = 1.0/10000.0;
    size_t expansion_size = 0;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "more precise resolution of "
			   "epi-alleles", "<mapped-reads> <regions>");
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_itr);
    opt_parse.add_opt("size", 's', "min AMR size", false, min_amr_size);
    opt_parse.add_opt("mean-size", 'm', "mean AMR size", false, mean_amr_size);
    opt_parse.add_opt("max-size", 'M', "max AMR size", false, max_amr_size);
    opt_parse.add_opt("expected", 'e', "expected number of AMRs", false, exp_amrs);
    opt_parse.add_opt("epiread", 'E', "reads in epiread format", false, EPIREAD_FORMAT);
    opt_parse.add_opt("expand", 'x', "bases to expand regions", false, expansion_size);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("progress", 'P', "print progress info", false, PROGRESS);
    opt_parse.add_opt("chrom", 'c', "dir of chroms (.fa extn)", false, chroms_dir);
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested() || leftover_args.size() != 2) {
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
    const string reads_file_name(leftover_args.front());
    const string regions_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE) {
      cerr << "min amr size:\t" << min_amr_size << endl;
      cerr << "mean amr size:\t" << mean_amr_size << endl;
      cerr << "max amr size:\t" << max_amr_size << endl;
    }

    // get the distribution for coverages
    vector<double> duration_probs;
    duration_probs.resize(max_amr_size);
    for (size_t i = 0; i < max_amr_size; ++i)
      duration_probs[i] = log(gsl_ran_geometric_pdf(i, 1.0/mean_amr_size));
    
    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!check_sorted(regions))
      throw SMITHLABException("regions not sorted: " + regions_file);
    expand_regions(expansion_size, regions);
    collapse(regions);
    
    vector<GenomicRegion> amrs;
    vector<epiread> the_reads;
    EpireadIO eio(reads_file_name, VERBOSE, EPIREAD_FORMAT, chroms_dir);
    
    for (size_t regions_idx = 0; regions_idx < regions.size(); ++regions_idx) {
      
      if (VERBOSE)
	cerr << "LOADING MAPPED READS" << endl;
      
      vector<epiread> reads;
      eio.load_reads(regions[regions_idx], reads);

      const size_t first_read_offset = adjust_read_offsets(reads);
      const size_t n_cpgs = get_n_cpgs(reads);
      if (VERBOSE)
	cerr << "READS:\t" << reads.size() << endl
	     << "OFFSET:\t" << first_read_offset << endl
	     << "TOTAL CPGS:\t" << n_cpgs << endl;
      
      // Declare the tables (essentially DP tables) that will be filled
      // during the computation
      vector<double> single_scores, pair_scores;
      vector<size_t> single_lookback, pair_lookback;
      
      dynamic_programming_segmentation(VERBOSE, PROGRESS, reads, max_itr, 
				       low_prob, high_prob, duration_probs, 
				       min_amr_size, mean_amr_size, max_amr_size, 
				       exp_amrs, 
				       single_scores, single_lookback, 
				       pair_scores, pair_lookback);
      
      // Trace-back through the DP tables to get the AMRs
      vector<pair<size_t, size_t> > segments;
      traceback(single_scores, pair_scores, single_lookback, pair_lookback, segments);
      
      vector<GenomicRegion> current_region_amrs;
      get_amrs(regions[regions_idx], first_read_offset, segments, 
	       current_region_amrs);
      eio.convert_coordinates(current_region_amrs);
      copy(current_region_amrs.begin(), current_region_amrs.end(),
	   back_inserter(amrs));
    }
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    copy(amrs.begin(), amrs.end(), 
	 std::ostream_iterator<GenomicRegion>(out, "\n"));
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
