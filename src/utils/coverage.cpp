/*    coverage: a program for counting the coverage
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith and Elena Harris and Song Qiang
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
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "QualityScore.hpp"

#include "FileIterator.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::ofstream;


static void
add_contribution(const size_t offset, const GenomicRegion &r,
		   size_t &count) {  
    const size_t position = offset - r.get_start();
    //    assert(position < r.seq.length());
    if (position >= r.get_width())
      return;//throw SMITHLABException("ERROR: Reads must be sorted by chromosome and end position.");

    ++count;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static bool
precedes(const GenomicRegion &r, const size_t offset) {
  return r.get_end() <= offset;
}


static bool
succeeds(const GenomicRegion &r, const size_t offset) {
  return r.get_start() > offset;
}


static void
advance(const size_t first, const size_t last,
		const GenomicRegion &chrom_region, 
		FileIterator<GenomicRegion> &regions,
        const size_t max_length) {
  while (regions.last_is_good() && 
		 chrom_region.same_chrom(*regions.get_last()) &&
		 !succeeds(*regions.get_last(), last+max_length)) {
	  	 	 regions.increment_last();
  	  }

  //   if (regions.last_is_good() != reads.last_is_good())
  //     throw SMITHLABException("read and map files seem out of sync");
  while (regions.first_is_good() && 
	 chrom_region.same_chrom(*regions.get_first()) &&
	 precedes(*regions.get_first(), first)) {
    regions.increment_first();
  }

  //   if (regions.first_is_good() != reads.first_is_good())
  //     throw SMITHLABException("read and map files seem out of sync");
}


static void
scan_chromosome(const size_t &chrom_len, const GenomicRegion &chrom_region,
				FileIterator<GenomicRegion> &regions, 
				std::ostream &out,
                const size_t max_length) {
  
  const string chrom_name(chrom_region.get_chrom());
  size_t i = 0; 
  for (i = 0; i < chrom_len - 1 && regions.first_is_good(); ++i) {
	  advance(i, i, chrom_region, regions, max_length);
	  size_t count = 0;
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
    		  j != regions.get_last(); ++j)
    	  add_contribution(i, *j, count);
      if(count>0)
    	  out << chrom_name << "\t" << i << "\t"  
    		  << i+1 << "\t" << count << endl;
  }
}


static void
read_chrom_lengths(const string chrom_file,
					vector<string> &chrom_names,
					vector<size_t> &chrom_lengths) {
	std::ifstream in(chrom_file.c_str());
	string line, name;
	size_t length;
	while (std::getline(in, line))
	{
		std::istringstream iss(line);
		iss >> name >> length;
		chrom_names.push_back(name);
		chrom_lengths.push_back(length);
	}
}


static void
advance_chromosome(const GenomicRegion &chrom_region, 
		   FileIterator<GenomicRegion> &regions) {

  while (regions.last_is_good() && 
	 (*regions.get_last() < chrom_region)) {
    assert(regions.last_is_good());
    regions.increment_last();
  }
  while (regions.first_is_good() && 
	 (*regions.get_first() < chrom_region))
    regions.increment_first();
}


static void
scan_chroms(const bool VERBOSE,
			const string &outfile, const vector<string> &chrom_names,
			const vector<size_t> &chrom_lengths,
			FileIterator<GenomicRegion> &regions,
            const size_t max_length) {

  // Get the output stream
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  for (size_t i = 0; i < chrom_names.size(); ++i) {
      if (VERBOSE) cerr << "[SCANNING=" << chrom_names[i] << "]";
      //TODO: WHAT HAPPENS IF A CHROM IS MISSING??
      const GenomicRegion chrom_region(chrom_names[i], 0, 0);
      advance_chromosome(chrom_region, regions);
     
      scan_chromosome(chrom_lengths[i], chrom_region,
			regions, out, max_length);
    if (VERBOSE) cerr << " [DONE]" << endl;
  }
}



int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;    
    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";
    
    size_t BUFFER_SIZE = 100000;
    size_t max_length = 10000;
    string out_stat;

    double cutoff = -std::numeric_limits<double>::max();
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Program to count coverages",
			   "<mapped-reads>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("chrom", 'c', "file of chrom lengths",
		      true , chrom_file);
    opt_parse.add_opt("max_length", 'L', "The maximal read length of the input file (default 10000)",
                      false , max_length);
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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/    
    vector<string> chrom_names;
    vector<size_t> chrom_lengths;
    read_chrom_lengths(chrom_file,chrom_names, chrom_lengths);
    
    FileIterator<GenomicRegion> regions(mapped_reads_file, BUFFER_SIZE);
    scan_chroms(VERBOSE, outfile, chrom_names, chrom_lengths, regions, max_length);
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
