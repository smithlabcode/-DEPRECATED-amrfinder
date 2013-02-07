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
#include "MappedRead.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::ofstream;


static void
add_contribution(const size_t offset, const MappedRead &r,
		   size_t &count) {  
    const size_t position = offset - r.r.get_start();
    //    assert(position < r.seq.length());
    if (position >= r.seq.length())
      return;//throw SMITHLABException("ERROR: Reads must be sorted by chromosome and end position.");

    ++count;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static bool
precedes(const MappedRead &r, const size_t offset) {
  return r.r.get_end() <= offset;
}


static bool
succeeds(const MappedRead &r, const size_t offset) {
  return r.r.get_start() > offset;
}


static void
advance(const size_t first, const size_t last,
		const GenomicRegion &chrom_region, 
		FileIterator<MappedRead> &regions,
        const size_t max_length) {
  while (regions.last_is_good() && 
		 chrom_region.same_chrom(regions.get_last()->r) &&
		 !succeeds(*regions.get_last(), last+max_length)) {
	  	 	 regions.increment_last();
  	  }

  //   if (regions.last_is_good() != reads.last_is_good())
  //     throw SMITHLABException("read and map files seem out of sync");
  while (regions.first_is_good() && 
	 chrom_region.same_chrom(regions.get_first()->r) &&
	 precedes(*regions.get_first(), first)) {
    regions.increment_first();
  }

  //   if (regions.first_is_good() != reads.first_is_good())
  //     throw SMITHLABException("read and map files seem out of sync");
}


static void
scan_chromosome(const string &chrom, const GenomicRegion &chrom_region,
				FileIterator<MappedRead> &regions, 
				std::ostream &out,
                const size_t max_length) {
  
  const string chrom_name(chrom_region.get_chrom());
  size_t i = 0; 
  size_t chrom_len = chrom.length();
  for (i = 0; i < chrom_len - 1 && regions.first_is_good(); ++i) {
	  advance(i, i, chrom_region, regions, max_length);
	  size_t count = 0;
      for (vector<MappedRead>::const_iterator j(regions.get_first());
    		  j != regions.get_last(); ++j)
    	  add_contribution(i, *j, count);
      out << chrom_name << "\t" << i << "\t+\t"  
    		  << i+1 << "\t" << count << endl;
  }
}


static void
identify_chromosomes(const bool VERBOSE, const string chrom_file,
		     const string fasta_suffix, vector<string> &chrom_files) {
  if (VERBOSE)
    cerr << "[IDENTIFYING CHROMS] ";
  if (isdir(chrom_file.c_str())) 
    read_dir(chrom_file, fasta_suffix, chrom_files);
  else chrom_files.push_back(chrom_file);
  if (VERBOSE) {
    cerr << "[DONE]" << endl 
	 << "chromosome files found (approx size):" << endl;
    for (vector<string>::const_iterator i = chrom_files.begin();
	 i != chrom_files.end(); ++i)
      cerr << *i << " (" << roundf(get_filesize(*i)/1e06) << "Mbp)" << endl;
    cerr << endl;
  }
}


static void
advance_chromosome(const GenomicRegion &chrom_region, 
		   FileIterator<MappedRead> &regions) {

  while (regions.last_is_good() && 
	 (regions.get_last()->r < chrom_region)) {
    assert(regions.last_is_good());
    regions.increment_last();
  }
  while (regions.first_is_good() && 
	 (regions.get_first()->r < chrom_region))
    regions.increment_first();
}


static void
fix_chrom_names(vector<string> &chrom_names)
{
  // make sure the chrom names don't have spaces

  for (size_t i = 0; i < chrom_names.size(); ++i) {
    const size_t chr_name_end = chrom_names[i].find_first_of(" \t");
    if (chr_name_end != string::npos)
      chrom_names[i].erase(chr_name_end);
  }
}


static void
scan_chroms(const bool VERBOSE,
			const string &outfile, const vector<string> &chrom_files, 
			FileIterator<MappedRead> &regions,
            const size_t max_length) {

  // Get the output stream
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  for (size_t i = 0; i < chrom_files.size(); ++i) {
    const string fn(strip_path_and_suffix(chrom_files[i]));
    if (VERBOSE)
      cerr << "[LOADING CHROM FILE=" << fn << "]";
    vector<string> chrom_names, chroms;
    read_fasta_file(chrom_files[i].c_str(), chrom_names, chroms);
    fix_chrom_names(chrom_names);
    for (size_t j = 0; j < chroms.size(); ++j) {
      if (VERBOSE) cerr << "[SCANNING=" << chrom_names[j] << "]";
      //TODO: WHAT HAPPENS IF A CHROM IS MISSING??
      const GenomicRegion chrom_region(chrom_names[j], 0, 0);
      advance_chromosome(chrom_region, regions);
     
      scan_chromosome(chroms[j], chrom_region,
			regions, out, max_length);
    }
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
    opt_parse.add_opt("chrom", 'c', "file or dir of chroms (FASTA format; .fa suffix)",
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
    vector<string> chrom_files;
    identify_chromosomes(VERBOSE, chrom_file, fasta_suffix, chrom_files);
    sort(chrom_files.begin(), chrom_files.end());
    
    FileIterator<MappedRead> regions(mapped_reads_file, BUFFER_SIZE);
    scan_chroms(VERBOSE, outfile, chrom_files, regions, max_length);
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
