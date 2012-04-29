/*    cpg2base: A utility program for converting genomic intervals
 *    specifying a range in terms of CpGs into one in terms of
 *    bases. IMPORTANT: the CpG interval is closed, since a half-open
 *    CpG interval would present semantic issues in converting back to
 *    bases.
 *
 *    Copyright (C) 2012 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
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

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}


static void
collect_cpgs(const string &s, unordered_map<size_t, size_t> &cpgs) {
  const size_t lim = s.length() - 1;
  size_t cpg_count = 0;
  for (size_t i = 0; i < lim; ++i)
    if (is_cpg(s, i)) {
      cpgs[cpg_count] = i;
      ++cpg_count;
    }
}


static void
identify_chromosomes(const string chrom_file, const string fasta_suffix, 
		     unordered_map<string, string> &chrom_files) {
  vector<string> the_files;
  if (isdir(chrom_file.c_str())) {
    read_dir(chrom_file, fasta_suffix, the_files);
    for (size_t i = 0; i < the_files.size(); ++i)
      chrom_files[strip_path_and_suffix(the_files[i])] = the_files[i];
  }
  else chrom_files[strip_path_and_suffix(chrom_file)] = chrom_file;
}



static void
convert_coordinates(const unordered_map<size_t, size_t> &cpgs, 
		    GenomicRegion &region)  {
  const unordered_map<size_t, size_t>::const_iterator 
    start_itr(cpgs.find(region.get_start()));
  const unordered_map<size_t, size_t>::const_iterator
    end_itr(cpgs.find(region.get_end()));
  if (start_itr == cpgs.end() || end_itr == cpgs.end())
    throw SMITHLABException("could not convert:\n" + region.tostring());
  region.set_start(start_itr->second);
  region.set_end(end_itr->second);
}


int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    
    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "A utility program for converting "
			   "genomic intervals in CpGs to genomic intervals in bases",
			   "<cpg-intervals>");
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("chrom", 'c', "file or dir of chroms (.fa extn)", true , chrom_file);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string interval_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    unordered_map<string, string> chrom_files;
    identify_chromosomes(chrom_file, fasta_suffix, chrom_files);
    if (VERBOSE)
      cerr << "CHROMS:\t" << chrom_files.size() << endl;
    
    std::ifstream in(interval_file.c_str());
    if (!in) 
      throw SMITHLABException("cannot open input file " + interval_file);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
    unordered_map<size_t, size_t> cpgs;
    vector<string> chrom_names, chroms;
    GenomicRegion region;
    GenomicRegion chrom_region("chr0", 0, 0);
    while (!in.eof() && in >> region) {
      // get the correct chrom if it has changed
      if (!region.same_chrom(chrom_region)) {
	const unordered_map<string, string>::const_iterator 
	  fn(chrom_files.find(region.get_chrom()));
	if (fn == chrom_files.end())
	  throw SMITHLABException("could not find chrom: " + region.get_chrom());
	chrom_names.clear();
	chroms.clear();
	read_fasta_file(fn->second.c_str(), chrom_names, chroms);
	if (chrom_names.size() > 1)
	  throw SMITHLABException("multiple chroms/file: " + fn->second);
	if (VERBOSE)
	  cerr << "PROCESSING: " << chrom_names.front() << endl;
	collect_cpgs(chroms.front(), cpgs);
	chrom_region.set_chrom(chrom_names.front());
      }
      convert_coordinates(cpgs, region);
      out << region << endl;
    }
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
