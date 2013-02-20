/*
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith and Fang Fang
 *
 *    Authors: Andrew D. Smith and Fang Fang
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

#include "BAMIO.hpp"

#include <smithlab_utils.hpp>
#include <smithlab_os.hpp>
#include <MappedRead.hpp>
#include <iostream>
#include <fstream>
#include <numeric>

using std::string;
using std::vector;
using std::ifstream;
using std::tr1::unordered_map;
using std::streampos;

#include "api/BamReader.h"
#include "api/BamAlignment.h"

using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefVector;
using BamTools::BamReader;
using BamTools::RefData;

static MappedRead
BamToMappedRead(const unordered_map<size_t, string> &chrom_lookup,
		   	   const BamAlignment &ba) {
  const unordered_map<size_t, string>::const_iterator
    the_chrom(chrom_lookup.find(ba.RefID));
  if (the_chrom == chrom_lookup.end())
    throw SMITHLABException("no chrom with id: " + toa(ba.RefID));

  const string chrom = the_chrom->second;
  const size_t start = ba.Position;
  const size_t end = start + ba.Length;
  const string name = ba.Name;
  const float score = ba.MapQuality;
  const char str='+';
  const string seq = ba.QueryBases;
  const string qual = ba.Qualities;

  MappedRead mr;
  mr.r = GenomicRegion(chrom,start,end,name,score,str);
  mr.seq = seq;
  mr.scr = qual;
  return mr;
}

size_t
BAMIO::max_fragment_size = 1000ul;

static void
identify_chromosomes(const string chrom_file, const string fasta_suffix, 
		     unordered_map<string, string> &chrom_files) {
  vector<string> the_files;
  if (isdir(chrom_file.c_str())) {
    read_dir(chrom_file, fasta_suffix, the_files);
    for (size_t i = 0; i < the_files.size(); ++i)
      chrom_files[strip_path_and_suffix(the_files[i])] = the_files[i];
  }
  else the_files.push_back(chrom_file);
}


BAMIO::BAMIO(const string &rfn, const bool v, 
		     const string &chrom_dir) :
  reads_file_name(rfn), position(0), VERBOSE(v){
  static const string fasta_suffix = "fa";
  identify_chromosomes(chrom_dir, fasta_suffix, chrom_seq_files);
}

static void
backup_to_start_of_current_record(std::ifstream &in) {
  if (in.tellg() > 0)
    do {
      in.seekg(-1, std::ios_base::cur);
    }
    while (in.tellg() > 0 && in.peek() != '\n' && in.peek() != '\r');
}

streampos
BAMIO::find_last_position_before(const string &query_chrom, 
				     const size_t query_pos,
				     const size_t range_low, 
				     const size_t range_high,
				     std::ifstream &in) {
  size_t high_pos = range_high;
  in.seekg(high_pos);
  assert(in.tellg() >= 0);
  backup_to_start_of_current_record(in);
  
  size_t low_pos = range_low;
  in.seekg(low_pos);
  assert(in.tellg() >= 0);
  backup_to_start_of_current_record(in);
  
  string chrom, seq;
  size_t start = 0ul;
  
  while ((high_pos>low_pos)&&(high_pos - low_pos > 1)) {
    const size_t mid_pos = (low_pos + high_pos)/2;
    in.seekg(mid_pos);
    assert(in.tellg() >= 0);
    backup_to_start_of_current_record(in);
    if (!(in >> chrom >> start >> seq))
      throw SMITHLABException("problem loading reads");
  
    if (chrom < query_chrom || (chrom == query_chrom && start + 
				seq.length() <= query_pos))
      low_pos = mid_pos;
    else
      high_pos = mid_pos;
  }
  
  while ((in >> chrom >> start >> seq) && 
	 (chrom < query_chrom || (chrom == query_chrom && start + 
				  seq.length() <= query_pos)))
    low_pos = in.tellg();
  
  return low_pos;
}


streampos
BAMIO::find_first_position_after(const string &query_chrom, 
				     const size_t query_pos,
				     const size_t range_low, 
				     const size_t range_high,
				     std::ifstream &in) {
  size_t high_pos = range_high;
  in.seekg(high_pos);
  assert(in.tellg() >= 0);
  backup_to_start_of_current_record(in);
  
  size_t low_pos = range_low;
  in.seekg(low_pos);
  assert(in.tellg() >= 0);
  backup_to_start_of_current_record(in);
  
  string chrom, seq;
  size_t start = 0ul;

  while (high_pos - low_pos > 1) {
    const size_t mid_pos = (low_pos + high_pos)/2;
    in.seekg(mid_pos);
    assert(in.tellg() >= 0);
    backup_to_start_of_current_record(in);

    if (!(in >> chrom >> start >> seq))
      throw SMITHLABException("problem loading reads");
    
    if (chrom < query_chrom || (chrom == query_chrom && start <= query_pos))
      low_pos = mid_pos;
    else
      high_pos = mid_pos;
  }
  
  return high_pos;
}


streampos
BAMIO::find_last_position_before_mr(const string &query_chrom, 
					const size_t query_pos,
					const size_t range_low, 
					const size_t range_high,
					std::ifstream &in) {
  
  assert(in.tellg() >= 0);
  size_t high_pos = range_high;
  in.seekg(high_pos);
  backup_to_start_of_current_record(in);
  
  size_t low_pos = range_low;
  in.seekg(low_pos);
  backup_to_start_of_current_record(in);
  
  MappedRead mr;
  in.seekg(low_pos);
  
  while (high_pos - low_pos > 1) {
    
    const size_t mid_pos = (low_pos + high_pos)/2;
    in.seekg(mid_pos);
    backup_to_start_of_current_record(in);
    
    if (!(in >> mr))
      throw SMITHLABException("problem loading epireads");
    
    if (mr.r.get_chrom() < query_chrom || (mr.r.get_chrom() == query_chrom && 
					   mr.r.get_start() + max_fragment_size <= 
					   query_pos))
      low_pos = mid_pos;
    else high_pos = mid_pos;
  
  }

  while ((in >> mr) && mr.r.get_chrom() == query_chrom && mr.r.get_end() < query_pos)
    low_pos = in.tellg();
  
  return low_pos;
}


streampos
BAMIO::find_first_position_after_mr(const string &query_chrom, 
					const size_t query_pos,
					const size_t range_low, 
					const size_t range_high,
					std::ifstream &in) {
  assert(in.tellg() >= 0);
  
  size_t high_pos = range_high;
  in.seekg(high_pos);
  backup_to_start_of_current_record(in);
  
  size_t low_pos = range_low;
  in.seekg(low_pos);
  backup_to_start_of_current_record(in);
  
  MappedRead mr;
  // string chrom, seq;
  // size_t start = 0ul;
  in.seekg(low_pos);

  while (high_pos - low_pos > 1) {
    
    const size_t mid_pos = (low_pos + high_pos)/2;
    in.seekg(mid_pos);
    backup_to_start_of_current_record(in);
    
    if (!(in >> mr))
      throw SMITHLABException("problem loading reads");
    
    if (mr.r.get_chrom() < query_chrom || (mr.r.get_chrom() == query_chrom && 
					   mr.r.get_start() <= query_pos))
      low_pos = mid_pos;
    else high_pos = mid_pos;
    
  }
  
  return high_pos;
}


inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}


static void
collect_cpgs(const string &s, unordered_map<size_t, size_t> &cpgs) {
  const size_t lim = s.length() - 1;
  size_t cpg_count = 0;
  for (size_t i = 0; i < lim; ++i) {
    if (is_cpg(s, i)) {
      cpgs[i] = cpg_count;
      ++cpg_count;
    }
  }
}


static bool
convert_meth_states_pos(const string &chrom, 
			const unordered_map<size_t, size_t> &cpgs, 
			const MappedRead &mr, epiread &r) {
  const size_t width = mr.r.get_width();
  const size_t offset = mr.r.get_start();

  size_t cpg_count = 0;
  string states;
  size_t first_cpg = std::numeric_limits<size_t>::max();
  size_t last_cpg = first_cpg;
  for (size_t i = 0; i < width; ++i) {
    if (is_cpg(chrom, offset + i)) {
      if (mr.seq[i] == 'C') {
	states += 'C';
	++cpg_count;
      }
      else if (mr.seq[i] == 'T') {
	states += 'T';
	++cpg_count;
      }
      else states += 'N';
      if (first_cpg == std::numeric_limits<size_t>::max()) {
	first_cpg = i;
      }
      last_cpg = i;
    }
  }
  if (first_cpg != std::numeric_limits<size_t>::max()) {
    const size_t the_cpg = cpgs.find(offset + first_cpg)->second;
    r.pos = the_cpg;
    r.seq = states;
  }
  return cpg_count > 0;
}


static bool
convert_meth_states_neg(const string &chrom, 
			const unordered_map<size_t, size_t> &cpgs, 
			MappedRead mr, epiread &r) {
  const size_t width = mr.r.get_width();
  const size_t offset = mr.r.get_start();
  
  // NEED TO TAKE REVERSE COMPLEMENT FOR THE NEGATIVE STRAND ONES!
  revcomp_inplace(mr.seq);
  
  size_t cpg_count = 0;
  string states;
  size_t first_cpg = std::numeric_limits<size_t>::max();
  size_t last_cpg = first_cpg;
  for (size_t i = 0; i < width; ++i) {
    if (is_cpg(chrom, offset + i - 1)) {
      if (mr.seq[i] == 'G') {
	states += 'C';
	++cpg_count;
      }
      else if (mr.seq[i] == 'A') {
	states += 'T';
	++cpg_count;
      }
      else states += 'N';
      if (first_cpg == std::numeric_limits<size_t>::max()) {
	first_cpg = i;
      }
      last_cpg = i;
    }
  }
  if (first_cpg != std::numeric_limits<size_t>::max()) {
    assert(cpgs.find(offset + first_cpg - 1) != cpgs.end());
    const size_t the_cpg = cpgs.find(offset + first_cpg - 1)->second;
    r.pos = the_cpg;
    r.seq = states;
  }
  return cpg_count > 0;
}


bool
BAMIO::load_reads_next_chrom(size_t &refID, string &chrom, vector<epiread> &the_reads) {
  if (position == std::numeric_limits<size_t>::max()) return false;
  
  string prev_chrom;
  the_reads.clear();
  
  // open and check the file  
  BamReader reader;
  reader.Open(reads_file_name);

  // Get header and reference
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  unordered_map<size_t, string> chrom_lookup;
  for (size_t i = 0; i < refs.size(); ++i)
    chrom_lookup[i] = refs[i].RefName;
  reader.Jump(refID,position);

  unordered_map<size_t, size_t> cpgs;
  vector<string> chrom_names, chrom_seqs;
  MappedRead mr;
  BamAlignment bam;
  while ((reader.GetNextAlignment(bam)) && 
		  (chrom_lookup[bam.RefID] == prev_chrom || prev_chrom.empty())) {
	  mr = BamToMappedRead(chrom_lookup, bam);
	if (prev_chrom.empty()) {
		if (VERBOSE)
	  	  std::cerr << "[STATUS] loading reads " 
	  		    << mr.r.get_chrom() << std::endl;
	  	const unordered_map<string, string>::const_iterator 
	  	  fn(chrom_seq_files.find(mr.r.get_chrom()));
	  	if (fn == chrom_seq_files.end())
	  	  throw SMITHLABException("could not find chrom: " + mr.r.get_chrom());
	  	chrom_names.clear();
	  	chrom_seqs.clear();
	  	read_fasta_file(fn->second.c_str(), chrom_names, chrom_seqs);
	  	assert(chrom_names.size() == 1);
	  	collect_cpgs(chrom_seqs.front(), cpgs);
	}
	epiread r;
	const bool has_cpgs = mr.r.pos_strand() ? 
	convert_meth_states_pos(chrom_seqs.front(), cpgs, mr, r) :
	convert_meth_states_neg(chrom_seqs.front(), cpgs, mr, r);
	if (has_cpgs)
		the_reads.push_back(r);
	prev_chrom = mr.r.get_chrom();
  }
  reader.Close(); 
  chrom = prev_chrom;
  if (prev_chrom.empty()) return false;
	      
  vector<size_t> tmp(cpgs.size(), 0ul);
  for (unordered_map<size_t, size_t>::iterator i(cpgs.begin()); 
		  i != cpgs.end(); ++i)
	tmp[i->second] = i->first;
	cpg_rev_lookup[chrom].swap(tmp);
	    
	if (the_reads.empty()) return false;
}


void
BAMIO::convert_coordinates(vector<GenomicRegion> &amrs) const {
    unordered_map<string, vector<size_t> >::const_iterator current_cpgs;
    for (size_t i = 0; i < amrs.size(); ++i) {
      if (i == 0 || !amrs[i].same_chrom(amrs[i-1])) {
	current_cpgs = cpg_rev_lookup.find(amrs[i].get_chrom());
	assert(current_cpgs != cpg_rev_lookup.end());
      }
      assert(current_cpgs->second.size() > amrs[i].get_end());
      amrs[i].set_start(current_cpgs->second[amrs[i].get_start()]);
      amrs[i].set_end(current_cpgs->second[amrs[i].get_end()]);
    }
}
