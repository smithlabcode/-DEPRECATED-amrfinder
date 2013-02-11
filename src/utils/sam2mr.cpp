/* 
 *    Copyright (C) 2009 University of Southern California and
 *                       Andrew D. Smith
 *                       Song Qiang
 *                       Fang Fang
 *
 *    Authors: Song Qiang, Fang Fang
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

// #define NDEBUG

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <cmath>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"


using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream_iterator;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::make_pair;

inline bool static 
is_pairend(const size_t flag) {return flag & 0x1;}

inline bool static 
is_singlend(const size_t flag) {return !is_pairend(flag);}

inline bool static
is_mapping_paired(const size_t flag) {return flag & 0x2;}

inline bool static 
is_unmapped(const size_t flag) {return flag & 0x8;}

inline bool static 
is_mapped(const size_t flag) {return !(is_unmapped(flag));}

inline bool static
is_revcomp(const size_t flag) {return flag & 0x10;}

inline bool static
is_mate1(const size_t flag) {return flag & 0x40;}

void static
apply_CIGAR(const string &seq, const string &qual,
            const string &CIGAR, string &new_seq, string &new_qual)
{
    assert(seq.size() == qual.size());
    assert(new_seq.size() == 0 && new_qual.size() == 0);
    size_t n;
    char op;
    size_t i = 0;

    std::istringstream iss(CIGAR);
    while (iss >> n >> op)
    {
        switch (op)
        {
        case 'M':
            new_seq += seq.substr(i, n);
            new_qual += qual.substr(i, n);
            i += n;
            break;
        case 'I':
            i += n;
            break;
        case 'D':
            new_seq += string(n, 'N');
            new_qual += string(n, 'B');
            break;
        case 'S':
            i += n;
            break;
        case 'H':
            ;
            break;
        case 'P':
            ;
            break;
        case '=':
            ;
            break;
        case 'X':
            ;
            break;
        }
    }
    
    assert(i == seq.length());
    assert(new_seq.size() == new_qual.size());
}

inline static void
get_mismatch(const string &align_score, size_t &mismatch)
{
    mismatch = atoi(align_score.substr(5).c_str());
}

inline static void
get_strand(const size_t &flag, string &strand)
{
	if(is_revcomp(flag))
		strand="-";
	else
		strand="+";
}


int 
main(int argc, const char **argv) 
{
    try 
    {
        string infile("/dev/stdin");
        string outfile;

        bool VERBOSE;
        
        /****************** COMMAND LINE OPTIONS ********************/
        OptionParser opt_parse(strip_path(argv[0]),
                               "convert SAM output of Bowtie to MR format of RMAP"
                               "");
        opt_parse.add_opt("outfile", 'o', "reads output file", 
                          OptionParser::REQUIRED, outfile);
        opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

        vector<string> leftover_args;
        opt_parse.parse(argc, argv, leftover_args);
        if (argc == 1 || opt_parse.help_requested()) 
        {
            cerr << opt_parse.help_message() << endl;
            return EXIT_SUCCESS;
        }
        if (opt_parse.about_requested()) 
        {
            cerr << opt_parse.about_message() << endl;
            return EXIT_SUCCESS;
        }
        if (opt_parse.option_missing()) 
        {
            cerr << opt_parse.option_missing_message() << endl;
            return EXIT_SUCCESS;
        }
        if (leftover_args.size() > 0) infile = leftover_args.front();
        /****************** END COMMAND LINE OPTIONS *****************/

        std::ifstream in(infile.c_str());
        std::ofstream out(outfile.c_str());

        
        string line;
        size_t nline = 1;
        while (std::getline(in, line))
        {
            string name, chrom, CIGAR, mate_name, seq,
                qual, MF, Aq, align_score;

            size_t flag, start, mapq_score, mate_start;
            int seg_len;
            std::istringstream iss(line);
            iss >> name;
            if(name.substr(0,1)!="@"){
            	if (iss >> flag >> chrom >> start >> mapq_score >> CIGAR
            			>> mate_name >> mate_start >> seg_len >> seq >> qual 
            			>> MF >> Aq >> align_score)
            	{
					if (is_pairend(flag))
					{
						assert(out.good());
						if(!(flag & 0x4)){ // only working with mapped reads
							--start; // SAM are 1-based
		
							size_t mismatch;
							get_mismatch(align_score, mismatch);
						
							string strand;
							get_strand(flag, strand);
	
							string new_seq, new_qual;
							apply_CIGAR(seq, qual, CIGAR, new_seq, new_qual);
		
							if (is_revcomp(flag))
							{
								revcomp_inplace(new_seq);
								std::reverse(new_qual.begin(), new_qual.end());
							}
							
							if (is_mate1(flag)){
								out << chrom << "\t" << start << "\t"
										<< start + new_seq.length() << "\t"
										<< (name + "/1") << "\t" << mismatch << "\t"
										<< strand << "\t" << new_seq << "\t"
										<< new_qual << "\t" << endl;
							}
							else{
								out << chrom << "\t" << start << "\t"
										<< start + new_seq.length() << "\t"
										<< (name + "/2") << "\t" << mismatch << "\t"
										<< strand << "\t" << new_seq << "\t"
										<< new_qual << "\t" << endl;
							}
								
						}
		
					}
					else 
					{
						if(!(flag & 0x4)){// only working with mapped reads
							--start; // SAM are 1-based
		
							size_t mismatch;
							get_mismatch(align_score, mismatch);
								
							string strand;
							get_strand(flag, strand);
						
							if (is_revcomp(flag))
							{
								revcomp_inplace(seq);
								std::reverse(qual.begin(), qual.end());
							}
						 
							string new_seq, new_qual;
							apply_CIGAR(seq, qual, CIGAR, new_seq, new_qual);
						
							out << chrom << "\t" << start << "\t"
									<< start + new_seq.length() << "\t"
									<< name << "\t" << mismatch << "\t"
									<< strand << "\t" << new_seq << "\t"
									<< new_qual << "\t" << endl;
						}
					}
				}
				else
				{
					cerr << "Line " << nline << ": " << line << endl;
				}
            }
            nline++;
        }

        in.close();
        out.close();
    }
    catch (const SMITHLABException &e) 
    {
        cerr << "ERROR:\t" << e.what() << endl;
        return EXIT_FAILURE;
    }
    catch (std::bad_alloc &ba) 
    {
        cerr << "ERROR: could not allocate memory" << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
