This is the README file for the first release of amrfinder version 1.
amrfinder is a program for identifying allele-specific DNA methylation
from bisulfite short-read sequencing technology
(such as Solexa/Illumina).


CONTACT INFORMATION:
========================================================================
Fang Fang: ffang@usc.edu
Andrew D Smith: andrewds@usc.edu

SYSTEM REQUIREMENTS:
========================================================================
The amrfinder software will only run on UNIX-type system with GNU Scientific
Library (GSL, http://www/gmi/prg/software/gsl) and GNU Compiler Collection 
(GCC, http://gcc.gnu.org/). Also, amrfinder will only run on 64-bit machines.


INSTALLATION:
========================================================================
Unpack the archive and change into the archive
directory. Then type 'make install'. A 'bin' directory will be created
in the current directory, and it will contain the program
binaries. These can be moved around, and also do not depend on any
dynamic libraries, so they should simply work when executed.


USAGE EXAMPLES:
========================================================================
Each program included in this software package will print a list of
options if executed without any command line arguments. Many of the
programs use similar options (for example, output files are specified
with '-o'). For the most basic usage of scanning allelically methylted 
regions (AMRs), use the command:

     amrfinder -o output.bed -c chroms_dir input_reads.mr

The output will appear in output.bed, and the output is in BED format
(see the UCSC Genome Browser Help documentation for details of this
format). Each line of the file indicates the location of an AMR in the 
CpG numbering format. The input read file in a format called 'MappedRead',
which consists of 8 columns with the first 6 the same as in BED format and
the last two columns being the sequence and quality score of the read.

LICENSE
========================================================================
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
