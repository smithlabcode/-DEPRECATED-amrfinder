#!/usr/bin/perl

=head1 NAME

mergemates.pl - a program for merging mates in paired-end sequencing data;

=head1 SYNOPSIS

perl mergemates.pl [OPTIONS]

 Options:
 -i               the input mapped read file (.mr)
 -o               the output merged file (.bed)
 -h,-help         brief help information
 -man             full documentation

=cut


package main;
use strict;
use warnings;
use IO::File;
use Getopt::Long ();    #   Resist name-space pollution!
use Pod::Usage ();      #   Ditto!
use File::Basename;

MAIN:
{
   #   Check arguments.
    my( $help, $man, $input, $output);
    $help=0;
    $man=0;

    Getopt::Long::GetOptions(
        'h|help' =>  \$help,
        'man' =>  \$man,
        'i=s' => \$input,  
        'o=s' => \$output,  
    );

    #   Check for requests for help or for man (full documentation):

    Pod::Usage::pod2usage( -exitstatus => 2 ) if ( $help );
    Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 2 ) if ( $man );

    #   Check for required variables.
        unless ( defined($input) && defined($output))
    {
        Pod::Usage::pod2usage( -exitstatus => 2 );
    }
    
    open(INPUT, "$input")||die"Can not open the input file.\n";
    open(OUTPUT,">$output")||die"Can not open the output file.\n";
    
    my $seq;
    my $pre_seq;
	my @items;
	my @pre_items;
	my @names;
	my @pre_names;
	my $diff_chrom_count=0;
	my $flag=1;
	
	my $frag=basename($output,".bed")."_fragLen.txt";
	open(FRAG,">$frag")||die"Can not open the fragment length file.\n";
	my $fragLen;
	
	$pre_seq = <INPUT>;
	chomp($pre_seq);
	@pre_items = split(/[\t ]+/, $pre_seq);
	@pre_names = split("/",$pre_items[3]);
	while($seq=<INPUT>){
		chomp($seq);
		@items = split(/[\t ]+/, $seq);
		@names = split("/",$items[3]);
		
		if($names[0] eq $pre_names[0]){
			if($items[0] eq $pre_items[0]){
				if($items[1]<$pre_items[1]){
					print OUTPUT "$items[0]\t$items[1]\t$pre_items[2]\t$names[0]\t0\t+\n";
					$fragLen = $pre_items[2]-$items[1];
					print FRAG "$fragLen\n";
				}
				else{
					print OUTPUT "$items[0]\t$pre_items[1]\t$items[2]\t$names[0]\t0\t+\n";
					$fragLen = $items[2]-$pre_items[1];
					print FRAG "$fragLen\n";
				}
			}
			else{
				$diff_chrom_count++;
			}
			$flag=2;
			if(defined($pre_seq = <INPUT>)){
				chomp($pre_seq);
				@pre_items = split(/[\t ]+/, $pre_seq);
				@pre_names = split("/",$pre_items[3]);
				$flag=1;
			}
		}
		else{
			print OUTPUT "$pre_items[0]\t$pre_items[1]\t$pre_items[2]\t$pre_items[3]\t$pre_items[4]\t$pre_items[5]\n";
			$fragLen = $pre_items[2]-$pre_items[1];
			print FRAG "$fragLen\n";
			$pre_seq=$seq;
			@pre_items = split(/[\t ]+/, $pre_seq);
			@pre_names = split("/",$pre_items[3]);
			$flag=1;
		}
	}
	if($flag==1){
		print OUTPUT "$pre_items[0]\t$pre_items[1]\t$pre_items[2]\t$pre_items[3]\t$pre_items[4]\t$pre_items[5]\n";
		$fragLen = $pre_items[2]-$pre_items[1];
		print FRAG "$fragLen\n";
	}
	print "Mismapped mates to different chromosomes: $diff_chrom_count\n";
	close(INPUT);
	close(OUTPUT); 
	close(FRAG);
}

