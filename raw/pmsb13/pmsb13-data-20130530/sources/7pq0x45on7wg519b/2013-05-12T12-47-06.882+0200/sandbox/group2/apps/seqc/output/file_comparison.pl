#!/usr/bin/env perl
use strict; 
use warnings;

use Text::Diff;
use File::Compare;

my $f1 = $ARGV[0];
my $f2 = $ARGV[1];
my $outfile = $ARGV[2];
open FILE1, "$f1" or die "Could not open file: $! \n";
open FILE2, "$f2" or die "Could not open file: $! \n";

if (compare($f1,$f2) == 0) {
	print "They're equal\n";
}
else  
{
	open (OUTFILE, ">$outfile") or die "Cannot open $outfile for writing \n";
	my $diffs = diff $f1 => $f2;

	print OUTFILE $diffs;
	close OUTFILE;
}
close(FILE1);
close(FILE2);

exit 0;
