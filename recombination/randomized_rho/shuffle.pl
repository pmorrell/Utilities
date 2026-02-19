#!/usr/bin/perl -w
# written by Peter Morrell, April 2011, St. Paul, MN
use List::Util qw(shuffle);

my ($fileName) = @ARGV;
open IN, "< $fileName" || die "Can't open $fileName";

while(<IN>) {
    chomp;
    push(@sample_name, $_);
	}

@random_name = shuffle(@sample_name);

 foreach (@random_name) {
 	print $_ . "\n";
 } 
