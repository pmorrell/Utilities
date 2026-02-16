#!/usr/bin/perl -w

use strict;

# getCol.pl metaData colFirst colLast file
# Extract a range of columns from a tab-delimited file.
# The program takes four arguments:
#   1. The number of columns of metaData (sample name, etc.)
#   2. The first column of data needed
#	3. The last column of data needed
#	4. The name of the file

my ($metaData, $colFirst, $colLast,  $fileName) = @ARGV;

# switch to Perl 0-based indexing
$metaData -= 1;
$colFirst -= 1;
$colLast -= 1;

open IN, "< $fileName" || die "Can't open $fileName";

while(<IN>) {
  chomp;
  my @line = split(/\t/, $_);  # split tab delimited line ($_) and make an array

my @slice = @line[0..$metaData,$colFirst..$colLast];	

@slice = join("\t", @slice);
print @slice, "\n";

}

close(IN);
exit;

