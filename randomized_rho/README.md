# Code for running Hudson's two site likelihood recombination rate estimator

This shell script calls two Perl scripts to subsample a large fasta file.

Recombination rate estimation in maxhap uses lookup tables of a fixed size. This code is designed to take a random set of individuals from a list, draw a fixed number individuals, and then calculate rho some number of times. This could all be made clearer in the code.

The file here include:
Horse_list.txt    List of individual samples, this list is randomized and then a fixed number of individuals is chose.    
shuffle.pl        A Perl script that reorders sample names in a flat text file

Dependencies:
selectSeqs.pl    From Naoki Takebashi's [sequence analysis scripts] \
                 (http://raven.iab.alaska.edu/~ntakebay/teaching/programming/perl-scripts/selectSeqs.pl)
tohudson2001     From Kevin Thornton [sequtils](https://github.com/molpopgen/sequtils)


