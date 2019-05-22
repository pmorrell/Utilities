#!/usr/bin/perl -wpl -i.bak 

# convert names in USDA manifest files to standard names that are sortable
# adding leading "0s" to name so things will sort into order

# remove a couple of columns with nothing but missing data
s/,,/,/g;

# remove a couple of markers names from the columns that had no data
s/PI290175,//g;
s/PI356742,//g;


# replaces sample prefix & tab with prefix only
# I don't remember why?
s/(CIho)\t/$1/g;
s/(PI)\t/$1/g;

# convert all CIho and PI names so they are the same length and can be sorted
s/CIho(\d\d\d\b)/CIho00$1/g;
s/CIho(\d\d\d\d\b)/CIho0$1/g;
# convert PI names with 4 digits
s/PI(\d\d\d\d\b)/PI00$1/g;
# convert PI names with 5 digits
s/PI(\d\d\d\d\d\b)/PI0$1/g;

# following lines will clean up "BOPA1_" or "BOPA2_" in SNP name
# applied to genotyping file, "TCAPbarley9K_data.csv_ABconverted.csv"
s/BOPA[1-2]_//g;
# changes dashes to underscores in SNP names; R and other tools don't like dashes in names
s/\w-\w/_/g;



# lines below are for sample list
#s/\b[0-9]{4}\b/00$&/g;
#s/\b[0-9]{5}\b/0$&/g;
#s/\b\s\b//g;
