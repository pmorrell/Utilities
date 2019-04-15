# Identifying mis-mapped raw_reads

- to identify improperly mapped reads, need to compare the original
fastq to a BAM file

- extract read names from a fastq file
- read a gzipped fastq file, strip out the read names, sort and remove '@'
`cd /home/morrellp/pmorrell/shared/Datasets/NGS/Barley_Exome/Barley_Genomic_Prediction/Rasmusson`
`zcat Rasmusson_R1.fastq.gz | awk '/^@HISEQ/ {print $1}' | sed -e s/\@//g > ~/scratch/raw_reads.txt`

- extract read names from the matching BAM file
`cd ~/shared/Projects/Selective_Sweeps/sequence_handling/Skylar/SAM_Processing/Picard`
`samtools view Rasmusson.bam | cut -f 1 > ~/scratch/mapped_reads.txt`

- my scratch is at /panfs/roc/scratch/pmorrell
- identify reads that in the mapped reads file but not from the original fastq
- the lscpu command will identify the number of available cores
`lscpu`
- lab nodes have 16 cores!
- GNU sort can use all the cores at once
`time sort -u --parallel=16 raw_reads.txt >raw_reads_sort.txt &`\
real	8m40.872s\
user	16m5.420s\
sys	0m12.326s\

`time sort -u --parallel=2 mapped_reads.txt >mapped_reads_sort.txt &`
real	4m29.568s
user	7m48.458s
sys	0m5.446s

- is sorting faster with more processors?
- no!
`time sort -u --parallel=16 mapped_reads.txt >mapped_reads_sort.txt &`
real	7m4.601s
user	13m38.326s
sys	0m12.079s

- don't bother with the parallel operation!
`time sort -u mapped_reads.txt >mapped_reads_sort.txt &`
real	4m21.660s
user	8m10.404s
sys	0m7.201s


- using comm with options below will return reads in mapped reads (in BAM) that aren't in the original fastq
`comm -23 mapped_reads.txt raw_reads.txt`

- for this comparison, there are no reads in the BAM that did not come from the appropriate fastq
