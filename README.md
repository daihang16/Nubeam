# Nucleotide be a matrix (Nubeam)
Nubeam is a reference-free approach to analyze short sequencing reads. 
It represents nucleotides by matrices, 
transforms a read into a product of matrices, 
and based on which assigns numbers to reads. 
A sequencing sample, which is a collection of reads, becomes a collection of numbers that form an empirical distribution. 
Then the genetic difference between samples is quantified by the distance between empirical distributions.
## Compiling:
### Dependency
`zlib` is required to compile. To install `zlib`, run the following commands:

`wget https://www.zlib.net/zlib1211.zip`

`unzip zlib1211.zip`

`cd zlib-1.2.11/`

`./configure`

`make`

`sudo make install`
### Compile nubeam-dedup
Run the following commands:

`wget --no-check-certificate --content-disposition https://github.com/daihang16/nubeam/archive/master.zip`

`unzip nubeamdedup-master.zip`

`cd nubeamdedup-master/`

`make`
## Usage:
```
./nubeam [qtf, rgc_beta, rgc_res, cad, cad2]

./nubeam qtf [-iodwSnfh]
compute quadriples for reads in fastq format.
produces prefix.quad.gz (gc content is within) and prefix.quad.log.
-i : input filename
-o : output prefix
-d : length of the reads (default d=75).
-w : sliding window size (default w=d).
-S : sliding window step (default S=w).
-n : number of missing nucleotide allowed.
-f : value, plus 33 is the PHRED quality value of fastq reads.
-h : print this help

./nubeam rgc_beta [-ioh]
perform regression on gc contents from read quantification and output regression coefficients.
produces prefix.beta.log.
-i : input file name.
-o : output prefix.
-h : print this help

./nubeam rgc_res [-ioh beta]
regress out gc contents from read quantification and output residuals.
produces prefix.nogc.gz and prefix.nogc.log.
-i : input file name.
-beta : beta file name.
-o : output prefix.
-h : print this help

./nubeam cad [-iombh bf]
compute pariwise distances of a set; the inputs are nubeam qtf outputs.
produces prefix.cad.log.
-i : specifies input file which is output of nubeam.
-o : output prefix (prefix.log contains pairwise distance matrix)
-m : choice of methods: h2 (Hellinger distance), cos (Cosine dissimilarity).
-b : designating the number of bins per column of scores
-bf : the file describing how to partition the bins
-h : print this help

./nubeam cad2 [-ijombh bf]
compute pariwise cross distances between two sets; the inputs are nubeam qtf outputs.
produces prefix.cad2.log.
-i : specifies input file of first set, which is output of nubeam.
-j : specifies input file of second set, which is output of nubeam.
-o : output prefix (prefix.log contains pairwise distance matrix)
-m : choice of methods: h2 (Hellinger distance), cos (Cosine dissimilarity).
-b : designating the number of bins per column of scores
-bf : the file describing how to partition the bins
-h : print this help
```
## Examples:
- Quantify reads

  `./nubeam qtf -i S1.fq -o S1.fq -d 75 -a 0 -n 0 -f 0`
  
  Quantify the reads in input file `S1.fq`, with the read length of 75, adaptor size of 0, `N` not allowed in read, the output file name will be `S1.fq.quad.gz`. The output file has six columns of numbers: first four columns are Nubeam quadruplets for reads, the last two columns are GC counts for reads.

- Regress out GC content
  - Obtain regression coeffients
  
    First combined all the output files produced by `qtf` together:
  
    `cat S1.fq.quad.gz S2.fq.quad.gz S3.fq.quad.gz > all.quad.gz`
  
    Then calculate the regression coefficients for GC count:

    `./nubeam rgc_beta -i all.quad.gz -o all.quad`
  
    The regression coeffients are in `all.quad.beta.log`.
  - Obtain residuals
  
    For each original output files produced by `qtf`:
    
    `./nubeam rgc_res -i S1.fq.quad.gz -beta all.quad.beta.log -o S1.fq.quad`
    
    The residuals will be written to `S1.fq.quad.nogc.gz`.
  
- Quantify pair-wise distance
  - Calculate within-group distances
  
    `./nubeam cad -o output -m h2 -b 10 -bf bin.txt -i S1.fq.quad.nogc.gz -i S2.fq.quad.nogc.gz -i S3.fq.quad.nogc.gz`
  
    For `n` samples, the command calculate `n(n-1)/2` Hellinger distances. The number of bins partitioned for R4 space is `10^4`, if the `bin.txt` exists, it will be used for partitioning; if not, the partitioning will be calculated and written to `bin.txt`. The distance matrix is at the end of `output.cad.log`.
  
  - Calculate between-group distances
  
    `./nubeam cad2 -o output -m h2 -b 10 -bf bin.txt -i S1.fq.quad.nogc.gz -i S2.fq.quad.nogc.gz -i S3.fq.quad.nogc.gz -j S4.fq.quad.gz -j S5.fq.quad.gz`
  
    For a group of `n` samples and a group of `m` samples, the command calculate `nm` Hellinger distances. The number of bins partitioned for R4 space is `10^4`, if the `bin.txt` exists, it will be used for partitioning; if not, the partitioning will be calculated and written to `bin.txt`. The distance matrix is at the end of `output.cad2.log`.
    
