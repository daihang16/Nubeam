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
```
## Examples:
- For single-end reads
  - Consider reads from complementary strand (default)
  
    `./nubeam-dedup -i read.fq`
    
    The command gives the following output on screen:
    
    `Output unique reads to /current/working/directory/read.uniq.fastq`
  - Do not consider reads from complementary strand
    
    `./nubeam-dedup -i read.fq -s 0`
  - Consider reads from complementary strand (default), output gzipped file with a compression level of 6, output removed duplicated reads 
  
    `./nubeam-dedup -i read.fq -z 6 -r 1`
    
    The command gives the following output on screen:
    
    `Output removed duplicated reads to /current/working/directory/read.removed.fastq.gz`
    
    `Output unique reads to /current/working/directory/read.uniq.fastq.gz`
