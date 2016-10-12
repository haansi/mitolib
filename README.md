# greenvariantcaller - Convert your BAM file to Variants

<a href="https://www.java.com/"><img src="http://img.shields.io/badge/language-java-brightgreen.svg" alt="Language" data-canonical-src="http://img.shields.io/badge/language-java-brightgreen.svg" style="max-width:100%;"></a></p>

This is a simple Java implementation to read BAM files based on the HGSJDK API (https://github.com/samtools/htsjdk) from small genomes like the mitochondrial genome into pileup files and performce a naive variant calling. Variant files of human mitochondrial genomes can be used for **mitochondrial haplogroup classification** with [HaploGrep2](http://haplogrep.uibk.ac.at).


## Usage Examples

### Default Command for bam2var
This command writes the variants and the raw pileup file in an output folder 

```bash
java -jar greenVC-0.1.jar bam2var --in data/HG01500.IBS.exome.MT.bam --out resultfolder  
                                  --ref data/rcrs.fasta  --VAF 0.2 --QUAL 20

```

### Default Command for haplochecker
This command writes from the variants generated with the naive **bam2var** to a haplogrep input file, by splitting it in major/ minor allele profiles in order to check for sample contamination 

```bash
java -jar greenVC-0.1.jar haplochecker --in variants.txt --out haplogrepinput.hsd   --VAF 0.05 

```
