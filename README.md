# greenvariantcaller - Convert your BAM file to Variants

<a href="https://www.java.com/"><img src="http://img.shields.io/badge/language-java-brightgreen.svg" alt="Language" data-canonical-src="http://img.shields.io/badge/language-java-brightgreen.svg" style="max-width:100%;"></a></p>

This is a simple Java implementation to read BAM files based on the HGSJDK API (https://github.com/samtools/htsjdk) from small genomes like the mitochondrial genome into pileup files and performce a naive variant calling. Variant files of human mitochondrial genomes can be used for **mitochondrial haplogroup classification** with [HaploGrep2](http://haplogrep.uibk.ac.at).


## Usage Examples

### Default Command
This command writes the variants and the raw pileup file in an output folder 

```bash
java -jar greenVC-0.1.jar bam2var --in data/HG01500.IBS.exome.MT.bam --out resultfolder  --ref data/rcrs.fasta  --VAF 0.2 --QUAL 20

```

