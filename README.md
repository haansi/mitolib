# greenvariantcaller - Convert your BAM file to Variants

<a href="https://www.java.com/"><img src="http://img.shields.io/badge/language-java-brightgreen.svg" alt="Language" data-canonical-src="http://img.shields.io/badge/language-java-brightgreen.svg" style="max-width:100%;"></a></p>

This is a simple Java implementation to read BAM files based on the HTSJDK API (https://github.com/samtools/htsjdk) from small genomes like the mitochondrial genome into pileup files and performce a naive variant calling. Variant files of human mitochondrial genomes can be used for **mitochondrial haplogroup classification** with [HaploGrep2](http://haplogrep.uibk.ac.at).


## Usage Examples

### bam2var
This command writes the variants and the raw pileup file in an output folder 

```bash
java -jar greenVC-0.1.jar bam2var --in data/HG01500.IBS.exome.MT.bam --out resultfolder  
                                  --ref data/rcrs.fasta  --VAF 0.2 --QUAL 20

```

### haplocheck
This command writes from the variants generated with the naive **bam2var** to a haplogrep input file, by splitting it in major/ minor allele profiles in order to check for sample contamination 

```bash
java -jar greenVC-0.1.jar haplocheck --in variants.txt --out haplogrepinput.hsd   --VAF 0.05 

```

### haplocheck-mtDNA-Server
This command writes from the heteroplasmies.txt file generated with the **mtDNA-Server** (https://mtDNA-Server.uibk.ac.at) to a haplogrep input file, by splitting it in major/ minor allele profiles in order to check for sample contamination. Example call: 

```bash
java -jar greenVC-0.1.jar haplocheck-mtDNA-Server --in heteroplasmies.txt --out haplogrepinput.hsd  --VAF 0.05 

```

### lofreq
This command reads the variants generated with **LoFreq** (http://csb5.github.io/lofreq/) and generates to a haplogrep input file, by splitting the heteroplasmic variants in major/ minor allele profiles in order to check for sample contamination 

```bash
java -jar greenVC-0.1.jar lofreq --in inputfile.vcf --out haplogrepinput.hsd 
```

### GUI
Not providing any parameter opens a simple GUI, which performs the ** bam2var ** naive variant calling.
