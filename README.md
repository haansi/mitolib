# mitolib - a library for mtDNA genome data processing

<a href="https://www.java.com/"><img src="http://img.shields.io/badge/language-java-brightgreen.svg" alt="Language" data-canonical-src="http://img.shields.io/badge/language-java-brightgreen.svg" style="max-width:100%;"></a></p>

⛔️ DEPRECATED!!!

- Please use <b>Mutserve</b> for variant calling https://github.com/seppinho/mutserve 
- or <b>HaploCheck</b> for contamination detection https://github.com/genepi/haplocheck

Be aware of the limitations when using the outdated Tools:

- splitter - Split mitochondrial variants and heteroplasmies from <a href="https://mtdna-server.uibk.ac.at/index.html">mtDNA-Server</a>
- contChecker - Compare mitochondrial profiles from extended report in <a href="http://haplogrep.uibk.ac.at/">HaploGrep 2</a>
- lofreq - Split mitochondrial variants according the VCF file generated with <a href="http://csb5.github.io/lofreq/">LoFreq</a> as Variant file (txt)
- bam2var - naive variant caller (input BAM file, output Variant file txt)
- haplochecker - check for contamination in mtDNA NGS data (BAM file) based on <a href="http://phylotree.org/">Phylotree 17</a> 

```bash
java -jar mitolib-0.1.2.jar haplochecker --in Bamfile.bam --out outputfolder --ref rCRS.fasta  
--QUAL 20 --MAPQ 200 --VAF 0.01
```

- haplochecker2 - check for contamination in <a href="https://mtdna-server.uibk.ac.at/index.html">mtDNA-Server</a>  raw data  (large txt.file) based on <a href="http://phylotree.org/">Phylotree 17</a>   

```bash
java -jar mitolib-0.1.2.jar haplochecker2 --in raw.txt --out outputfolder --VAF 0.01 
```
