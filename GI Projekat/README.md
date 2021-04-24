
# Variant calling using binomial distribution

Variant calling algoritm using binomial distribution done as a research project for course Computational Genomics.

------------
Text of the project is given in "Projektni zadatak.pdf"
Presentation about this project is "Variant calling using binomial distribution.pptx"
Youtube video for this presenation: [Variant calling using binomial distribution](https://www.youtube.com/watch?v=J-Ce-A3Jj40)
Pseudo code of the solution is in pseudo_code.ipynb

- Notebook variant_call_binom.ipynb contains python scripts that parse pileup file and generate .vcf file using binomial distribution for variant calling
- Notebook metrics.ipynb contains python scripts for evaluation of generated .vcf files by comparing them to the reference .vcf file generated using bcftools-call tool.
- Notebook metric_graphs.ipynb has graphic representation of evaluation metrics

As pileup file and .vcf files are too large to store on github, they are being pulled from CGC servers.
For demo purposes, in folder short_pileup_file_and_vcf, can be found fileup file of ~10000 positions and vcf with the same amount of positions.
Variant call results with short pileup file can be found in folder short_variant_call_results

Metric .csv sheets for both regular and short variant call evaluation can be found in folder "metric sheets"

------------

#### Resourses used:

Pileup format:

[Wikipedia - Pileup format](https://en.wikipedia.org/wiki/Pileup_format)
[Sourceforge - Samtools - Pileup format](http://samtools.sourceforge.net/pileup.shtml)
[Samtools-mpileup Manual page](http://www.htslib.org/doc/samtools-mpileup.html)

Variant calling:

[Melbourne Bioinformatics - Tutorials - Introduction to Variant detection](https://www.melbournebioinformatics.org.au/tutorials/tutorials/var_detect_advanced/var_detect_advanced_background/)
[Variant identification and analysis](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/)

VCF:
[VCF Version 4.2 Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
[Wikipedia - Variant Call Format](https://en.wikipedia.org/wiki/Variant_Call_Format)
[Youtube: Understanding VCF file | Variant Call Format Part 1/3](https://www.youtube.com/watch?v=xHPm0DKAS7c)
[Youtube: Understanding VCF file | Variant Call Format Part 2/3](https://www.youtube.com/watch?v=WV3Pls1_z_4)
[Youtube: Understanding VCF file | Variant Call Format Part 3/3](https://www.youtube.com/watch?v=2P4EItXCtFI)
[Understanding VCF format](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/)
