Variant calling algoritm using binomial distribution done as a research project for course Computational Genomics.

Text of the project is given in "Projektni zadatak.pdf"
Presentation about this project is "Variant calling using binomial distribution.pptx"
Youtube video for this presenation: https://www.youtube.com/watch?v=J-Ce-A3Jj40
Pseudo code of the solution is in pseudo_code.ipynb

Notebook variant_call_binom.ipynb contains python scripts that parse pileup file and generate .vcf file using binomial distribution for variant calling
Notebook metrics.ipynb contains python scripts for evaluation of generated .vcf files by comparing them to the reference .vcf file generated using bcftools-call tool.
Notebook metric_graphs.ipynb has graphic representation of evaluation metrics

As pileup file and .vcf files are too large to store on github, they are being pulled from CGC servers.
For demo purposes, in folder short_pileup_file_and_vcf, can be found fileup file of ~10000 positions and vcf with the same amount of positions.
Variant call results with short pileup file can be found in folder short_variant_call_results

Metric .csv sheets for both regular and short variant call evaluation can be found in folder "metric sheets"

Resourses used:

Pileup format:
https://en.wikipedia.org/wiki/Pileup_format
http://samtools.sourceforge.net/pileup.shtml
http://www.htslib.org/doc/samtools-mpileup.html

Variant calling:
https://www.melbournebioinformatics.org.au/tutorials/tutorials/var_detect_advanced/var_detect_advanced_background/

VCF:
https://samtools.github.io/hts-specs/VCFv4.2.pdf
https://en.wikipedia.org/wiki/Variant_Call_Format
https://www.youtube.com/watch?v=xHPm0DKAS7c
https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/

and many more...




