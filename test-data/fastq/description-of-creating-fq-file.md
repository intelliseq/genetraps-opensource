# Files description

capn3.1.fq and capn3.2.fq files were created by extracting reads from chr15:42377802-42397802 from sample 267 aligned to GRCh38.
The sequence contains rs80338800; rs762960425 CAPN3 on chr15:42387802

# Files preparation

The code below was used to create fastq files
``` sh
samtools view -b 267.markdup.realigned.recalibrated.bam chr15:42377802-42397802 > test.bam
samtools sort -n -o test-sorted.bam test.bam
bedtools bamtofastq -i test-sorted.bam -fq capn3.1.fq -fq2 capn3.2.fq
```
