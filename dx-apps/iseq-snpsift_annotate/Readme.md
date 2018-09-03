**Input:**
* **vcf** - file to be annotated, can be compressed,
* **output_name** - output file will be named output_name.vcf,
* **database** - compressed file to annotate with (index is produced within the app).


**Available databses (originally used in this order):**

| Dataset | vcf name | runtime* | final cost* |
|----|----|----:|----:|
| clinvar | clinvar.v2.vcf.gz | 0:00:48 | $0.0039 |
| mitomap diseases | mitomap-diseases.vcf.gz | 0:00:49 | $0.0051 |
| gnomad-exomes | gnomad-exomes-all-populations.vcf.gz | 0:02:39 | $0.0127 |
| gnomad-genomes | gnomad-genomes-all-populations.vcf.gz |0:27:26 |$0.1308|
| 1000 genomes | kg.vcf.gz |  0:02:52 | $0.0179|
| mitomap | mitomap.vcf.gz |0:00:50|$0.0039|
| sift | sift.vcf.gz |0:03:14| $0.0155|
| cadd | cadd.vcf.gz | 1:11:18 | $0.4421 |
| m-cap | m-cap.vcf.gz |0:01:36|$0.0100|
| iseq | iseq-population-frequencies.vcf.gz |0:01:02| $0.0065|


*on `test_FBN1.vcf.gz`


**dx run command**

`dx run iseq_snpsift_annotate -ivcf='test_FBN1.vcf.gz' -ioutput_name='test_FBN1_annotated_clinvar' -idatabase='genetraps-resources:vcf/clinvar.v2.vcf.gz' `

