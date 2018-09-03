**Input:**
* **vcf** - vcf or vcf.gz file to filter (after joint_genotyping),
* **output_name** - output file will be named output_name.vcf,
* **command** - snpsift command to filter, without "", e.g. `( QUAL >= 30 )` or `(isHom( GEN[0] ) & isRef( GEN[0] ))` or `(ANN[ANY].IMPACT='HIGH' | ANN[ANY].IMPACT='MODERATE)`. For more information look here: [http://snpeff.sourceforge.net/SnpSift.html#filter]

**warning**

Using isHom() and countHom(): a no call "./." is considered as a homozygote. To avoid this try using `!(GEN[?].GT='./.')`.

**dx run command**

`dx run iseq_snpsift_filter -ivcf='test_FBN1.vcf.gz' -ioutput_name='test_FBN1_filtered' -icommand='(isHom( GEN[0] ) & isRef( GEN[0] ))'`
