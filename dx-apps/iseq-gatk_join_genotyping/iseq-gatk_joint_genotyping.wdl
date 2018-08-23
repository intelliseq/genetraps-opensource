task iseq_gatk_joint_genotyping {

Array[File] GVCFs
String output_name
File reference

command <<<
tar -xvf ${reference}

gatk CombineGVCFs -R `pwd`"/GRCh38.no_alt_analysis_set.fa" -V ${sep=" -V " GVCFs} -O cohort.g.vcf.gz
gatk --java-options "-Xmx4g" GenotypeGVCFs -R `pwd`"/GRCh38.no_alt_analysis_set.fa" -V cohort.g.vcf.gz -O ${output_name}.vcf.gz

>>>

runtime {
        docker: "broadinstitute/gatk"
	}

output {
     File result = "${output_name}.vcf.gz"
	}

}
