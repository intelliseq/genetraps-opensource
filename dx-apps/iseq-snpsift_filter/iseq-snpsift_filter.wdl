task iseq_snpsift_filter {

File vcf
String command
String output_name

command <<<

zcat -f ${vcf} | java -jar /tmp/snpEff/SnpSift.jar filter "${command}" > ${output_name}.vcf

>>>

runtime {
        docker: "alexcoppe/snpsift"
	}

output {
     File result = "${output_name}.vcf"
	}

}

