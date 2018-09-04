task iseq_snpsift_annotate{

File vcf
File database
String output_name

command <<<

input_dir=`dirname ${vcf}`
dbn=`basename ${database}`
echo '### DOCKER TABIX LOG ###' > ${output_name}_annotate_log.txt
docker run --rm -v `pwd`:`pwd` -v $input_dir:$input_dir -w $input_dir bionode/htslib tabix -f ${database} > $input_dir/"$dbn".tbi 2>>${output_name}_annotate_log.txt
echo '### DOCKER SNPSIFT LOG ###' >> ${output_name}_annotate_log.txt
docker run --rm -v `pwd`:`pwd` -v $input_dir:$input_dir -w $input_dir alexcoppe/snpsift annotate ${database} ${vcf} > ${output_name}.vcf 2>>${output_name}_annotate_log.txt

>>>

runtime {
	dx_instance_type: "mem1_ssd2_x2"
}

output {
	File result = "${output_name}.vcf"
	File log = "${output_name}_annotate_log.txt"
	}

}

