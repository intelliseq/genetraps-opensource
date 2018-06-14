task BwaMem {

  File fastq_1
  File fastq_2

  String? sample_id = "emptyid"

  Int? index = "0"

  String? base_file_name = "basename"

  String? RG_PL = "ILLUMINA"

  File? ref_fasta = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa"
  File? ref_fasta_index = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.fai"
  File? ref_dict = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.dict"
  File? ref_alt
  File? ref_sa = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.sa"
  File? ref_amb = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.amb"
  File? ref_bwt = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.bwt"
  File? ref_ann = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.ann"
  File? ref_pac = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.pac"

  #String? bwa_commandline = "-K 100000000 -v 3 -Y"

  String docker_image
  Int? num_cpu = "2"

  command <<<

    echo "--- LOG: Setting environmental variables ---"
    echo "--- LOG: Setting environmental variables ---" | cat
    echo "--- LOG: Setting environmental variables ---" >&2
    RG_ID=$(gzip -cd ${fastq_1} | head -1 | cut -d ':' -f 3,4 | sed 's/:/\./g')
    RG_PU="$RG_ID"".""${sample_id}"
    RG_LB="${sample_id}"".library"
    RG_SM="${sample_id}"

    echo "--- LOG: Getting bwa version ---" >&2
    mkdir /tmp/results
    dx-docker run -v /tmp/results:/tmp/results --entrypoint "/bin/sh" intelliseq/bwa:latest -c "bwa 2>&1 | grep -e '^Version' | sed 's/Version: //' > /tmp/results/version.txt"


    #set -o -e -x pipefail
    echo "--- LOG: bwa ---" >&2
    #set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    dx-docker run -v /tmp/results:/tmp/results -v /home/dnanexus/execution/inputs:/home/dnanexus/execution/inputs --entrypoint "/bin/sh" intelliseq/bwa:latest -c "
    bwa mem -t ${num_cpu} -R "@RG\tID:""$RG_ID""\tPU:""$RG_PU""\tPL:${RG_PL}\tLB:""$RG_LB""\tSM:""$RG_SM" -K 100000000 -p -v 3 $bash_ref_fasta ${fastq_1} ${fastq_2}  
    | samblaster 2> /tmp/results/${base_file_name}_'${index}'.samblaster.stderr.log
    | samtools sort -@ ${num_cpu} - > /tmp/results/${base_file_name}_'${index}'.aln.bam
    "
    #2> {base_file_name}_"{index}".samblaster.stderr.log \
    #| samtools sort -@ ${num_cpu} - > {base_file_name}_"{index}".aln.bam 2> {base_file_name}_{index}.samtools-sort.stderr.log
    echo "--- LOG: printing results ---" >&2
    ls /tmp/results | cat >&2
    echo "--- LOG: finito ---" >&2

    #bwa mem -t ${num_cpu} -R "@RG\tID:""$RG_ID""\tPU:""$RG_PU""\tPL:${RG_PL}\tLB:""$RG_LB""\tSM:""$RG_SM" ${fastq_1} ${fastq_2} 2> ${base_file_name}_${index}.bwa.stderr.log \
    #| samblaster 2> ${base_file_name}_"${index}".samblaster.stderr.log \
    #| samtools sort -@ ${num_cpu} - > ${base_file_name}_"${index}".aln.bam 2> ${base_file_name}_${index}.samtools-sort.stderr.log

  >>>

  runtime {
    #docker: docker_image
    #dx_instance_type: "mem1_hdd1_x1"
  }


  output {
    #File bwalog = "/tmp/results/stdout.log"
    #File bwaerr = "/tmp/results/stderr.log"
    #File output_bam = "/tmp/results/my.bam"
    #File output_bam = "${base_file_name}_${index}.aln.bam"
    #File bwa_stderr_log = "${base_file_name}_${index}.bwa.stderr.log"
    File deduplication_stderr_log = "/tmp/results/${base_file_name}_${index}.samblaster.stderr.log"
    #File samtools_stderr_log = "${base_file_name}_${index}.samtools-sort.stderr.log"
    Array[String] version = read_lines("/tmp/results/version.txt")
  }
}

