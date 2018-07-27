workflow IseqBwaMem {

  String sample_id = "001"
  String ref_name = "grch38"

  String RG_PL = "Illumina"

  Array[File] fastq_1
  Array[File] fastq_2

  File ref_fasta = "dx://project-FBb60yQ0FYJ8p0bY32XX4Zz6:file-FBb62X00FYJ6GY7ZJZv05GJy"
  File? ref_fasta_index = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.fai"
  File? ref_dict = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.dict"
  File? ref_alt
  File? ref_sa = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.sa"
  File? ref_amb = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.amb"
  File? ref_bwt = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.bwt"
  File? ref_ann = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.ann"
  File? ref_pac = "dx://Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.pac"

  Int compression_level = 5

  String bwa_docker = "intelliseq/bwa:latest"
  String sambamba_docker = "intelliseq/bwa:latest"

  String BwaMem_num_cpu = "8"

  String MergeBams_num_cpu = "8"

  String base_file_name = sample_id + "." + ref_name

  # Map reads to reference, mark duplicates ans sort
  #index = 0
  scatter (index in range(length(fastq_1))) {

    call BwaMem {
      input:

        fastq_1 = fastq_1[index],
        fastq_2 = fastq_2[index],

        sample_id = sample_id,
        index = index,

        RG_PL = RG_PL,

        base_file_name = base_file_name,

        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_alt = ref_alt,
        ref_sa = ref_sa,
        ref_amb = ref_amb,
        ref_bwt = ref_bwt,
        ref_ann = ref_ann,
        ref_pac = ref_pac,

        #compression_level = compression_level,

        num_cpu = BwaMem_num_cpu,

        docker_image = bwa_docker
     }
   }

   call MergeBams {

     input:

      aln_bams = BwaMem.output_bam,

      base_file_name = base_file_name,
      compression_level = compression_level,

      num_cpu = MergeBams_num_cpu,

      docker_image = sambamba_docker
   }

   String version = BwaMem.version[0][0]

   output {

     File output_bam = MergeBams.output_bam
     File output_bam_bai = MergeBams.output_bam_bai

     Array[File] bwa_stderr_log = BwaMem.bwa_stderr_log
     Array[File] deduplication_stderr_log = BwaMem.deduplication_stderr_log
     Array[File] samtools_stderr_log = BwaMem.samtools_stderr_log
     File sambamba_stderr_log = MergeBams.sambamba_stderr_log
     String bwa_version = version

   }

}


task BwaMem {

  File fastq_1
  File fastq_2

  String? sample_id = "emptyid"

  Int? index = "0"

  String? base_file_name = "basename"

  String? RG_PL = "ILLUMINA"

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File? ref_alt
  File ref_sa
  File ref_amb
  File ref_bwt
  File ref_ann
  File ref_pac

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
    dx-docker run -v /tmp/results:/tmp/results -v /home/dnanexus/execution/inputs:/home/dnanexus/execution/inputs --entrypoint "/bin/sh" intelliseq/bwa:latest -c "bwa mem -t ${num_cpu} -R '@RG\tID:''$RG_ID''\tPU:''$RG_PU''\tPL:${RG_PL}\tLB:''$RG_LB''\tSM:''$RG_SM' -K 100000000 -p -v 3 $bash_ref_fasta ${fastq_1} ${fastq_2} 2> /tmp/results/${base_file_name}_${index}.bwa.stderr.log | samblaster 2> /tmp/results/${base_file_name}_'${index}'.samblaster.stderr.log | samtools sort -@ ${num_cpu} - > /tmp/results/${base_file_name}_'${index}'.aln.bam 2> /tmp/results/${base_file_name}_${index}.samtools-sort.stderr.log"
    echo "--- LOG: printing results ---" >&2
    ls /tmp/results | cat >&2
    echo "--- LOG: finito ---" >&2

  >>>

  runtime {
    #docker: docker_image
    dx_instance_type: "mem2_ssd1_x4"
  }


  output {
    File output_bam = "/tmp/results/${base_file_name}_${index}.aln.bam"
    File bwa_stderr_log = "/tmp/results/${base_file_name}_${index}.bwa.stderr.log"
    File deduplication_stderr_log = "/tmp/results/${base_file_name}_${index}.samblaster.stderr.log"
    File samtools_stderr_log = "/tmp/results/${base_file_name}_${index}.samtools-sort.stderr.log"
    Array[String] version = read_lines("/tmp/results/version.txt")
  }
}

task MergeBams {

  Array[File] aln_bams

  String base_file_name

  Int compression_level
  String num_cpu

  String docker_image

  Int number_of_bam_pieces = length(aln_bams)
  File first_aln_bam = aln_bams[0]

  command <<<

    mkdir /tmp/results
    echo "--- base file name: " $base_file_name
    if [ "${number_of_bam_pieces}" -eq "1" ]
    then
      mv  ${sep=' ' aln_bams} /home/dnanexus/execution/inputs/${base_file_name}.aln.bam
      echo "--- indexing..."
      dx-docker run hello-world
      dx-docker run -v /tmp/results:/tmp/results -v /home/dnanexus/execution/inputs:/home/dnanexus/execution/inputs --entrypoint "/bin/sh" intelliseq/sambamba:latest -c "sambamba index -t ${num_cpu} /home/dnanexus/execution/inputs/${base_file_name}.aln.bam 2> /tmp/results/${base_file_name}.sambamba.stderr.log"
      echo '{"error": {"type": "AppInternalError", "message": "Error while running sambamba index"}}' > job_error.json
    else
      dx-docker run -v /tmp/results:/tmp/results -v /home/dnanexus/execution/inputs:/home/dnanexus/execution/inputs --entrypoint "/bin/sh" intelliseq/sambamba:latest -c "sambamba merge -t ${num_cpu} -l ${compression_level} /home/dnanexus/execution/inputs/${base_file_name}.aln.bam ${sep=' ' aln_bams} 2>  /tmp/results/${base_file_name}.sambamba.stderr.log"
      echo '{"error": {"type": "AppInternalError", "message": "Error while running sambamba merge"}}' > job_error.json
    fi
    ls /home/dnanexus/execution/inputs/
  >>>

  runtime {
    dx_instance_type: "mem1_ssd1_x8"
  }

  output {
    File output_bam = "/home/dnanexus/execution/inputs/${base_file_name}.aln.bam"
    File output_bam_bai = "/home/dnanexus/execution/inputs/${base_file_name}.aln.bam.bai"
    File sambamba_stderr_log = "/tmp/results/${base_file_name}.sambamba.stderr.log"
  }

}
