#   ----------
#   IseqBwaMem
#   ----------
#
#   Katarzyna Kolanek
#   06/06/2018
#
#   This workflow takes an array of left and an array of right Illumina paired-end reads in FASTQ format,
#   belonging to the same sample.
#   The workflow performs the following actions:
#     + maps reads to a reference genome with the BWA-MEM algorithm
#     + adds read groups information to resulting bam
#     + sorts reads by coordinate
#     + marks duplicated reads
#
#
#   Read Groups mini-guide:
#   https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
#



#  WORKFLOW DEFINITION
#  *******************

workflow IseqBwaMem {

  String sample_id = "001"
  String ref_name = "grch38"

  String RG_PL = "Illumina"

  Array[File] input_fastqs_1
  Array[File] input_fastqs_2

  File? ref_fasta = "Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa"
  File? ref_fasta_index = "Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.fai"
  File? ref_dict = "Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.dict"
  File? ref_alt = "Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.alt"
  File? ref_sa = "Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.sa"
  File? ref_amb = "Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.amb"
  File? ref_bwt = "Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.bwt"
  File? ref_ann = "Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.ann"
  File? ref_pac = "Intelliseq - Resources:Reference genomes/GRCh38-no-alt-analysis-set/GRCh38.no_alt_analysis_set.fa.pac"

  String bwa_commandline = "-K 100000000 -v 3 -Y $bash_ref_fasta"
  Int compression_level = 5

  String bwa_docker = "lltw/bwa:latest"
  String sambamba_docker = "lltw/sambamba"

  String BwaMem_num_cpu = "8"

  String MergeBams_num_cpu = "8"

  String base_file_name = sample_id + "." + ref_name

  # Map reads to reference, mark duplicates ans sort

  scatter (index in range(length(input_fastqs_1))) {

    call BwaMem {
      input:

        fastq_1 = input_fastqs_1[index],
        fastq_2 = input_fastqs_2[index],

        sample_id = sample_id,
        index = index,

        RG_PL = RG_PL,

        base_file_name = base_file_name,

        bwa_commandline = bwa_commandline,

        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_alt = ref_alt,
        ref_sa = ref_sa,
        ref_amb = ref_amb,
        ref_bwt = ref_bwt,
        ref_ann = ref_ann,
        ref_pac = ref_pac,

        compression_level = compression_level,

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


# TASKS DEFINITIONS
# *****************

task BwaMem {

  File fastq_1
  File fastq_2

  String sample_id

  Int index

  String bwa_commandline
  String base_file_name

  String RG_PL

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File? ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  Int compression_level
  String num_cpu

  String docker_image

  command <<<

    RG_ID=$(zcat ${fastq_1} | head -1 | cut -d ':' -f 3,4 | sed 's/:/\./g')
    RG_PU="$RG_ID"".""${sample_id}"
    RG_LB="${sample_id}"".library"
    RG_SM="${sample_id}"

    bwa 2>&1 | grep -e '^Version' | sed 's/Version: //' > version.txt

    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    bwa mem -t ${num_cpu} -R "@RG\tID:""$RG_ID""\tPU:""$RG_PU""\tPL:${RG_PL}\tLB:""$RG_LB""\tSM:""$RG_SM" ${bwa_commandline} ${fastq_1} ${fastq_2} 2> ${base_file_name}_${index}.bwa.stderr.log \
    | samblaster 2> ${base_file_name}_"${index}".samblaster.stderr.log \
    | samtools sort -@ ${num_cpu} - > ${base_file_name}_"${index}".aln.bam 2> ${base_file_name}_${index}.samtools-sort.stderr.log

  >>>

  runtime {
    docker: docker_image
    dx_instance_type: "mem1_ssd1_x8"
  }


  output {
    File output_bam = "${base_file_name}_${index}.aln.bam"
    File bwa_stderr_log = "${base_file_name}_${index}.bwa.stderr.log"
    File deduplication_stderr_log = "${base_file_name}_${index}.samblaster.stderr.log"
    File samtools_stderr_log = "${base_file_name}_${index}.samtools-sort.stderr.log"
    Array[String] version = read_lines("version.txt")
  }
}

task MergeBams {

  Array[File] aln_bams

  String base_file_name

  Int compression_level
  String disk_size
  String mem_size
  String num_cpu

  String docker_image

  Int number_of_bam_pieces = length(aln_bams)
  File first_aln_bam = aln_bams[0]

  command <<<

    if [ "${number_of_bam_pieces}" -eq "1" ]
    then
    mv  ${sep=' ' aln_bams} ${base_file_name}.aln.bam
    sambamba index -t ${num_cpu} ${base_file_name}.aln.bam 2>  ${base_file_name}.sambamba.stderr.log
    else

      sambamba merge -t ${num_cpu} -l ${compression_level} ${base_file_name}.aln.bam ${sep=' ' aln_bams} 2>  ${base_file_name}.sambamba.stderr.log
    fi

  >>>

  runtime {
    docker: docker_image
    dx_instance_type: "mem1_ssd1_x8"
  }

  output {
    File output_bam = "${base_file_name}.aln.bam"
    File output_bam_bai = "${base_file_name}.aln.bam.bai"
    File sambamba_stderr_log = "${base_file_name}.sambamba.stderr.log"
  }

}
