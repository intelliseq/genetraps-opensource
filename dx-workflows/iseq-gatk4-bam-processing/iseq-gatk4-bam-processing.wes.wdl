#   --------------------------------------------
#   PostProcessingForVariantDiscovery_GATK4_WES
#   --------------------------------------------
#
#   Katarzyna Kolanek
#   15/06/2018
#
#   This workflow performs post-processing of aligned, sorted and deduped BAM
#   for variant discovery with GATK4. The workflow is a modified version of
#   the followinng GATK4 workfow:
#
#   https://github.com/gatk-workflows/gatk4-data-processing
#
#   The initial steps from original workflow were moved to another workflow
#   (iseq-bwa-mem, IseqBwaMem) and this workflow performs only Base (Quality
#   Score) Recalibration,
#
#   The workflow is optimized to process WES data on DNAnexus.
#


#  WORKFLOW DEFINITION
#  *******************

workflow PostProcessingForVariantDiscovery_GATK4_WES {

  File input_bam

  String sample_id
  String ref_name

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File sequence_grouping_file
  File sequence_grouping_with_unmapped_file

  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  String gotc_docker = "broadinstitute/genomes-in-the-cloud:2.3.0-1501082129"
  String picard_docker = "broadinstitute/genomes-in-the-cloud:2.3.0-1501082129"
  String gatk_docker = "broadinstitute/gatk:4.0.0.0"
  String python_docker = "python:2.7"

  String gotc_path = "/usr/gitc/"
  String picard_path =  "/usr/gitc/"
  String gatk_path = "/gatk/gatk"

  Int compression_level = 5

  String BaseRecalibrator_java_opt = "-Xms4000m"
  String BaseRecalibrator_dnanexus_instance_type = "mem1_ssd1_x4"
  # "-Xms4000m", "6 GB", 200
  #  Memory: 7.5, Cores: 4, Storage: 80


  String GatherBqsrReports_java_opt = "-Xms3000m"
  String GatherBqsrReports_dnanexus_instance_type = "mem1_ssd1_x2"
  # "-Xms3000m", "3500 MB", 100
  #  Memory: 3.8, Cores: 1, Storage: 32

  String ApplyBQSR_java_opt = "-Xms3000m"
  String ApplyBQSR_dnanexus_instance_type = "mem1_ssd2_x2"
  # "-Xms3000m", "3500 MB", 200
  #  Memory: 3.8, Cores: 2, Storage: 80

  String GatherBamFiles_java_opt = "-Xms2000m"
  String GatherBamFiles_dnanexus_instance_type = "mem2_hdd2_x1"
  # "-Xms2000m", "3 GB", 400
  #  Memory: 3.8, Cores: 1, Storage: 410

  String base_file_name = sample_id + "." + ref_name


  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel

  Array[Array[String]] sequence_grouping = read_tsv(sequence_grouping_file)
  Array[Array[String]] sequence_grouping_with_unmapped = read_tsv(sequence_grouping_with_unmapped_file)

  scatter (index in range(length(sequence_grouping))) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:

        index = index,

        input_bam = input_bam,

        recalibration_report_filename = base_file_name + ".recal_data.csv",

        sequence_group_interval = sequence_grouping[index],

        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,

        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,

        docker_image = gatk_docker,
        gatk_path = gatk_path,

        BaseRecalibrator_dnanexus_instance_type = BaseRecalibrator_dnanexus_instance_type,
        base_file_name = base_file_name,

        java_opt = BaseRecalibrator_java_opt
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports {
    input:

      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = base_file_name + ".recal_data.csv",

      base_file_name = base_file_name,

      docker_image = gatk_docker,
      gatk_path = gatk_path,
      GatherBqsrReports_dnanexus_instance_type = GatherBqsrReports_dnanexus_instance_type,

      java_opt = GatherBqsrReports_java_opt
  }

  scatter (index in range(length(sequence_grouping_with_unmapped))) {

    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:

      index = index,

      input_bam = input_bam,

      output_bam_basename = base_file_name + ".recalibrated",

      base_file_name = base_file_name,

      recalibration_report = GatherBqsrReports.output_bqsr_report,
      sequence_group_interval = sequence_grouping_with_unmapped[index],

      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,

      docker_image = gatk_docker,
      gatk_path = gatk_path,

      ApplyBQSR_dnanexus_instance_type = ApplyBQSR_dnanexus_instance_type,

      java_opt = ApplyBQSR_java_opt
  }
}

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:

      input_bams = ApplyBQSR.recalibrated_bam,

      output_bam_basename = base_file_name,

      docker_image = picard_docker,
      picard_path = picard_path,

      GatherBamFiles_dnanexus_instance_type = GatherBamFiles_dnanexus_instance_type,

      compression_level = compression_level,
      java_opt = GatherBamFiles_java_opt
  }

  # Outputs that will be retained when execution is complete
  output {

    Array[File] BaseRecalibrator_stdout_stderr_logs  = BaseRecalibrator.BaseRecalibrator_stdout_stderr_log
    File GatherBqsrReports_stdout_stderr_log  = GatherBqsrReports.GatherBqsrReports_stdout_stderr_log
    Array[File] ApplyBQSR_stdout_stderr_logs = ApplyBQSR.ApplyBQSR_stdout_stderr_log

    File bqsr_report = GatherBqsrReports.output_bqsr_report
    File analysis_ready_bam = GatherBamFiles.output_bam
    File analysis_ready_bam_index = GatherBamFiles.output_bam_index
    File analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5
  }

}

# TASK DEFINITIONS

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {

  Int index

  File input_bam

  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String base_file_name

  String java_opt

  String docker_image
  String gatk_path
  String BaseRecalibrator_dnanexus_instance_type

  command {

    samtools index ${input_bam}

    ${gatk_path} --java-options "${java_opt}" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${recalibration_report_filename} \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
      -L ${sep=" -L " sequence_group_interval} &> ${base_file_name}_${index}.BaseRecalibrator.stdout.stderr.log
  }
  runtime {
    docker: docker_image
    dx_instance_type: "${BaseRecalibrator_dnanexus_instance_type}"
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
    File BaseRecalibrator_stdout_stderr_log = "${base_file_name}_${index}.BaseRecalibrator.stdout.stderr.log"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
# Note that when run from GATK 3.x the tool is not a walker and is invoked differently.
task GatherBqsrReports {

  Array[File] input_bqsr_reports
  String output_report_filename

  String base_file_name

  String docker_image
  String gatk_path

  String GatherBqsrReports_dnanexus_instance_type

  String java_opt

  command {
    ${gatk_path} --java-options "${java_opt}" \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename} &> ${base_file_name}.GatherBqsrReports.stdout.stderr.log
    }
  runtime {
    docker: docker_image
    dx_instance_type: "${GatherBqsrReports_dnanexus_instance_type}"
  }
  output {
    File output_bqsr_report = "${output_report_filename}"
    File GatherBqsrReports_stdout_stderr_log = "${base_file_name}.GatherBqsrReports.stdout.stderr.log"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {

  Int index

  File input_bam

  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String base_file_name

  String java_opt

  String docker_image
  String gatk_path
  String ApplyBQSR_dnanexus_instance_type


  command {

    samtools index ${input_bam}

    ${gatk_path} --java-options "${java_opt}" \
      ApplyBQSR \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${output_bam_basename}.bam \
      -L ${sep=" -L " sequence_group_interval} \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities &> ${base_file_name}_${index}.ApplyBQSR.stdout.stderr.log
  }
  runtime {
    docker: docker_image
    dx_instance_type: "${ApplyBQSR_dnanexus_instance_type}"
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
    File ApplyBQSR_stdout_stderr_log= "${base_file_name}_${index}.ApplyBQSR.stdout.stderr.log"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {

  Array[File] input_bams
  String output_bam_basename

  Int compression_level

  String java_opt

  String GatherBamFiles_dnanexus_instance_type

  String docker_image
  String picard_path

  command {
    java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar ${picard_path}picard.jar \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true
    }
  runtime {
    docker: docker_image
    dx_instance_type: "${GatherBamFiles_dnanexus_instance_type}"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}
