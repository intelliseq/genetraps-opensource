task iseq_bwa {

    File fastq_file_1
    File fastq_file_2
    File reference
    String? RG_PL = "ILLUMINA"
    String? sample_id = "emptyid"

    command <<<

    set -x

    RG_ID=$(gzip -cd ${fastq_file_1} | head -1 | cut -d ':' -f 3,4 | sed 's/:/\./g')
    RG_PU="$RG_ID"".""${sample_id}"
    RG_LB="${sample_id}"".library"
    RG_SM="${sample_id}"

    REFERENCE_ARCHIVE=${reference}
    REFERENCE_DIRECTORY=$(dirname $REFERENCE_ARCHIVE)

    tar -xvf ${reference}
    mkdir /tmp/results

    echo "==== REFERENCE DIR ===="
    ls -la $REFERENCE_DIRECTORY
    echo "==== PWD ===="
    ls -la `pwd`

    dx-docker run -v /tmp/results:/tmp/results -v `pwd`:`pwd` -v /home/dnanexus/execution/inputs:/home/dnanexus/execution/inputs --entrypoint "/bin/sh" intelliseq/bwa:latest -c "bwa mem -t 4 -R '@RG\tID:''$RG_ID''\tPU:''$RG_PU''\tPL:${RG_PL}\tLB:''$RG_LB''\tSM:''$RG_SM' -K 100000000 -p -v 3 `pwd`"/GRCh38.no_alt_analysis_set.fa" ${fastq_file_1} ${fastq_file_2} 2> /tmp/results/${sample_id}.bwa.stderr.log | samblaster 2> /tmp/results/${sample_id}.samblaster.stderr.log | samtools sort -@ 4 - > /tmp/results/${sample_id}.aln.bam 2> /tmp/results/${sample_id}.samtools-sort.stderr.log"
  
    >>>
    runtime {
      dx_instance_type: "mem2_ssd1_x4"
    }
    output {
      File output_bam = "/tmp/results/${sample_id}.aln.bam"
      File bwa_log = "/tmp/results/${sample_id}.bwa.stderr.log"
    }
}
