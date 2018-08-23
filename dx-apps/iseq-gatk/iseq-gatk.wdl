task iseq_gatk_haplotype_caller {

    File bam
    File reference
    Array[File] vcfs
    String? sample_id = "emptyid"
    String? interval
    File? intervalFile

    command <<<

    INTERVAL=""
    if [[ -n "${interval}" ]]; then INTERVAL="-L ${interval}" ; fi
    if [[ -n "${intervalFile}" ]]; then INTERVAL="-L ${intervalFile}" ; fi

    set -x
    KNOWN_SITES=`echo ${sep=' ' vcfs} | sed 's/ / --known-sites /g' | sed -e 's/^/--known-sites /'`
    REFERENCE_ARCHIVE=${reference}
    REFERENCE_DIRECTORY=$(dirname $REFERENCE_ARCHIVE)
    INPUTS="/home/dnanexus/execution/inputs"
    tar -xvf ${reference}

    mkdir /tmp/results

    for file in ${sep=' ' vcfs}  ; do
        dx-docker run -v $INPUTS:$INPUTS -v `pwd`:/genome --entrypoint /bin/bash broadinstitute/gatk:4.0.8.0 -c "gatk IndexFeatureFile -F $file" 
    done
    dx-docker run -v $INPUTS:$INPUTS -v `pwd`:/genome --entrypoint /bin/bash broadinstitute/gatk:4.0.8.0 -c "gatk BaseRecalibrator $KNOWN_SITES -R /genome/GRCh38.no_alt_analysis_set.fa -I ${bam} -O /genome/req.table"
    dx-docker run -v $INPUTS:$INPUTS -v `pwd`:/genome --entrypoint /bin/bash broadinstitute/gatk:4.0.8.0 -c "gatk ApplyBQSR -R /genome/GRCh38.no_alt_analysis_set.fa -I ${bam} --bqsr-recal-file /genome/req.table -O /genome/recalibrated.bam"
    dx-docker run -v $INPUTS:$INPUTS -v `pwd`:/genome --entrypoint /bin/bash broadinstitute/gatk:4.0.8.0 -c "gatk --java-options '-Xmx16G' HaplotypeCaller $INTERVAL -R /genome/GRCh38.no_alt_analysis_set.fa -I /genome/recalibrated.bam -ERC GVCF -O /genome/emptyid.g.vcf"
  
    >>>
    runtime {
      dx_instance_type: "mem2_ssd1_x4"
    }
    output {
      File output_gatk = "emptyid.g.vcf"
    }
}
