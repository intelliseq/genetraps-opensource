task iseq_fastqc {
    File fastq_file
    command <<<
    mkdir /tmp/fastqc
    cp ${fastq_file} fastq_file.fq.gz
    /opt/FastQC/fastqc fastq_file.fq.gz --outdir="/tmp/fastqc/"
    cp /tmp/fastqc/*.html summary_html
    unzip -p /tmp/fastqc/*.zip */fastqc_data.txt > summary_txt
    >>>
    runtime {
        docker: "intelliseq/fastqc"
    }
    output {
        File summary_html = "summary_html"
        File summary_txt = "summary_txt"
    }
}
