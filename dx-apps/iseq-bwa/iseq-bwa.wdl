task iseq_fastqc {
    File fastq_file_1
    File fastq_file_2
    File reference
    command <<<
      tar reference
      
    >>>
    runtime {
    }
    output {
    }
}
