task iseq_tar {
    String? name = "archive"
    Array[File] files
    command <<<
      LIST_OF_FILES=`echo ${sep=' ' files}`
      dx-docker run -v /home/dnanexus/execution/inputs:/home/dnanexus/execution/inputs intelliseq/tar:latest -cf /home/dnanexus/execution/inputs/${name}.tar $LIST_OF_FILES
    >>>
    output {
      File archive = "/home/dnanexus/execution/inputs/${name}.tar"
    }
}
