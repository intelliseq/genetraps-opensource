task iseq_tar {
    String? name = "archive"
    Array[File] files
    command <<<
      LIST_OF_FILES=`echo ${sep=' ' files} | sed "s/ .*\// /g" | sed "s/.*\///g"`
      echo $LIST_OF_FILES
      ls -la /home/dnanexus/execution/inputs
      cd /home/dnanexus/execution/inputs
      tar -cf ${name}.tar *
      #dx-docker run -v /home/dnanexus/execution/inputs:/home/dnanexus/execution/inputs intelliseq/tar:latest -cf /home/dnanexus/execution/inputs/${name}.tar 
    >>>
    output {
      File archive = "/home/dnanexus/execution/inputs/${name}.tar"
    }
}
