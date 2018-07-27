#   --------------------------
#   IseqCreateSequenceGrouping
#   --------------------------
#
#   Katarzyna Kolanek
#   30/05/2018
#
#   This is a workflow containg a task "CreateSequenceGroupingTSV" taken from https://github.com/gatk-workflows/gatk4-data-processing
#   The workfow generates {ref_name}.sequence_grouping.txt and {ref_name}.sequence_grouping_with_unmapped.txt files for a specified
#   reference genome ({ref_name}). Those files are reqired by iseq-gatk4-bam-processing workflow.


# WORKFLOW DEFINITION
# *******************

workflow IseqCreateSequenceGrouping {

  String ref_name
  File ref_dict

  String python_docker = "python:2.7"

  # Create list of sequences for scatter-gather parallelization
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict,
      ref_name = ref_name,
      docker_image = python_docker
  }

  output {
    File sequence_grouping_file = CreateSequenceGroupingTSV.sequence_grouping_file
    File sequence_grouping_with_unmapped_file = CreateSequenceGroupingTSV.sequence_grouping_with_unmapped_file
  }

}


# TASKS DEFINITIONS
# *****************

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {

  File ref_dict
  String ref_name

  String docker_image

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("${ref_name}.sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("${ref_name}.sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    docker: docker_image
    dx_instance_type: "mem1_ssd1_x2"
  }
  output {
    File sequence_grouping_file = "${ref_name}.sequence_grouping.txt"
    File sequence_grouping_with_unmapped_file = "${ref_name}.sequence_grouping_with_unmapped.txt"
  }
}
