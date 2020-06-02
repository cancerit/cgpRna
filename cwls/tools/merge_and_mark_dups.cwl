#!/usr/bin/env cwl-runner

class: CommandLineTool

id: "merge-mark-dups"

label: "merge BAMs and mark dupliates"

cwlVersion: v1.0

doc: |
  ![build_status](https://quay.io/repository/wtsicgp/cgprna/status)
  A Docker container for the cgpRna mapping flow. See the [cgpRna](https://github.com/cancerit/cgpRna) website for more information.

  Please read the relevant [changes](https://github.com/cancerit/cgpRna/blob/dev/CHANGES.md) when upgrading.

  Parameters for a CWL definition are generally described in a json file, but parameters can be provided on the command line.

  To see the parameters descriptions please run: cwltool --tool-help path_to.cwl

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/wtsicgp/cgprna:2.5.0"

hints:
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 8000

inputs:
  sorted_bams:
    doc: "BAM files to be merged."
    type:
      type: array
      items: File
      inputBinding:
        prefix: I=
        separate: false

  threads:
    type: int?
    doc: "Number of threads to use."
    inputBinding:
      prefix: markthreads=
      separate: false

  out_bam_name:
    type: string?
    default: 'out_merged.bam'
    inputBinding:
      prefix: O=
      separate: false

  out_bam_index_name:
    type: string?
    doc: "if specified, make sure it matches '{out_bam_name}.bai'."
    default: 'out_merged.bam.bai'
    inputBinding:
      prefix: indexfilename=
      separate: false

  md5_file_name:
    type: string?
    default: 'out_merged.bam.md5'
    inputBinding:
      prefix: md5filename=
      separate: false

  dup_met_file_name:
    type: string?
    default: 'out_merged.bam.met'
    inputBinding:
      prefix: M=
      separate: false

outputs:
  dup_marked_merged_bam:
    type: File
    outputBinding:
      glob: $(inputs.out_bam_name)
    secondaryFiles:
    - .bai

  dup_marked_bam_dup_met:
    type: File
    outputBinding:
      glob: $(inputs.md5_file_name)

  dup_marked_bam_md5:
    type: File
    outputBinding:
      glob: $(inputs.dup_met_file_name)

baseCommand: ["bammarkduplicates2", "md5=1", "index=1"]

$schemas:
  - https://schema.org/version/latest/schema.rdf

$namespaces:
  s: http://schema.org/

s:codeRepository: https://github.com/cancerit/cgpRna
s:license: https://spdx.org/licenses/AGPL-3.0-only

s:author:
  - class: s:Person
    s:email: mailto:cgphelp@sanger.ac.uk
    s:name: Yaobo Xu
