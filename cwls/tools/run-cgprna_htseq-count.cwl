#!/usr/bin/env cwl-runner

class: CommandLineTool

id: "run-cgprna_htseq-count"

label: "cgpRna htseq-count flow"

cwlVersion: v1.0

doc: |
  ![build_status](https://quay.io/repository/wtsicgp/cgprna/status)
  A Docker container for the cgpRna mapping flow. See the [cgpRna](https://github.com/cancerit/cgpRna) website for more information.

  Please read the relevant [changes](https://github.com/cancerit/cgpRna/blob/dev/CHANGES.md) when upgrading.

  Parameters for a CWL definition are generally described in a json file, but parameters can be provided on the command line.

  To see the parameters descriptions please run: cwltool --tool-help path_to.cwl

requirements:
  - class: DockerRequirement
    dockerPull: "cgprna:2.3.4"

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 4000

inputs:
  sample_bam:
    type: File
    doc: "Input BAM file, in which reads are mapped to a reference genome (NOT transcriptome)."
    inputBinding:
      prefix: --input
      separate: true

  reference:
    type: File
    doc: "A reference GTF file."
    inputBinding:
      prefix: --reference
      separate: true

# outputs:
#   star_sample_bam:
#     type: File
#     outputBinding:
#       glob: $(inputs.sample_name).star.Aligned.out.bam

#   star_transcriptome_bam:
#     type: File
#     outputBinding:
#       glob: $(inputs.sample_name).star.AlignedtoTranscriptome.out.bam

#   star_transcriptome_bam_index:
#     type: File
#     outputBinding:
#       glob: $(inputs.sample_name).star.AlignedtoTranscriptome.out.bam.bai

#   mark_dup_bam:
#     type: File
#     outputBinding:
#       glob: $(inputs.sample_name).bam

#   mark_dup_bam_index:
#     type: File
#     outputBinding:
#       glob: $(inputs.sample_name).bam.bai

#   mark_dup_bam_dup_met:
#     type: File
#     outputBinding:
#       glob: $(inputs.sample_name).bam.met

#   mark_dup_bam_md5:
#     type: File
#     outputBinding:
#       glob: $(inputs.sample_name).bam.md5

baseCommand: ["run-cgprna", "count"]

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

$namespaces:
  s: http://schema.org/

s:codeRepository: https://github.com/cancerit/cgpRna
s:license: https://spdx.org/licenses/AGPL-3.0-only

s:author:
  - class: s:Person
    s:email: mailto:cgphelp@sanger.ac.uk
    s:name: Yaobo Xu
