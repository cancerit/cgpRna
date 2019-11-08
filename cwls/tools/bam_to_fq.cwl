#!/usr/bin/env cwl-runner

class: CommandLineTool

id: "bam-to-fq"

label: "convert bam to unzipped fq"

cwlVersion: v1.0

doc: |
  ![build_status](https://quay.io/repository/wtsicgp/cgprna/status)
  A Docker container for the cgpRna mapping flow. See the [cgpRna](https://github.com/cancerit/cgpRna) website for more information.

  Please read the relevant [changes](https://github.com/cancerit/cgpRna/blob/dev/CHANGES.md) when upgrading.

  Parameters for a CWL definition are generally described in a json file, but parameters can be provided on the command line.

  To see the parameters descriptions please run: cwltool --tool-help path_to.cwl

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/wtsicgp/cgprna:2.4.1"
  - class: InlineJavascriptRequirement

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 8000

inputs:
  in_bam:
    doc: "input BAM."
    type: File
    inputBinding:
      prefix: filename=
      separate: false

  output_prefix:
    type: string?
    default: 'bamtofastq_converted'

  matched_out_name:
    type: string
    doc: "output name of the fastq file containing first mate of matched pairs in the BAM."
    default: ""
    inputBinding:
      valueFrom: |
        ${
          var extension = "matched.1_1.fq";
          var eles;
          if (self) {
            eles = [inputs.output_prefix, self, extension]
          } else {
            eles = [inputs.output_prefix, extension]
          }
          return eles.join(".")
        }
      prefix: F=
      separate: false

  matched_2_out_name:
    type: string
    doc: "output name of the fastq file containing second mate of matched pairs in the BAM."
    default: ""
    inputBinding:
      valueFrom: |
        ${
          var extension = "matched.1_2.fq";
          var eles;
          if (self) {
            eles = [inputs.output_prefix, self, extension]
          } else {
            eles = [inputs.output_prefix, extension]
          }
          return eles.join(".")
        }
      prefix: F2=
      separate: false

  single_end_out_name:
    type: string
    doc: "output name of the fastq file containing single-ended reads in the BAM."
    default: ""
    inputBinding:
      valueFrom: |
        ${
          var extension = "singe_ended.fq";
          var eles;
          if (self) {
            eles = [inputs.output_prefix, self, extension]
          } else {
            eles = [inputs.output_prefix, extension]
          }
          return eles.join(".")
        }
      prefix: S=
      separate: false

  unmatcted_out_name:
    type: string
    doc: "output name of the fastq file containing first mate of unmatched pairs in the BAM."
    default: ""
    inputBinding:
      valueFrom: |
        ${
          var extension = "unmatched.1_1.fq";
          var eles;
          if (self) {
            eles = [inputs.output_prefix, self, extension]
          } else {
            eles = [inputs.output_prefix, extension]
          }
          return eles.join(".")
        }
      prefix: O=
      separate: false

  unmatcted_2_out_name:
    type: string
    doc: "output name of the fastq file containing second mate of unmatched pairs in the BAM."
    default: ""
    inputBinding:
      valueFrom: |
        ${
          var extension = "unmatched.1_2.fq";
          var eles;
          if (self) {
            eles = [inputs.output_prefix, self, extension]
          } else {
            eles = [inputs.output_prefix, extension]
          }
          return eles.join(".")
        }
      prefix: O2=
      separate: false

outputs:
  output_fqs:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.matched.1_[12].fq"

baseCommand: ["bamtofastq", "exclude=SECONDARY,SUPPLEMENTARY"]

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
