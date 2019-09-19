#!/usr/bin/env cwl-runner

class: CommandLineTool

id: "run-cgprna_map"

label: "cgpRna mapping"

cwlVersion: v1.0

doc: |
  ![build_status](https://quay.io/repository/wtsicgp/cgprna/status)
  A Docker container for the cgpRna mapping flow. See the [cgpRna](https://github.com/cancerit/cgpRna) website for more information.

  Please read the relevant [changes](https://github.com/cancerit/cgpRna/blob/dev/CHANGES.md) when upgrading.

  Parameters for a CWL definition are generally described in a json file, but parameters can be provided on the command line.

  To see the parameters descriptions please run: cwltool --tool-help path_to.cwl

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/wtsicgp/cgprna:2.3.4"

hints:
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 30000
    outdirMin: 20000

inputs:
  raw_reads:
    doc: "RAW read input, can be a bam file, or one of a pair of FastQ file. (optionally gzip compressed)."
    type:
      type: array
      items: File
    inputBinding:
      prefix: --input
      separate: true
      itemSeparator: ' '

  reference:
    type: File
    doc: "The core STAR reference and a GTF file bundled in a tar.gz."
    inputBinding:
      prefix: --reference
      separate: true

  sample_name:
    type: string
    doc: "Sample name, which will used to prefix output file names and SM tag in the BAM file header."
    default: ''
    inputBinding:
      prefix: --sample-name
      separate: true
      shellQuote: true

  threads:
    type: int?
    doc: "Number of threads to use."
    inputBinding:
      prefix: --threads
      separate: true
      shellQuote: true

  rg_id_tag:
    type: string?
    doc: "Readgroup ID tag value in the output BAM. Default: 1 or taken from the input raw BAM file."
    inputBinding:
      prefix: --rg-id-tag
      separate: true
      shellQuote: true

  lb_tag:
    type: string?
    doc: "Sequencing library tag value in the output BAM header. Default: None or taken from the input raw BAM file."
    inputBinding:
      prefix: --lb-tag
      separate: true
      shellQuote: true

  ds_tag:
    type: string?
    doc: "Description tag value in the output BAM header. Default: None or taken from the input raw BAM file."
    inputBinding:
      prefix: --ds-tag
      separate: true
      shellQuote: true

  pl_tag:
    type: string?
    doc: "Platform tag value in the output BAM header. Default: None or taken from the input raw BAM file."
    inputBinding:
      prefix: --pl-tag
      separate: true
      shellQuote: true

  pu_tag:
    type: string?
    doc: "Platform unit tag value in the output BAM header. Default: None or taken from the input raw BAM file."
    inputBinding:
      prefix: --pu-tag
      separate: true
      shellQuote: true

outputs:
  star_transcriptome_bam:
    type: File
    outputBinding:
     glob: $(inputs.sample_name).star.AlignedtoTranscriptome.out.bam
    secondaryFiles:
    - .bai

  dup_marked_bam:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).bam
    secondaryFiles:
    - .bai

  dup_marked_bam_dup_met:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).bam.met

  dup_marked_bam_md5:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).bam.md5

baseCommand: ["run-cgprna", "map"]

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
