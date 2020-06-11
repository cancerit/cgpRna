#!/usr/bin/env cwl-runner

class: CommandLineTool

id: "run-cgprna_defuse"

label: "cgpRna defuse"

cwlVersion: v1.0

doc: |
  ![build_status](https://quay.io/repository/wtsicgp/cgprna/status)
  A Docker container for the cgpRna mapping flow. See the [cgpRna](https://github.com/cancerit/cgpRna) website for more information.

  Please read the relevant [changes](https://github.com/cancerit/cgpRna/blob/dev/CHANGES.md) when upgrading.

  Parameters for a CWL definition are generally described in a json file, but parameters can be provided on the command line.

  To see the parameters descriptions please run: cwltool --tool-help path_to.cwl

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/wtsicgp/cgprna:2.6.1"
  - class: InlineJavascriptRequirement

hints:
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 10000
    outdirMin: 10000

inputs:
  reads:
    doc: "RAW read input, can be bam files and FastQ files (optionally gzip compressed)."
    type: 
      type: array
      items: File
      inputBinding:
        itemSeparator: ' '
    inputBinding:
      prefix: --input
      separate: true

  reference:
    type: File
    doc: "The Tophat fusion reference files bundled in a tar.gz."
    inputBinding:
      prefix: --reference
      separate: true

  sample_name:
    type: string
    doc: "Sample name, which will used to prefix output file names and SM tag in the BAM file header."
    inputBinding:
      prefix: --sample-name
      separate: true

  gene_build:
    type: string?
    doc: "gene build folder name, if specified it'll be used to locate the folder within the reference bundle."
    inputBinding:
      prefix: --gene-build
      separate: true

  threads:
    type: int?
    doc: "Number of threads to use."
    inputBinding:
      prefix: --threads
      separate: true

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).defuse-fusion.normals.ext.filtered.txt

baseCommand: ["run-cgprna", "defuse"]

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
